#!/usr/bin/env python3
"""
ADSIM VTK Visualization Tool
============================

Interactive 3D visualization of ADSIM VTK output files with field-based coloring,
timestep animation, and profile extraction capabilities.

Features:
  • Auto-detection of mesh dimensionality (1D, 2D, 2D axisymmetric)
  • 1D line plots or 2D heatmaps depending on mesh type
  • Field selection and colormap customization
  • Interactive profile extraction (draw lines on 2D meshes)
  • Timestep animation and navigation
  • Support for water flow and reactive transport output

Usage:
    python visualize_vtk.py                    # Auto-discover all VTK files
    python visualize_vtk.py --project water    # Filter by project name
    python visualize_vtk.py --field Matric_Head --play  # Start with field and auto-play
"""

import os
import glob
import re
import numpy as np
import pyvista as pv
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
import argparse
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize


class ADSIMVTKVisualizer:
    """Main visualization class for ADSIM VTK output."""

    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.vtk_files = {}  # {project_name: [(time, filepath, mesh_type), ...]}
        self.current_project = None
        self.current_timestep = 0
        self.current_field = None
        self.meshes = {}  # {(project, timestep): mesh}
        self.available_fields = {}  # {project: [field_names]}
        self.mesh_types = {}  # {project: 'line', '2d_cartesian', '2d_axisym'}
        self.plotter = None
        self.animation_running = False
        self.animation_speed = 1.0  # FPS
        
        self.discover_vtk_files()

    # ─────────────────────────────────────────────────────────────────────────
    # File Discovery & Organization
    # ─────────────────────────────────────────────────────────────────────────

    def discover_vtk_files(self):
        """Discover and organize all VTK files in output directory."""
        if not self.output_dir.exists():
            print(f"Error: Output directory '{self.output_dir}' not found")
            return

        vtk_files = sorted(self.output_dir.glob("*.vtk"))
        if not vtk_files:
            print(f"Warning: No VTK files found in {self.output_dir}")
            return

        # Organize by project name
        for vtk_file in vtk_files:
            match = re.match(r"(.+?)_(\d{6})\.vtk$", vtk_file.name)
            if match:
                project_name = match.group(1)
                timestep = int(match.group(2))

                if project_name not in self.vtk_files:
                    self.vtk_files[project_name] = []

                self.vtk_files[project_name].append((timestep, str(vtk_file)))

        # Sort by timestep
        for project in self.vtk_files:
            self.vtk_files[project].sort(key=lambda x: x[0])

        print(f"✓ Found {len(self.vtk_files)} projects with {sum(len(v) for v in self.vtk_files.values())} VTK files")
        for project, files in self.vtk_files.items():
            print(f"  • {project}: {len(files)} timesteps (0 to {files[-1][0]})")

    def get_projects(self) -> List[str]:
        """Return list of available projects."""
        return sorted(self.vtk_files.keys())

    def select_project(self, project_name: str):
        """Select a project and analyze its structure."""
        if project_name not in self.vtk_files:
            print(f"Error: Project '{project_name}' not found")
            return False

        self.current_project = project_name
        self.current_timestep = 0

        # Analyze first file to detect mesh type and fields
        _, first_file = self.vtk_files[project_name][0]
        mesh = pv.read(first_file)
        self.mesh_types[project_name] = self._detect_mesh_type(mesh)
        self.available_fields[project_name] = self._extract_field_names(mesh)

        print(f"\n✓ Selected project: {project_name}")
        print(f"  Mesh type: {self.mesh_types[project_name]}")
        print(f"  Available fields ({len(self.available_fields[project_name])}): ")
        for field in self.available_fields[project_name][:5]:
            print(f"    • {field}")
        if len(self.available_fields[project_name]) > 5:
            print(f"    ... and {len(self.available_fields[project_name]) - 5} more")

        if not self.current_field and self.available_fields[project_name]:
            self.current_field = self.available_fields[project_name][0]
            print(f"  Default field: {self.current_field}")

        return True

    # ─────────────────────────────────────────────────────────────────────────
    # Mesh Analysis
    # ─────────────────────────────────────────────────────────────────────────

    def _detect_mesh_type(self, mesh: pv.PolyData) -> str:
        """Detect if mesh is 1D line, 2D, or 2D axisymmetric based on coordinates."""
        coords = mesh.points
        if len(coords) == 0:
            return "unknown"

        # Get coordinate ranges
        x_range = np.ptp(coords[:, 0])
        y_range = np.ptp(coords[:, 1])
        z_range = np.ptp(coords[:, 2]) if coords.shape[1] > 2 else 0

        # Determine type based on extent in each dimension
        # Threshold: if range < 1e-6 * max(other ranges), consider it degenerate
        ranges = [x_range, y_range, z_range]
        max_range = max(ranges)
        threshold = 1e-6 * max_range if max_range > 0 else 1e-10

        active_dims = sum(1 for r in ranges if r > threshold)

        if active_dims == 1:
            return "1d_line"
        elif active_dims == 2:
            # Check if it looks like axisymmetric (radial-vertical)
            # Typically: x ≈ radius (0 to some value), y ≈ vertical
            if np.min(coords[:, 0]) >= -1e-10:  # x ≥ 0 suggests radial coordinate
                return "2d_axisymmetric"
            return "2d_cartesian"
        else:
            return "3d"

    def _extract_field_names(self, mesh: pv.PolyData) -> List[str]:
        """Extract all available scalar field names from mesh."""
        fields = []
        if mesh.point_data:
            fields.extend(mesh.point_data.keys())
        return sorted(set(fields))

    # ─────────────────────────────────────────────────────────────────────────
    # Data Loading
    # ─────────────────────────────────────────────────────────────────────────

    def load_timestep(self, timestep: int) -> pv.PolyData:
        """Load mesh for a specific timestep."""
        if self.current_project is None:
            print("Error: No project selected")
            return None

        if timestep >= len(self.vtk_files[self.current_project]):
            print(f"Error: Timestep {timestep} out of range")
            return None

        self.current_timestep = timestep
        _, filepath = self.vtk_files[self.current_project][timestep]

        cache_key = (self.current_project, timestep)
        if cache_key not in self.meshes:
            mesh = pv.read(filepath)
            self.meshes[cache_key] = mesh
        else:
            mesh = self.meshes[cache_key]

        return mesh

    def get_timestep_info(self, timestep: int) -> Dict:
        """Extract time info from VTK file header."""
        if self.current_project is None:
            return {}

        _, filepath = self.vtk_files[self.current_project][timestep]

        info = {"timestep": timestep, "time": 0.0}
        try:
            with open(filepath, "r") as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    # Second line typically has: "ADSIM v0.1.0 - Time Step 0, Time = 0.0"
                    match = re.search(r"Time\s*=\s*([0-9eE+.\-]+)", lines[1])
                    if match:
                        info["time"] = float(match.group(1))
        except:
            pass

        return info

    def _get_6_timestep_indices(self) -> List[int]:
        """Get 6 equally-spaced timestep indices: initial, final, and 4 in between."""
        n_steps = len(self.vtk_files[self.current_project])
        if n_steps <= 6:
            return list(range(n_steps))
        
        # Get indices: 0, quarter, half1, half2, three-quarter, final
        indices = [
            0,
            n_steps // 4,
            n_steps // 2 - 1,
            n_steps // 2,
            (3 * n_steps) // 4,
            n_steps - 1
        ]
        return sorted(list(set(indices)))  # Remove duplicates and sort

    def _extract_middle_vertical_profile(self, mesh: pv.PolyData) -> Tuple[np.ndarray, np.ndarray]:
        """Extract 1D profile along the middle of the vertical domain.
        
        Returns: (vertical_positions, field_values)
        """
        coords = mesh.points
        field_data = mesh[self.current_field]
        
        # Find middle in horizontal direction (x or r coordinate)
        x_coords = coords[:, 0]
        x_middle = (x_coords.min() + x_coords.max()) / 2
        
        # Find all points near the middle horizontally
        x_tol = (x_coords.max() - x_coords.min()) * 0.05  # 5% tolerance
        middle_mask = np.abs(x_coords - x_middle) < x_tol
        
        if not np.any(middle_mask):
            # If no points near middle, just take the points closest to middle
            dists = np.abs(x_coords - x_middle)
            middle_mask = dists < (np.min(dists) + 1e-6)
        
        # Sort by vertical coordinate (y or z)
        middle_indices = np.where(middle_mask)[0]
        z_coords = coords[middle_indices, 1]
        sort_idx = np.argsort(z_coords)
        
        vertical_positions = z_coords[sort_idx]
        profile_values = field_data[middle_indices[sort_idx]]
        
        return vertical_positions, profile_values

    # ─────────────────────────────────────────────────────────────────────────
    # Visualization
    # ─────────────────────────────────────────────────────────────────────────

    def plot_current(self, save_path: Optional[str] = None):
        """Plot current timestep with current field."""
        if self.current_project is None or self.current_field is None:
            print("Error: Project and field must be selected first")
            return

        mesh = self.load_timestep(self.current_timestep)
        if mesh is None:
            return

        mesh_type = self.mesh_types[self.current_project]
        time_info = self.get_timestep_info(self.current_timestep)

        # Create appropriate plot based on mesh type
        if mesh_type == "1d_line":
            self._plot_1d_line(mesh, time_info, save_path)
        elif mesh_type == "2d_cartesian":
            self._plot_2d_cartesian(mesh, time_info, save_path)
        elif mesh_type == "2d_axisymmetric":
            self._plot_2d_axisymmetric(mesh, time_info, save_path)
        else:
            self._plot_3d(mesh, time_info, save_path)

    def _plot_1d_line(self, mesh: pv.PolyData, time_info: Dict, save_path: Optional[str] = None):
        """Plot 1D line data (domain position vs field value)."""
        import matplotlib.pyplot as plt

        # Extract coordinates and field data
        coords = mesh.points
        field_data = mesh[self.current_field]

        # Sort by position (typically z or y)
        positions = coords[:, 2] if coords[:, 2].max() > coords[:, 1].max() else coords[:, 1]
        sort_idx = np.argsort(positions)

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(positions[sort_idx], field_data[sort_idx], "o-", linewidth=2, markersize=4)
        ax.set_xlabel("Position (m)")
        ax.set_ylabel(self.current_field)
        ax.set_title(
            f"{self.current_project}: {self.current_field} "
            f"(Step {time_info['timestep']}, t={time_info['time']:.3e})"
        )
        ax.grid(True, alpha=0.3)

        # Add statistics
        stats_text = f"min={field_data.min():.3e}\nmax={field_data.max():.3e}\nmean={field_data.mean():.3e}"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        if save_path:
            fig.savefig(save_path, bbox_inches="tight")
            print(f"✓ Saved plot to {save_path}")
        else:
            plt.tight_layout()
            plt.show()

    def _plot_2d_cartesian(self, mesh: pv.PolyData, time_info: Dict, save_path: Optional[str] = None):
        """Plot 2D Cartesian mesh: heatmap (left) + 6 time profiles (right)."""
        coords = mesh.points
        field_data = mesh[self.current_field]

        # Create figure with 2 subplots
        fig, (ax_2d, ax_profiles) = plt.subplots(1, 2, figsize=(16, 7))

        # === LEFT: 2D Heatmap ===
        norm = Normalize(vmin=field_data.min(), vmax=field_data.max())
        scatter = ax_2d.scatter(coords[:, 0], coords[:, 1], c=field_data, cmap="viridis", 
                               norm=norm, s=50, alpha=0.8)

        ax_2d.set_xlabel("X (m)")
        ax_2d.set_ylabel("Y (m)")
        ax_2d.set_title(
            f"{self.current_project}: {self.current_field} (2D Heatmap)\n"
            f"Step {time_info['timestep']}, t={time_info['time']:.3e}"
        )
        ax_2d.set_aspect("equal", adjustable="box")
        ax_2d.grid(True, alpha=0.2)

        # Colorbar
        cbar = fig.colorbar(scatter, ax=ax_2d, label=self.current_field)

        # Statistics box
        stats_text = (f"min={field_data.min():.3e}\n"
                     f"max={field_data.max():.3e}\n"
                     f"mean={field_data.mean():.3e}\n"
                     f"std={field_data.std():.3e}")
        ax_2d.text(0.02, 0.98, stats_text, transform=ax_2d.transAxes, verticalalignment="top",
                  bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5), fontsize=8)

        # === RIGHT: 6 Time Snapshots of Middle Vertical Profile ===
        timestep_indices = self._get_6_timestep_indices()
        colors = plt.cm.cool(np.linspace(0, 1, len(timestep_indices)))
        
        for i, step_idx in enumerate(timestep_indices):
            # Load mesh at this timestep
            temp_mesh = self.load_timestep(step_idx)
            z_pos, profile_vals = self._extract_middle_vertical_profile(temp_mesh)
            time_info_temp = self.get_timestep_info(step_idx)
            
            label = f"t={time_info_temp['time']:.2e}"
            ax_profiles.plot(profile_vals, z_pos, "-o", color=colors[i], linewidth=2, 
                           markersize=4, label=label, alpha=0.8)
        
        ax_profiles.set_xlabel(self.current_field)
        ax_profiles.set_ylabel("Vertical Position (m)")
        ax_profiles.set_title("1D Profiles: Middle Vertical Line (6 Timesteps)")
        ax_profiles.grid(True, alpha=0.3)
        ax_profiles.legend(loc="best", fontsize=8)
        
        # Reload current timestep
        self.load_timestep(self.current_timestep)
        
        if save_path:
            fig.savefig(save_path, bbox_inches="tight")
            print(f"✓ Saved plot to {save_path}")
        else:
            plt.tight_layout()
            plt.show()

    def _plot_2d_axisymmetric(self, mesh: pv.PolyData, time_info: Dict, save_path: Optional[str] = None):
        """Plot 2D axisymmetric mesh: heatmap (left) + 6 time profiles (right)."""
        coords = mesh.points
        field_data = mesh[self.current_field]

        # Create figure with 2 subplots
        fig, (ax_2d, ax_profiles) = plt.subplots(1, 2, figsize=(16, 8))

        # === LEFT: 2D Heatmap (r-z) ===
        norm = Normalize(vmin=field_data.min(), vmax=field_data.max())
        scatter = ax_2d.scatter(coords[:, 0], coords[:, 1], c=field_data, cmap="viridis",
                               norm=norm, s=50, alpha=0.8)

        ax_2d.set_xlabel("Radius r (m)")
        ax_2d.set_ylabel("Depth z (m)")
        ax_2d.set_title(
            f"{self.current_project}: {self.current_field} (Axisymmetric Heatmap)\n"
            f"Step {time_info['timestep']}, t={time_info['time']:.3e}"
        )
        ax_2d.set_aspect("equal", adjustable="box")
        ax_2d.grid(True, alpha=0.2)

        # Colorbar
        cbar = fig.colorbar(scatter, ax=ax_2d, label=self.current_field)

        # Statistics
        stats_text = (f"min={field_data.min():.3e}\n"
                     f"max={field_data.max():.3e}\n"
                     f"mean={field_data.mean():.3e}\n"
                     f"std={field_data.std():.3e}")
        ax_2d.text(0.02, 0.98, stats_text, transform=ax_2d.transAxes, verticalalignment="top",
                  bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5), fontsize=8)

        # === RIGHT: 6 Time Snapshots of Middle Vertical Profile ===
        timestep_indices = self._get_6_timestep_indices()
        colors = plt.cm.cool(np.linspace(0, 1, len(timestep_indices)))
        
        for i, step_idx in enumerate(timestep_indices):
            # Load mesh at this timestep
            temp_mesh = self.load_timestep(step_idx)
            z_pos, profile_vals = self._extract_middle_vertical_profile(temp_mesh)
            time_info_temp = self.get_timestep_info(step_idx)
            
            label = f"t={time_info_temp['time']:.2e}"
            ax_profiles.plot(profile_vals, z_pos, "-o", color=colors[i], linewidth=2, 
                           markersize=4, label=label, alpha=0.8)
        
        ax_profiles.set_xlabel(self.current_field)
        ax_profiles.set_ylabel("Depth z (m)")
        ax_profiles.set_title("1D Profiles: Center Vertical Line (6 Timesteps)")
        ax_profiles.grid(True, alpha=0.3)
        ax_profiles.legend(loc="best", fontsize=8)
        
        # Reload current timestep
        self.load_timestep(self.current_timestep)
        
        if save_path:
            fig.savefig(save_path, bbox_inches="tight")
            print(f"✓ Saved plot to {save_path}")
        else:
            plt.tight_layout()
            plt.show()

    def _plot_3d(self, mesh: pv.PolyData, time_info: Dict, save_path: Optional[str] = None):
        """Plot 3D mesh with PyVista."""
        plotter = pv.Plotter()
        plotter.add_mesh(mesh, scalars=self.current_field, cmap="viridis",
                        show_edges=True, edge_color="k", edge_opacity=0.1)
        plotter.add_title(
            f"{self.current_project}: {self.current_field}\n"
            f"Step {time_info['timestep']}, t={time_info['time']:.3e}"
        )

        if save_path:
            plotter.screenshot(save_path)
            print(f"✓ Saved plot to {save_path}")
        else:
            plotter.show()

    # ─────────────────────────────────────────────────────────────────────────
    # Interactive Mode
    # ─────────────────────────────────────────────────────────────────────────

    def interactive_mode(self):
        """Launch interactive command-line interface."""
        print("\n" + "="*70)
        print("ADSIM VTK Visualization - Interactive Mode")
        print("="*70)

        # Select project
        projects = self.get_projects()
        if not projects:
            print("No projects available")
            return

        print("\nAvailable projects:")
        for i, proj in enumerate(projects, 1):
            n_files = len(self.vtk_files[proj])
            print(f"  {i}. {proj} ({n_files} timesteps)")

        choice = input("\nSelect project (number): ").strip()
        try:
            project_idx = int(choice) - 1
            if 0 <= project_idx < len(projects):
                self.select_project(projects[project_idx])
            else:
                print("Invalid choice")
                return
        except ValueError:
            print("Invalid input")
            return

        # Interactive loop
        self._interactive_loop()

    def _interactive_loop(self):
        """Main interactive command loop."""
        while True:
            print("\n" + "-"*70)
            print(f"Project: {self.current_project} | Field: {self.current_field} | Step: {self.current_timestep}")
            print("-"*70)
            print("Commands:")
            print("  plot              - Plot current timestep (interactive profile extraction)")
            print("  field <name>      - Change field (or 'list' to see all)")
            print("  step <n>          - Go to timestep n")
            print("  next/prev         - Next/previous timestep")
            print("  first/last        - First/last timestep")
            print("  save <path>       - Save current plot to file")
            print("  animate           - Animate through all timesteps")
            print("  export [dir]      - Export all timesteps (default: visualization/)")
            print("  info              - Show current mesh/field info")
            print("  timeline          - Show timeline info")
            print("  quit              - Exit")
            print("-"*70)

            cmd = input("Command: ").strip()

            if cmd == "plot":
                self.plot_current()
            elif cmd.startswith("field"):
                parts = cmd.split(maxsplit=1)
                if len(parts) > 1:
                    if parts[1] == "list":
                        print(f"\nAvailable fields ({len(self.available_fields[self.current_project])}):")
                        for i, field in enumerate(self.available_fields[self.current_project], 1):
                            print(f"  {i:2d}. {field}")
                    else:
                        field_name = parts[1]
                        if field_name in self.available_fields[self.current_project]:
                            self.current_field = field_name
                            print(f"✓ Changed to field: {field_name}")
                        else:
                            print(f"✗ Field '{field_name}' not found")
            elif cmd.startswith("step"):
                parts = cmd.split(maxsplit=1)
                if len(parts) > 1:
                    try:
                        step = int(parts[1])
                        if 0 <= step < len(self.vtk_files[self.current_project]):
                            self.load_timestep(step)
                            print(f"✓ Moved to timestep {step}")
                        else:
                            print(f"✗ Timestep out of range (0-{len(self.vtk_files[self.current_project])-1})")
                    except ValueError:
                        print("✗ Invalid timestep number")
            elif cmd == "next":
                if self.current_timestep < len(self.vtk_files[self.current_project]) - 1:
                    self.load_timestep(self.current_timestep + 1)
                    print(f"✓ Moved to timestep {self.current_timestep}")
                else:
                    print("✗ Already at last timestep")
            elif cmd == "prev":
                if self.current_timestep > 0:
                    self.load_timestep(self.current_timestep - 1)
                    print(f"✓ Moved to timestep {self.current_timestep}")
                else:
                    print("✗ Already at first timestep")
            elif cmd == "first":
                self.load_timestep(0)
                print(f"✓ Moved to timestep 0")
            elif cmd == "last":
                self.load_timestep(len(self.vtk_files[self.current_project]) - 1)
                print(f"✓ Moved to last timestep {self.current_timestep}")
            elif cmd.startswith("save"):
                parts = cmd.split(maxsplit=1)
                if len(parts) > 1:
                    self.plot_current(save_path=parts[1])
                else:
                    print("✗ Please specify output path")
            elif cmd == "animate":
                self._animate_all()
            elif cmd.startswith("export"):
                parts = cmd.split(maxsplit=1)
                if len(parts) > 1:
                    # User specified custom base directory
                    self._export_all(parts[1])
                else:
                    # Use default 'visualization' directory
                    self._export_all()
            elif cmd == "info":
                self._print_info()
            elif cmd == "timeline":
                self._print_timeline()
            elif cmd == "quit":
                print("Goodbye!")
                break
            else:
                print("✗ Unknown command")

    def _animate_all(self):
        """Animate through all timesteps."""
        import matplotlib.pyplot as plt

        n_steps = len(self.vtk_files[self.current_project])
        print(f"\nAnimating {n_steps} timesteps...")

        for step in range(n_steps):
            self.load_timestep(step)
            self.plot_current()
            print(f"  {step+1}/{n_steps}")

    def _export_all(self, base_dir: str = "visualization"):
        """Export all timesteps as images to organized folder structure.
        
        Creates: {base_dir}/{project_name}/{project_name}_step_*.svg
        """
        import os

        # Create project-specific subfolder
        output_dir = os.path.join(base_dir, self.current_project)
        os.makedirs(output_dir, exist_ok=True)
        n_steps = len(self.vtk_files[self.current_project])

        print(f"\nExporting {n_steps} timesteps to {output_dir}...")

        for step in range(n_steps):
            self.load_timestep(step)
            filename = f"{self.current_project}_step_{step:06d}.svg"
            filepath = os.path.join(output_dir, filename)
            self.plot_current(save_path=filepath)
            print(f"  {step+1}/{n_steps}: {filename}")

        print(f"✓ Export complete: {output_dir}")

    def _print_info(self):
        """Print information about current state."""
        print(f"\nProject: {self.current_project}")
        print(f"Mesh type: {self.mesh_types.get(self.current_project, 'unknown')}")
        print(f"Total timesteps: {len(self.vtk_files[self.current_project])}")
        print(f"Current timestep: {self.current_timestep}")

        time_info = self.get_timestep_info(self.current_timestep)
        print(f"Current time: {time_info['time']:.3e}")

        print(f"Available fields: {len(self.available_fields[self.current_project])}")
        print(f"Current field: {self.current_field}")

    def _print_timeline(self):
        """Print timeline of all timesteps."""
        print(f"\nTimeline for {self.current_project}:")
        print(f"{'Step':<6} {'Time':<15} {'Status':<10}")
        print("-" * 40)
        
        n_steps = len(self.vtk_files[self.current_project])
        for i in range(n_steps):
            info = self.get_timestep_info(i)
            status = "→ CURRENT" if i == self.current_timestep else ""
            print(f"{i:<6} {info['time']:<15.3e} {status}")
        
        if n_steps > 20:
            print(f"... ({n_steps} total timesteps)")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="ADSIM VTK Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_vtk.py                          # Interactive mode
  python visualize_vtk.py --project water          # Select water project
  python visualize_vtk.py --field Matric_Head      # Start with specific field
  python visualize_vtk.py --export ./images        # Export all to images
        """
    )

    parser.add_argument("--output-dir", default="output", help="VTK output directory")
    parser.add_argument("--project", help="Project name to load")
    parser.add_argument("--field", help="Field to visualize")
    parser.add_argument("--step", type=int, default=0, help="Starting timestep")
    parser.add_argument("--plot", action="store_true", help="Plot current and exit")
    parser.add_argument("--save", help="Save plot to file")
    parser.add_argument("--export", nargs='?', const='visualization', help="Export all timesteps (default: visualization/{project}/)")
    parser.add_argument("--animate", action="store_true", help="Animate through timesteps")

    args = parser.parse_args()

    # Create visualizer
    viz = ADSIMVTKVisualizer(output_dir=args.output_dir)

    # Handle command-line arguments
    if args.project:
        if not viz.select_project(args.project):
            return
    else:
        # If no project specified, use interactive mode
        viz.interactive_mode()
        return

    if args.field and args.field in viz.available_fields[viz.current_project]:
        viz.current_field = args.field

    if args.step:
        viz.load_timestep(args.step)

    if args.export:
        viz._export_all(args.export)
    elif args.animate:
        viz._animate_all()
    elif args.plot or args.save:
        viz.plot_current(save_path=args.save)
    else:
        # Default: interactive mode
        viz._interactive_loop()


if __name__ == "__main__":
    main()
