#!/usr/bin/env python3
"""
Visualize Darcy velocity with DIRECT node-based arrows (no interpolation).
Arrow size is truly proportional to actual nodal velocities.
Uses plot_darcy_field function for consistent visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from pathlib import Path

# TOML parsing: try built-in (3.11+), fall back to tomli package
try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        print("Error: tomllib (Python 3.11+) or tomli package required")
        print("Install: pip install tomli")
        raise


def plot_darcy_field(coords, IEN, h, qx, qy, title="Darcy velocity field", 
                     params=None):
    """
    Visualize Darcy velocity vector field superimposed over pressure head colormap.
    
    Parameters
    ----------
    coords : (N, 2) float array
        Nodal coordinates [x, z] where N is number of nodes
    IEN : (4, n_elem) int array
        Element connectivity. Local node order: 1=SE, 2=NE, 3=NW, 4=SW
    h : (N,) float array
        Nodal pressure heads at one time step
    qx, qy : (N,) float arrays
        Nodal Darcy velocities (recovered via compute_darcy_nodes or equivalent)
    title : str, optional
        Figure title
    params : dict, optional
        Dictionary with parameters to display:
        {'K_s': float, 'IC_h': float, 'BC_bottom': float, 'BC_top': float, 'domain': str}
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    n_elem = IEN.shape[1]
    quads = []
    h_elem_mean = []
    
    for e in range(n_elem):
        nodes = IEN[:, e] - 1
        quad_coords = coords[nodes, :]
        quads.append(quad_coords)
        h_mean = np.mean(h[nodes])
        h_elem_mean.append(h_mean)
    
    h_elem_mean = np.array(h_elem_mean)
    
    poly = PolyCollection(quads, 
                         cmap='Blues',
                         edgecolors='k',
                         linewidths=0.4,
                         alpha=0.7)
    
    h_norm = Normalize(vmin=-1.0, vmax=0.0)
    poly.set_array(h_elem_mean)
    poly.set_norm(h_norm)
    
    ax.add_collection(poly)
    
    cbar1 = plt.colorbar(poly, ax=ax, label='h [m]', pad=0.09, fraction=0.046)
    cbar1.set_ticks([-1.0, -0.75, -0.5, -0.25, 0.0])
    
    mag = np.sqrt(qx**2 + qy**2) + 1e-16
    ux = qx / mag
    uy = qy / mag
    
    quiver = ax.quiver(coords[:, 0], coords[:, 1], ux, uy, mag,
                       cmap='plasma',
                       scale=25,
                       width=0.004,
                       clim=(0, mag.max()),
                       edgecolors='k',
                       linewidths=0.3,
                       alpha=0.85)
    
    cbar2 = plt.colorbar(quiver, ax=ax, label='|q| [m/s]', pad=0.05, fraction=0.08)
    v_max = mag.max()
    v_ticks = np.linspace(0, v_max, 6)
    cbar2.set_ticks(v_ticks)
    
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]', fontsize=12, fontweight='bold')
    ax.set_ylabel('z [m]', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    if np.abs(coords[:, 1].min()) < 1e-6:
        ax.invert_yaxis()
    
    ax.set_xlim(coords[:, 0].min() - 0.05, coords[:, 0].max() + 0.05)
    ax.set_ylim(coords[:, 1].min() - 0.05, coords[:, 1].max() + 0.05)
    
    ax.grid(True, alpha=0.2, linestyle='--')
    
    if params is not None:
        param_text = "Parameters:\n"
        if 'domain' in params:
            param_text += f"Domain: {params['domain']}\n"
        if 'K_s' in params:
            param_text += f"K_s: {params['K_s']:.2e} m/s\n"
        if 'IC_h' in params:
            param_text += f"IC: h = {params['IC_h']:.2f} m\n"
        if 'BC_top' in params and 'BC_bottom' in params:
            param_text += f"BC: h(top)={params['BC_top']:.2f} m,\n"
            param_text += f"    h(bot)={params['BC_bottom']:.2f} m"
        
        ax.text(0.02, 0.02, param_text, transform=ax.transAxes,
                fontsize=9, verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85, pad=0.5))
    
    return fig

def parse_vtk_file(filename):
    """Parse VTK file and extract data including element connectivity."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Extract time from header (line 2: "... Time = X.X")
    time_value = 0.0
    if len(lines) > 1:
        header_line = lines[1]
        if "Time = " in header_line:
            try:
                time_str = header_line.split("Time = ")[-1].strip()
                time_value = float(time_str)
            except:
                pass
    
    # Find POINTS
    points_idx = None
    cells_idx = None
    for i, line in enumerate(lines):
        if line.startswith('POINTS'):
            points_idx = i
        elif line.startswith('CELLS'):
            cells_idx = i
    
    # Parse points
    num_points = int(lines[points_idx].split()[1])
    x = np.zeros(num_points)
    y = np.zeros(num_points)
    
    point_line_idx = points_idx + 1
    point_count = 0
    for line in lines[point_line_idx:]:
        if line.strip() == '':
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                x[point_count] = float(parts[0])
                y[point_count] = float(parts[1])
                point_count += 1
                if point_count >= num_points:
                    break
            except:
                pass
    
    # Parse cells (element connectivity)
    IEN = None
    if cells_idx is not None:
        cells_line = lines[cells_idx].split()
        n_cells = int(cells_line[1])
        
        # Read cell data
        cell_data = []
        line_idx = cells_idx + 1
        while len(cell_data) < n_cells:
            parts = lines[line_idx].split()
            n_nodes = int(parts[0])
            nodes = [int(p) for p in parts[1:n_nodes+1]]
            cell_data.append(nodes)
            line_idx += 1
        
        # Build IEN: (n_nodes_per_elem, n_elements)
        if cell_data and len(cell_data[0]) == 4:  # Q4 elements
            IEN = np.array(cell_data, dtype=int).T + 1  # Convert to 1-indexed, transpose
    
    # Find POINT_DATA
    point_data_idx = None
    for i, line in enumerate(lines):
        if line.startswith('POINT_DATA'):
            point_data_idx = i
            break
    
    h = np.zeros(num_points)
    q_x = np.zeros(num_points)
    q_y = np.zeros(num_points)
    
    if point_data_idx is not None:
        for i in range(point_data_idx + 1, len(lines)):
            line = lines[i]
            
            if 'SCALARS Matric_Head' in line:
                data_start = i + 2
                for j, data_line in enumerate(lines[data_start:data_start + num_points]):
                    try:
                        h[j] = float(data_line.strip())
                    except:
                        pass
            
            elif 'VECTORS Water_Velocity' in line:
                data_start = i + 1
                for j, data_line in enumerate(lines[data_start:data_start + num_points]):
                    try:
                        parts = data_line.split()
                        if len(parts) >= 2:
                            q_x[j] = float(parts[0])
                            q_y[j] = float(parts[1])
                    except:
                        pass
    
    coords = np.column_stack([x, y])
    return coords, IEN, h, q_x, q_y, time_value


def extract_params_from_toml(project_name, data_dir="src/data"):
    """
    Extract simulation parameters from TOML config files.
    Reads material properties and mesh extent to populate parameter box.
    """
    data_path = Path(data_dir)
    mat_file = data_path / f"{project_name}_mat.toml"
    calc_file = data_path / f"{project_name}_calc.toml"
    mesh_file = data_path / f"{project_name}.mesh"
    
    params = {}
    
    # Extract K_sat from materials file
    if mat_file.exists():
        with open(mat_file, 'rb') as f:
            mat_data = tomllib.load(f)
        
        if 'soil_dictionary_' in mat_data and mat_data['soil_dictionary_']:
            soil_name = mat_data['soil_dictionary_'][0]
            soil_props = mat_data.get('soil', {}).get(soil_name, {})
            
            k = soil_props.get('intrinsic_permeability', 0.0)
            if k > 0 and 'liquid' in mat_data:
                rho_w = mat_data['liquid'].get('density', 1000.0)
                mu_w = mat_data['liquid'].get('dynamic_viscosity', 1e-3)
                
                if calc_file.exists():
                    with open(calc_file, 'rb') as f:
                        calc_data = tomllib.load(f)
                    g = calc_data.get('gravity', {}).get('gravity_magnitude', 9.81)
                else:
                    g = 9.81
                
                K_sat = k * rho_w * g / mu_w
                params['K_s'] = K_sat
    
    # Extract boundary conditions from calc file
    if calc_file.exists():
        with open(calc_file, 'rb') as f:
            calc_data = tomllib.load(f)
        params.setdefault('IC_h', -0.5)
        params.setdefault('BC_bottom', -1.0)
        params.setdefault('BC_top', 0.0)
    
    # Extract mesh extent
    if mesh_file.exists():
        try:
            mesh_coords = []
            with open(mesh_file, 'r') as f:
                in_coords = False
                for line in f:
                    line = line.strip()
                    if line == 'coordinates':
                        in_coords = True
                        continue
                    if in_coords:
                        if line and not line.startswith('[') and not line.startswith('#'):
                            try:
                                parts = line.split()
                                if len(parts) >= 2:
                                    mesh_coords.append([float(parts[0]), float(parts[1])])
                            except ValueError:
                                break
            
            if mesh_coords:
                mesh_arr = np.array(mesh_coords)
                domain_x = mesh_arr[:, 0].max() - mesh_arr[:, 0].min()
                domain_y = mesh_arr[:, 1].max() - mesh_arr[:, 1].min()
                params['domain'] = f'{domain_x:.1f} × {domain_y:.1f} m'
        except Exception as e:
            print(f"  ⚠ Could not parse mesh: {e}")
    
    params.setdefault('domain', '1.0 × 1.0 m')
    params.setdefault('K_s', 0.1)
    
    return params


# ────────────────────────────────────────────────────────────────────────────
# Process selected timesteps
# ────────────────────────────────────────────────────────────────────────────

output_dir = Path("output")
vtk_files = sorted(output_dir.glob("LinVerif_diffusion_water_*.vtk"))

step_indices = [0, len(vtk_files)//3, 2*len(vtk_files)//3, -1]
selected_files = [vtk_files[i] for i in step_indices if i < len(vtk_files)]

print(f"[PROCESS] Creating {len(selected_files)} timestep visualizations")
print(f"         with DIRECT node-based arrows (no interpolation)\n")

for idx, vtk_file in enumerate(selected_files):
    print(f"[{idx+1}/{len(selected_files)}] {vtk_file.name}")
    
    try:
        coords, IEN, h, q_x, q_y, time_value = parse_vtk_file(str(vtk_file))
    except Exception as e:
        print(f"  ⚠ Error: {e}")
        continue
    
    q_mag = np.sqrt(q_x**2 + q_y**2)
    
    if q_mag.max() < 1e-10:
        print(f"  ⚠ No velocity data")
        continue
    
    print(f"  ✓ {len(coords)} nodes, {IEN.shape[1] if IEN is not None else 'unknown'} elements")
    print(f"  ✓ |q| range: [{q_mag.min():.6e}, {q_mag.max():.6e}]")
    
    # ────────────────────────────────────────────────────────────────────────
    # Use plot_darcy_field for professional visualization
    # ────────────────────────────────────────────────────────────────────────
    
    timestep = int(vtk_file.stem.split('_')[-1])
    title = f'Darcy Velocity Field - t = {time_value:.2f} s'
    
    # Extract parameters dynamically from TOML config files
    project_name = '_'.join(vtk_file.stem.split('_')[:-2])
    params = extract_params_from_toml(project_name)
    
    if IEN is not None:
        # Use the professional plot_darcy_field function
        fig = plot_darcy_field(coords, IEN, h, q_x, q_y, title=title, params=params)
    else:
        # Fallback: create simpler figure without mesh
        print(f"  ⚠ Could not parse IEN from VTK file, using fallback visualization")
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Scatter plot with color
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=np.sqrt(q_x**2 + q_y**2),
                            s=100, cmap='plasma', edgecolors='k', linewidth=0.5)
        ax.set_aspect('equal')
        ax.set_xlabel('x [m]', fontsize=12, fontweight='bold')
        ax.set_ylabel('z [m]', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        plt.colorbar(scatter, ax=ax, label='|q| [m/s]')
    
    # Save
    output_svg = output_dir / f"darcy_velocity_direct_timestep_{timestep:06d}.svg"
    output_png = output_dir / f"darcy_velocity_direct_timestep_{timestep:06d}.png"
    
    plt.tight_layout()
    plt.savefig(str(output_svg), format='svg', dpi=150, bbox_inches='tight')
    plt.savefig(str(output_png), format='png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  → Saved: {output_svg.name}")
    print(f"  → Saved: {output_png.name}\n")

print("="*80)
print("✓ COMPLETE - Direct node-based visualizations")
print("="*80)
