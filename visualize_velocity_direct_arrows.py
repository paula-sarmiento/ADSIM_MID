#!/usr/bin/env python3
"""
Visualize Darcy velocity with DIRECT node-based arrows (no interpolation).
Arrow size is truly proportional to actual nodal velocities.
Uses plot_darcy_field function for consistent visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from plot_darcy_field import plot_darcy_field

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
    
    # Define parameters for LinVerif_diffusion verification case
    params = {
        'domain': '1.0 × 1.0 m',
        'K_s': 1.0e-6,  # ConstantSoil parameter
        'IC_h': -0.5,   # Initial condition: h = -0.5 m
        'BC_bottom': -1.0,  # Boundary condition at z=0 (bottom)
        'BC_top': 0.0       # Boundary condition at z=1 (top)
    }
    
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
