# ADSIM VTK Visualization Tool

A powerful interactive Python tool for visualizing ADSIM simulation output files (VTK format) with support for 1D, 2D, and axisymmetric meshes.

## Features

✓ **Auto-discovery** of all VTK files in the output directory  
✓ **Mesh type detection** (1D line, 2D Cartesian, 2D axisymmetric)  
✓ **Field visualization** with color mapping and statistics  
✓ **Automatic 6-timestep profile extraction** – shows vertical profiles at initial, final, and 4 intermediate times  
✓ **Timestep navigation** – step through simulation results  
✓ **Animation support** – playback through all timesteps  
✓ **Batch export** – save all timesteps as images  
✓ **Command-line interface** – easy access to all features  

## Installation

### Prerequisites

- Python 3.7+
- pip

### Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Or install manually
pip install pyvista numpy matplotlib
```

## Quick Start

### Interactive Mode (Recommended)

```bash
python visualize_vtk.py
```

This launches an interactive menu where you can:
- Select a project (water flow, reactive transport, etc.)
- Choose which field to visualize
- Navigate through timesteps
- Extract 1D profiles from 2D plots
- Export visualizations

### Command-Line Mode

Visualize a specific project:
```bash
python visualize_vtk.py --project Consolidation1D --field Total_Concentration --plot
```

Export all timesteps to SVG files:
```bash
python visualize_vtk.py --project LinVerif_diffusion_water --export ./output_images
```

## Usage Guide

### Interactive Commands (15 total)

Once in interactive mode, available commands are:

| Command | Description | Example |
|---------|-------------|---------|
| `plot` | Display current timestep with field coloring | `plot` |
| `field list` | Show all available fields | `field list` |
| `field <name>` | Switch to a different field | `field Matric_Head` |
| `step <n>` | Go to timestep n | `step 5` |
| `next` / `prev` | Move forward/backward one timestep | `next` |
| `first` / `last` | Jump to first/last timestep | `first` |
| `timeline` | Show all timesteps and their times | `timeline` |
| `save <path>` | Save current plot to SVG file | `save output.svg` |
| `animate` | Show animation through all timesteps | `animate` |
| `export [dir]` | Export all timesteps (default: `visualization/`) | `export` or `export ./custom` |
| `info` | Show mesh and field information | `info` |
| `quit` | Exit the program | `quit` |

### Profile Extraction (2D Meshes)

For 2D meshes, each visualization automatically displays:

1. **Left panel** – Full 2D domain as a heatmap scatter plot with colorbar
   - X-axis: Horizontal position (x for Cartesian, r for axisymmetric)
   - Y-axis: Vertical position (y for Cartesian, z for axisymmetric)
   - Colors: Field values (min to max range)
   - Statistics box: min, max, mean, std dev

2. **Right panel** – Automatic vertical profile extraction showing 6 timesteps
   - **6 overlaid profiles** from initial state to final state
   - Each profile extracted from the **vertical center line** of the domain
   - **Color gradient** indicates time progression (cyan = early, magenta = late)
   - **Legend** shows time values for each profile
   - X-axis: Field value
   - Y-axis: Vertical position (same as left plot)

**Key advantages:**
- No manual interaction needed – profiles are extracted automatically
- Shows field evolution over time in a single comparison plot
- Visualizes infiltration, diffusion, and other transient processes
- Includes representative samples (initial, final, 4 intermediate times)

**Example profiles:**
```
For 16-step simulation:
  Timestep 0   → Initial condition
  Timestep 4   → 25% through simulation  
  Timestep 7   → Mid-lower state
  Timestep 8   → Mid-upper state
  Timestep 12  → 75% through simulation
  Timestep 15  → Final state
```

## 2D Visualization Layout

For 2D meshes (Cartesian or axisymmetric), each figure contains two subplots:

**Layout Diagram:**
```
┌──────────────────────────────────────────────────────────────┐
│  TestDif: Matric_Head (Axisymmetric)  │  1D Profiles: Center │
│          Step 0, t=1.000e+03         │   Vertical Line      │
│                                      │   (6 Timesteps)      │
│  ┌──────────────────────┐  colorbar  │  ┌─────────────────┐ │
│  │                      │ ┌─────┐    │  │    ▲            │ │
│  │   2D Heatmap         │ │ max │    │  │    │ t=1.5e+03 │ │
│  │   Scatter Plot       │ │     │    │  │  y │ ╱╲       │ │
│  │                      │ │     │    │  │    │╱  ╲      │ │
│  │   Statistics Box:    │ │ 0.0 │    │  │    │     ╲    │ │
│  │   min=-1.0e+08       │ │     │    │  │    │      ╲   │ │
│  │   max=+1.0e+08       │ │ min │    │  │    │       ▼  │ │
│  │   mean, std dev      │ └─────┘    │  │    └─────────────┘ │
│  │                      │            │  │ Legend: 6 timesteps│
│  └──────────────────────┘            │  └─────────────────┘ │
│  ↓ Vertical domain (m)               │  Field value (x-axis)│
│  ← Horizontal domain (r/x in m)      │                      │
└──────────────────────────────────────────────────────────────┘
```

**Left Subplot (Heatmap):**
- Colored scatter plot showing all mesh nodes
- Colors represent field values (viridis colormap)
- Colorbar shows value scale and range
- Statistics box displays min, max, mean, and std dev
- Axes labeled with appropriate units and coordinates

**Right Subplot (Time Series Profiles):**
- 6 overlaid line profiles extracted from domain center
- Color gradient from cyan (early) to magenta (late)
- Each profile shows vertical distribution of field values
- Legend identifies each timestep with its time value
- Grid lines for easy value reading

## Visualization Types

### 1D Line Plots
For 1D mesh problems (vertical columns), shows:
- X-axis: position along the line
- Y-axis: field value
- Includes min/max/mean statistics

### 2D Heatmaps  
For 2D Cartesian problems, displays:
- Colored scatter plot of all nodes
- Field value encoded as color
- Colorbar with min/max scale
- Option to extract 1D profiles along custom lines

### 2D Axisymmetric
For axisymmetric (r-z) problems:
- Radius on x-axis, depth on y-axis
- Same coloring and profile extraction as Cartesian
- Good for radial diffusion/consolidation problems

## Command-Line Options

```bash
python visualize_vtk.py [OPTIONS]

Options:
  --output-dir DIR      Directory containing VTK files (default: output)
  --project NAME        Select project to load (skips interactive selection)
  --field NAME          Select field to visualize
  --step N              Starting timestep (default: 0)
  --plot                Plot and exit (non-interactive mode)
  --save PATH           Save plot to file
  --export [DIR]        Export all timesteps (default: visualization/{project}/)
  --animate             Animate through all timesteps
```

Examples:
```bash
# View Consolidation1D at timestep 10
python visualize_vtk.py --project Consolidation1D --step 10 --plot

# Export water flow simulation to default visualization/ folder
python visualize_vtk.py --project LinVerif_diffusion_water --export
# Creates: visualization/LinVerif_diffusion_water/*.svg

# Export to custom base directory
python visualize_vtk.py --project Consolidation1D --export ./my_results
# Creates: my_results/Consolidation1D/*.svg

# Quick look at a specific field
python visualize_vtk.py --project testBC --field Degree_of_Carbonation --plot --save carbonation.svg
```

## Output Files

### Single Plot Save
- **Format**: SVG (vector graphics - scalable without quality loss)
- **Location**: User-specified path (e.g., `output.svg`)

### Batch Export (Organized by Project)
When exporting all timesteps, files are organized in project-specific folders:

```
visualization/                          # Default base directory (or custom --export path)
├── Consolidation1D/
│   ├── Consolidation1D_step_000000.svg
│   ├── Consolidation1D_step_000001.svg
│   └── ... (21 timesteps)
├── LinVerif_diffusion_water/
│   ├── LinVerif_diffusion_water_step_000000.svg
│   └── ... (51 timesteps)
└── TestDif/
    └── ... (16 timesteps)
```

**File Naming Pattern**: `{project_name}_step_{XXXXXX}.svg`
- `XXXXXX`: Zero-padded 6-digit timestep number
- **Size**: Typically 50-200 KB per image (varies with mesh complexity)

## Data Points

### Water Flow Fields (Consolidation/Diffusion)
- `Matric_Head` – Pressure head (m)
- `Water_Content` – Volumetric water content (-)
- `Saturation` – Water saturation (-)
- `Water_Pressure` – Water pore pressure (Pa)
- `Water_Velocity` – Water flow velocity (m/s)

### Reactive Transport Fields
- `Total_Concentration` – Total gas concentration (mol/m³)
- `Absolute_Pressure` – Gas pressure (Pa)
- `Gas1_Concentration` – Individual gas species (mol/m³)
- `Gas1_Concentration_Rate` – Gas reaction rates
- `Degree_of_Carbonation` – CaCO₃ formation (-)
- `CaCO3_Concentration` – Carbonate concentration (mol/m³)
- `Temperature` – Temperature (K)

## Project Structure

```
ADSIM_MID/
├── visualize_vtk.py           # Main visualization script
├── requirements.txt            # Python dependencies
├── README_VISUALIZATION.md     # This file
├── output/                     # VTK simulation output
│   ├── Consolidation1D_000000.vtk
│   ├── Consolidation1D_000001.vtk
│   ├── ...
│   ├── LinVerif_diffusion_water_000000.vtk
│   └── ...
└── visualization/              # Organized visualization exports (auto-created)
    ├── Consolidation1D/
    ├── LinVerif_diffusion_water/
    └── ...
```

## Troubleshooting

### "No VTK files found"
- Make sure `output/` directory exists in current working directory
- Run from the ADSIM_MID root directory
- Check file naming: should be `{project}_{6-digit-number}.vtk`

### "Field not found"
- Use `field list` to see available fields for your project
- Field names are case-sensitive (e.g., `Matric_Head` not `matric_head`)

### Slow performance with large meshes
- Use `--step` to skip to specific timesteps instead of exporting all
- Profile extraction may be slow for meshes with 10k+ nodes
- Try a different field with faster rendering

### Plot shows but looks wrong
- Check mesh type with `info` command
- Verify coordinates are in expected range (use `plot` to visually inspect)
- Some meshes may need coordinate flipping (not yet implemented)

## Technical Details

### Mesh Type Detection
The script automatically detects mesh dimensionality by analyzing coordinate ranges:
- **1D**: Active range in only one dimension
- **2D Cartesian**: Active ranges in x and y, z ≈ 0
- **2D Axisymmetric**: Active ranges in r (x) and z (y), with r ≥ 0
- **3D**: Active ranges in all three dimensions

### Profile Extraction Algorithm
For each point along the drawn line:
1. Calculate position (r, z) at that distance
2. Find nearest mesh node
3. Use that node's field value

This provides accurate visualization without full interpolation.

### File Caching
Loaded VTK files are cached in memory for performance. To clear:
- Restart the script
- Or select a different project and return

## Advanced Usage

### Batch Processing Multiple Projects
Export all projects to organized folders:
```bash
# Creates: visualization/{project}/ for each
for project in Consolidation1D LinVerif_diffusion_water TestDif testBC; do
    echo "Exporting $project..."
    python visualize_vtk.py --project $project --export
done
```

### Export to Custom Base Directory
```bash
# Export all projects to custom location
python visualize_vtk.py --project Consolidation1D --export /media/external_drive/simulations
# Creates: /media/external_drive/simulations/Consolidation1D/*.svg
```

## Performance Notes

- **Startup time**: 2-5 seconds (VTK file discovery and parsing)
- **Plot time**: <1 second for 2D plots with <1000 nodes
- **Memory**: ~10-50 MB for typical simulations
- **Export time**: ~10 seconds per project (20 timesteps)

## Related Files

- [src/write_vtk.jl](src/write_vtk.jl) – VTK file generation (Julia)
- [plot_celia_kernel.jl](plot_celia_kernel.jl) – Existing manual VTK parsing
- [plot_src_output_linear.jl](plot_src_output_linear.jl) – Linear diffusion plotting

## Future Enhancements

Potential additions (not yet implemented):
- Vector field visualization (arrows for velocity)
- Time series comparison (overlaid multiple timesteps)
- Cross-section extraction at fixed X or Y coordinates
- Video export (MP4/GIF animation)
- Coordinate system customization
- Custom colormap upload
- Statistical analysis (trend over time)

## License

Same as ADSIM project (MIT)

## Questions?

Refer to the ADSIM README or contact the development team.
