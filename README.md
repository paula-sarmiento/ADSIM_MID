<div align="center">
  <img src="Problemtype/ADSIM_2025.gid/images/Logo_ADSIM.svg" alt="ADSIM Logo" width="400"/>
  
  # ADSIM
  
  **Advection-Diffusion for Soil Improvement and Modification**
  
  [![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/luisez1988/ADSIM/releases)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![Julia](https://img.shields.io/badge/julia-1.8+-purple.svg)](https://julialang.org/)
  [![GiD](https://img.shields.io/badge/GiD-Problem%20Type-blue.svg)](https://www.gidsimulation.com/)
  
  *A powerful simulation tool for modeling advection-diffusion processes in soil systems with multi-gas transport capabilities*
  
</div>

---

## 🌟 Overview

🚧 **WORK IN PROGRESS!** 🚧

ADSIM is a comprehensive finite element simulation framework designed for modeling complex advection-diffusion phenomena in soil environments. ADSIM is currently being built in julia and already has a GiD problem type that provides an intuitive interface for simulating multi-gas transport, concentration fields, reaction kinetics, and flow dynamics in porous media.

### Key Features

- 🔬 **Multi-Gas Transport Modeling** - Dynamic handling of multiple gas species with independent properties
- 📊 **Plain and axisymmetric** - Flexible solver configurations for planar axisymmetric problems
- ⚙️ **Advanced Boundary Conditions** - Comprehensive control over concentration, flow, and pressure boundaries
- 🎯 **Dynamic Field Generation** - Automatic adaptation of input fields based on material definitions
- 📝 **TOML-Based Configuration** - Human-readable parameter files for easy setup
- 🖥️ **GiD Integration** - Full pre/post-processing capabilities through GiD interface

---

## 🚀 Quick Start

### Download Pre-built Executable (Recommended)

**Latest Release: v0.1.0** 🎉

1. Download the latest Windows executable from the [Releases page](https://github.com/luisez1988/ADSIM/releases/latest)
2. Extract the ZIP file to your desired location
3. Run `ADSIM_app\bin\ADSIM.exe <project_name>` from the command line
4. Example: `ADSIM.exe --version` to check installation

### Building from Source

#### Prerequisites

- **Julia** 1.8 or higher ([Download Julia](https://julialang.org/downloads/))
- **GiD** (Pre/Post-processor) for mesh generation
- **Git** for cloning the repository

#### Installation Steps

1. **Clone this repository:**
   ```bash
   git clone https://github.com/luisez1988/ADSIM.git
   cd ADSIM
   ```

2. **Install Julia dependencies:**
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. **Run from source:**
   ```bash
   cd src
   julia -t8 kernel.jl <project_name>
   ```

4. **(Optional) Build standalone executable:**
   ```bash
   julia build_app.jl
   ```
   This creates a self-contained executable in the `ADSIM_app/` folder.

### GiD Problem Type Installation

1. Copy the `Problemtype/ADSIM_2025.gid` folder to your GiD problem types directory:
   ```
   <GiD_Installation>/problemtypes/
   ```

2. Restart GiD and select "ADSIM_2025" as your problem type

---

## 📖 Documentation

### Core Components

| Component | Description |
|-----------|-------------|
| **Materials** | Define soil properties and gas characteristics |
| **Boundary Conditions** | Set up concentration, flow, and pressure constraints |
| **Initial Conditions** | Configure starting concentrations and temperatures |
| **Calculation Parameters** | Specify solver settings, time steps, and units |

### Configuration Files

- 📄 [**Materials TOML Instructions**](Problemtype/Materials%20TOML%20writting%20instructions.md) - Define material properties
- 📄 [**Calculation Parameters Guide**](Problemtype/Calculation_Parameters_TOML_Instructions.md) - Configure simulation settings
- 📄 [**Mesh File Writer**](Problemtype/MESH_FILE_WRITER_README.md) - Mesh generation documentation
- 📄 [**Dynamic Gas Fields**](Problemtype/DYNAMIC_GAS_FIELDS_IMPLEMENTATION.md) - Multi-gas field implementation

---

## 🔧 Usage

### 1. Setup Your Model

Define your materials in the GiD interface, including:
- Soil properties (porosity, permeability, etc.)
- Gas species and their diffusion coefficients
- Temperature-dependent parameters

### 2. Apply Boundary Conditions

Configure boundary conditions for:
- **Concentration BC** - Fixed concentration values at boundaries
- **Flow BC** - Mass flux specifications
- **Total Pressure BC** - Pressure constraints

### 3. Set Initial Conditions

Specify initial states:
- Gas concentrations for each species
- Temperature distribution
- Initial pressure fields

### 4. Configure Calculation Parameters

Set simulation parameters using TOML format:
```toml
[units]
length = "m"
time = "s"
temperature = "K"

[solver]
solver_type = "2D-Plane"

[time]
total_simulation_time = 1.0
time_step = 0.01
```

### 5. Generate Mesh and Run

Use GiD's meshing tools to create your computational grid, then execute the simulation.

---

## 🎨 Features in Detail

### Dynamic Gas Field Generation

ADSIM automatically adapts input fields based on your material definitions. When you add or remove gas species, the interface dynamically updates to provide appropriate concentration and flow fields for each gas.

### Multi-Physics Coupling

The solver handles coupled processes including:
- Advective transport
- Diffusive transport
- Temperature effects
- Pressure-driven flow

### Flexible Solver Options

- 2D Plane strain/stress analysis
- Axisymmetric configurations
- 3D volumetric simulations

---

## 📁 Project Structure

```
ADSIM/
├── src/                          # Julia source code
│   ├── kernel.jl                 # Main execution kernel
│   ├── ADSIM.jl                  # Module definition
│   ├── version.jl                # Version management
│   ├── fully_explicit_solver.jl  # Solver implementation
│   ├── read_*.jl                 # Input file readers
│   ├── write_vtk.jl              # VTK output writer
│   ├── data/                     # Test cases and examples
│   └── output/                   # Simulation results
├── Problemtype/
│   ├── ADSIM_2025.gid/          # GiD problem type files
│   │   ├── images/               # Icons and logos
│   │   ├── scripts/              # TCL processing scripts
│   │   ├── xml/                  # XML configuration files
│   │   └── *.tcl, *.bat, *.spd  # Problem type definitions
│   └── *.md                      # Documentation files
├── scripts/
│   └── update_version.jl         # Version update utility
├── .github/
│   └── workflows/                # CI/CD automation
├── build_app.jl                  # Executable build script
├── VERSION                       # Version file
├── Project.toml                  # Julia package metadata
├── CHANGELOG.md                  # Version history
├── LICENSE                       # MIT License
└── README.md                     # This file
```

---

## 🔄 Releases and Versioning

ADSIM follows [Semantic Versioning](https://semver.org/). Each release is tagged and includes:

- **Compiled Windows executable** (standalone, no Julia installation required)
- **Source code** (for building on other platforms)
- **GiD problem type package**
- **Release notes** with changelog

### Current Release

**Version 0.1.0** (Alpha) - First public release with core functionality

See [CHANGELOG.md](CHANGELOG.md) for detailed version history and [Releases](https://github.com/luisez1988/ADSIM/releases) for downloads.

### Creating a New Release

For maintainers, follow the [Release Template](.github/RELEASE_TEMPLATE.md) for the complete release process

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit pull requests, report bugs, or suggest new features.

### Development Guidelines

1. Follow existing code structure and naming conventions
2. Update documentation for new features
3. Test thoroughly with various gas configurations
4. Add examples for new functionality

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
Copyright (c) 2025 luisez1988
```

---

## 👥 Authors

- **Luis E. Zambrano-Cruzatty** - [@luisez1988](https://github.com/luisez1988)

---

## 🔗 Links

- [GiD Website](https://www.gidsimulation.com/)
- [Project Repository](https://github.com/luisez1988/ADSIM)
- [Issue Tracker](https://github.com/luisez1988/ADSIM/issues)

---

## 📧 Contact

For questions, suggestions, or collaboration opportunities, please open an issue on GitHub.

---

<div align="center">
  
  **Made with ❤️ for the geotechnical and environmental engineering community**
  
  <img src="Problemtype/ADSIM_2025.gid/images/g6.png" alt="ADSIM Icon" width="64"/>
  
</div>
