# TOML writting instructions for material data

This document shows the structure of the TOML file for the material properties.

## Header

The file contains "commented" info about the version of ADSIM used.
```toml
# ADSIM material file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"
```

## material properties
The following is the structure of the gas properties.
```toml
gas_dictionary_=[gas_1, gas_2, ... gas_n]

#gas properties
[gas.gas_1]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

[gas.gas_2]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

.
.
.

[gas.gas_n]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

#liquid properties
[liquid]
dynamic_viscosity= val
density= val
specific_heat= val

#soil dictionary
soil_dictionary= [soil_1, soil_2, ... soil_n]

#soil properties
[soil.soil_1]
#physical properties
specific_gravity= val
porosity= val
saturation= val
residual_water_content= val
granular_tortuosity= val
intrinsic_permeability= val
lime_content= val
residual_lime= val
reaction_rate= val
#thermal properties
specific_heat_solids= val
#SWRC properties (optional, defaults provided)
swrc_model= "None"  #options: "None", "Van Genuchten", "Cavalcante", "Brook and Corey"
swrc_max_anw= val  #optional, defaults to 0.0
swrc_saturation_max_anw= val  #optional, defaults to 0.0
swrc_vg_alpha= val  #optional, Van Genuchten parameter, defaults to 0.0
swrc_vg_n= val  #optional, Van Genuchten parameter, defaults to 0.0
swrc_cav_delta= val  #optional, Cavalcante parameter, defaults to 0.0
swrc_bc_air_entry_pressure= val  #optional, Brooks-Corey parameter, defaults to 0.0
swrc_bc_lambda= val  #optional, Brooks-Corey parameter, defaults to 0.0

[soil.soil_2]
#physical properties
specific_gravity= val
porosity= val
saturation= val
residual_water_content= val
granular_tortuosity= val
intrinsic_permeability= val
lime_content= val
residual_lime= val
reaction_rate= val
#thermal properties
specific_heat_solids= val
#SWRC properties (optional)
swrc_model= "None"
swrc_max_anw= val
swrc_saturation_max_anw= val
swrc_vg_alpha= val
swrc_vg_n= val
swrc_cav_delta= val
swrc_bc_air_entry_pressure= val
swrc_bc_lambda= val
.
.
.
[soil.soil_n]
#physical properties
specific_gravity= val
porosity= val
saturation= val
residual_water_content= val
granular_tortuosity= val
intrinsic_permeability= val
lime_content= val
residual_lime= val
lime_reaction_rate= val  #optional
#thermal properties
specific_heat_solids= val
#SWRC properties (optional)
swrc_model= "None"
swrc_max_anw= val
swrc_saturation_max_anw= val
swrc_vg_alpha= val
swrc_vg_n= val
swrc_cav_delta= val
swrc_bc_air_entry_pressure= val
swrc_bc_lambda= val

```