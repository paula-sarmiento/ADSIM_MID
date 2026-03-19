#===============================================================================
# ADSIM Material Data Writer (TOML format)
# This script writes material properties for gases, liquids, and soils
# to a .toml file for use in the FEM code
#===============================================================================

proc ADSIM::WriteMaterialData { filename } {
    # Initialize the file
    GiD_WriteCalculationFile init $filename
    
    # Get root of XML tree
    set root [$::gid_groups_conds::doc documentElement]
    
    # Extract liquid properties ONCE (reuse for both [liquid] section and K_sat computation)
    set xp_liquid {//container[@n="materials"]/container[@n="m_liquid"]}
    set liquid_container [$root selectNodes $xp_liquid]
    
    set liquid_props(dyn_visc) 1.0e-3
    set liquid_props(density) 1000.0
    
    if {$liquid_container != ""} {
        set temp_visc [$liquid_container selectNodes {string(value[@n="Dynamic_viscosity_"]/@v)}]
        set temp_dens [$liquid_container selectNodes {string(value[@n="density_"]/@v)}]
        
        if {$temp_visc != "" && $temp_visc != "0.0"} { 
            set liquid_props(dyn_visc) $temp_visc
        }
        if {$temp_dens != "" && $temp_dens != "0.0"} {
            set liquid_props(density) $temp_dens
        }
    }
    
    # Write header
    ADSIM::WriteMaterialHeader
    
    # Write all dictionaries first (TOML requirement: all top-level keys before tables)
    ADSIM::WriteGasDictionary $root
    ADSIM::WriteSoilDictionary $root
    GiD_WriteCalculationFile puts ""
    
    # Write gas properties
    ADSIM::WriteGasProperties $root
    
    # Write liquid properties (pass extracted values)
    ADSIM::WriteLiquidProperties $root liquid_props
    
    # Write soil properties (pass extracted liquid properties for K_sat computation)
    ADSIM::WriteSoilProperties $root liquid_props
    
    # Close the file
    GiD_WriteCalculationFile end
}

#===============================================================================
# Write material file header
#===============================================================================
proc ADSIM::WriteMaterialHeader { } {
    GiD_WriteCalculationFile puts "# ADSIM material file header (need to add disclaimer and license)"
    GiD_WriteCalculationFile puts "# ADSIM_version = \"2025 v0.1.0\""
    GiD_WriteCalculationFile puts "# File_created = \"[clock format [clock seconds] -format %Y-%m-%d]\""
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write gas dictionary (must be at top level before any tables)
#===============================================================================
proc ADSIM::WriteGasDictionary { root } {
    # Get all gas materials
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    
    if {[llength $gas_blocks] == 0} {
        return
    }
    
    # Build gas dictionary array
    set gas_names [list]
    foreach gas_block $gas_blocks {
        set gas_name [$gas_block @name]
        lappend gas_names "\"$gas_name\""
    }
    
    # Write gas dictionary as array
    GiD_WriteCalculationFile puts "gas_dictionary_ = \[[join $gas_names ", "]\]"
}

#===============================================================================
# Write gas properties
#===============================================================================
proc ADSIM::WriteGasProperties { root } {
    # Get all gas materials
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    
    if {[llength $gas_blocks] == 0} {
        return
    }
    
    # Write gas properties
    GiD_WriteCalculationFile puts "# Gas properties"
    foreach gas_block $gas_blocks {
        set gas_name [$gas_block @name]
        GiD_WriteCalculationFile puts "\[gas.\"$gas_name\"\]"
        
        # Dynamic viscosity
        set dyn_visc [$gas_block selectNodes {string(value[@n="Dynamic_viscosity_"]/@v)}]
        GiD_WriteCalculationFile puts "dynamic_viscosity = $dyn_visc"
        
        # Molar mass
        set molar_mass [$gas_block selectNodes {string(value[@n="molar_mass_"]/@v)}]
        GiD_WriteCalculationFile puts "molar_mass = $molar_mass"
        
        # Diffusion coefficient
        set diff_coef [$gas_block selectNodes {string(value[@n="diffusion_coefficient_1"]/@v)}]
        GiD_WriteCalculationFile puts "diff_coefficient = $diff_coef"
        
        GiD_WriteCalculationFile puts ""
    }
}

#===============================================================================
# Write liquid properties (using pre-extracted values)
#===============================================================================
proc ADSIM::WriteLiquidProperties { root liquid_props_array } {
    upvar $liquid_props_array liquid_props
    
    # Get liquid container to fetch specific heat (only this is not pre-extracted)
    set xp_liquid {//container[@n="materials"]/container[@n="m_liquid"]}
    set liquid_container [$root selectNodes $xp_liquid]
    
    GiD_WriteCalculationFile puts "# Liquid properties"
    GiD_WriteCalculationFile puts "\[liquid\]"
    
    # Use pre-extracted properties
    GiD_WriteCalculationFile puts "dynamic_viscosity = $liquid_props(dyn_visc)"
    GiD_WriteCalculationFile puts "density = $liquid_props(density)"
    
    # Specific heat
    if {$liquid_container != ""} {
        set spec_heat [$liquid_container selectNodes {string(value[@n="specific_heat_water_"]/@v)}]
        GiD_WriteCalculationFile puts "specific_heat = $spec_heat"
    } else {
        GiD_WriteCalculationFile puts "specific_heat = 0.0"
    }
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write soil dictionary (must be at top level before any tables)
#===============================================================================
proc ADSIM::WriteSoilDictionary { root } {
    # Get all soil materials
    set xp_soils {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    if {[llength $soil_blocks] == 0} {
        return
    }
    
    # Build soil dictionary array
    set soil_names [list]
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        lappend soil_names "\"$soil_name\""
    }
    
    # Write soil dictionary as array
    GiD_WriteCalculationFile puts "soil_dictionary_ = \[[join $soil_names ", "]\]"
}

#===============================================================================
# Write soil properties (using pre-extracted liquid properties for K_sat)
#===============================================================================
proc ADSIM::WriteSoilProperties { root liquid_props_array } {
    upvar $liquid_props_array liquid_props
    
    # Get all soil materials
    set xp_soils {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    if {[llength $soil_blocks] == 0} {
        return
    }
    
    # Write soil properties
    GiD_WriteCalculationFile puts "# Soil properties"
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        GiD_WriteCalculationFile puts "\[soil.\"$soil_name\"\]"
        
        # Get the first group for this soil (assumes all groups have same values)
        set groups [$soil_block selectNodes {condition[@n="soil_basic"]/group}]
        
        if {[llength $groups] > 0} {
            set first_group [lindex $groups 0]
            
            # Physical properties
            GiD_WriteCalculationFile puts "# Physical properties"
            
            set spec_grav [$first_group selectNodes {string(.//value[@n="specific_gravity"]/@v)}]
            GiD_WriteCalculationFile puts "specific_gravity = $spec_grav"
            
            set porosity [$first_group selectNodes {string(.//value[@n="initial_porosity_"]/@v)}]
            GiD_WriteCalculationFile puts "porosity = $porosity"
            
            set saturation [$first_group selectNodes {string(.//value[@n="saturation_"]/@v)}]
            GiD_WriteCalculationFile puts "saturation = $saturation"
            
            set tort [$first_group selectNodes {string(.//value[@n="granular_tortuosity_"]/@v)}]
            GiD_WriteCalculationFile puts "granular_tortuosity = $tort"
            
            set perm [$first_group selectNodes {string(.//value[@n="intrinsic_permeabiliy_"]/@v)}]
            GiD_WriteCalculationFile puts "intrinsic_permeability = $perm"
            
            set res_water [$first_group selectNodes {string(.//value[@n="residual_water_content_"]/@v)}]
            GiD_WriteCalculationFile puts "residual_water_content = $res_water"
            
            set lime [$first_group selectNodes {string(.//value[@n="lime_content_"]/@v)}]
            GiD_WriteCalculationFile puts "lime_content = $lime"

            set reaction_rate [$first_group selectNodes {string(.//value[@n="reaction_rate_"]/@v)}]
            GiD_WriteCalculationFile puts "lime_reaction_rate = $reaction_rate"
            
            set res_lime [$first_group selectNodes {string(.//value[@n="lime_impure_"]/@v)}]
            GiD_WriteCalculationFile puts "residual_lime = $res_lime"
            
            # Thermal properties
            GiD_WriteCalculationFile puts "# Thermal properties"
            
            set spec_heat_solid [$first_group selectNodes {string(.//value[@n="specific_heat_solid_"]/@v)}]
            GiD_WriteCalculationFile puts "specific_heat_solids = $spec_heat_solid"
            
            # SWRC (Soil Water Retention Curve) properties
            GiD_WriteCalculationFile puts "# SWRC properties"
            
            set swrc_model [$first_group selectNodes {string(.//value[@n="swrc_model"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_model = \"$swrc_model\""
            
            set max_anw [$first_group selectNodes {string(.//value[@n="max_anw"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_max_anw = $max_anw"
            
            set saturation_max_anw [$first_group selectNodes {string(.//value[@n="saturation_max_anw"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_saturation_max_anw = $saturation_max_anw"
            
            # Van Genuchten parameters
            set vg_alpha [$first_group selectNodes {string(.//value[@n="vg_alpha"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_vg_alpha = $vg_alpha"
            
            set vg_n [$first_group selectNodes {string(.//value[@n="vg_n"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_vg_n = $vg_n"
            
            # Cavalcante parameters
            set cav_delta [$first_group selectNodes {string(.//value[@n="cav_delta"]/@v)}]
            GiD_WriteCalculationFile puts "swrc_cav_delta = $cav_delta"
            
            # Write hydraulic_properties subsection for SWRC model
            GiD_WriteCalculationFile puts ""
            GiD_WriteCalculationFile puts "\[soil.\"$soil_name\".hydraulic_properties\]"
            
            # Map physical properties to SWRC parameters
            GiD_WriteCalculationFile puts "swrc_model = \"$swrc_model\""
            GiD_WriteCalculationFile puts "theta_s = $porosity"
            GiD_WriteCalculationFile puts "theta_r = $res_water"
            
            # Compute K_sat from intrinsic permeability
            # K_sat = (k * rho * g) / mu
            # where: k=intrinsic_permeability [m^2], rho=density [kg/m^3], g=9.81 [m/s^2], mu=dynamic_viscosity [Pa·s]
            set g 9.81
            if {$perm != "" && $perm != "0.0"} {
                set k_sat [expr {($perm * $liquid_props(density) * $g) / $liquid_props(dyn_visc)}]
            } else {
                set k_sat 1.0e-5
            }
            GiD_WriteCalculationFile puts "K_sat = $k_sat"
            
            # Write SWRC fitting parameters
            if {$swrc_model eq "Van_Genuchten"} {
                GiD_WriteCalculationFile puts "alpha = $vg_alpha"
                GiD_WriteCalculationFile puts "n_param = $vg_n"
            } elseif {$swrc_model eq "Cavalcante"} {
                GiD_WriteCalculationFile puts "delta = $cav_delta"
            }
        }
        
        GiD_WriteCalculationFile puts ""
    }
}
