# Dynamic Gas Concentration Fields Implementation

## Overview
The concentration fields for gases are now dynamically generated based on the number of gases defined in the Materials section. This means that when you add or remove gas materials, the corresponding concentration input fields will automatically update.

## What Changed

### 1. **InitialConcentrations.xml**
- Replaced static `Concentration_gas_1` field with dynamic field generation
- Added `update_proc` to trigger updates when materials change
- Added `<dynamicnode>` element that calls TCL procedure to generate fields

### 2. **Concentration_BC.xml** (Concentration Boundary Conditions)
- Updated to dynamically generate concentration fields for each gas
- Fields automatically match the gases defined in Materials

### 3. **Flow_BC.xml** (Flow Boundary Conditions)
- Updated to dynamically generate flow fields for each gas
- Fields automatically match the gases defined in Materials

### 4. **ADSIM_2025.tcl** (TCL Procedures)
Added three new dynamic field generator procedures:

#### `ADSIM_2025::GetGasConcentrationFields`
- Reads the number of gas blockdata items from Materials
- Generates concentration value fields for Initial Concentrations
- Each field is named `Concentration_gas_N` and labeled with the gas name

#### `ADSIM_2025::GetGasConcentrationBCFields`
- Generates concentration boundary condition fields
- Fields are named `concentration_gas_N_` for compatibility with existing scripts

#### `ADSIM_2025::GetGasFlowBCFields`
- Generates flow boundary condition fields
- Fields are named `flow_gas_N` for compatibility with existing scripts

## How It Works

1. **User defines gases in Materials:**
   - Navigate to Materials → Gas properties
   - Add gas blockdata items (Gas 1, Gas 2, etc.)
   - Each gas has properties like dynamic viscosity, molar mass, diffusion coefficient

2. **System automatically generates concentration fields:**
   - When you open Initial Concentrations condition
   - System counts the number of gas definitions in Materials
   - Generates one concentration field per gas
   - Field labels include the gas name from Materials

3. **Example:**
   ```
   If Materials has:
   - Gas 1 (CO2)
   - Gas 2 (O2)
   - Gas 3 (N2)
   
   Then Initial Concentrations will show:
   - Concentration CO2 [mol/L^3]
   - Concentration O2 [mol/L^3]
   - Concentration N2 [mol/L^3]
   ```

## Benefits

1. **Automatic synchronization:** No need to manually update concentration fields
2. **Less error-prone:** Field count always matches number of gases
3. **Better labeling:** Fields use actual gas names for clarity
4. **Scalable:** Works for any number of gases (1 to N)
5. **Consistent:** Same approach used for Initial Concentrations, Concentration BC, and Flow BC

## Compatibility

- Existing calculation scripts in `Generate_geom_dictionaries.tcl` will continue to work
- Field names follow the same convention (`concentration_gas_1_`, `flow_gas_1`, etc.)
- Default behavior: If no gases are defined, shows one default concentration field

## Usage Notes

- **Important:** Always define your gas materials first before setting concentrations
- If you add/remove gases in Materials, the concentration fields will update when you:
  - Reopen the tree
  - Navigate to the Initial Concentrations condition
  - Refresh the GiD interface
- The dynamic generation happens at runtime, not when files are saved

## Technical Details

### XML Dynamic Node Syntax
```xml
<dynamicnode command="ADSIM_2025::GetGasConcentrationFields" args="args"/>
```

### TCL Procedure Pattern
```tcl
proc ADSIM_2025::GetGasConcentrationFields { domNode args } {
    # Get root XML document
    set root [$domNode ownerDocument documentElement]
    
    # Query for all gas blockdata items
    set xp {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp]
    set num_gases [llength $gas_blocks]
    
    # Generate value fields dynamically
    set result ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        # Get gas name and create field
        ...
        append result "<value n=\"...\" pn=\"...\" v=\"0.0\"/>\n"
    }
    
    return $result
}
```

## Future Enhancements

Potential improvements:
1. Auto-update fields when Materials change without needing to reopen
2. Add validation to ensure concentrations sum correctly
3. Link concentration units to gas properties
4. Add tooltips showing which gas each field corresponds to

---

**Implementation Date:** November 14, 2025  
**Developer:** ADSIM Development Team  
**Version:** ADSIM v0.x.x
