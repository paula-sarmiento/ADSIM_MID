# SWRC Parameter Management - All Implementation Options

**Date:** January 28, 2026  
**Status:** Evaluation Phase  
**Purpose:** Compare all possible approaches for SWRC parameter visibility management

---

## Problem Summary

SWRC (Soil Water Retention Curve) parameters need to be displayed based on the model selected for **each individual soil**, but the current `state=` attribute mechanism searches globally across all soils, causing cross-contamination.

**Current Bug:** When Soil 1 uses "Van Genuchten" and Soil 2 uses "Cavalcante", both sets of parameters appear for both soils.

---

## All Available Options

### Option 1: Show All Parameters Always (Remove State Checks)

**Concept:** Delete `state=` attributes, display all 9 SWRC parameters for every soil regardless of selection.

**Implementation:**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" editable="0"/>
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  
  <!-- All parameters always visible -->
  <value n="vg_alpha" pn="Van Genuchten α [L^-1]" v="0.0"/>
  <value n="vg_n" pn="Van Genuchten n [-]" v="0.0"/>
  <value n="vg_m" pn="Van Genuchten m [-]" v="0.0"/>
  <value n="cav_alpha" pn="Cavalcante α [L^-1]" v="0.0"/>
  <value n="cav_n" pn="Cavalcante n [-]" v="0.0"/>
  <value n="cav_m" pn="Cavalcante m [-]" v="0.0"/>
  <value n="cav_lambda" pn="Cavalcante λ [-]" v="0.0"/>
</container>
```

**User Experience:**
```
Soil 1 - Van Genuchten selected:
├─ SWRC Model: [Van Genuchten ▼]
├─ Anw: 0.0
├─ Van Genuchten α: 0.036      ← User fills these
├─ Van Genuchten n: 1.56       ← User fills these
├─ Van Genuchten m: 0.359      ← User fills these
├─ Cavalcante α: 0.0           ← Visible but ignored
├─ Cavalcante n: 0.0           ← Visible but ignored
├─ Cavalcante m: 0.0           ← Visible but ignored
└─ Cavalcante λ: 0.0           ← Visible but ignored
```

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐⭐⭐ | Delete 9 lines of XML |
| **User Experience** | ⭐⭐⭐ | Long list, but clear labels |
| **Data Handling** | ⭐⭐⭐⭐⭐ | Trivial - all fields static |
| **Maintainability** | ⭐⭐⭐⭐⭐ | Very simple |
| **Robustness** | ⭐⭐⭐⭐⭐ | No bugs possible |
| **Follows Code Pattern** | ⭐⭐⭐⭐⭐ | Standard static fields |

**Pros:**
- ✅ Ultra simple - 5 minute implementation
- ✅ Zero bugs
- ✅ All data persists automatically
- ✅ WriteMaterialData.tcl works as-is (just add SWRC writing)
- ✅ Parameter names clearly indicate which model

**Cons:**
- ❌ Long parameter list (11 fields total)
- ❌ User might fill wrong parameters
- ❌ Not elegant visually

**Best For:** Quick fix, simple projects, users familiar with SWRC models

---

### Option 2: Accept Current Buggy Behavior

**Concept:** Keep `state=` checks as-is, document that parameters from all models appear.

*(Omitted - not recommended)*

---

### Option 3: Restructure as Separate SWRC Materials

**Concept:** SWRC configs become separate blockdata like Gas materials.

*(Omitted - breaks conceptual model where SWRC is soil property)*

---

### Option 4: Dynamic Node Generation ⭐ (Original Favorite)

**Concept:** Use `<dynamicnode>` to generate only relevant parameters for each soil's selected model.

**Implementation:**

**XML:**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" 
         editable="0" actualize_tree="1"/>
  
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  
  <!-- Dynamic generation based on selection -->
  <dynamicnode command="ADSIM_2025::GenerateSWRCParameters" args="args"/>
</container>
```

**TCL Procedure:**
```tcl
proc ADSIM_2025::GenerateSWRCParameters { domNode args } {
    set root [$domNode ownerDocument documentElement]
    set parent_blockdata [$domNode selectNodes {ancestor::blockdata}]
    
    if {$parent_blockdata == ""} { return "" }
    
    # Get SWRC model for THIS specific soil
    set xp {.//container[@n="soil_swrc"]/value[@n="swrc_model"]/@v}
    set model [$parent_blockdata selectNodes "string($xp)"]
    
    set result ""
    
    if {$model eq "Van Genuchten"} {
        # Try to preserve existing values if present
        set existing_alpha [$parent_blockdata selectNodes {string(.//value[@n="vg_alpha"]/@v)}]
        if {$existing_alpha eq ""} { set existing_alpha "0.0" }
        
        append result "<value n=\"vg_alpha\" pn=\"Van Genuchten α \[L^-1\]\" v=\"$existing_alpha\" "
        append result "help=\"Van Genuchten alpha parameter.\"/>\n"
        
        # Similar for vg_n, vg_m...
        
    } elseif {$model eq "Cavalcante"} {
        # Generate Cavalcante parameters...
    }
    
    return $result
}
```

**User Experience:**
```
Soil 1 - Van Genuchten:
├─ SWRC Model: [Van Genuchten ▼]
├─ Anw: 0.0
├─ Van Genuchten α: 0.036      ← ONLY these 3 appear
├─ Van Genuchten n: 1.56
└─ Van Genuchten m: 0.359

Soil 2 - Cavalcante:
├─ SWRC Model: [Cavalcante ▼]
├─ Anw: 0.0
├─ Cavalcante α: 0.02          ← ONLY these 4 appear
├─ Cavalcante n: 2.0
├─ Cavalcante m: 0.5
└─ Cavalcante λ: 1.5
```

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐ | TCL procedure + XML changes + WriteMaterialData updates |
| **User Experience** | ⭐⭐⭐⭐⭐ | Perfect - only relevant parameters |
| **Data Handling** | ⭐⭐⭐ | Needs conditional reading |
| **Maintainability** | ⭐⭐⭐⭐ | Clear logic |
| **Robustness** | ⭐⭐⭐⭐ | Good if tested properly |
| **Follows Code Pattern** | ⭐⭐⭐ | New pattern (no current examples) |

**Pros:**
- ✅ Perfect per-soil visibility
- ✅ Clean, minimal interface
- ✅ Cleaner saved XML files
- ✅ Scalable to more models

**Cons:**
- ⚠️ Medium development effort
- ⚠️ Data persistence needs testing
- ⚠️ WriteMaterialData.tcl needs conditional logic
- ⚠️ Value preservation when switching models uncertain
- ⚠️ No existing example in codebase

**Best For:** Professional interface, multiple soil types, proper software architecture

---

### Option 5: Collapsible Containers with Visual Grouping 🆕

**Concept:** Keep all parameters but organize into collapsible containers per model. Use clear visual hierarchy.

**Implementation:**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" 
         editable="0" icon="icon_select.png"/>
  
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  
  <!-- Van Genuchten Container -->
  <container n="vg_params" pn="▼ Van Genuchten Parameters" 
             icon="folder.png" open="0"
             help="Fill these ONLY if Van Genuchten is selected above.">
    <value n="vg_alpha" pn="α [L^-1]" v="0.0" 
           help="ONLY for Van Genuchten model. Leave 0.0 if not using."/>
    <value n="vg_n" pn="n [-]" v="0.0"
           help="ONLY for Van Genuchten model. Leave 0.0 if not using."/>
    <value n="vg_m" pn="m [-]" v="0.0"
           help="ONLY for Van Genuchten model. Leave 0.0 if not using."/>
  </container>
  
  <!-- Cavalcante Container -->
  <container n="cav_params" pn="▼ Cavalcante Parameters" 
             icon="folder.png" open="0"
             help="Fill these ONLY if Cavalcante is selected above.">
    <value n="cav_alpha" pn="α [L^-1]" v="0.0"
           help="ONLY for Cavalcante model. Leave 0.0 if not using."/>
    <value n="cav_n" pn="n [-]" v="0.0"
           help="ONLY for Cavalcante model. Leave 0.0 if not using."/>
    <value n="cav_m" pn="m [-]" v="0.0"
           help="ONLY for Cavalcante model. Leave 0.0 if not using."/>
    <value n="cav_lambda" pn="λ [-]" v="0.0"
           help="ONLY for Cavalcante model. Leave 0.0 if not using."/>
  </container>
</container>
```

**User Experience:**
```
Soil 1 - Van Genuchten selected:
├─ SWRC Model: [Van Genuchten ▼]
├─ Anw: 0.0
├─ ▼ Van Genuchten Parameters      ← User expands this
│  ├─ α: 0.036
│  ├─ n: 1.56
│  └─ m: 0.359
├─ ▶ Cavalcante Parameters          ← Collapsed, ignored

```

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐⭐ | Just XML restructuring |
| **User Experience** | ⭐⭐⭐⭐ | Clear visual grouping |
| **Data Handling** | ⭐⭐⭐⭐⭐ | Static fields - trivial |
| **Maintainability** | ⭐⭐⭐⭐⭐ | Very simple |
| **Robustness** | ⭐⭐⭐⭐⭐ | No bugs possible |
| **Follows Code Pattern** | ⭐⭐⭐⭐⭐ | Standard containers |

**Pros:**
- ✅ Simple implementation (just XML)
- ✅ Clear visual organization
- ✅ Containers can be collapsed
- ✅ All data persists
- ✅ No TCL code needed
- ✅ WriteMaterialData.tcl simple

**Cons:**
- ⚠️ Still shows all parameters (but organized)
- ⚠️ User could expand wrong container
- ⚠️ Slightly longer than Option 4

**Best For:** Quick improvement over Option 1, visual clarity without complexity

---

### Option 6: Conditional Containers (State on Container Level) 🆕

**Concept:** Apply `state=` check at container level instead of individual values. May work better if GiD evaluates container state differently.

**Implementation:**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" 
         editable="0" actualize_tree="1"/>
  
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  
  <!-- Entire container visibility controlled -->
  <container n="vg_params" pn="Van Genuchten Parameters" 
             state="[ADSIM_2025::CheckSWRCModel {Van Genuchten}]">
    <value n="vg_alpha" pn="α [L^-1]" v="0.0"/>
    <value n="vg_n" pn="n [-]" v="0.0"/>
    <value n="vg_m" pn="m [-]" v="0.0"/>
  </container>
  
  <container n="cav_params" pn="Cavalcante Parameters"
             state="[ADSIM_2025::CheckSWRCModel {Cavalcante}]">
    <value n="cav_alpha" pn="α [L^-1]" v="0.0"/>
    <value n="cav_n" pn="n [-]" v="0.0"/>
    <value n="cav_m" pn="m [-]" v="0.0"/>
    <value n="cav_lambda" pn="λ [-]" v="0.0"/>
  </container>
</container>
```

**User Experience:** *Same as current but organized in containers*

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐⭐ | Minimal XML changes |
| **User Experience** | ⭐⭐⭐ | Same bug as current |
| **Data Handling** | ⭐⭐⭐⭐⭐ | Static fields |
| **Maintainability** | ⭐⭐⭐⭐ | Organized structure |
| **Robustness** | ⭐⭐ | **Same global bug** |
| **Follows Code Pattern** | ⭐⭐⭐⭐⭐ | Standard pattern |

**Pros:**
- ✅ Better organization than flat list
- ✅ Simple implementation
- ✅ Static fields persist

**Cons:**
- ❌ **Still has the global visibility bug**
- ❌ Doesn't solve the core problem

**Best For:** Not recommended - doesn't fix the actual issue

---

### Option 7: Hybrid - Static Fields + Smart TCL Warning 🆕

**Concept:** Show all parameters (Option 1) but add TCL procedure that validates and warns user if they filled parameters for wrong model.

**Implementation:**

**XML (same as Option 1):**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" 
         editable="0" actualize_tree="1"/>
  
  <!-- All parameters visible -->
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  <value n="vg_alpha" pn="Van Genuchten α [L^-1]" v="0.0"/>
  <!-- ... all other parameters ... -->
</container>
```

**TCL Validation (called when generating files):**
```tcl
proc ADSIM::ValidateSWRCParameters { root } {
    set xp_soils {//container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    set warnings [list]
    
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        set groups [$soil_block selectNodes {condition[@n="soil_basic"]/group}]
        
        if {[llength $groups] > 0} {
            set first_group [lindex $groups 0]
            set model [$first_group selectNodes {string(.//value[@n="swrc_model"]/@v)}]
            
            # Check for values in wrong parameters
            if {$model eq "Van Genuchten"} {
                # Check if Cavalcante parameters are non-zero
                set cav_alpha [$first_group selectNodes {string(.//value[@n="cav_alpha"]/@v)}]
                
                if {$cav_alpha != "" && $cav_alpha != "0.0"} {
                    lappend warnings "$soil_name: Using Van Genuchten but Cavalcante α is non-zero ($cav_alpha)"
                }
            }
            # Similar checks for other models...
        }
    }
    
    # Show warnings to user
    if {[llength $warnings] > 0} {
        set msg "⚠️ SWRC Parameter Warnings:\n\n"
        foreach warn $warnings {
            append msg "• $warn\n"
        }
        append msg "\nThese values will be ignored. Continue anyway?"
        
        set answer [tk_messageBox -title "SWRC Validation" -message $msg \
                    -icon warning -type okcancel]
        
        if {$answer eq "cancel"} {
            return 0
        }
    }
    
    return 1
}
```

**User Experience:**
1. User fills parameters freely
2. When generating files, validation runs
3. Popup warns about parameters filled for wrong model
4. User can continue or go back to fix

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐ | XML simple + TCL validation |
| **User Experience** | ⭐⭐⭐⭐ | Feedback prevents errors |
| **Data Handling** | ⭐⭐⭐⭐⭐ | Static fields |
| **Maintainability** | ⭐⭐⭐⭐ | Clear validation logic |
| **Robustness** | ⭐⭐⭐⭐ | Catches user errors |
| **Follows Code Pattern** | ⭐⭐⭐⭐⭐ | Uses existing Calculate pattern |

**Pros:**
- ✅ Simple to implement
- ✅ Catches user mistakes
- ✅ No data handling complexity
- ✅ Helpful feedback

**Cons:**
- ⚠️ Validation only happens at file generation
- ⚠️ Still shows all parameters
- ⚠️ User could ignore warnings

**Best For:** Safety layer on top of Option 1, preventing common errors

---

### Option 8: Multiple Conditions (One Per Model) 🆕

**Concept:** Create separate condition for each SWRC model. User assigns the appropriate condition.

**Implementation:**
```xml
<blockdata n="soil_material" name="Soil 1">
  <condition n="soil_basic" pn="Physical and thermal properties">
    <!-- Basic soil properties -->
  </condition>
  
  <!-- User chooses ONE of these conditions -->
  <condition n="swrc_none" pn="SWRC: None" ov="surface" ovm="element"
             help="No SWRC model applied.">
    <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  </condition>
  
  <condition n="swrc_van_genuchten" pn="SWRC: Van Genuchten" ov="surface" ovm="element"
             help="Van Genuchten SWRC model.">
    <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
    <value n="vg_alpha" pn="α [L^-1]" v="0.0"/>
    <value n="vg_n" pn="n [-]" v="0.0"/>
    <value n="vg_m" pn="m [-]" v="0.0"/>
  </condition>
  
  <condition n="swrc_cavalcante" pn="SWRC: Cavalcante" ov="surface" ovm="element"
             help="Cavalcante SWRC model.">
    <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
    <value n="cav_alpha" pn="α [L^-1]" v="0.0"/>
    <value n="cav_n" pn="n [-]" v="0.0"/>
    <value n="cav_m" pn="m [-]" v="0.0"/>
    <value n="cav_lambda" pn="λ [-]" v="0.0"/>
  </condition>
</blockdata>
```

**User Experience:**
```
Soil 1
├─ Physical and thermal properties (assigned to elements)
└─ SWRC Properties
   ├─ SWRC: None
   ├─ SWRC: Van Genuchten           ← User assigns this condition
   ├─ SWRC: Cavalcante

```

User assigns appropriate SWRC condition to soil elements, only relevant parameters appear.

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐ | Restructure XML + update reading |
| **User Experience** | ⭐⭐⭐ | Extra assignment step |
| **Data Handling** | ⭐⭐⭐ | Need to detect which condition |
| **Maintainability** | ⭐⭐⭐ | Multiple conditions to maintain |
| **Robustness** | ⭐⭐⭐⭐ | Clean separation |
| **Follows Code Pattern** | ⭐⭐⭐⭐ | Standard GiD pattern |

**Pros:**
- ✅ Perfect per-soil visibility
- ✅ Standard GiD mechanism
- ✅ Clear separation

**Cons:**
- ⚠️ Extra assignment step for user
- ⚠️ Could forget to assign SWRC condition
- ⚠️ WriteMaterialData needs to check which condition is assigned
- ⚠️ More XML to maintain (4 conditions)

**Best For:** Projects where different soils truly use different models

---

### Option 9: XML Include Per Model 🆕

**Concept:** Use GiD's `<include>` mechanism to load parameter definitions dynamically.

**Implementation:**

**Materials.xml:**
```xml
<container n="soil_swrc" pn="Partially saturated properties">
  <value n="swrc_model" pn="SWRC Model" v="None" 
         values="None,Van Genuchten,Cavalcante" 
         editable="0" actualize_tree="1"/>
  
  <value n="anw" pn="Air-water interfacial area Anw" v="0.0"/>
  
  <!-- Conditional includes based on some mechanism -->
  <include path="xml/swrc_vangenuchten.xml" 
           if="[ADSIM_2025::CheckSWRCModel {Van Genuchten}]"/>
  <include path="xml/swrc_cavalcante.xml"
           if="[ADSIM_2025::CheckSWRCModel {Cavalcante}]"/>
</container>
```

**swrc_vangenuchten.xml:**
```xml
<value n="vg_alpha" pn="Van Genuchten α [L^-1]" v="0.0"/>
<value n="vg_n" pn="Van Genuchten n [-]" v="0.0"/>
<value n="vg_m" pn="Van Genuchten m [-]" v="0.0"/>
```

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐ | If conditional includes work |
| **User Experience** | ⭐⭐⭐⭐⭐ | Perfect if it works |
| **Data Handling** | ⭐⭐⭐ | Unknown behavior |
| **Maintainability** | ⭐⭐⭐ | Multiple files |
| **Robustness** | ⭐ | **May not work - includes might not support conditions** |
| **Follows Code Pattern** | ⭐⭐⭐ | Similar to BC includes |

**Pros:**
- ✅ Modular structure
- ✅ Clean separation
- ✅ Perfect visibility IF it works

**Cons:**
- ❌ **GiD includes likely don't support conditional loading**
- ❌ **Same global state issue**
- ❌ Unproven approach
- ⚠️ Multiple files to maintain

**Best For:** Experimental - needs testing to see if GiD supports this

---

### Option 10: Pre-Generation Script with User Prompts 🆕

**Concept:** Before generating files, show dialog asking user which models are used, then only write those parameters.

**Implementation:**

**TCL Procedure:**
```tcl
proc ADSIM::PromptSWRCModels { } {
    # Get all soils and their models
    set root [customlib::GetBaseRoot]
    set xp_soils {//container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    set soil_models [dict create]
    
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        set groups [$soil_block selectNodes {condition[@n="soil_basic"]/group}]
        
        if {[llength $groups] > 0} {
            set first_group [lindex $groups 0]
            set model [$first_group selectNodes {string(.//value[@n="swrc_model"]/@v)}]
            dict set soil_models $soil_name $model
        }
    }
    
    # Show confirmation dialog
    set msg "SWRC Models Selected:\n\n"
    dict for {soil model} $soil_models {
        append msg "• $soil: $model\n"
    }
    append msg "\nParameters will be written only for selected models. Continue?"
    
    set answer [tk_messageBox -title "SWRC Models" -message $msg \
                -icon info -type okcancel]
    
    return [expr {$answer eq "ok"}]
}

proc ADSIM::WriteSoilSWRC { soil_block } {
    # Get model type
    set model [$soil_block selectNodes {string(.//value[@n="swrc_model"]/@v)}]
    
    # Write ONLY the parameters for selected model
    if {$model eq "Van Genuchten"} {
        set alpha [$soil_block selectNodes {string(.//value[@n="vg_alpha"]/@v)}]
        # ... write VG params only
    } elseif {$model eq "Cavalcante"} {
        set alpha [$soil_block selectNodes {string(.//value[@n="cav_alpha"]/@v)}]
        # ... write Cavalcante params only
    }
    # Ignore parameters for other models
}
```

| Criteria | Rating | Notes |
|----------|--------|-------|
| **Development Effort** | ⭐⭐⭐ | TCL dialog + selective writing |
| **User Experience** | ⭐⭐⭐⭐ | Clear feedback |
| **Data Handling** | ⭐⭐⭐⭐ | Clean selective output |
| **Maintainability** | ⭐⭐⭐⭐ | Clear logic |
| **Robustness** | ⭐⭐⭐⭐ | Explicit confirmation |
| **Follows Code Pattern** | ⭐⭐⭐⭐ | Similar to existing Calculate |

**Pros:**
- ✅ Can use Option 1 (show all params)
- ✅ Output only writes relevant params
- ✅ User confirmation reduces errors
- ✅ Simple implementation

**Cons:**
- ⚠️ Still shows all parameters in GUI
- ⚠️ Extra dialog step

**Best For:** Combination with Option 1 for safety + clean output

---

## Comparison Matrix

| Option | Dev Effort | UX Quality | Data Simple | Robustness | Best Use Case |
|--------|-----------|-----------|-------------|------------|---------------|
| **1. Show All** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | Quick fix |
| **4. Dynamic Node** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | Professional interface |
| **5. Collapsible Containers** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | Visual organization |
| **6. Container State** | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ❌ Still buggy |
| **7. Smart Validation** | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Error prevention |
| **8. Multiple Conditions** | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | Strict separation |
| **9. XML Includes** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐ | ❓ Experimental |
| **10. Pre-Gen Prompt** | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Safety layer |

---

## Recommended Combinations

### Approach A: Quick & Simple
**Option 5 (Collapsible Containers) + Option 7 (Validation)**
- Easy to implement
- Good visual organization
- Catches user errors
- No complex data handling

### Approach B: Professional & Clean
**Option 4 (Dynamic Node)**
- Best user experience
- Proper per-soil behavior
- Requires thorough testing
- More development time

### Approach C: Ultra Simple
**Option 1 (Show All) + Option 10 (Pre-Gen Prompt)**
- Minimal effort
- Works immediately
- Confirmation before output
- Simple data handling

---

## Decision Criteria

Choose based on:

1. **Time Available:**
   - Limited: Option 1 or 5
   - Moderate: Option 4
   - Plenty: Option 4 + 7

2. **User Expertise:**
   - Beginners: Option 5 or 7 (helpful guidance)
   - Experts: Option 1 (fast, no hand-holding)

3. **Project Complexity:**
   - Few soils, one model: Option 1
   - Many soils, mixed models: Option 4

4. **Maintenance Priority:**
   - Want simple code: Option 1 or 5
   - Want best UX: Option 4

---

## My Top 3 Recommendations

### 🥇 **Best Overall: Option 4 (Dynamic Node)**
Perfect per-soil visibility, clean interface, scalable

### 🥈 **Best Balance: Option 5 (Collapsible Containers)**
Simple implementation, good UX, no complexity

### 🥉 **Best Quick Fix: Option 5 + Option 7 (Containers + Validation)**
Organized display + error prevention, minimal effort

---

**What would you like to evaluate further or implement?**
