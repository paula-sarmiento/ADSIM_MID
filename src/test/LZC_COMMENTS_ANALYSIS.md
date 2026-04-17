"""
LCZ COMMENTS ANALYSIS — implicit_richards_solver.jl

Three architectural issues to address:

1. LINE 279: Neumann BC Assembly
   ─────────────────────────────
   Comment: "Assemble flow BC a priori. See how flow is assembled for 
            carbonation equation and then use directly."
   
   Issue: Flow BCs (Neumann) are handled inline in assemble_richards!()
          by just adding q_boundary_water to residual. The function
          apply_neumann_edge_richards!() exists but is never called.
   
   Status: Currently working but not following ADSIM pattern (inconsistent)
   
   Solution: Check how carbonation equation pre-assembles flow BCs, 
             then apply same pattern to Richards solver.

───────────────────────────────────────────────────────────────────────

2. LINE 412: Richards Cache Building
   ────────────────────────────────
   Comment: "I think this should be done when shape functions are 
            initialized and after the mesh is read once, then shared 
            across solvers. This way we avoid building the same cache 
            multiple times if we run multiple solvers sequentially."
   
   Issue: Cache is built INSIDE implicit_richards_solver() at line 412-413,
          called EVERY TIME the solver runs. But kernel.jl already calls
          initialize_shape_functions!(mesh) at Step 6/8.
   
   Current flow:
     kernel.jl [Step 6] → initialize_shape_functions!(mesh)
     kernel.jl [Step 8] → implicit_richards_solver()
                           └→ build_richards_cache(mesh) ← REDUNDANT
   
   Solution: Move `build_richards_cache()` to kernel.jl AFTER Step 6,
             pass `cache` as argument to solver instead of building it inside.
             This:
             - Builds cache only once
             - Shares cache if multiple solvers run sequentially
             - Follows ADSIM's "initialize once in kernel" pattern

───────────────────────────────────────────────────────────────────────

3. LINE 417: Element Properties Precomputation
   ──────────────────────────────────────────
   Comment: "This should go when we initialize and validate materials. 
            That will make it consistent with ADSIM's workflow"
   
   Issue: elem_props is computed INSIDE implicit_richards_solver() at line 417,
          called EVERY TIME the solver runs. But kernel.jl already validates
          materials at Step 3/8.
   
   Current flow:
     kernel.jl [Step 3] → validate_swrc_parameters(materials)
     kernel.jl [Step 3.1] → compute_K_sat_runtime!(materials, calc_params)
     kernel.jl [Step 8] → implicit_richards_solver()
                           └→ precompute_element_water_props() ← REDUNDANT
   
   Solution: Move `precompute_element_water_props()` to kernel.jl AFTER Step 3.1,
             pass `elem_props` as argument to solver instead of computing it inside.
             This:
             - Computes properties only once
             - Happens immediately after K_sat is computed
             - Follows ADSIM's "initialize once in kernel" pattern

───────────────────────────────────────────────────────────────────────

ARCHITECTURAL PATTERN IN ADSIM:

Kernel.jl executes physics initialization SEQUENTIALLY:
  Step 1: Read mesh                    → mesh object
  Step 2: Read materials              → materials object  
  Step 3.1: Validate & compute K_sat  → updated materials
  Step 3.2-3.3: Validate SWRC params  → verified materials
  [New] Step 3.4: Precompute elem water props → elem_props object
  Step 4: Initialize variables        → global arrays
  Step 5: Apply BCs and ICs           → populated arrays
  Step 6: Initialize shape functions  → shape cache
  [New] Step 6.5: Build Richards cache → cache object
  Step 8: Call solver(mesh, cache, elem_props, ...) → solution

Benefits:
- Each component computed ONCE, not repeatedly
- Solvers receive precomputed data, not raw data
- Consistent with ADSIM workflow
- Easier to profile and optimize
- Easy to extend with new solvers

───────────────────────────────────────────────────────────────────────

IMPLEMENTATION PLAN:

□ 1. Move build_richards_cache() to kernel.jl
     - Call after initialize_shape_functions!() [Step 6]
     - Define as `cache = build_richards_cache(mesh)`
     - Pass `cache` to implicit_richards_solver() 

□ 2. Move precompute_element_water_props() to kernel.jl
     - Call after compute_K_sat_runtime!() [Step 3.1]
     - Define as `elem_props = precompute_element_water_props(mesh, materials)`
     - Pass `elem_props` to implicit_richards_solver()

□ 3. Update implicit_richards_solver() signature
     - Add `cache::RichardsCache` argument
     - Add `elem_props::Vector{ElementWaterProps}` argument
     - Remove the cache-building and elem_props lines
     - Remove the log_print lines for "Building Richards cache" and 
       "Precomputing element water properties"
     - Update docstring

□ 4. Check Neumann BC assembly pattern
     - Compare with carbonation equation implementation
     - Consider consolidating flow BC assembly

"""
