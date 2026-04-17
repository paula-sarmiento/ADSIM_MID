# ADSIM Water Solver - Issues & Required Changes

## COMPREHENSIVE ANALYSIS OF PROBLEMS & FIXES NEEDED

Generated: April 17, 2026
Status: Code Review Complete

---

## 1. FUNCTION SIGNATURE MISMATCHES

### Issue 1.1: picard_richards!() Signature Mismatch
**File:** `src/implicit_richards_solver.jl` (line 317)
**Problem:** The function has different signature than what tests are trying to call

**Current Signature:**
```julia
function picard_richards!(
    h_curr     :: Vector{Float64},
    h_prev     :: Vector{Float64},
    mesh,
    elem_props :: Vector{ElementWaterProps},
    Δt         :: Float64,
    e_g        :: Vector{Float64},
    A          :: SparseMatrixCSC{Float64, Int},
    cache      :: RichardsCache;
    tol        :: Float64 = 1e-8,
    max_iter   :: Int     = 100,
    ω          :: Float64 = 1.0
)
```

**What Tests Are Calling:**
```julia
picard_result = picard_richards!(mesh, materials, cache, elem_props, calc_params, dt)
```

**Fix Required:**
- Either provide a wrapper function that matches the test interface, OR
- Update tests to match the actual picard_richards! signature
- The full signature needs: h_curr, h_prev, mesh, elem_props, Δt, e_g (gravity), A (sparse matrix), cache
- Missing: gravity vector construction and sparse matrix creation before Picard call

---

### Issue 1.2: initialize_shape_functions!() - Unclear Input
**File:** `src/shape_functions.jl`
**Problem:** Shape functions module interface not clear - appears to pre-allocate data during kernel initialization but not clear what gets stored/how.

**Needed Changes:**
- Document what mesh data is modified by `initialize_shape_functions!(mesh)`
- Ensure it can be called once and reused across multiple solvers

---

## 2. GLOBAL VARIABLE MANAGEMENT ISSUES

### Issue 2.1: Global Variable Pollution
**File:** `src/initialize_variables.jl` (multiple locations)
**Problem:** Heavy reliance on global variables (h, theta_w, P_boundary_water, q_boundary_water, S_r, P_water, v_water) makes testing and solver composition difficult

**Affected Variables:**
```julia
global h, theta_w, S_r, P_water, v_water  # Water state
global P_boundary_water                     # BC masking
global q_boundary_water                     # Neumann flux
```

**Issues:**
1. Must be initialized via `zero_variables!()`, not obvious from solver context
2. No clean way to save/restore state between solver calls
3. Hard to run multiple solvers in sequence without state contamination
4. Makes testing impossible without complex setup

**Fix Required:**
- Option A: Pass these as arguments (large refactor)
- Option B: Package into a struct (cleaner, recommended)
- Option C: Provide explicit initialization function with clear documentation

**Proposed Solution (Option B):**
```julia
struct WaterState
    h          :: Vector{Float64}
    theta_w    :: Vector{Float64}
    S_r        :: Vector{Float64}
    P_water    :: Vector{Float64}
    v_water    :: Matrix{Float64}
end

struct WaterBoundaryConditions
    P_boundary_water   :: Matrix{Int}
    q_boundary_water   :: Vector{Float64}
end
```

---

## 3. BOUNDARY CONDITION ENFORCEMENT ISSUES

### Issue 3.1: BC Enforcement in Assembly is Inefficient
**File:** `src/implicit_richards_solver.jl` (line 237-258 in assemble_richards!)
**Problem:** The BC masking uses `findnz()` in a loop, which is O(n²) or worse

**Current Code:**
```julia
for i in 1:mesh.num_nodes
    if P_boundary_water[i, 1] == 0  # BC node
        rows, cols, vals = findnz(A)  # ← Called for EVERY node!
        for k in eachindex(rows)
            if rows[k] == i
                A.nzval[k] = rows[k] == cols[k] ? 1.0 : 0.0
            end
        end
        R[i] = 0.0
    end
end
```

**Fix Required:**
- Replace with efficient sparse matrix row operations:
```julia
# Build row indices once
row_indices = [A.rowval[k] for k in A.colptr[1]:A.colptr[end]-1]

# Mark BC node rows once
bc_nodes = findall(==(0), P_boundary_water[:, 1])
for i in bc_nodes
    # Clear row i in sparse matrix
    for k in A.colptr[i]:(A.colptr[i+1]-1)
        if A.rowval[k] == i
            A.nzval[k] = 1.0  # Diagonal
        else
            A.nzval[k] = 0.0  # Off-diagonal
        end
    end
    R[i] = 0.0
end
```

**Or Better:**  Use SparseMatrixCSC droptol! or rebuild selectively.

---

### Issue 3.2: BC Re-enforcement in Time Loop
**File:** `test_phase1_infiltration_proper.jl` (inside time loop)
**Problem:** Manually re-applying BCs each time step is wrong - should be enforced only by assembly masking at each Picard iteration

**Current Test Code:**
```julia
for step in 1:n_steps
    # Re-apply Dirichlet BC at TOP (θ_w = 0.26)  ← WRONG LOCATION!
    for (node_id, theta_bc) in initial_bc_theta
        theta_w[node_id] = theta_bc
    end
    
    # Then solve...
    picard_result = picard_richards!(...)
end
```

**Fix Required:**
- Remove manual BC re-application from time loop
- Picard iteration should handle BC enforcement via assembly masking
- Only initialization needs to set BC values once: `apply_water_volumetric_content_bc!()`
- Assembly in Picard should then ENFORCE them via P_boundary_water masking

**Corrected Approach:**
```julia
for step in 1:n_steps
    # NO manual re-application - assembly handles it
    picard_result = picard_richards!(h, h_prev, mesh, elem_props, Δt, e_g, A, cache)
    h_prev .= h  # Update for next time step
end
```

---

## 4. INITIALIZATION PRIORITY ISSUES

### Issue 4.1: apply_all_initial_conditions!() Not Properly Documented
**File:** `src/initialize_variables.jl`
**Problem:** The function exists but its priority order vs. BC priority is unclear

**Questions:**
- Does `apply_all_initial_conditions!()` call both `apply_water_initial_conditions!()` AND `apply_water_boundary_conditions!()`?
- What's the correct sequence: IC first, then BC? Or are they separate?
- Looking at kernel.jl line 212-217, there are separate calls to apply_all_initial_conditions! and initialize_all_flows!, but water BC application happens separately

**Investigation Result:**
- Yes, `apply_all_initial_conditions!()` applies ICs for all phases
- Then separate `apply_water_boundary_conditions!()` should be called
- But `kernel.jl` doesn't show explicit water BC application - it's hidden in `apply_water_volumetric_content_bc!()` called within test

**Fix Required:**
- Make kernel.jl explicit about water BC application order
- Document that BC should override IC at prescribed nodes
- Ensure `apply_all_initial_conditions!()` doesn't implicitly apply BCs

---

## 5. MISSING HELPER FUNCTIONS FOR SPARSE MATRIX OPERATIONS

### Issue 5.1: BC Masking Needs Sparse Matrix Utility
**File:** `src/implicit_richards_solver.jl`
**Problem:** Matrix row zeroing for Dirichlet BCs is ad-hoc, error-prone

**Fix Required:**
Create utility function:
```julia
function enforce_dirichlet_row!(A::SparseMatrixCSC, i::Int, bc_value::Float64, R::Vector)
    """
    Zero row i of sparse matrix A, set diagonal to 1, clear residual.
    Maintains sparse structure for efficient Picard iteration.
    """
    for k in A.colptr[i]:(A.colptr[i+1]-1)
        if A.rowval[k] == i
            A.nzval[k] = 1.0  # Diagonal = 1
        else
            A.nzval[k] = 0.0  # Off-diagonal = 0
        end
    end
    R[i] = 0.0
end
```

---

## 6. MISSING CACHE MANAGEMENT FOR EFFICIENCY

### Issue 6.1: Cache Built Twice
**File:** `src/implicit_richards_solver.jl` (line 407-409)
**Problem:** Both test and kernel build cache separately, then element properties again

**Current Code in kernel:**
```julia
# kernel.jl line 207
initialize_shape_functions!(mesh)

# Later in implicit_richards_solver (line 407)
cache = build_richards_cache(mesh)
elem_props = precompute_element_water_props(mesh, materials)
```

**Fix Required:**
- Cache should be built ONCE after shape functions initialization
- Pass the same cache object to all solvers
- Avoid recomputing element properties multiple times

**Recommended Kernel Modification:**
```julia
# kernel.jl,  after shape function init
initialize_shape_functions!(mesh)
richards_cache = build_richards_cache(mesh)
elem_props = precompute_element_water_props(mesh, materials)

# Pass to solver
implicit_richards_solver(mesh, materials, ..., 
    prebuilt_cache = (cache=richards_cache, elem_props=elem_props))
```

---

## 7. NUMERICAL/PHYSICS ISSUES

### Issue 7.1: Backward Euler Stability Not Documented
**File:** `src/implicit_richards_solver.jl`
**Problem:** As implemented with lumped mass, backward Euler should be unconditionally stable BUT: no documentation or validation of this property

**Missing:**
- Energy stability analysis documentation
- Test for mass conservation verification
- Validation that theta_w remains in [0, porosity] after each step

**Fix Required:**
Add validation after each Picard solve:
```julia
function validate_water_solution!(mesh, materials)
    """Check that solution is physically reasonable"""
    for node_id in 1:mesh.num_nodes
        soil_id = mesh.node_to_soil[node_id]
        n = materials.soils[soil_id].porosity
        theta_res = materials.soils[soil_id].water.residual_water_content
        
        if theta_w[node_id] < theta_res || theta_w[node_id] > n
            @warn "Water content out of bounds at node $node_id: $(theta_w[node_id])"
            return false
        end
    end
    return true
end
```

---

## 8. TEST STRUCTURE ISSUES

### Issue 8.1: test_phase1_infiltration_proper.jl Has Undefined Variables
**File:** `test_phase1_infiltration_proper.jl` (line 150+)
**Problem:** References global variables that won't be defined in the scope:

- `interior_nodes` - needs to be computed
- `h, theta_w, P_water, v_water` - global variables that need explicit access
- `mesh.node_to_soil` - may not exist in mesh structure
- Return value of `picard_richards!()` called as tuple unpacking

**Fix Required:**
```julia
# After apply_water_boundary_conditions!
interior_nodes = [i for i in 1:mesh.num_nodes 
                  if !haskey(mesh.volumetric_content_bc, i)]

# Check returned iterations, not tuple
picard_iters = picard_richards!(h, h_prev, mesh, elem_props, dt, 
                                [gx, gy], A, cache)
```

---

## 9. ASSEMBLY CORRECTNESS ISSUES

### Issue 9.1: P_boundary_water Row Dimension
**File:** `src/initialize_variables.jl` (line 661)
**Problem:** P_boundary_water allocated as `ones(Int, Nnodes, 1)` (matrix) but accessed as vector

**Current Code:**
```julia
P_boundary_water = ones(Int, Nnodes, 1)  # Matrix [Nnodes × 1]

# Later in assembly:
if P_boundary_water[i, 1] == 0  # Accessing [i, 1] ✓ (correct for matrix)
```

**But** when applying BCs:
```julia
for (node_id, theta_bc) in mesh.volumetric_content_bc
    P_boundary_water[node_id, 1] = 0  # Setting [i, 1] ✓ (okay)
```

**Issue:** Mixing vector and matrix access is confusing

**Fix Required:**
```julia
P_boundary_water = ones(Int, Nnodes)  # Use vector, not matrix!

# Then everywhere:
if P_boundary_water[i] == 0  # Simpler
    P_boundary_water[node_id] = 0
end
```

---

## SUMMARY TABLE: Issues & Priority

| # | Issue | Severity | Category | Fix Time |
|----|-------|----------|----------|----------|
| 1.1 | picard_richards!() signature mismatch | **CRITICAL** | Interface | 1 hour |
| 3.1 | BC masking inefficient with findnz() | HIGH | Performance | 30 min |
| 3.2 | Manual BC re-application in time loop | **CRITICAL** | Correctness | 15 min |
| 2.1 | Global variable pollution | MEDIUM | Design | 2-3 hours |
| 4.1 | Unclear initialization sequence | MEDIUM | Documentation | 30 min |
| 9.1 | P_boundary_water matrix vs vector | MEDIUM | Clarity | 15 min |
| 5.1 | Sparse matrix utility missing | LOW | Code Quality | 30 min |
| 6.1 | Cache built multiple times | LOW | Efficiency | 30 min |
| 7.1 | Stability validation missing | LOW | Testing | 1 hour |
| 8.1 | Test has undefined variables | HIGH | Testing | 30 min |

---

## RECOMMENDED FIX ORDER

1. **PHASE 1 (CRITICAL - Required to Run):**
   - Fix Issue 1.1: Create wrapper or adapt test signatures
   - Fix Issue 3.2: Remove manual BC re-application from time loop
   - Fix Issue 8.1: Define all variables in test scope
   - Fix Issue 9.1: Make P_boundary_water a vector

2. **PHASE 2 (IMPORTANT - Correctness):**
   - Fix Issue 3.1: Optimize BC masking implementation
   - Fix Issue 4.1: Document initialization sequence clearly
   - Fix Issue 7.1: Add physical solution validation

3. **PHASE 3 (NICE-TO-HAVE - Design):**
   - Fix Issue 2.1: Refactor globals into structs (big refactor)
   - Fix Issue 5.1: Add sparse matrix utilities
   - Fix Issue 6.1: Build cache once, pass through

---

## FILES NEEDING EDITS

**CRITICAL:**
- `src/implicit_richards_solver.jl` (BC assembly, add wrapper, fix structure)
- `src/initialize_variables.jl` (P_boundary_water as vector, document BC sequence)
- `test_phase1_infiltration_proper.jl` (remove manual BC application, fix scope)

**IMPORTANT:**
- `src/kernel.jl` (document exact initialization sequence, cache management)
- `src/shape_functions.jl` (document what gets stored)

**NICE-TO-HAVE:**
- Create new file: `src/water_solver_utilities.jl` (sparse matrix ops, validation)
- Documentation: `WATER_SOLVER_README.md`

