# Matrices and vectors in ADSIM

This is documentation to create the matrices and vectors variables for ADSIM. 

## Lumped mass system

The formulation of ADSIM is based on fully explicit lumped mass system formulation. Thus the lumped mass matrix can be stored in variable of the form `M_L(Nnodes)`. This mass matrix is the same for all gases.

## Concentration flows

There are five concentration flows vectors that need to be calculated, namely:
1. Advection flow
2. Gravitational flow
3. Diffusion flow
4. Boundary flows
5. Sink/source flows

Each flow is particular for each gas species. Thus the flows should be specified as `q(Nnodes, Ngases)`. The only flow that is only defined for one gas is the sink/source term. Which will be only activated for a unique gas at this moment, and that is CO2.

## Initialization of boundary flows

For boundary conditions where flow is assigned, these must be initialized from data written in the GID interface.

Note that the GID interfae only permits the assignment of unirform steady normal flows. Thus, the parameter assigned to the nodes, must consider the length of influence associated with that node as in the example below.

.------ --------.--------
 h/2       h/2     h/2

## Finding flows in pressure nodes

In the future, it would be important to determine the flow that passes trought areas defined by the user.

# Writting vtk output

ADSIM uses vtk templates to write the output for visualization. 

The output consists of mesh data with nodal and element information. 

At this stage, only the advection and diffusion flows could be measured at the Gauss points, but this will be implemented later.

Currently all data are mapped to the nodes including:
1. Concentrations of all gases (scalar)
2. Total concentration (scalar)
3. Absolute pressure (scalar)
4. Rate of concentrations (scalar)
5. Rate of reactions (scalar)
6. Lime concentration (scalar)
7. CO2 concentration (scalar)
8. CaCO3 concentration (scalar)
9. Degree of carbonation ($DoC$)
10. Volumetric binder content (scalar)
11. Gas seepage velocity (vector)
12. Temperature (scalar)
13. Rate of temperature (scalar)

The output files should be correctly serialized so that Paraview can understand the 

# Time step calculation

The time step is calculated using the minor of three time scales defined below:
$$    
    \Delta t= C_N \  \text{min} \left\{ \cfrac{h_{min}^2 \tau}{\theta_g D_{max}}, h_{min}^2 \left(\cfrac{ \mu_{g}}{C_g^i K T R } \right)_{min} , \cfrac{1}{ \kappa \theta_w C_{CO_2,max} }  \right\}
$$

Where $0 < C_N \leq 1$ represents the Courant number and $h$ indicates the characteristic mesh size. Numbers in the numerator or denominator are maximized or minimized to ensure the safe calculation of the critical time step.

$\kappa=$ Tortuosity factor <br>
$\theta_g=n-\theta_w $ volumetric gas concentration <br>
$D_{max}=$ maximum diffusion gas coefficient <br>
$\mu_g=$ Minimum gas dynamic viscosity <br>
$C_g^i=$ Maximum initial gas concentration <br>
$C_{CO_2,max}=$ maximum concentration of $CO_2$

The code loops all elements to get the minimum characteristic length and gets the maximum or minimum parameters needed to calculate the critical time step.

Finally, the code multiplies this minimum by the Courant number $C_N$ given in the calculation parameters.