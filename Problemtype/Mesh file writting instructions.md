# Mesh file writing instructions      

The mesh in ADSIM is written in plain text files with extension `.mesh`. The structure follows the following order:

## Header
This indicates the version of the ADSIM used.
## Counters
The counters are written in the following format:

``` 
MESH numNodes Ntotal numElements NelmTotal
```
`Ntotal` is the total number of nodes

`NelmTotal` is the total number of elements

## Nodal coordinates

The file must be structured in the following manner:
```
coordinates
0.0 0.0
1.0 0.0
1.0 1.0
0.0 1.0
end coordinates
```
The coordinates must be wrapped on the `coordinates` and `end coordinates` text. The coordinates follow and the number of coordinates must be exactly `Ntotal`

## Elements

The file must be written in the following manner:
```
elements
1 2 3 4
4 5 6 7
2 5 6 3
end elements
```

The nodes ids corresponding that define the element conductivities are wrapped around `elements` and `end elements`. The order of the information implicitly defines the element id.

## Concentration boundary conditions

Use the following structure:
```
concentration_bc
counter
.
.
.
node_id value_gas1 value_gas2 ... value_gasn
.
.
.
end concentration_bc
```

The structure is wrapped as shown. Each line contains the node id and gas concentrations assigned to that node id. Eg. `2 0.5 0.25` indicates that node 2 has a concentration of 0.5 $\text{mol}/\text{L}^3$ of gas 1 and 0.25 $\text{mol}/\text{L}^3$ of gas 2. The `counter` is needed to initialize the arrays that will contain the data.

## Uniform flow boundary conditions
Use the following structure:
```
uniform_flow_bc
counter
.
.
.
node_id value_gas1 value_gas2 ... value_gasn
.
.
.
end uniform_flow_bc
```

Each line contains the node id within a line where the condition was applied. After this, the value of normal flow is set, but it will be recalculated in real time to consider length of elements. Eg. `2 1.0 0.0` indicates that node 2 has a flow of 0.5 $\text{mol}/\text{L}^2 \text{T}$ of gas 1 and 0.25 $\text{mol}/\text{L}^2 \text{T}$ of gas 2.

## Absolute pressure boundary condition

Use the following structure:
```
absolute_pressure
counter
.
.
.
NodeId value
.
.
.
ebd absolute_pressure
```
Here the value of pressure is assigned to the `NodeId`. Eg, `1 250` indicates node 1 has an absolute pressure of 250  $\text{F} /\text{L}^2$

## Initial gas concentrations

Use the following structure:
```
initial_concentrations
counter
.
.
.
ElmId value_gas1 value_gas2 ... value_nasn
.
.
.
end initial_concentrations
```
E.g., `5 0.5 0.2` indicates element 5 has an initial concentration on its nodes of 0.5 and 0.2 $\text{mol}/\text{L}^3$ for gas 1 and 2 respectively.

## Initial temperature

use the following structure:
```
initial_temperature
counter
.
.
.
ElmId value
.
.
.
end initial_temperature
```
E.g., `5 295` indicates element 5 has an initial temperature on its nodes of 295 $\text{T}_p$. Later, a subroutine will assign temperatures to nodes.

## Material assignation
Use the following structure:
```
material
.
.
.
ElmID MatIndex
.
.
.
end material
```
The structure consists of a list of elements with the material assigned to each. E.g., `12 1` indicates material index 1 is assigned to element 12.