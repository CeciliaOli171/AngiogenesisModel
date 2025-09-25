# Angiogenesis Model for Endometriotic Lesions 

## Introduction 

Our overall goal is to develop an agent-based model to describe endometriotic lesions' mechanics. We start by establishing an angiogenesis model for endometriotic lesions. Indeed, angiogenesis is one of the main phenomenon contributing to the implantation and proliferation of lesions. By establishing a vascular network, they have access to hormones and nutrients promoting their growth.	We choose the overlapping spheres method for our agent-based model. We describe the molecular concentrations around the lesion using reaction-diffusion PDEs. We use the open source simulation package Chaste [1] in C++ to perform simulations and the open source software Paraview [2] to view the results. 

## Agent-Based Model

In our model, we do not consider individual cells but sections of the vessel; however, for purposes of notation, we refer to them as cells. We represent the system using an Overlapping Spheres (OS) model where each cell is a sphere characterised by its radius and the coordinates of its centre. We consider two types of cells, vessel tips (VTs) and vessel segments (VSs). VTs start and end the vessel; they migrate in response to signals sent by the extra-cellular matrix and the chemical environment. SCs form the body of the vessel. Depending on the type of cells, different forces are applied. According to Newton’s second law and assuming negligible inertia, the motion of cells are described by the following system:
```math
\frac{dx}{dt} = \sum F_{ext}  
		= F_{R, c} + F_{M, c} + F_{H, c} + F_{P, c} + F_{A, c}, \; c \in \mathcal{C} 
```

### Cell Types

We distinguish three types of cell which have different mechanical and division properties: vessel tips (VTs), vessel segments (VSs), branching cells. 

VTs are a section of vasculature. They are composed of several type of cells and are located at the end of blood vessel branches. Activated by VEGF released by the lesion, they lead the movement of the vessel toward it.

VSs are located inside the blood vessels. They also are a part of the blood vessels and so composed by different type of cells. They follow the vessel tips movement, however, unlike the latter, they are not activated by any exterior factor.

Branching cells are VSs located at branching points 3.2. The difference between them and regular VSs resides in the following assumption: the mechanical force exerted by their three neighbours is high enough to neglect the angular force exerted by each branch. Their interaction with their micro-environment is described by the random force.

Additionally, we apply boundary conditions to the first cell starting the blood vessel - that we can consider as the centre of mass of the system. This cell is fully constrained, thus, it is not moving.

### Cell Cycle and Division 

VTs and VSs differ also by their dividing property: only VTs divide. According to a probability of sprouting ($P_{sprout}$), they divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division). Literature shows that VEGF influences endothelial cells proliferation [Cristofanilli2002,Costa2016,Becit2001]. VEGF concentration needs to reach a certain threshold in order to activate Ang-2. The combination of VEGF and Ang-2 leads to the creation of new sprouts [RamsauerDAmore2002]. Therefore the sprouting probability $P_{sprout}$ is a function of Ang-2 concentration $a$ at a time $t$ and position $x$: 
```math
		P_{sprout}(x, t) = \lambda \frac{c(x, t)^n}{c(x, t)^n+ K^n} 
```
with $\lambda$ the maximum sprouting rate.


In our model, VTs can divide into either a VS or another VT. When the daughter cell is a VS, it is added in the vessel branch. When it divides into a VT, the two daughter cells are positioned on the same plane, with the same distance separating them from the closest VS. This VS becomes a branching point; we say that the vessel is sprouting or that it creates a new sprout. In order to do that, we follow the following steps:

1. we start by calculating the vector $x_{branch}$ describing the branch direction using the coordinates of the parent before division cell $x_{parent}$ and of its closest neighbour in the branch $x_{neighbour}$ which gives us the following expression:
  ```math
   x_{branch} = x_{neighbour} - x_{parent} = (a, b, c)^T, \; a, b, c \in \mathbb{R}
  ```
   
3. from the $x_{branch}$ vector, we can calculate a perpendicular vector:
  ```math
  x_{perpendicular} = \lambda (-b, a, 0)^T + \mu (0, -c, b)^T + \nu (-c, 0, a)^T
  ```
with $\lambda, \mu, \nu$ chosen randomly in the interval $\[-1, 1\]$.

3. we obtain finally the daughter cells coordinates:
   ```math
   x_{daughter}^{1} = x_{parent} - 0.5 x_{perpendicular} + 0.5 x_{branch}
   ```
   ```math
   x_{daughter}^2 = x_{parent} + 0.5 x_{perpendicular} + 0.5 x_{branch}
   ```


When the VT is ready to divide, if the distance between it and its closest neighbour is less than a fixed length $l_{ref}$, then it does not divide and connects with the closest neighbour instead, creating a loop and a new branching point. This reflects a behaviour known as anastomosis, which enables the flow of blood through the vascular network. 


The cell cycle duration is stochastic, its value is selected using a uniform distribution $U(a,b)$. In our model, we choose to use $a= 12.0h$ and $b= 14.0h$ independently of any time observed experimentally.


### Forces 

We define:

• $\mathcal{C}$ the set containing all the cells of the system 

• $\mathcal{P}$ the subset containing all the vessel tips

• $\mathcal{S}$ the subset containing all the vessel segments


#### Mechanical Force

Cells interact with each other only if their centres are separated by less than a cut-off distance $l_c$. Only cells in the same branch experience this force; branching segments belong to all the branches they are connected to. For all cells, we use the following mechanical force [Meineke2001,Germano2023,Osborne2017]: 


```math 
\sum_{i \in \mathcal{N}} S_d l_d \log(1 + \frac{x_{i,c}}{l_d}) \frac{x_i - x_c}{|x_i - x_c|} if l_{i,c} - l_d < |x_i - x_c| < l_{i,c}
```

```math
F_{M, c} = \sum_{i \in \mathcal{N}} S_d x_{i, c} \exp(- \rho \frac{x_{i, c}}{l_d}) \frac{x_i - x_c}{|x_i - x_c|} if l_{i,c} < |x_i - x_c| < l_c
```

```math
\textbf{0} if |x_i - x_c| > l_c 
```

with $c \in \mathcal{C}$ and $x_{n,c} = |x_n - x_c| - l_{n, c}$. $\mathcal{N}$ is the set of the node's neighbours in the same branch, $S_d$ is the spring stiffness, $l_d$ is the cut-off distance for the connectivity of the nodes and $l_{n,c}$ is the rest length of the spring between the cells $n$ and $c$. 


#### Random Force

All the cells are interacting with their environment and those interactions are represented by a random force:

```math
	F_{R, c} =  \sigma \mathbf{\xi}_c, \; c \in \mathcal{C}
```

where $\sigma$ is the sensitivity of the cell to random fluctuations and $\xi_c$ is a random unit vector chosen from a uniform distribution. 


#### Chemotactic Force

VT movement is subject to the molecular factors present around them. Among them, VEGF is responsible for their migration toward the lesion [Donnez1998, Takehara2004]. For VT, we use the following chemotactic force:
```math
	F_{H, p} = \chi H(\nabla c, c, h), \; p \in \mathcal{P}
```
where $\chi$ is the chemotactic sensitivity, $c$ the VEGF concentration that is time and space dependent and $h$ is a fixed gradient of other chemoattractant factors. 


#### Persistence Force

VT move along the ECM fibres [DaubMerks2013,YouStallcup2017,McLaren2000]. For VT, we use a persistence force from Perfhal et al. [Perfhal2017]:
```math
	F_{P, p} = \omega_p \frac{x_p(t) - x_p(t-\tau)}{|x_p(t) - x_p(t-\tau)|}, \; p \in \mathcal{P}
```
where $\omega_p$ is the persistence coefficient and $\tau$ is the directional persistence time.


#### Angular Force

The angular force describes the interactions between the VS and its micro-environment. VS follow the VT movement and therefore align with the \acrshort{ecm} fibres [DaubMerks2013,YouStallcup2017,McLaren2000]. This force is not applied to branching segments. For VS, we use an angular force from Perfhal et al. [Perfhal2017]: 
```math
	F_{A, s} = \omega_a (\alpha_{angular} - \pi) \frac{(x_b - x_s) + (x_c - x_s)}{|(x_b - x_s) + (x_c - x_s)|}, \; s \in \mathcal{S}, \; b,c \in \mathcal{C}-\{s\} 
```
where $\omega_a$ is the angular spring constant and $\alpha_{angular}$ the angle made by the VS considered and its two closest neighbours in the same branch.


## Molecular Concentrations Model

### Vascular Endothelial Growth Factor (VEGF) PDE
VEGF is released by the lesion and activates the VTs so they migrate toward the lesion. We describe the VEGF concentration $c(x,t)$ at position $x$ and time $t$ by the following PDE:
```math
  \frac{\partial c}{\partial t}(x,t) = D_c \Delta c(x,t) - \epsilon_c n(x,t) c(x,t) - (M_c-A_c) c(x,t)
```
where $D_c$ is the diffusion coefficient, $\epsilon_c$ the consumption rate by \acrshort{vt}s, $A_c$ the creation rate, $M_c$ the decay rate and $n(x,t)$ the density of \acrshort{vt} at position $x$ and time $t$. 


### Initial and Boundary Conditions
we consider that initially the concentration of \acrshort{vegf} in the \acrshort{ecm} is equal to the baseline concentration $c_0$ since the lesion has not started to diffuse \acrshort{vegf} yet. We assume that the \acrshort{vegf} concentration is maximal equal to $c_{\max}$ at the lesion border $\Gamma$ at all times therefore we use a Dirichlet boundary condition. We summarise those conditions in the following system:
```math
    c(x, 0) &= c_0, \forall \x \notin \Gamma, (IC) \\ 
    c(x \in \Gamma, t) &= c_{\max}, \forall t \geq 0, (BC)
```


### Steady-State Model
If we consider only the VEGF concentration and that the diffusion is the same in the whole $y$ and $z$ axis, then we can consider a steady state for the PDE in one dimension. As in Daub and Merks [DaubMerks2013], we neglect the consumption term by the VTs. We obtain a one-dimensional ODE, similar to Daub and Merks [DaubMerks2013], that can be solved analytically when decay is greater than creation (i.e. $M_c > A_c$):
```math
\forall x \in \Omega, c(x) = (c_{\max}-c_0) e^{-K_c x} + c_0 with \; \; K_c = \sqrt{\frac{M_c-A_c}{D_c}}
```
    

## Blood Flow Model


## Coupling Models

### Cell Division 

Experimental results show that there is an increase in VEGF and Ang-2 concentration and number of branches close to the lesion [Donnez1998,Hur2006,Jingting2008]. We note $c_{\min}$ the minimal VEGF concentration value necessary in order to activate Ang-2 and to observe the first branches in the blood network. We consider for the function $F$ in $P_{sprout}$ equation, a Hill function with a Hill coefficient $n$ and half saturation constant $K$: 
```math
    P_{sprout}{(x, t)} = F(c, x, t) = \lambda \frac{a{(x, t)}^n}{a{(x, t)}^n+K_a^n} = \lambda \frac{c{(x, t)}^n}{c{(x, t)}^n+K^n}
```
with $K = \frac{K_a}{A}$.


As we do not include the consumption of VEGF in our system, we account for the expected reduced level of VEGF local to a sprouting event by introducing a minimal length $l_B$ between two consecutive sprouts. Similarly, anastomosis can happen only if the branch the vessel tip leading has also a minimal length $l_B$.


### Chemotactic Force

For the lesion hypothesis, VT are moving along the VEGF gradient, therefore when we couple the two models, the equation for the chemotactic force becomes:
```math
F_{H, p} = \chi (\nabla c+{(-h, 0, 0)}^T), p \in \mathcal{P}
```


## Code 

### Image Analysis Pipeline

To be modified and completed for the comparison between the mathematical model and the image analysis.

Quantitative comparison: 
* number of cells in one plane 
* neighbouring distance 

### Model Analysis 

ParametersSensitivityRunner: class written for the model analysis. 

Model analysis files:
* Sprouting probability tests
* Forces coefficients tests

### Source 

Source files:
* Forces: angular, chemotactic, directional persistence, linear mechanical force, persistence, random
* Cell cycle: sprouting, anastomosis (comprised inside daughter cell file)
* Cell mutation: branching cells, daughter cell, pinned cell (BC), tip cell, vessel segment
* Writers (used for the analysis of the model): birth time cells, branches number, consecutive branches, tortuosity

### Test

Test files:
* Test of all the forces in the model individually
* Test of the model in 2D and 3D 

## Future Goals 
 
1. couple this angiogenesis model with blood flow model 

## Bibliography

[1] G. R. Mirams, C. J. Arthurs, M. O. Bernabeu, R. Bordas, J. Cooper, A. Corrias, Y. Davit, S.-J. Dunn, A. G. FleVTher, D. G. Harvey, M. E. Marsh, J. M. Osborne, P. Pathmanathan, J. Pitt-Francis, J. Southern, N. Zemzemi, and D. J Gavaghan. Chaste: An open source c++ library for computational physiology and biology. PLOS Computational Biology, 9(3)(e1002970), 2013.

[2] J. Ahrens, B. Geveci, and C. Law. Paraview: An end-user tool for large-data visualization. in visualization handbook. Elsevier, pages 717–731, 2005.

[35] F. A. Meineke, C. S. Potten, and M. Loeffler. Cell migration and organization in the intestinal crypt using a lattice-free model. Cell Proliferation, 34(4):253–266, 2001.

[42] H. Perfahl, B. D. Hughes, T. Alarco ́n, P. K. Maini, M. C. Lloyd, M. Reuss, and H. M Byrne. 3d hybrid modelling of vascular network formation. Journal of Theoretical Biology, 414:254–268, 2017.

[43] J. Pitt-Francis, P. Pathmanathan, M. O. Bernabeu, R. Bordas, J. Cooper, A. G. Fletcher, G. R. Mirams, P. Murray, J. M. Osborne, A. Walter, S. J. Chapman, A. Garny, I. M. M. Van Leeuwen, P. K. Maini, B. Rodr ́ıguez, S. L. Waters, J. P. Whiteley, H. M. Byrne, and D. J Gavaghan. Chaste: A test-driven approach to software development for biological modelling. Computer Physics Communica- tions, 180(12):2452–2471, 2009.
