# Angiogenesis Model for Endometriotic Lesions 

## Introduction 

Our overall goal is to develop an agent-based model to describe endometriotic lesions' mechanics. We start by establishing an angiogenesis model for endometriotic lesions. Indeed, angiogenesis is one of the main phenomenon contributing to the implantation and proliferation of lesions. By establishing a vascular network, they have access to hormones and nutrients promoting their growth.	We choose the overlapping spheres method for our agent-based model. We describe the molecular concentrations around the lesion using reaction-diffusion PDEs. We use the open source simulation package Chaste [37] in C++ to perform simulations and the open source software Paraview [2] to view the results. 

## Methodology

In our model, we do not consider individual cells but sections of the vessel; however, for purposes of notation, we refer to them as cells. We represent the system using an Overlapping Spheres (OS) model where each cell is a sphere characterised by its radius and the coordinates of its centre. We consider two types of cells, Tip Cells (TCs) and Stalk Cells (SCs). TCs start and end the vessel; they migrate in response to signals sent by the extra-cellular matrix and the chemical environment. SCs form the body of the vessel. Depending on the type of cells, different forces are applied. According to Newton’s second law and assuming negligible inertia, the motion of cells are described by the following system:
```math
\frac{dx}{dt} = \sum F_{ext}  
		= F_{R, c} + F_{M, c} + F_{H, c} + F_{P, c} + F_{A, c}, \; c \in \mathcal{C} 
```

### Cell Types

We distinguish three types of cell which have different mechanical and division properties: Tip Cells (TCs), Stalk Cells (SCs), branching cells. 

TCs are a section of vasculature. They are composed of several type of cells and are located at the end of blood vessel branches. Activated by VEGF released by the lesion, they lead the movement of the vessel toward it.

SCs are located inside the blood vessels. They also are a part of the blood vessels and so composed by different type of cells. They follow the tip cells movement, however, unlike the latter, they are not activated by any exterior factor.

Branching cells are SCs located at branching points 3.2. The difference between them and regular SCs resides in the following assumption: the mechanical force exerted by their three neighbours is high enough to neglect the angular force exerted by each branch. Their interaction with their micro-environment is described by the random force.

Additionally, we apply boundary conditions to the first cell starting the blood vessel - that we can consider as the centre of mass of the system. This cell is fully constrained, thus, it is not moving.

### Cell Cycle and Division 

We use a stochastic cell cycle implemented within Chaste as “Uniform Cell Cycle” [37], [43]. The cell cycle duration is stochastic, its value is selected using a uniform distribution on an interval [a,b], where a and b are respectively the minimum and maximum length of cell cycle. In our model, a = 12.0 and b = 14.0. If the cell does not divide, then the cell cycle duration is set to be infinite [37], [43].

In our model, TCs can divide into either a SC or another TC. When the daughter cell is a SC, it is added in the vessel branch. When it divides into a TC, the two daughter cells are positioned on the same plane, with the same distance separating them from the closest SC. This SC becomes a branching point; we say that the vessel is sprouting or that it creates a new sprout. In order to do that, we follow the following steps:

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

We set up the following division property: the leading cell of the branch divides. To decide if the resulting daughter cell discretises into a TC or a SC, we set a probability to create a new sprout $=P_{sprout}$. The probability of sprouting depends on the VEGF concentration. The closer the tip cells get to the lesion, the more they will divide into new branches. For a position $\textbf{x}$ and a time $t$, we consider four different formulas:
```math
		P_{sprout}(\textbf{x}, t) &= \lambda c(\textbf{x}, t) \; \; \text{linear function}
		P_{sprout}(\textbf{x}, t) &= \lambda \frac{c(\textbf{x}, t)}{c(\textbf{x}, t)+ K} \; \; \text{hill function ($n=1$)}
		P_{sprout}(\textbf{x}, t) &= \lambda \frac{c(\textbf{x}, t)^2}{c(\textbf{x}, t)^2+ K^2} \; \; \text{hill function ($n=2$)}
		P_{sprout}(\textbf{x}, t) &= \lambda \frac{c(\textbf{x}, t)^3}{c(\textbf{x}, t)^3+ K^3} \; \; \text{hill function ($n=3$)}
```
with $\lambda$ the maximum sprouting rate.

We consider in our model the apoptosis event, necessary to implement the blood flow in our blood network. If the mature TC has a distance smaller than a cut-off length $l_{merge}$, then two scenarios can happen. If the neighbour cell is another TC, then the two branches merge and create a loop; the TCs become SCs. If the neighbour cell is a SC, then the cells connect; the TC becomes a SC and the neighbour cell becomes a branching cell (if it was not already the case).


We consider in our model the anastomosis event, necessary to implement the blood flow in our blood network. If the mature TC has a distance smaller than a cut-off length $l_{merge}$, then two scenarios can hapen. If the neighbour cell is another TC, then the two branches merge and create a loop; the TCs become SCs. If the neighbour cell is a SC, then the cells connect; the TC becomes a SC and the neighbour cell becomes a branching cell if it was not already the case. 


### Molecular Concentrations 

We write $n(\textbf{x},t)$ the density of TCs at position $\textbf{x}$ and time $t$.

#### Vascular Endothelial Growth Factor (VEGF)
VEGF is released by the lesion and activates the TCs so they migrate toward the lesion. We describe the VEGF concentration $c(\textbf{x},t)$ at position $\textbf{x}$ and time $t$ by the following PDE:
```math
  \frac{\partial c}{\partial t}(\textbf{x},t) &= D_c \Delta c(\textbf{x},t) - \epsilon_c n(\textbf{x},t) c(\textbf{x},t) + (A_c-M_c) c(\textbf{x},t) \\
  c(\textbf{x} \in \Gamma, t) &= c_{max}, \; t>0 && \text{(BC)}  \\
  c(\textbf{x}, 0) &= c_0(\x) && \text{(IC)} 
```
with $D_c$ corresponds to the diffusion coefficient, $A_c$ is the creation coefficient, $M_c$ is the decay coefficient and $\epsilon_c$ is the consumption by TCs coefficient. We consider that the VEGF concentration is equal to $c_0$ at the border $\Gamma$ of the endometriotic lesion and equal to $0$ in the other borders of the mesh.

If we consider only the VEGF concentration and that the diffusion is the same in the whole y and z axis, then we can consider a steady state of the PDE in 1D and neglect the consumption term by the tip cells (which is low compared to the other term of the equation). We obtain a simple 1D ODE that can be solved analytically: 
```math 
c(x) = c_{max} \frac{1}{1 - e^{-2 \sqrt{\frac{M_c - A_c}{D_c}} l_{max}}} (e^{- \sqrt{\frac{M_c - A_c}{D_c}} x} - e^{\sqrt{\frac{M_c - A_c}{D_c}} (x-2 l_{max})})
```


#### Matrix-Metallo Proteinase (MMP)
To be continued. 


### Forces 

We define:

• $\mathcal{C}$ the set containing all the cells of the system 

• $\mathcal{P}$ the subset containing all the tip cells

• $\mathcal{S}$ the subset containing all the stalk cells


#### Mechanical Force

We suppose that the cell connectivity can be represented by a spring. Cells interact with each other only if their centre are separated by less than a cut-off distance $l_d$. We use the  mechanical force already implemented in Chaste. We add the following condition:  if two cells belong to the same branch or at least one of them is a branching cell then they are connected and exert a force on each other. For all cells, we use Meineke et al. [35] mechanical force: 

$F_{M, c} =$

```math 
\sum_{n \in \mathcal{N}} S_d l_d \log(1 + \frac{|x_n - x_c| - l_{n, c}}{l_{n, c}}) if |x_n - x_c| < l_{n,c}
```

```math
\sum_{n \in \mathcal{N}} S_d (|x_n - x_c| - l_{n, c}) \exp(-\frac{|x_n - x_c| - l_{n, c}}{2 l_d}) if l_{n,c} < |x_n - x_c| < l_d
```

```math
\textbf{0} if |x_n - x_c| > l_d
```

with $c \in \mathcal{C}$. $\mathcal{N}$ is the set of the node's neighbours in the same branch, $S_d$ is the spring stiffness, $l_d$ is the cut-off distance for the connectivity of the nodes and $l_{n,c}$ is the rest length of the spring between the cells $n$ and $c$. 

#### Random Force

All the cells are interacting with their environment and those interactions are represented by a random force. For all cells, the random force is:

```math
	F_{R, c} =  \sigma \mathbf{\xi}_c, \; c \in \mathcal{C}
```

with $\sigma$ the sensitivity of the cell to random fluctuations and $\mathbf{\xi}_c$ a random unit vector. 

#### Chemotactic Force

TC's movement is subject to the chemicals present around them. In the case of endometriosis, the lesions are releasing VEGF that induces the blood vessel growth towards it. For TCs, the chemotactic force is:
```math
	F_{H, p} = \chi (c_x, 0, 0)^T, \; p \in \mathcal{P}
```
with $\chi$ the chemotactic sensitivity and $c_x$ the VEGF gradient.  

#### Persistence Force

TCs tend to move toward the same direction through the extra-cellular matrix fibres. The persistence force describe this motion. For TCs, the persistence force is:
```math
	F_{P, p} = \omega_p \frac{x_p(t) - x_p(t-\tau)}{|x_p(t) - x_p(t-\tau)|}, \; p \in \mathcal{P}
```
with $\omega_p$ the persistence coefficient and $\tau$ the timescale.

#### Angular Force

The angular force describes the interactions between the \acrshort{sc} and its micro-environment which makes it stabilise and align with its neighbours. For SCs, the angular force is: 
```math
	F_{A, s} = \omega_a (\alpha_{angular} - \pi) \frac{(x_b - x_s) + (x_c - x_s)}{|(x_b - x_s) + (x_c - x_s)|}, \; s \in \mathcal{S}, \; b,c \in \mathcal{C}-\{s\} 
```
with $\omega_a$ the angular spring constant. $\alpha_{angular}$ is the angle made by the SC considered and its two closest neighbours in the same branch. 

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
 
1. continue the quantitative analysis to establish the parameters.
2. adapt the model to endometriotic lesions.
3. couple this angiogenesis model with blood flow model for resulting nutrient transport to the lesion.

## Bibliography

[2] J. Ahrens, B. Geveci, and C. Law. Paraview: An end-user tool for large-data visualization. in visualization handbook. Elsevier, pages 717–731, 2005.

[35] F. A. Meineke, C. S. Potten, and M. Loeffler. Cell migration and organization in the intestinal crypt using a lattice-free model. Cell Proliferation, 34(4):253–266, 2001.

[37] G. R. Mirams, C. J. Arthurs, M. O. Bernabeu, R. Bordas, J. Cooper, A. Corrias, Y. Davit, S.-J. Dunn, A. G. Fletcher, D. G. Harvey, M. E. Marsh, J. M. Osborne, P. Pathmanathan, J. Pitt-Francis, J. Southern, N. Zemzemi, and D. J Gavaghan. Chaste: An open source c++ library for computational physiology and biology. PLOS Computational Biology, 9(3)(e1002970), 2013.

[42] H. Perfahl, B. D. Hughes, T. Alarco ́n, P. K. Maini, M. C. Lloyd, M. Reuss, and H. M Byrne. 3d hybrid modelling of vascular network formation. Journal of Theoretical Biology, 414:254–268, 2017.

[43] J. Pitt-Francis, P. Pathmanathan, M. O. Bernabeu, R. Bordas, J. Cooper, A. G. Fletcher, G. R. Mirams, P. Murray, J. M. Osborne, A. Walter, S. J. Chapman, A. Garny, I. M. M. Van Leeuwen, P. K. Maini, B. Rodr ́ıguez, S. L. Waters, J. P. Whiteley, H. M. Byrne, and D. J Gavaghan. Chaste: A test-driven approach to software development for biological modelling. Computer Physics Communica- tions, 180(12):2452–2471, 2009.
