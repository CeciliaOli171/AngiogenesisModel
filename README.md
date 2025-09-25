# Angiogenesis Model for Endometriotic Lesions 

## Introduction 

We develop an agent-based model to describe endometriotic lesions' mechanics. We start by establishing an angiogenesis model for endometriotic lesions. Indeed, angiogenesis is one of the main phenomenon contributing to the implantation and proliferation of lesions. By establishing a vascular network, they have access to hormones and nutrients promoting their growth.	We choose the overlapping spheres method for our agent-based model. We describe the molecular concentrations around the lesion using reaction-diffusion PDEs. We use the open source simulation package Chaste [1,2] in C++ to perform simulations and the open source software Paraview [3] to view the results. 

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

Branching cells are VSs located at branching points. The difference between them and regular VSs resides in the following assumption: the mechanical force exerted by their three neighbours is high enough to neglect the angular force exerted by each branch. Their interaction with their micro-environment is described by the random force.

Additionally, we apply boundary conditions to the first cell starting the blood vessel - that we can consider as the centre of mass of the system. This cell is fully constrained, thus, it is not moving.

### Cell Cycle and Division 

VTs and VSs differ also by their dividing property: only VTs divide. According to a probability of sprouting ($P_{sprout}$), they divide into either a VS and a VT, lengthening the vessel (asymmetric division); or into two VTs, creating a new sprout (symmetric division). Literature shows that VEGF influences endothelial cells proliferation [4,5,6]. VEGF concentration needs to reach a certain threshold in order to activate Ang-2. The combination of VEGF and Ang-2 leads to the creation of new sprouts [7]. Therefore the sprouting probability $P_{sprout}$ is a function of Ang-2 concentration $a$ at a time $t$ and position $x$: 
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

Cells interact with each other only if their centres are separated by less than a cut-off distance $l_c$. Only cells in the same branch experience this force; branching segments belong to all the branches they are connected to. For all cells, we use the following mechanical force [8,9,10]: 


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

VT movement is subject to the molecular factors present around them. Among them, VEGF is responsible for their migration toward the lesion [11,12]. For VT, we use the following chemotactic force:
```math
	F_{H, p} = \chi H(\nabla c, c, h), \; p \in \mathcal{P}
```
where $\chi$ is the chemotactic sensitivity, $c$ the VEGF concentration that is time and space dependent and $h$ is a fixed gradient of other chemoattractant factors. 


#### Persistence Force

VT move along the ECM fibres [13,14,15]. For VT, we use a persistence force from Perfhal et al. [16]:
```math
	F_{P, p} = \omega_p \frac{x_p(t) - x_p(t-\tau)}{|x_p(t) - x_p(t-\tau)|}, \; p \in \mathcal{P}
```
where $\omega_p$ is the persistence coefficient and $\tau$ is the directional persistence time.


#### Angular Force

The angular force describes the interactions between the VS and its micro-environment. VS follow the VT movement and therefore align with the ECM fibres [13,14,15]. This force is not applied to branching segments. For VS, we use an angular force from Perfhal et al. [16]: 
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
where $D_c$ is the diffusion coefficient, $\epsilon_c$ the consumption rate by VTs, $A_c$ the creation rate, $M_c$ the decay rate and $n(x,t)$ the density of VT at position $x$ and time $t$. 


### Initial and Boundary Conditions
we consider that initially the concentration of VEGF in the ECM is equal to the baseline concentration $c_0$ since the lesion has not started to diffuse VEGF yet. We assume that the VEGF concentration is maximal equal to $c_{\max}$ at the lesion border $\Gamma$ at all times therefore we use a Dirichlet boundary condition. We summarise those conditions in the following system:
```math
    c(x, 0) = c_0, \forall \x \notin \Gamma, (IC) 
```
```math
    c(x \in \Gamma, t) = c_{\max}, \forall t \geq 0, (BC)
```


### Steady-State Model
If we consider only the VEGF concentration and that the diffusion is the same in the whole $y$ and $z$ axis, then we can consider a steady state for the PDE in one dimension. As in Daub and Merks [13], we neglect the consumption term by the VTs. We obtain a one-dimensional ODE, similar to Daub and Merks [13], that can be solved analytically when decay is greater than creation (i.e. $M_c > A_c$):
```math
\forall x \in \Omega, c(x) = (c_{\max}-c_0) e^{-K_c x} + c_0 with \; \; K_c = \sqrt{\frac{M_c-A_c}{D_c}}
```
    

## Blood Flow Model


## Coupling Models

### Cell Division 

Experimental results show that there is an increase in VEGF and Ang-2 concentration and number of branches close to the lesion [11,17,18]. We note $c_{\min}$ the minimal VEGF concentration value necessary in order to activate Ang-2 and to observe the first branches in the blood network. We consider for the function $F$ in $P_{sprout}$ equation, a Hill function with a Hill coefficient $n$ and half saturation constant $K$: 
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

### Model Analysis 

* Number of branches
* Time to reach lesion 
* Percentage of vascularisation
* Cell density
* Network velocity
* Expected number of branches 

### Source 

* Forces: angular, chemotactic, directional persistence, linear mechanical force, persistence, random
* Cell cycle: sprouting, anastomosis (comprised inside daughter cell file)
* Cell mutation: branching cells, daughter cell, pinned cell (BC), tip cell, vessel segment
* Writers (used for the analysis of the model): birth time cells, branches number, consecutive branches, tortuosity

### Test

* Test of all the forces in the models individually (TestForcesModel.hpp)
* Test of the angiogenesis model in 2D and 3D (TestAngiogenesisModel.hpp)
* Test of the molecular model (TestSolvingMolecularConcentrationsPdes.hpp)
* Test of the hybrid model with a constant VEGF concentration (TestAngiogenesisWithVegfConcentrationConstant.hpp)
* Test of the hybrid model in the steady-state case of VEGF concentration PDE (TestAngiogenesisWithVegfConcentrationAnalyticalApproximationOfPde.hpp)
* Test of the hybrid model with a VEGF concentration PDE(TestAngiogenesisWithVegfConcentrationConstant.hpp)

### Image Analysis Pipeline

To be modified and completed for the comparison between the mathematical model and the image analysis.

Quantitative comparison: 
* number of cells in one plane 
* neighbouring distance 


## Future Goals 
 
1. couple this angiogenesis model with blood flow model 
2. write a pipeline between image analysis and model to validate it 
3. adapt the model for endometrial angiogenesis 

## Bibliography

[1] Mirams, G. R., Arthurs, C. J., Bernabeu, M. O., Bordas, R., Cooper, J., Corrias, A., Davit, Y., Dunn, S.-J., Fletcher, A. G., Harvey, D. G., Marsh, M. E., Osborne, J. M., Pathmanathan, P., Pitt-Francis, J., Southern, J., Zemzemi, N., & Gavaghan, D. J. (2013). Chaste: An Open Source C++ Library for Computational Physiology and Biology. PLOS Computational Biology, 9(3), e1002970. https://doi.org/10.1371/journal.pcbi.1002970

[2] Pitt-Francis, J., Pathmanathan, P., Bernabeu, M. O., Bordas, R., Cooper, J., Fletcher, A. G., Mirams, G. R., Murray, P., Osborne, J. M., Walter, A., Chapman, S. J., Garny, A., Van Leeuwen, I. M. M., Maini, P. K., Rodríguez, B., Waters, S. L., Whiteley, J. P., Byrne, H. M., & Gavaghan, D. J. (2009). Chaste: A test-driven approach to software development for biological modelling. Computer Physics Communications, 180(12), 2452–2471. https://doi.org/10.1016/j.cpc.2009.07.019

[3] Ahrens, J., Geveci, B., & Law, C. (2005). ParaView: An End-User Tool for Large-Data Visualization. In Visualization Handbook (pp. 717–731). Elsevier. https://doi.org/10.1016/B978-012387582-2/50038-1

[4] Cristofanilli, M., Charnsangavej, C., & Hortobagyi, G. N. (2002). Angiogenesis modulation in cancer research: Novel clinical approaches. Nature Reviews Drug Discovery, 1(6), 415–426. https://doi.org/10.1038/nrd819

[5] Costa, G., Harrington, K. I., Lovegrove, H. E., Page, D. J., Chakravartula, S., Bentley, K., & Herbert, S. P. (2016). Asymmetric division coordinates collective cell migration in angiogenesis. Nature Cell Biology, 18(12), 1292–1301. https://doi.org/10.1038/ncb3443

[6] Becit, N., Ceviz, M., Koçak, H., Yekeler, İ., Ünlü, Y., Çelenk, Ç., & Akın, Y. (2001). The Effect of Vascular Endothelial Growth Factor on Angiogenesis.An Experimental Study. European Journal of Vascular and Endovascular Surgery, 22(4), 310–316. https://doi.org/10.1053/ejvs.2001.1468

[7] Ramsauer, M., & D’Amore, P. A. (2002). Getting Tie(2)d up in angiogenesis. Journal of Clinical Investigation, 110(11), 1615–1617. https://doi.org/10.1172/JCI0217326

[8] Meineke, F. A., Potten, C. S., & Loeffler, M. (2001). Cell migration and organization in the intestinal crypt using a lattice‐free model. Cell Proliferation, 34(4), 253–266. https://doi.org/10.1046/j.0960-7722.2001.00216.x

[9] Germano, D. P. J., Zanca, A., Johnston, S. T., Flegg, J. A., & Osborne, J. M. (2022). Free and interfacial boundaries in individual-based models of multicellular biological systems [Preprint]. Biophysics. https://doi.org/10.1101/2022.12.13.520331

[10] Osborne, J. M., Fletcher, A. G., Pitt-Francis, J. M., Maini, P. K., & Gavaghan, D. J. (2017). Comparing individual-based approaches to modelling the self-organization of multicellular tissues. PLOS Computational Biology, 13(2), e1005387. https://doi.org/10.1371/journal.pcbi.1005387

[11] Donnez, J., Smoes, P., Gillerot, S., Casanas-Roux, F., & Nisolle, M. (1998). Vascular endothelial growth factor (VEGF) in endometriosis. Human Reproduction, 13(6), 1686–1690. https://doi.org/10.1093/humrep/13.6.1686

[12] Takehara, M., Ueda, M., Yamashita, Y., Terai, Y., Hung, Y.-C., & Ueki, M. (2004). Vascular endothelial growth factor A and C gene expression in endometriosis. Human Pathology, 35(11), 1369–1375. https://doi.org/10.1016/j.humpath.2004.07.020

[13] Daub, J. T., & Merks, R. M. H. (2013). A Cell-Based Model of Extracellular-Matrix-Guided Endothelial Cell Migration During Angiogenesis. Bulletin of Mathematical Biology, 75(8), 1377–1399. https://doi.org/10.1007/s11538-013-9826-5

[14] You, W.-K., & Stallcup, W. (2017). Localization of VEGF to Vascular ECM Is an Important Aspect of Tumor Angiogenesis. Cancers, 9(8), 97. https://doi.org/10.3390/cancers9080097 

[15] McLaren, J. (2000). Vascular endothelial growth factor and endometriotic angiogenesis. Human Reproduction Update, 6(1), 45–55. https://doi.org/10.1093/humupd/6.1.45 

[16] Perfahl, H., Hughes, B. D., Alarcón, T., Maini, P. K., Lloyd, M. C., Reuss, M., & Byrne, H. M. (2017). 3D hybrid modelling of vascular network formation. Journal of Theoretical Biology, 414, 254–268. https://doi.org/10.1016/j.jtbi.2016.11.013

[17] Hur, S. E. (2006). Angiopoietin-1, angiopoietin-2 and Tie-2 expression in eutopic endometrium in advanced endometriosis. Molecular Human Reproduction, 12(7), 421–426. https://doi.org/10.1093/molehr/gal049

[18] Jingting, C., Yangde, Z., Yi, Z., Mengxiong, L., Rong, Y., Yu, Z., Guoqing, P., & Lixiu, P. (2008). Expression of heparanase and angiopoietin-2 in patients with endometriosis. European Journal of Obstetrics & Gynecology and Reproductive Biology, 136(2), 199–209. https://doi.org/10.1016/j.ejogrb.2006.09.018