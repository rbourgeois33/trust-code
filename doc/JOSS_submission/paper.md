---
title: 'TRUST: An HPC platform for thermohydraulic applications'
tags:
  - thermohydraulic
  - cfd
  - c++
  - hpc
  - gpu
  - nuclear
authors:
  - name: Adrien Bruneton
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Pierre Ledac
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Elie Saikali
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Anida Khizar
    equal-contrib: true
    affiliation: 1
  - name: Luc Lecointre
    equal-contrib: true
    affiliation: 1
  - name: Guillaume Jomée
    equal-contrib: true
    affiliation: 1
  - name: Rémi Bourgeois
    equal-contrib: true
    affiliation: 1
  - name: Paul Zehner
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: LCAN, Research engineer, CEA, France
   index: 1

date: 25 September 2025
bibliography: paper.bib

---

# TRUST: An Open Source Thermohydraulic Platform from CEA

## 1. Introduction & Historical Background

**TRUST** is an open-source numerical platform developed at the **French Alternative Energies and Atomic Energy Commission (CEA)**.  
It originates from decades of development of in-house codes at CEA for **nuclear reactor thermal-hydraulics** and **computational fluid dynamics (CFD)** applications. In order to consolidate expertise and provide a modern, extensible framework, CEA launched TRUST as a unified, open, and modular platform.

The project was opened to the community to foster collaboration between research institutes, academia, and industry. Today, TRUST is used both as a **production-grade simulation tool** and as a **research platform** to test new models, discretization schemes, and parallel algorithms.

---

## 2. Objectives & Philosophy

The guiding principles behind TRUST are:

- **Generality:** support for a wide range of flow regimes, from incompressible to compressible, single-phase to multiphase.
- **Modularity:** new physics, discretization methods, and numerical schemes can be plugged in without rewriting the entire code.
- **Performance:** scalable execution on HPC architectures through MPI and, more recently, GPU acceleration.
- **Openness:** released under a BSD license, making it accessible for academic, industrial, and collaborative projects.
- **Validation:** a large set of benchmarks and reference cases ensure physical fidelity and numerical robustness.

---

## 3. Main Features & Physics Capabilities

### 3.1 Flow Regimes
- **Incompressible flows:** robust solvers for low-Mach flows, suitable for hydrodynamic studies where density is nearly constant.
- **Compressible flows:** capabilities extended to handle moderate to high Mach numbers, enabling aerodynamics and shock-related applications.
- **Multiphase flows:** models for interfacial flows, stratified regimes, and phase change phenomena, with ongoing developments for nuclear safety studies.

### 3.2 Thermohydraulics & Heat Transfer
- Conjugate heat transfer between fluids and solids.
- Energy equations with various source terms (radiation, conduction, convection).
- Interfaces to external thermodynamic property libraries such as **CoolProp** or internal CEA databases.

### 3.3 Numerical Methods
- Support for several discretization families:
  - Finite Volume methods on structured or unstructured meshes.
  - Finite Element formulations for more complex geometries.
  - Hybrid discretizations and discontinuous Galerkin in development.
- Linear and nonlinear solvers with preconditioning strategies tailored for CFD applications.
- Temporal discretization: implicit and semi-implicit time advancement schemes.

### 3.4 Boundary Conditions & Source Terms
- Wide variety of boundary conditions: inflow/outflow, periodic, wall functions, thermal boundaries.
- Source terms representing pumps, resistances, porous media, or physical models relevant for reactor systems.

---

## Notion of BALTIK

**TRUST** is used as the base code of several other CEA codes. You can build your code on TRUST, taking advantage of the data structure and HPC framework.  

## 5. Example Applications

TRUST has been applied in a variety of contexts:

- **Nuclear thermal-hydraulics:** safety studies of pressurized water reactors, cooling system modeling, decay heat removal.
- **General CFD problems:** turbulent flow simulations in pipes, channels, and industrial geometries.
- **Multiphysics coupling:** combining TRUST with external codes for structural mechanics, neutronics, or chemical kinetics.
- **Research & teaching:** a platform for testing new turbulence models, studying discretization effects, or training students in HPC CFD.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
