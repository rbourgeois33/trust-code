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


## 0. Size

- Paragraph 1 and 2 = 1 page
- Another page for the rest, small parts: 3 and 5, more room for 4 
  
## 1. Introduction & Historical Background

**TRUST** is an open-source numerical platform developed at the **French Alternative Energies and Atomic Energy Commission (CEA)**.  
It originates from decades of development of in-house codes at CEA for **nuclear reactor thermal-hydraulics** and **computational fluid dynamics (CFD)** applications. In order to consolidate expertise and provide a modern, extensible framework, CEA launched TRUST as a unified, open, and modular platform.

The project was opened to the community to foster collaboration between research institutes, academia, and industry. Today, TRUST is used both as a **production-grade simulation tool** and as a **research platform** to test new models, discretization schemes, and parallel algorithms.

---

## 2. Physical and numerical capabilities

- Conduction 
- Navier-Stokes incompressible -> DNS / QC / WC
- Multiphase && PolyMAC
- Coupling

## 3. Code structure 

- C++ OO
- Modularity
- CI
- Solvers : Petsc, ...
- IO : cgns
  
## 4. HPC capabilities

- MPI
- Kokkos / GPU
- PDI
- Scalling and big runs

## 5. Derived applications

- Idea 
- Pub for TrioCFD

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
