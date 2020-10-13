# ABMEv.jl
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vboussange.github.io/ABMEv.jl/stable) -->
<!-- For now we only direct to dev documentation. In the future, one will need to deploy a ssh key to and use TagBot. -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://vboussange.github.io/ABMEv.jl/dev)
[![pipeline status](https://gitlab.ethz.ch/bvictor/abmev/badges/master/pipeline.svg)](https://gitlab.ethz.ch/bvictor/abmev/-/commits/master)
[![coverage report](https://gitlab.ethz.ch/bvictor/abmev/badges/master/coverage.svg)](https://gitlab.ethz.ch/bvictor/abmev/-/commits/master)

<div align="center"> <img
src="docs/src/assets/logo.png"
alt="ABMEv.jl logo" width="400"></img> </div>

This is a suite for simulating an Agent Based Model that captures the evolutionary dynamics of a population in a multidimensional space.

## Installation
```julia
using Pkg;
Pkg.add("https://gitlab.ethz.ch/bvictor/abmev.git")
```
This will download latest version from git repo and download all dependencies.
To check out from an other branch than master, one has to do the trick
```julia
using Pkg;
Pkg.add("ABMEv#no_C_matrix")
```
## Getting started
```julia
using ABMEv
```
