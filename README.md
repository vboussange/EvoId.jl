# ABMEv.jl
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vboussange.github.io/ABMEv.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://vboussange.github.io/ABMEv.jl/dev)

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
