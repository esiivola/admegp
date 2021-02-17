# Info
Source code for the paper:

Siivola E, Weber S, Vehtari A. Qualifying drug dosing regimens in pediatrics using
Gaussian processes. Statistics in Medicine. 2021;1–18. https://doi.org/10.1002/sim.8907

Link to the online version of the paper [https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8907](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8907)

NOTE: Stan >= 2.11 is required

## Quick Start
Run oral_1cmt_allo.R in src folder

## Introduction
These files are an attempt to model maturation in clearance with the
help of Gaussian Processes. The code for PKPD is mostly based on Sebastian
Webers' work on building PKPD models in Stan. Pharmacokinetic (PK)
models describe the relationship of drug concentration in a given
subject over time as a function of the dosing history. PK models
facilitate mass action kinetics laws to describe the drug absorption,
distribution and elimination process with a compartmental
approximation.

## NONMEM Data Sets
A widely known data format to describe the longitudinal data of each
patient is the NONMEM format. As such convenient functions are
provided to map a NONMEM data set to the internal format.

# Files
- `src/stan/*.stan` different stan model definitions
- `src/stan/utils/*` helper functions for Stan
- `src/utils.R` functions for generating data, setting defaults and plotting
- `src/oral_1ct_allo.R` the main file
