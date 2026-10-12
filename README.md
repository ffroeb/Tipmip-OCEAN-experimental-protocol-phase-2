This repository contains details of the experimental protocol of the second phase of the TIPMIP-OCEAN project. See http://tipmip.org/ocean-domain/ and Swingedouw et al., 2026

The experimental protocol is described in detail in Fröb et al, (submitted), TIPMIP-OCEAN experimental protocol phase 2: Southern Ocean and Combined Dual-Hemispheric Freshwater Forcing. GMD [Full reference to come]

The attached files are:

create_SO_hose_mask.py - python code to create netcdf files and plot the spatial mask for applying idealised freshwater fluxes in the Southern Ocean (i.e., a total of 0.3Sv). Freshwater fluxes decrease exponentially with increasing distance from the coast off Antarctica, and are set to zero beyond 300 km. Follows definition of Swingedouw et al., 2026 and Jackson et al., 2023

grid.nc - grid file for BLOM, NorESM2, CMIP6


References
Swingedouw, D., et al.: TIPMIP-OCEAN experimental protocol phase 1: Tipping dynamics of the AMOC, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2026-1698, 2026.
Jackson, L. C., et al.: Understanding AMOC stability: the North Atlantic Hosing Model Intercomparison Project, Geosci. Model Dev., 16, 1975–1995, https://doi.org/10.5194/gmd-16-1975-2023, 2023.

