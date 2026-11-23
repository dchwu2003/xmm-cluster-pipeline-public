"""Initial-model prototype: predict the starting ICM temperature of a region from band fluxes measured
before any fit, trained on the fit log.

* ``data``: load the bundled per-spectrum table (``data/initmodel/spectra.csv``);
* ``features``: the feature matrix (hardness ratios, sky fractions, surface brightness, z, N_H, ...);
* ``cv``: leave-one-cluster-out evaluation of extra trees against the catalogue-kT baseline;
* ``replay``: summary of refits started from different initial values;
* ``band_fluxes``: how the band fluxes are measured from the spectra (reference; needs the spectra).
"""
