This is a FABM-based adaptation of the ECOSMO marine biogeochemical model, originally developed for the Baltic/North Sea (Neumann, 2000; Daewel and Schrum 2013),
and later adapted by NERSC for applications to the Nordic Seas and Arctic Ocean (Yumruktepe et al., 2022), and used within the TOPAZ-ECOSMO operational CMEMS product.

The development of the niva-ecosmo code was originally funded by the Horizon Europe project EU-INTERCHANGE as part of an effort to develop efficient and well-performing
biogeochemical components of high-resolution digital twins for the Atlantic/Arctic (within NorHAPS model) and the Norwegian coastal region (Norkyst model).

The current niva-ecosmo version was based on the version of the NERSC operational code within the FABM during August 2025.
It has since diverged from the NERSC code in terms of formulation details and coding style. 
All these changes are documented and explained within dated comments at the top of the files concerned.

The code herein is intended to provide a "lean", computationally-efficient biogeochemical module of medium complexity, that can be suitable
for incorporation into high-resolution, 3D, coupled physical-biogeochemical models.
A key metric of efficient is the Run Time Ratio (RTR) between coupled physical-biogeochemical and physics-only simulations.
We aim for a target RTR value less than 3 such that long (e.g. decadal) simulations over large domains and at high-resolution
do not become prohibitive in terms of time demands and computational cost.

The source code is provided in subfolder 'src', while 'example_configs' contains example fabm.yaml files and figures showing comparisons of
model output time series with observations and skill assessments based on 4D point-to-point matchups to in situ observational data.




References:

Daewel and Schrum (2013) (D13), doi:10.1016/j.jmarsys.2013.03.008

Neumann (2000), Journal of Marine Systems 25 (2000) 405–419

Yumruktepe et al. (2022) (Y22), doi:10.5194/gmd-15-3901-2022
