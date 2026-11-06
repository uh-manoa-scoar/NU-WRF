**Implementation of 1D ocean model in NU-WRF**

NU-WRF currently has the following 1-D OML modeling capability available:  
- A slab (bulk) mixed-layer model by Pollard, Rhines, and Thompson (1972) that captures near-inertial currents and mixed-layer deepening (and the corresponding entrainment) by wind stress.    
- A 3-D Price–Weller–Pinkel (1D-PWP) model (Lee and Chen 2012) is a 3-D extension of 1-D PWP model (Price 1981; Price et al. 1994), which is also a bulk mixed-layer model. It represents diurnal cycling and surface forcing without a prognostic turbulence budget. 

We are integrating a 1-D OML model based on Kantha & Clayson's layer model, which incorporates a prognostic turbulence state, into NU-WRF. Two model versions are under development:      
- Kantha & Clayson (1994, KC94): prognoses turbulent kinetic energy (TKE) and turbulence length scale, offering a computationally efficient representation of upper-ocean mixing      
- Kantha & Clayson (2004, KC04) add wave-driven processes (e.g., wave breaking and Langmuir turbulence) to capture sea-state effects on upper-ocean turbulence.


The KC94 model has been implemented as a new module, module_sf_kc1994oml.F, within the surface layer of NU-WRF. It is called through the module_surface_driver.F and the module_sf_ocean_driver.F. Similar to the already present 1-D ocean mixed layer models in NU-WRF, its use is controlled by the namelist keyword “sf_ocean_physics”, which, when set to “3” will activate the KC94 model. Ocean initial conditions are located in the wrfinput file and can be either uniformly defined by giving ocean temperature, salinity, and depth profile values in the namelist (ocean_t, ocean_z, ocean_z) or by directly replacing the ocean fields (OM_TMP, OM_S, OM_DEPTH) in the wrfinput file for a more accurate definition of the 3-D ocean fields. Additionally, values for ocean currents can be given by the user, modifying directly the wrfinput file (OM_U, OM_V) or simply left to 0 by default. 

Main modifications to the implementation of the KC94 in NU-WRF compared to the original KC94 model (Kantha & Clayson, 1994) have been the addition of bathymetry and land-sea mask information, making it suitable not only for open ocean applications but also for coastal and shelf applications in NU-WRF. This information is added by the user by modifying the OM_DEPTH field in the wrfinput (i.e., setting land values = -999). Over these points, a bottom friction option is turned on using a typical logarithmic function to determine the bottom drag coefficient. After each ocean model integration, the SST field is updated through the TSK variable and sent back to the surface layer scheme.

In the NU-WRF output file, the ocean variables OM_MLD, OM_TMP, and OM_S are updated, respectively, the ocean mixed layer depth, the ocean temperature, and salinity. In addition, several ocean mixing quantities can be output by the model, such as twice the turbulent kinetic energy, the turbulent length scale, and the eddy diffusivity and viscosity (respectively OM_Q2, OM_L, OM_KH, and OM_KM).

Similarly, KC04 is implemented as a new module (i.e., module_sf_kc2004oml.F) and activated by setting “sf_ocean_physics=4” in the namelist. This more recent scheme involves additional mixing terms related to wave breaking and Langmuir turbulence; thus, initial conditions of surface wave fields will be needed. The latter can eventually be taken from the wave variables newly implemented for the wave-based momentum and sea-spray parameterization (i.e., significant wave height, wave period), or if necessary, new wave fields will be made available in either the wrfinput or wrflowinp file. 


**References**

Pollard, R. T., P. B. Rhines, R. O R. Y. Thompson, 2006: *The deepening of the wind-mixed layer*. Geophysical Fluid Dynamics, 4, 381-404      
Price, J. F., 1981: *Upper Ocean Response to a Hurricane*. J. Phys. Oceanogr., 11, 153–175       
Price, J. F., T. B. Sanford, and G. Z. Forristall, 1994: *Forced Stage Response to a Moving Hurricane.* J. Phys. Oceanogr., 24, 233–260      
