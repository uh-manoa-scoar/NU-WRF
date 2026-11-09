# **Technical Description of Model Code, User Guidelines, and Guidance for Further Development**

## **1. Model Code**

The implementation of sea spray and wave-affected surface flux physics into the NU-WRF model is available in the branch `seolab_spray`, created from the `develop` branch of the `nu-wrf-dev` project in the `ASTG/NUWRF` group on GitLab. All source code modifications are viewable on GitLab, with the major changes summarized below.

- **WRF/dyn_em/module_first_rk_step_part1.F**: In subroutine `first_rk_step_part1`, passed new arguments from memory (wave input variables and spray output variables) to the call to the subroutine `surface_driver`.
- **WRF/phys/module_surface_driver.F**: 
  - In the subroutine `surface_driver`, passed new arguments to MYNN surface layer calls, i.e., the main subroutine `SFCLAY_mynn` and the sea ice wrapper subroutine `mynn_seaice_wrapper`, which itself calls `SFCLAY_mynn`.
  - In subroutine `mynn_seaice_wrapper`, managed passing of new arguments to calls of `SFCLAY_mynn`.
- **WRF/phys/module_sf_mynn.F**: Most changes to NU-WRF source code are in this module.
  - In subroutine `SFCLAY_mynn` (a call on 2D surface fields), managed passing of new arguments to call of subroutine `SFCLAY1D_mynn` (a call on 1D subsets of 2D fields).
  - In subroutine `SFCLAY1D_mynn`:
    - Added new options to `IF` statement that calculates momentum roughness length `ZNT`.
    - Added new options to `IF` statement that calculates scalar roughness lengths `Z_T` and `Z_Q`.
    - After calculations for stability and friction velocity (`UST`), added a call to the top-level spray heat flux subroutine `sprayHFs`.  This involved adding a large `IF` tree to manage spray options.
    - After calling `sprayHFs`, updated existing surface layer variables `MOL`, `qstar`, `QFX`, `LH`, and `HFX`.
    - Added calculations and new variables for 2-m conditions with spray.
  - Added new subroutines `edson_etal_2013_wave1`, `edson_etal_2013_wave2`, `sauvage_etal_2023_wave1`, and `sauvage_etal_2023_wave2` to calculate new `ZNT` options.
- **WRF/phys/module_sf_sprayHFs.F**: Added this file, which contains the top-level spray subroutine `sprayHFs` and additional subroutines that are used by it. This is external code provided by the project team.
- **WRF/Registry/Registry.EM_COMMON**: Added wave input variables and spray output variables.
- Updated **WRF/main/depend.common** and **WRF/phys/Makefile** to include new file **module_sf_sprayHFs.F**.


## **2. User Guidelines**

The new sea spray and wave physics are available in the MYNN surface layer scheme in NU-WRF.  The MYNN surface layer scheme is selected by setting `sf_sfclay_physics = 5` in the `&physics` section of the WRF input file `namelist.input`.  

The new sea spray and wave physics options are controlled using the `isftcflx` parameter in the `&physics` section of `namelist.input`.  New options are structured using integers in the format `5XXXX`.  The leading `5` indicates that the option belongs to the family of new wave and spray physics options provided by this work.  The second and third digit from the left specify the option to use for momentum roughness length, and these options select between the wave-based roughness options implemented in this work.  The available options are:

- `isftcflx = 500XX`: uses the COARE 3.5 wind-based model (Edson et al. 2013).  This is the baseline option for the new family of parameterizations.
- `isftcflx = 501XX`: uses the COARE 3.5 wave-based model based on peak wave age (Edson et al. 2013).  This option requires the wave input variable `WP_WAVE` (wave period), which should contain the *peak* wave period.
- `isftcflx = 502XX`: uses the COARE 3.5 wave-based model based on peak wave age and significant wave height (Edson et al. 2013).  This option requires the wave input variables `SWH_WAVE` (significant wave height) and `WP_WAVE`, where `WP_WAVE` should contain the *peak* wave period.
- `isftcflx = 503XX`: uses the modified COARE 3.5 wave-based model based on mean wave age and significant wave height (Sauvage et al. 2023).  This option requires the wave input variables `SWH_WAVE` and `WP_WAVE`, where `WP_WAVE` should contain the *mean* wave period.
- `isftcflx = 504XX`: uses the modified COARE 3.5 wave-based model based on wind-wave misalignment, peak wave age, and significant wave height (Sauvage et al. 2023).  This option requires the wave input variables `SWH_WAVE`, `WP_WAVE`, and `WDIR_WAVE` (peak wave direction), where `WP_WAVE` should contain the *peak* wave period.

The first and second digit from the right specify the option to use for heat fluxes, and these options select between the sea spray heat flux options implemented in this work.  The available options are:

- `isftcflx = 5XX00`: uses COARE 3.6 scalar roughness lengths, with no sea spray physics (i.e., bulk physics only).  This is the baseline option for the new family of parameterizations.
- `isftcflx = 5XX01`: uses COARE 3.6 scalar roughness lengths, with spray.  Uses original Barr et al. (2023) spray model coefficients and wave physics parameterized by wind per Barr et al. (2027, in prep).
- `isftcflx = 5XX02`: same as `5XX01` except that wave physics use input wave properties `SWH_WAVE`, `EPS_WAVE` (wave energy dissipation flux), `DCP_WAVE` (dominant/peak phase speed), and `MSS_WAVE` (mean square wave slope), as in Barr et al. (2023).

The keys above can be combined to provide a wide range of options for users to explore wave and sea spray physics.  For example, [give a couple of examples].

For options that depend on wave variables, the user must provide these on the mass point grid (i.e., unstaggered in x, y, and z directions) during preprocessing in wrflowinput files.  They will also need to activate these fields in the wrflowinput, output, and restart files as necessary by modifying the relevant lines in [WRF Registry file].
The sea spray calculations modify the WRF heat flux variables HFX, LH, and QFX within the surface layer call, and these changes then propagate through model calculations.  In this way, sea spray does not produce new state or forcing fields that interact with WRF outside of the surface layer call.  The sea spray code does, however, provide eleven [?] diagnostic output fields that characterize the module’s behavior and aid in scientific interpretation of its effects.  These can be turned on or off in WRF output via the Registry file at [path].  Five [?] of these fields are spray mass and heat fluxes, described as follows:

[List spray flux output fields]

The remaining six [?] fields are the changes to 2-m temperature, potential temperature, and specific humidity caused by spray and the actual 2-m values of these variables incorporating these changes.  These are:

[List spray surface layer variable output fields]

All calculations involving spray occur within the MYNN surface layer call.  In the WRF code, a surface layer diagnostics module [name] is called downstream of the MYNN surface layer call, and this diagnostics module overrides the calculations made within the MYNN surface layer call and writes the results to [variable names].  The diagnostics module operates on the total heat fluxes coming from the surface layer module and computes 2-m diagnostics assuming that the fluxes are due to bulk physics only.  When spray is present, this diagnostics module produces incorrect results due to the incorrect attribution of spray heat fluxes to bulk air-sea differences.  Thus, when using any of the isftxflx = 5XXXX options, the variables [variable names] should not be used but rather the new variables [variable names] provided by this work.  We note that when sea spray is not active, the variables [variable names] still produce the correct 2-m results for the bulk-only system.  Finally, we note that sea spray does not produce a direct affect on the wind profile in this work, so the diagnostic variables U10 and V10 are still valid.

Although the `mynn_seaice_wrapper` subroutine in **WRF/phys/module_surface_driver.F** has been updated to ingest the new sea spray input parameters, the interaction of sea spray with sea ice is not covered by the new spray model physics.  So, the sea spray physics should not be used in areas where sea ice is present.

**3. Description of Wave and Spray Implementation and Guidance for Future Updates**

A full record of code changes made to implement wave and sea spray physics into NU-WRF can be found in the Git repo for the NU-WRF code.  The spray physics calculations are performed by a top-level subroutine called sprayHFs; this subroutine and supporting subroutines are included in a new file module_sf_sprayHFs.F in the directory [path].  Most modifications to existing NU-WRF code to implement the wave and spray physics were made to the MYNN surface layer module file, module_sf_mynn.F.  Outside of this file, changes were also made to the surface driver, Registry, and other files as needed to bring wave variables to the surface layer module and to handle input and output fields in memory.  Within module_sf_mynn.F, the major code modifications (excluding variable definitions, input/output field passing, and other code management-related items) are as follows:

The IF loop defining the momentum roughness length (variable ZNT) was expanded to include the isftcflx = 5XXXX options.

The IF loop defining the thermal and moisture roughness lengths (variables ZT and ZQ) was expanded to include the isftcflx = 5XXXX options.

A large IF loop was added to define inputs and switches needed by sprayHFs, and to call sprayHFs to perform the spray calculations.
The turbulent flux scales MOL and QSTAR were updated for the effects of spray.
The heat fluxes HFX, QFX, and LH were updated for the effects of spray.
2-m reference fields for temperature, potential temperature, and specific humidity were updated for the effects of spray.
Subroutines were added or updated to include the new wave-dependent momentum roughness length formulations.

The final published configuration of the COARE-SPRAY model was not available at the time of writing of this report.  Should a future user or developer wish to update the NU-WRF code to align with this final configuration, the anticipated required changes are as follows:

[add these here after getting access to the code on Gitlab again]


**References**

- Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: *Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models.* J. Atmos. Sci., 80, 933 – 960, https://doi.org/10.1175/JAS-D-22-0126.1.
- Edson, J. B., V. Jampana, R. A. Weller, S. P. Bigorre, A. J. Plueddemann, C. W. Fairall, S. D. Miller, L. Mahrt, D. Vickers, and H. Hersbach, 2013: *On the exchange of momentum over the open ocean.* J. Phys. Oceanogr., 43, 1589 – 1610, https://doi.org/10.1175/JPO-D-12-0173.1.
- Sauvage, C., H. Seo, C. A. Clayson, and J. B. Edson, 2023: *Improving wave-based air-sea momentum flux parameterization in mixed seas.* J. Geophys. Res. Oceans, 128, e2022JC019277.
