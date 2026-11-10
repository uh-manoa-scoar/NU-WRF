# **Technical Description of Model Code, User Guidelines, and Guidance for Further Development**

## **1. Model Code**

The implementation of sea spray and wave-affected surface flux physics into the NU-WRF model is available in the branch `seolab_spray`, created from the `develop` branch of the `nu-wrf-dev` project in the `ASTG/NUWRF` group on GitLab. All source code modifications are viewable on GitLab, with the code changes (7 files modified/added) summarized below.

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

The folder **NU-WRF/wave_sea_spray/WRF** in this repo contains copies of the seven modified/added files that were put onto the NASA NU-WRF GitLab repo.  The directory tree in this repo mirrors that of the NASA GitLab repo.

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
- `isftcflx = 5XX03`: uses COARE 3.6 scalar roughness lengths, with spray.  Uses Control spray model coefficients (C3.6_CtlSpr) and wave physics parameterized by wind from Barr et al. (2027, in prep).
- `isftcflx = 5XX04`: same as `5XX03` except that wave physics use input wave properties `SWH_WAVE`, `EPS_WAVE`, `DCP_WAVE`, and `MSS_WAVE`, as in Barr et al. (2023).
- `isftcflx = 5XX05`: uses COARE 3.6 scalar roughness lengths, with spray.  Uses spray model coefficients optimized for the COARE 3.6 scalar roughness lengths (C3.6_OptSpr) and wave physics parameterized by wind from Barr et al. (2027, in prep).
- `isftcflx = 5XX06`: same as `5XX05` except that wave physics use input wave properties `SWH_WAVE`, `EPS_WAVE`, `DCP_WAVE`, and `MSS_WAVE`, as in Barr et al. (2023).
- `isftcflx = 5XX07`: uses COARE-SPRAY model with scalar roughness lengths and spray model coefficients optimized to fit direct covariance data, per Barr et al. (2027, in prep).  Wave physics parameterized by wind from Barr et al. (2027, in prep).
- `isftcflx = 5XX08`: same as `5XX07` except that wave physics use input wave properties `SWH_WAVE`, `EPS_WAVE`, `DCP_WAVE`, and `MSS_WAVE`, as in Barr et al. (2023).

The keys above can be combined to provide a wide range of options for users to explore wave and sea spray physics.  For example, `isftcflx = 50407` accounts for wind-wave misalignment in computing stress, requiring input wave properties `SWH_WAVE`, `WP_WAVE`, and `WDIR_WAVE`, and uses wind-based COARE-SPRAY for heat fluxes with spray, which does not require additional inputs of `EPS_WAVE`, `DCP_WAVE`, and `MSS_WAVE`.  Setting `isftcflx = 50000` uses wind-based COARE 3.5 for stress and COARE 3.6 for scalar roughness (no spray); this is considered the baseline or default for the new options introduced in this work.

For options that require input wave variables, the user must provide these on the mass point grid (i.e., unstaggered in x, y, and z directions) during preprocessing in `wrflowinput` files.  They will also need to activate these fields in the `wrflowinput`, `output`, and `restart` files as necessary by modifying the relevant lines in **WRF/Registry/Registry.EM_COMMON**.

The sea spray flux calculations modify the WRF heat flux variables `HFX`, `LH`, and `QFX` within the surface layer call, and these changes then propagate through model calculations.  In this way, sea spray does not produce new state or forcing fields that interact with WRF outside of the surface layer call.  The sea spray code does, however, provide eleven diagnostic output fields that characterize the module’s behavior and aid in scientific interpretation of its effects.  These can be turned on or off in WRF output via the Registry file at **WRF/Registry/Registry.EM_COMMON**.  Five of these fields are related to spray mass and heat fluxes, described as follows:

- `dHS1_SPR` [W m-2]: Change to surface sensible heat flux (`HFX`) due to spray, accounting for subgrid feedback between bulk and spray heat fluxes
- `dHL1_SPR` [W m-2]: Change to surface latent heat flux (`LH`) due to spray, accounting for subgrid feedback between bulk and spray heat fluxes
- `M_SPR` [kg m-2 s-1]: Generated spray mass flux
- `H_SNSPR` [W m-2]: Net spray sensible heat flux, accounting for direct heating and evaporative cooling
- `H_LSPR` [W m-2]: Spray latent heat flux from evaporation

The remaining six fields are the changes to 2-m air properties caused by subgrid spray-layer feedback and the resulting 2-m values of these variables incorporating these changes.  These are:

- `dTREF_SPR` [K]: Change to 2-m air temperature due to subgrid spray-layer feedback
- `dTHREF_SPR` [K]: Change to 2-m air potential temperature due to subgrid spray-layer feedback
- `dQREF_SPR` [kg kg-1]: Change to 2-m air specific humidity due to subgrid spray-layer feedback.  The spray module assumes that specific humidity and water vapor mixing ratio are equivalent, so this is also the change to 2-m water vapor mixing ratio due to subgrid spray-layer feedback.
- `T2_SPR` [K]: 2-m air temperature including effects of subgrid spray-layer feedback.  This is the sum of `dTREF_SPR` and `SFCLAY1D_mynn`'s internal calculation of the WRF output variable `T2`.
- `TH2_SPR` [K]: 2-m air potential temperature including effects of subgrid spray-layer feedback.  This is the sum of `dTHREF_SPR` and `SFCLAY1D_mynn`'s internal calculation of the WRF output variable `TH2`.
- `Q2_SPR` [kg kg-1]: 2-m air specific humidity including effects of subgrid spray-layer feedback.  This is the sum of `dQREF_SPR` and `SFCLAY1D_mynn`'s internal calculation of the WRF output variable `Q2`.  The spray module assumes that specific humidity and water vapor mixing ratio are equivalent, so this is also the 2-m water vapor mixing ratio including effects of subgrid spray-layer feedback.

All calculations involving spray occur within the MYNN surface layer call (`SFCLAY_mynn`), which is called within the surface driver subroutine (`surface_driver` in **WRF/phys/module_surface_driver.F**).  In `surface_driver`, a surface layer diagnostics module `SFCDIAGS` (**WRF/phys/module_sf_sfcdiags.F**) can be called downstream of `SFCLAY_mynn` (calling `SFCDIAGS` appears to be governed by options set for the land surface scheme).  `SFCDIAGS` overrides the calculations for `T2`, `TH2`, and `Q2` made within `SFCLAY_mynn`.  `SFCDIAGS` operates on the total heat fluxes coming from the surface layer module and computes 2-m diagnostics assuming that the fluxes are due to bulk physics only.  When spray is present, this produces incorrect results due to the incorrect attribution of spray heat fluxes to bulk air-sea differences.  Thus, when spray is present, the variables `T2`, `TH2`, and `Q2` in WRF output should not be used but rather the new variables `T2_SPR`, `TH2_SPR`, and `Q2_SPR` provided by this work.  We note that sea spray does not produce a direct affect on the wind profile in this work, so the diagnosed variables `U10` and `V10` are still valid.

The `mynn_seaice_wrapper` subroutine in **WRF/phys/module_surface_driver.F** has been updated to ingest the new sea spray input and output parameters, but the interaction of sea spray with sea ice is not covered by the new spray model physics.  So, the sea spray physics should not be used in areas where sea ice is present.

Sea spray effects are only applied in this work to the physical variables `HFX`, `QFX`, `LH`, `MOL`, and `qstar`.  The transfer coefficients `FLQC`, `FLHC`, `CHS`, `ch`, `CQS2`, `CHS2`, `Ck`, and `Cka`, which are diagnosed within `SFCLAY1D_mynn`, are not updated.  These may be updated in the future by/for a downstream user as necessary.

## **3. Guidance for Future Development**

The final published configuration of the COARE-SPRAY model was not available at the time of writing of this report.  Changes to update the current NU-WRF code to the final configuration are expected to be small, since the finalized COARE-SPRAY module will retain the overall structure and function of what already exists in NU-WRF.  Should a future user or developer wish to update the NU-WRF code to align with the final COARE-SPRAY configuration, the anticipated required changes are as follows:

- Bring in an updated version of **WRF/phys/module_sf_sprayHFs.F**.  This module is maintained by Ben Barr at https://github.com/bwbarr/sprayHFs under the branch `nu-wrf`.  Please coordinate with Ben Barr to ensure that the file is available and updated.
- Update code comments and keys to the `whichSSGF` input variable in **WRF/phys/module_sf_mynn.F** to reflect the final description and terminology for the COARE-SPRAY model configuration.

Please direct any questions regarding future changes to the spray code to Ben Barr at bbarr@ucar.edu.

## **References**

- Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models. *J. Atmos. Sci.*, **80**, 933 – 960, https://doi.org/10.1175/JAS-D-22-0126.1.
- Barr, B. W., H. Seo, J. B. Edson, C. Sauvage, C. A. Clayson, S. S. Chen, S. F. Zippel, B.-M. Kim, E. J. Thompson, and C. W. Fairall, *In Preparation* (anticipated 2027): Understanding and parameterizing bulk and sea spray heat fluxes in high winds using direct covariance air-sea heat flux observations.
- Edson, J. B., V. Jampana, R. A. Weller, S. P. Bigorre, A. J. Plueddemann, C. W. Fairall, S. D. Miller, L. Mahrt, D. Vickers, and H. Hersbach, 2013: On the exchange of momentum over the open ocean. *J. Phys. Oceanogr.*, **43**, 1589 – 1610, https://doi.org/10.1175/JPO-D-12-0173.1.
- Sauvage, C., H. Seo, C. A. Clayson, and J. B. Edson, 2023: Improving wave-based air-sea momentum flux parameterization in mixed seas. *J. Geophys. Res. Oceans*, **128**, e2022JC019277.
