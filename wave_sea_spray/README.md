**3. Technical Description of Model Codes, Implementation, and User Guidelines**

**3.1 Model Code**
The implementation of sea spray and wave physics into the NU-WRF model is available at [website] in the branch seolab_spray.

**3.2 User Guidelines for New Sea Spray and Wave Physics**

The new sea spray and wave physics are implemented in the MYNN surface layer scheme in NU-WRF.  The MYNN surface layer scheme is selected by setting sf_sfclay_physics = 5 in the &physics section of the WRF input file namelist.input.  

The new sea spray and wave physics options are controlled using the isftcflx parameter in the &physics section of namelist.input.  New options are structured using integers in the format 5XXXX.  The leading 5 indicates that the option belongs to the family of new wave and spray physics options provided by this work.  The second and third digit from the left specify the option to use for momentum roughness length, and these options select between the wave-based roughness options implemented in this work.  The available options are:

**isftcflx = 500XX**: uses the COARE 3.5 wind-based model.  This is the baseline option for the new family of parameterizations.

**isftcflx = 501XX:** uses the model based on peak wave age.  This option requires the wave input variables [].

**isftcflx = 502XX:** uses the model based on peak wave age and significant wave height.  This option requires the wave input variables [].

**isftcflx = 503XX:** uses the model based on mean wave age.  This option requires the wave input variables [].

**isftcflx = 504XX:** uses the model based on wind-wave misalignment angle.  This option requires the wave input variables [].

The first and second digit from the right specify the option to use for heat fluxes, and these options select between the sea spray heat flux options implemented in this work.  The available options are:

isftcflx = 5XX00: uses COARE 3.6 scalar roughness lengths, with no sea spray physics (i.e., bulk physics only).  This is the baseline option for the new family of parameterizations.
[add the other spray options here… don’t remember the order]

The keys above can be combined to provide a wide range of options for users to explore wave and sea spray physics.  For example, [give a couple of examples].

For options that depend on wave variables, the user must provide these on the mass point grid (i.e., unstaggered in x, y, and z directions) during preprocessing in wrflowinput files.  They will also need to activate these fields in the wrflowinput, output, and restart files as necessary by modifying the relevant lines in [WRF Registry file].
The sea spray calculations modify the WRF heat flux variables HFX, LH, and QFX within the surface layer call, and these changes then propagate through model calculations.  In this way, sea spray does not produce new state or forcing fields that interact with WRF outside of the surface layer call.  The sea spray code does, however, provide eleven [?] diagnostic output fields that characterize the module’s behavior and aid in scientific interpretation of its effects.  These can be turned on or off in WRF output via the Registry file at [path].  Five [?] of these fields are spray mass and heat fluxes, described as follows:

[List spray flux output fields]

The remaining six [?] fields are the changes to 2-m temperature, potential temperature, and specific humidity caused by spray and the actual 2-m values of these variables incorporating these changes.  These are:

[List spray surface layer variable output fields]

All calculations involving spray occur within the MYNN surface layer call.  In the WRF code, a surface layer diagnostics module [name] is called downstream of the MYNN surface layer call, and this diagnostics module overrides the calculations made within the MYNN surface layer call and writes the results to [variable names].  The diagnostics module operates on the total heat fluxes coming from the surface layer module and computes 2-m diagnostics assuming that the fluxes are due to bulk physics only.  When spray is present, this diagnostics module produces incorrect results due to the incorrect attribution of spray heat fluxes to bulk air-sea differences.  Thus, when using any of the isftxflx = 5XXXX options, the variables [variable names] should not be used but rather the new variables [variable names] provided by this work.  We note that when sea spray is not active, the variables [variable names] still produce the correct 2-m results for the bulk-only system.  Finally, we note that sea spray does not produce a direct affect on the wind profile in this work, so the diagnostic variables U10 and V10 are still valid.

**3.3 Description of Wave and Spray Implementation and Guidance for Future Updates**
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

Barr, B. W., S. S. Chen, and C. W. Fairall, 2023: *Sea-state-dependent sea spray and air-sea heat fluxes in tropical cyclones: A new parameterization for fully coupled atmosphere-wave-ocean models.* J. Atmos. Sci., 80, 933 – 960, https://doi.org/10.1175/JAS-D-22-0126.1.
Edson, J. B., V. Jampana, R. A. Weller, S. P. Bigorre, A. J. Plueddemann, C. W. Fairall, S. D. Miller, L. Mahrt, D. Vickers, and H. Hersbach, 2013: *On the exchange of momentum over the open ocean.* J. Phys. Oceanogr., 43, 1589 – 1610, https://doi.org/10.1175/JPO-D-12-0173.1.
Sauvage, C., H. Seo, C. A. Clayson, and J. B. Edson, 2023: *Improving wave-based air-sea momentum flux parameterization in mixed seas.* J. Geophys. Res. Oceans, 128, e2022JC019277.
