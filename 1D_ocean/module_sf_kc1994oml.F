!WRF:MODEL_LAYER:OCEAN PHYSICS
!
MODULE module_sf_kc1994oml
        
CONTAINS

!---------------------------------------------------------------------------          
    SUBROUTINE KC1994OML(ims,ime, jms,jme, kms,kme, &
                         its,ite, jts,jte, kts,kte, &
                         ids,ide, jds,jde, kds,kde, &
                         om_tmp,om_s,om_u,om_v, &
                         om_tmp_b,om_s_b,om_u_b,om_v_b, &
                         om_depth,om_ml_temp,om_ml_sal,om_ml_turb,om_ml_dens,om_lat, &
                         HFX,QFX,GSW,GLW,UST,U_PHY,V_PHY,TSK,LH,XLAND,RAINCV, &
                         RAINNCV,DT, &
                         om_q2,om_q2l,om_q2_b,om_q2l_b, &
                         om_l,om_km,om_kh,om_kq,om_kmb,om_khb, &
                         xtime,id,omdt,KB_tot)

!----------------------------------------------------------------------------
!
!    1D OCEAN MIXED LAYER MODEL BASED ON SECOND MOMENT CLOSURE.
!     
!    BASED ON: KANTHA, L.H. AND C.A. CLAYSON (1994)
!              AN IMPROVED MIXED LAYER MODEL FOR GEOPHYSICAL
!              APPLICATIONS, J. GEOPHYS. RES.,99,25235.
!
!    COPYRIGHT BY KANTHA AND CLAYSON
!    kantha@boulder.colorado.edu, cclayson@fsu.edu
!
!    IMPLEMENTED IN WRF OCEAN PHYSICS OCTOBER 2023 BY
!    SEOLAB: CESAR SAUVAGE, BEN BARR, AND HYODAE SEO
!    csauvage@whoi.edu, benjamin.barr@whoi.edu, hseo@whoi.edu
!
!         om_tmp        3D ocean temperature (k)
!         om_s	        3D ocean salinity (psu) 
!         om_u          3D ocean current in x-direction (m s-1)
!         om_v          3D ocean current in y-direction (m s-1)
!         om_tmp_b      3D ocean temperature (k), previous timestep
!         om_s_b	3D ocean salinity (psu), previous timestep
!         om_u_b        3D ocean current in x-direction (m s-1), previous timestep
!         om_v_b        3D ocean current in y-direction (m s-1), previous timestep
!         om_depth      3D ocean depth (m)
!         om_ml_temp    Temperature-based ocean mixed layer depth (m)
!         om_ml_sal     Salinity-based ocean mixed layer depth (m)
!         om_ml_turb    Turbulence-based ocean mixed layer depth (m)
!         om_ml_dens    Density-based ocean mixed layer depth (m)
!         om_lat        Latitude of ocean grid (deg), same as WRF model
!         HFX           WRF surface sensible heat flux (W m-2), (+) OCN --> ATM
!         QFX           WRF surface moisture flux (kg m-2 s-1), (+) OCN --> ATM
!         GSW           Net shortwave flux at surface (W m-2), (+) ATM --> OCN
!         GLW           Downward longwave flux at surface (W m-2), (+) ATM --> OCN
!         UST           u* from WRF model (air-side) (m s-1)
!         U_PHY         x-component of wind on WRF mass grid (m s-1)
!         V_PHY         y-component of wind on WRF mass grid (m s-1)
!         TSK           WRF surface temperature (K)
!         LH            WRF surface latent heat flux (W m-2), (+) OCN --> ATM
!         XLAND         WRF landmask (1 for land, 2 for water)
!         RAINCV        TIME-STEP CUMULUS PRECIPITATION (mm)
!         RAINNCV       TIME-STEP NONCONVECTIVE PRECIPITATION (mm)
!         DT            WRF timestep
!         om_q2         TKE times 2
!         om_q2l        q2 times turbulent length scale
!         om_q2_b       TKE times 2, previous timestep
!         om_q2l_b      q2 times turbulent length scale, previous timestep
!         om_l          Turbulent length scale (m)
!         om_km         Vertical kinematic viscosity (m2 s-1)
!         om_kh         Vertical diffusivity for scalars (temp, salinity, etc)
!         om_kq         Vertical diffusivity for turbulence quantities
!         om_kmb        Background values of om_km
!         om_khb        Background values of om_kh
!         xtime         Counter of model timesteps (I think). Used to determine if at IC.
!         id            Grid id number
!         omdt          Ocean timestep (min)
!         KB_tot        Number of total vertical gridpoints
!         KB            Number of vertical gridpoints depending on bathymetr (actual ocean points)
!
!    This subroutine uses nested subroutines that allow its declared variables
!    to be global with respect to nested subroutine calls.  This was done to 
!    retain the structure of the original KC 1994 code, which used global 
!    declarations via a Fortran 77 common block.  Since nested subroutines
!    do not require argument passing, returned (modified) fields for each
!    nested subroutine call are stated throughout in comments for clarity.
!
!-------------------------------------------------------------------------------   

        IMPLICIT NONE

        INTEGER,INTENT(IN) :: id,KB_tot
        INTEGER,INTENT(IN) :: ids,ide, jds,jde, kds,kde,  &
                              ims,ime, jms,jme, kms,kme,  &
                              its,ite, jts,jte, kts,kte
        REAL,INTENT(IN) :: xtime,omdt,DT
        REAL,DIMENSION(ims:ime,jms:jme),INTENT(IN) :: HFX,QFX,LH,GSW,GLW,UST,XLAND,RAINCV,RAINNCV
        REAL,DIMENSION(ims:ime,jms:jme),INTENT(INOUT) :: TSK
        REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: U_PHY,V_PHY
        real,dimension(ims:ime,KB_tot,jms:jme),INTENT(INOUT) :: OM_TMP,OM_S,OM_U,OM_V,OM_Q2,OM_Q2L
        real,dimension(ims:ime,KB_tot,jms:jme),INTENT(INOUT) :: OM_TMP_B,OM_S_B,OM_U_B,OM_V_B,OM_Q2_B,OM_Q2L_B
        real,dimension(ims:ime,KB_tot,jms:jme),INTENT(INOUT) :: OM_L,OM_KM,OM_KH,OM_KQ,OM_KMB,OM_KHB
        real,dimension(ims:ime,KB_tot,jms:jme),INTENT(IN) :: OM_DEPTH
        real,dimension(ims:ime,jms:jme),INTENT(IN) :: OM_LAT
        real,dimension(ims:ime,jms:jme),INTENT(OUT) :: OM_ML_TEMP,OM_ML_SAL,OM_ML_TURB,OM_ML_DENS

        INTEGER :: i,j,k,KBM1,KBM2,KI,LOOPDO,KB
        real,dimension(KB_tot) :: TF,T,TB,SF,S,SB,UF,U,UB,VF,V,VB,Q2F,Q2,Q2B,Q2LF,Q2L,Q2LB
        real,dimension(KB_tot) :: L,KM,KH,KQ,KMB,KHB
        real,dimension(KB_tot) :: Z,ZZ,DZ,DZZ
        real,dimension(KB_tot) :: A,C,VH,VHP,DM,SWRAD
        real,dimension(KB_tot) :: RHO,GM,GH,SM,SH,KN,SPROD,BPROD,PROD,DTEF
        REAL,DIMENSION(KB_tot) :: QLI,POD,RII,RIF,GDRDZ
        REAL :: WUSURF,WVSURF,WUBOT,WVBOT,WTSURF,WSSURF,WTSOLR,TERM,WFTOTL,WFSURF,WFSOLR
        REAL :: H,D,N,DH,PBAR,CC,ALAT,COR,EPS
        REAL :: DT1,DT2,DP1,DP2,DP3,DP4,BOTTOM,DELZ
        REAL :: TAUX,TAUY,EV0,PR0,QSOL,QBR,UG,VG,ZZD
        REAL :: USTAR,ZIN,UMOLPR,DUDZ,DVDZ,CBC,theta

        REAL,PARAMETER :: PR=1., KAPPA=0.40, PI=3.14159
        REAL,PARAMETER :: CEEDEE=1.25E-03, RHORAT=1.201E-03
        REAL,PARAMETER :: A1=0.92,B1=16.6,A2=0.74,B2=10.1,C1=0.08
        !!REAL,PARAMETER :: E1=1.8,E2=1.33,E3=1.0,E4=1.0,SQ=0.2,SQ1=0.41
        REAL,PARAMETER :: E1=1.8,E2=8.37,E3=1.0,E4=6.9,SQ=0.41,SQ1=1.0
        REAL,PARAMETER :: SEF=1.0,CONST1=6.51*SEF
        REAL,PARAMETER :: C2=0.7,C3=0.2,E5=0.0,E6=0.0,EM=0.0,FILT=0.0
        REAL,PARAMETER :: GEE=9.806,SMALL=1.E-10
        REAL,PARAMETER :: FCT=0.25
        REAL,PARAMETER :: DELT=0.2, DELS=0.1, DELR=1.E-5
        REAL,PARAMETER :: GRAV=9.806, SMOTH=0.10, UMOL=8.8E-7

        ! LEVEL 2 1/4 parameters
        REAL, PARAMETER :: CF1=A2*(1.-6.*A1/B1), &
                           CF2=3.*A2*(6.*A1+B2*(1.-C3)), &
                           CF3=A1*(1.-3.*C1-6.*A1/B1), &
                           CF4=9.*A1*(2.*A1+A2*(1.-C2)), &
                           CF5=9.*A1*A2, &
                           GHMIN=-0.25, &
                           GHMAX=0.98/(A2*(B1+12.*A1+3.*B2*(1.-C3)))
 
        ! Solar Absorption:
        ! Currently planning to use 2-band (Paulson and Simpson, which has most 
        ! available turbidity types) with a user setting to select type.  
        ! For now, hard-coding TYPE IB that is used in KC2004 code.
        INTEGER,PARAMETER :: SURABS=0,NBAND=2
        ! -------- Select water type ---------
        ! *** CLEAR ***
        !REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.60,  0.40/)
        !REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/ 40.00,  5.00/)
        ! *** TYPE I ***
        !REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.42,  0.58/)
        !REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/ 23.00,  0.35/)
        ! *** TYPE IA ***
        !REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.38,  0.62/)
        !REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/ 20.00,  0.60/)
        ! *** TYPE IB ***
        REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.33,  0.67/)
        REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/ 17.00,  1.00/)
        ! *** TYPE II ***
        !REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.23,  0.77/)
        !REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/ 14.00,  1.50/)
        ! *** TYPE III ***
        !REAL,PARAMETER,DIMENSION(NBAND) :: ALPH = (/  0.22,  0.78/)
        !REAL,PARAMETER,DIMENSION(NBAND) :: SCL  = (/  7.90,  1.40/)
        ! ------------------------------------
        REAL,DIMENSION(NBAND) :: SCLI
        DO N=1,NBAND
            SCLI(N)=1./SCL(N)
        ENDDO


        print*,'Running KC 1994 OML model: domain ',id


        ! Set up a few variables
        EPS=1.E-10
        DT1=omdt*60.    ! Convert timestep to [s]
        DT2=2.*DT1

        ! Loop over grid points and perform 1D calculations
        DO i = its, ite
            DO j = jts, jte

                ! ==========================================================
                ! 1. Set up vertical grid and coriolis parameter.
                ! ==========================================================

                ! Define KB depending on bathymetry
                ! check ocean depth field in wrfinput file set such as land is -999.9
                KB=KB_tot !1st guess
                DO K = 1, KB_tot
                        if (om_depth(i,K,j) .LT. -900.) THEN
                                KB = K-1
                                EXIT ! exit loop whenever KB value is found
                        ENDIF
                ENDDO       
                KBM1=KB-1
                KBM2=KB-2
                ! Vertical grid items. Must have KB = (BOTTOM/DELZ)+1
                BOTTOM = om_depth(i,KB,j)
                DELZ = om_depth(i,2,j) - om_depth(i,1,j)    ! Grid spacing, uniform required
                H=BOTTOM

                ! Set vertical grid and spacing using DEPTH().  
                ! Returns: Z, ZZ, DZ, and DZZ
                CALL DEPTH()

                ! Calculate coriolis parameter
                ALAT = om_lat(i,j)
                COR=1.46E-4*SIN(PI*ALAT/180.)

                ! ==========================================================
                ! 2. Initialize or import T, S, U, V, Q2, Q2L, L, KM, KH, KQ, KMB, and KHB.
                ! ==========================================================

                if (xtime .EQ. 0.) THEN    ! Initialize if xtime == 0
                    DO K=1,KB
                        ! Import from IC
                        TB(K) = om_tmp(i,K,j) - 273.15    ! [C]
                        SB(K) = om_s(i,K,j)    ! [psu]
                        UB(K) = om_u(i,K,j)    ! [m s-1]
                        VB(K) = om_v(i,K,j)    ! [m s-1]

                        ! Assign same value to current fields
                        T(K)=TB(K)
                        S(K)=SB(K)
                        U(K)=UB(K)
                        V(K)=VB(K)
 
                        ! For turbulence and diffusivities, initialize with nominal small values
                        Q2B(K)=EPS
                        Q2(K)=Q2B(K)
                        Q2LB(K)=SQRT(Q2(K))*UMOL*0.1
                        Q2L(K)=Q2LB(K)
                        L(K)=Q2LB(K)/Q2B(K)
                        KM(K)=0.1*UMOL
                        KH(K)=0.1*UMOL
                        KQ(K)=0.1*UMOL
                        KMB(K)=0.1*UMOL
                        KHB(K)=0.1*UMOL
                    ENDDO
                    KQ(1)=0.
                    KQ(KB)=0.
                
                else    ! Take fields from memory
                    DO K=1,KB
                        TB(K) = om_tmp_b(i,K,j) - 273.15
                        T(K)  = om_tmp(i,K,j) - 273.15
                        SB(K) = om_s_b(i,K,j)
                        S(K)  = om_s(i,K,j)
                        UB(K) = om_u_b(i,K,j)
                        U(K)  = om_u(i,K,j)
                        VB(K) = om_v_b(i,K,j)
                        V(K)  = om_v(i,K,j)
                        Q2B(K)  = om_q2_b(i,K,j)
                        Q2(K)   = om_q2(i,K,j)
                        Q2LB(K) = om_q2l_b(i,K,j)
                        Q2L(K)  = om_q2l(i,K,j)
                        L(K)   = om_l(i,K,j)
                        KM(K)  = om_km(i,K,j)
                        KH(K)  = om_kh(i,K,j)
                        KQ(K)  = om_kq(i,K,j)
                        KMB(K) = om_kmb(i,K,j)
                        KHB(K) = om_khb(i,K,j)
                    ENDDO
                endif

                ! ==========================================================
                ! 3. Determine surface fluxes of momentum, heat, solar radiation, and salt.
                ! ==========================================================

                !***Change to pass only lowest level of U_PHY and V_PHY***
                ! 10m wind speed
                !wsp10 = SQRT(U_PHY(i,kts,j)**2+V_PHY(i,kts,j)**2)
                !wsp10m = wsp10*log(10/1e-4)/log(zu/1e-4)
                theta = atan2(V_PHY(i,kts,j),U_PHY(i,kts,j))    ! Assume surf stress and wind vectors aligned
                TAUX = (1.225*UST(i,j)**2)*cos(theta)    ! X-momentum flux [N m-2], (+) EASTWARD ATM --> OCN
                TAUY = (1.225*UST(i,j)**2)*sin(theta)    ! Y-momentum flux [N m-2], (+) NORTHWARD ATM --> OCN
                ! QBR = Net turbulent + longwave surface flux [W m-2], (+) OCN --> ATM
                QBR = (0.985 * 5.67e-8 * om_tmp(i,1,j)**4) - GLW(i,j) + LH(i,j) + HFX(i,j)
                WTSURF = QBR * 0.246e-6    ! QBR/(rhosw*cpsw) [C m s-1], (+) OCN --> ATM
                QSOL = GSW(i,j)    ! Shortwave surface flux [W m-2], (+) ATM --> OCN
                WTSOLR = -QSOL * 0.246e-6    ! QSOL/(rhosw*cpsw) [psu m s-1], (+) OCN --> ATM
                ! EVO = Surface evaporation rate = LH/(rhosw*Lv) [m s-1], (+) OCN --> ATM
                EV0=LH(i,j)/((2.5008e6-2.3e3*(om_tmp(i,1,j)-273.15))*1.025e3)
                ! PR0=0.0    ! PR0 = Surface precip flux [m s-1], (+) ATM --> OCN
                PR0 = (RAINNCV(i,j) + RAINCV(i,j)) / (1000*DT) ! convert from mm to m/s
                WSSURF = -S(1)*(EV0-PR0)    ! Surface salinity flux [psu m s-1], (+) OCN --> ATM

                ! ==========================================================
                ! 4. Advance temperature and salinity equations, filter, and
                ! calculate density.
                ! ==========================================================

                ! Advance temperature equation using PROFT().
                ! Returns: TF
                CALL PROFT()

                ! Advance salinity equation using PROFS().
                ! Returns: SF
                CALL PROFS()

                ! Apply filter to remove time-splitting effects (Asselin 1972).
                DO K=1,KB
                    T(K)=T(K)+.5*SMOTH*(TF(K)+TB(K)-2.*T(K))
                    TB(K)=T(K)
                    T(K)=TF(K)
                    S(K)=S(K)+.5*SMOTH*(SF(K)+SB(K)-2.*S(K))
                    SB(K)=S(K)
                    S(K)=SF(K)
                ENDDO

                ! Calculate density based on T and S using DENS().
                ! Returns: RHO
                CALL DENS()

                ! ==========================================================
                ! 5. Determine surface momentum flux forcing and calculate currents.
                ! ==========================================================

                WUSURF = -TAUX*1.0e-3    ! TAUX/rhosw [m2 s-2], (+) EASTWARD OCN --> ATM
                WVSURF = -TAUY*1.0e-3    ! TAUY/rhosw [m2 s-2], (+) NORTHWARD OCN --> ATM

                ! For generalized use and after consulting with CAC, we are
                ! assuming zero geostrophic currents and no barotropic tides.
                ! This means that there is no uniquely zonal or meridional
                ! behavior that needs to be included in the U and V prognostic
                ! equations, and all velocities can remain in WRF grid-relative
                ! coordinates (i.e., no rotation is necessary).
                ! ***CONFIRM THIS WITH HYODAE AND CESAR***
                UG=0.
                VG=0.

                DO K=1,KBM1
                    UF(K)=UB(K)+DT2*COR*(V(K)-VG)
                    VF(K)=VB(K)-DT2*COR*(U(K)-UG)
                ENDDO

                ! ==========================================================
                ! 6. Advance x and y mean momentum equations and filter.
                ! ==========================================================

                ! Advance x-momentum equation using PROFU().
                ! Returns: UF, WUBOT
                CALL PROFU()

                ! Advance y-momentum equation using PROFV().
                ! Returns: VF, WVBOT
                CALL PROFV()

                ! Apply filter
                DO K=1,KB
                    U(K)=U(K)+.5*SMOTH*(UF(K)+UB(K)-2.*U(K))
                    V(K)=V(K)+.5*SMOTH*(VF(K)+VB(K)-2.*V(K))
                    UB(K)=U(K)
                    U(K)=UF(K)
                    VB(K)=V(K)
                    V(K)=VF(K)
                ENDDO

                ! ==========================================================
                ! 7. Advance Q2 and Q2L equations to solve for new TKE, 
                !    turbulence length scale, and eddy diffusivities.  Filter.
                ! ==========================================================

                ! Advance Q2 and Q2L equations using PROFQ().
                ! Returns: Q2F, Q2LF, KM, KH, KQ, KMB, KHB
                CALL PROFQ()

                ! Apply filter
                DO K=1,KB
                    Q2(K)=Q2(K)+.5*SMOTH*(Q2F(K)+Q2B(K)-2.*Q2(K))
                    Q2B(K)=Q2(K)
                    Q2(K)=Q2F(K)
                    Q2L(K)=Q2L(K)+.5*SMOTH*(Q2LF(K)+Q2LB(K)-2.*Q2L(K))
                    Q2LB(K)=Q2L(K)
                    Q2L(K)=Q2LF(K)
                ENDDO

                ! ==========================================================
                ! 8. Calculate mixed layer depths and send results to memory.
                ! ==========================================================

                ! Calculate four kinds of mixed layer depths using MLDPTH().
                ! Returns: DP1, DP2, DP3, DP4
                CALL MLDPTH()

                ! Send 2-D fields to memory
                if ((XLAND(i,j)*OM_TMP(i,1,j)).GE.1.5) then
                    TSK(i,j) = T(1)+273.15    ! Update sea surface temperature
                    om_ml_temp(i,j) = DP1    ! Temperature-based MLD
                    om_ml_sal(i,j)  = DP2    ! Salinity-based MLD
                    om_ml_turb(i,j) = DP3    ! Turbulence-based MLD
                    om_ml_dens(i,j) = DP4    ! Density-based MLD
                endif

                ! Send 3-D fields to memory
                if ((XLAND(i,j)*OM_TMP(i,1,j)).GE.1.5) then
                    DO K=1,KB
                        om_tmp(i,K,j)=T(K)+273.15
                        om_s(i,K,j)=S(K)
                        om_u(i,K,j)=U(K)
                        om_v(i,K,j)=V(K)
                        om_tmp_b(i,K,j)=TB(K)+273.15
                        om_s_b(i,K,j)=SB(K)
                        om_u_b(i,K,j)=UB(K)
                        om_v_b(i,K,j)=VB(K)

                        om_q2(i,K,j)=Q2(K)
                        om_q2l(i,K,j)=Q2L(K)
                        om_q2_b(i,K,j)=Q2B(K)
                        om_q2l_b(i,K,j)=Q2LB(K)

                        om_l(i,K,j)=L(K)
                        om_km(i,K,j)=KM(K)
                        om_kh(i,K,j)=KH(K)
                        om_kq(i,K,j)=KQ(K)
                        om_kmb(i,K,j)=KMB(K)
                        om_khb(i,K,j)=KHB(K)
                    ENDDO
                ENDIF

            ENDDO    ! End loop over grid j index
        ENDDO    ! End loop over grid i index

CONTAINS

!-------------------------------------------------------------
!------------------ NESTED SUBROUTINES -----------------------
!-------------------------------------------------------------

    ! ======================================================
    ! DEPTH() - Set up vertical grid items.
    ! ======================================================
    !   Reused intermediate variables (not used by top-level subroutine):
    !       None
    !   Returned variables (changed for use by top-level subroutine):
    !       Z, ZZ, DZ, DZZ
    ! ======================================================

    SUBROUTINE DEPTH()

        DO K=1,KB
            Z(K)=-DELZ*FLOAT(K-1)/BOTTOM
        ENDDO
        DO K=1,KBM1
            ZZ(K)=(Z(K)+Z(K+1))*0.5
        ENDDO
        DO K=1,KBM1
            DZ(K)=Z(K)-Z(K+1)
            DZZ(K)=ZZ(K)-ZZ(K+1)
        ENDDO
        ZZ(KB)=-1.01
        DZ(KB)=0.01
        DZZ(KB)=0.01
        DZZ(KB-1)=0.01

    END SUBROUTINE DEPTH

!-------------------------------------------------------------

    ! ======================================================
    ! DENS() - Compute density minus 1 using T and S.
    ! ======================================================
    !   Reused intermediate variables:
    !       None
    !   Returned variables:
    !       RHO
    ! ======================================================

    SUBROUTINE DENS()

        DO K=1,KBM1
            RHO(K)=S(K)*S(K)*S(K)*6.76786136E-6-S(K)*S(K)* &
                   4.8249614E-4+S(K)*8.14876577E-1-0.22584586E0
        ENDDO
        DO K=1,KBM1
            RHO(K)=RHO(K) &
                  *(T(K)*T(K)*T(K)*1.667E-8-T(K)*T(K)*8.164E-7 &
                  +T(K)*1.803E-5)
        ENDDO
        DO K=1,KBM1
            RHO(K)=RHO(K)+1.0 &
                  -T(K)*T(K)*T(K)*1.0843E-6+T(K)*T(K)*9.8185E-5 &
                  -T(K)*4.786E-3
        ENDDO
        DO K=1,KBM1
            RHO(K)=RHO(K) &
                  *(S(K)*S(K)*S(K)*6.76786136E-6-S(K)*S(K) &
                  *4.8249614E-4+S(K)*8.14876577E-1+3.895414E-2)
        ENDDO
        DO K=1,KBM1
            RHO(K)=RHO(K) &
                  -(T(K)-3.98)**2 *(T(K)+283.)/(503.57*(T(K)+67.26))
            RHO(K)=RHO(K)*1.E-3
        ENDDO

    END SUBROUTINE DENS

!-------------------------------------------------------------

    ! ======================================================
    ! PROFT() - Advance temperature equation.
    ! ======================================================
    !   Reused intermediate variables:
    !       DH, DM, A, C, VH, VHP, TERM
    !   Returned variables:
    !       TF
    ! ======================================================

    SUBROUTINE PROFT()

        ! SETS NUMBER OF BANDS TO CALCULATE SOLAR PENETRATION OVER
        IF (SURABS.EQ.0) THEN    ! Don't absorb any solar radiation at surface
            LOOPDO=NBAND
            WFTOTL=WTSURF
        ELSE    ! Absorb last component of solar radiation at surface
            LOOPDO=NBAND-1
            WFTOTL=WTSURF+ALPH(NBAND)*WTSOLR
        ENDIF

        DH=H

        ! COMPUTE SOLAR RADIATION AVAILABLE AT EACH MODEL LEVEL
        DO K=1,KB
            SWRAD(K)=0.
            DO N=1,LOOPDO
                SWRAD(K)=SWRAD(K)+ALPH(N)*EXP(Z(K)*DH*SCLI(N))
            ENDDO
            SWRAD(K)=SWRAD(K)*WTSOLR
        ENDDO

        ! CALCULATE TEMPERATURE CHANGE DUE TO SOLAR RADIATION AT EACH MODEL LEVEL
        DO K=1,KBM1
            DM(K)=TB(K)-DT2*(SWRAD(K)-SWRAD(K+1))/(DZ(K)*DH)
        ENDDO

        DO K=2,KBM1
            A(K-1)=-DT2*(KH(K)+KHB(K))/(DZ(K-1)*DZZ(K-1)*DH*DH)
            C(K)=-DT2*(KH(K)+KHB(K))/(DZ(K)*DZZ(K-1)*DH*DH)
        ENDDO

        ! Upper BC: Prescribed surface heat flux
        VH(1)=A(1)/(A(1)-1.)
        VHP(1)=(-DT2*WFTOTL/(-DZ(1)*DH)-DM(1))/(A(1)-1.)

        DO K=2,KBM2
            TERM=1.0/(A(K)+C(K)*(1.0-VH(K-1))-1.0)
            VH(K)=A(K)*TERM
            VHP(K)=(C(K)*VHP(K-1)-DM(K))*TERM
        ENDDO

       ! IF (KB .LT. KB_tot) THEN 
       !    ! IF RESIDUAL SW RADIATION IS ABSORBED AT THE BOTTOM: 
       !    TF(KBM1)=(C(KBM1)*VHP(KBM2)-DM(KBM1)+(DT2*SWRAD(KB)/DZ(KB))) &
       !            /(C(KBM1)*(1.-VH(KBM2))-1.)
       ! ELSE
       !    ! Lower BC: Zero heat flux at bottom
       !    TF(KBM1)=(C(KBM1)*VHP(KBM2)-DM(KBM1)) &
       !            /(C(KBM1)*(1.-VH(KBM2))-1.)
       ! ENDIF
       
       ! Lower BC: Zero heat flux at bottom
       TF(KBM1)=(C(KBM1)*VHP(KBM2)-DM(KBM1)) &
               /(C(KBM1)*(1.-VH(KBM2))-1.)

        DO K=2,KBM1
            KI=KB-K
            TF(KI)=VH(KI)*TF(KI+1)+VHP(KI)
        ENDDO

        ! Reset reused intermediate variables
        DH   = 1./0.
        DM   = 1./0.
        A    = 1./0.
        C    = 1./0.
        VH   = 1./0.
        VHP  = 1./0.
        TERM = 1./0.

    END SUBROUTINE PROFT

!-------------------------------------------------------------

    ! ======================================================
    ! PROFS() - Advance salinity equation.
    ! ======================================================
    !   Reused intermediate variables:
    !       DH, DM, A, C, VH, VHP, TERM
    !   Returned variables:
    !       SF
    ! ======================================================

    SUBROUTINE PROFS()

        DH=H

        DO K=1,KB
            DM(K)=SB(K)
        ENDDO

        DO K=2,KBM1
            A(K-1)=-DT2*(KH(K)+KHB(K))/(DZ(K-1)*DZZ(K-1)*DH*DH)
            C(K)=-DT2*(KH(K)+KHB(K))/(DZ(K)*DZZ(K-1)*DH*DH)
        ENDDO

        ! Upper BC: Prescribed surface salt flux
        VH(1)=A(1)/(A(1)-1.)
        VHP(1)=(-DT2*WSSURF/(-DZ(1)*DH)-DM(1))/(A(1)-1.)

        DO K=2,KBM2
            TERM=1./(A(K)+C(K)*(1.-VH(K-1))-1.)
            VH(K)=A(K)*TERM
            VHP(K)=(C(K)*VHP(K-1)-DM(K))*TERM
        ENDDO

        ! Lower BC: Zero salt flux at bottom
        SF(KBM1)=(C(KBM1)*VHP(KBM2)-SB(KBM1)) &
                /(C(KBM1)*(1.-VH(KBM2))-1.)

        DO K=2,KBM1
            KI=KB-K
            SF(KI)=VH(KI)*SF(KI+1)+VHP(KI)
        ENDDO

        ! Reset reused intermediate variables
        DH   = 1./0.
        DM   = 1./0.
        A    = 1./0.
        C    = 1./0.
        VH   = 1./0.
        VHP  = 1./0.
        TERM = 1./0.

    END SUBROUTINE PROFS

!-------------------------------------------------------------

    ! ======================================================
    ! PROFU() - Advance x-momentum equation.
    ! ======================================================
    !   Reused intermediate variables:
    !       DH, A, C, VH, VHP, CBC
    !   Returned variables:
    !       UF, WUBOT
    ! ======================================================

    SUBROUTINE PROFU()

        DH=H
        DO K=1,KB
            C(K)=0.5*(KM(K)+KM(K)+KMB(K)+KMB(K))
        ENDDO
        DO K=2,KBM1
            A(K-1)=-DT2*C(K)/(DZ(K-1)*DZZ(K-1)*DH*DH)
            C(K)=-DT2*C(K)/(DZ(K)*DZZ(K-1)*DH*DH)
        ENDDO

        ! Upper BC: Prescribed momentum flux
        VH(1)=A(1)/(A(1)-1.)
        VHP(1)=(-DT2*WUSURF/(-DZ(1)*DH)-UF(1))/(A(1)-1.)

        DO K=2,KBM2
            VHP(K)=1./(A(K)+C(K)*(1.-VH(K-1))-1.)
            VH(K)=A(K)*VHP(K)
            VHP(K)=(C(K)*VHP(K-1)-UF(K))*VHP(K)
        ENDDO

        ! Bottom BC: Calculate CBC to include bottom friction or
        ! set CBC = 0 for no bottom friction.
        ! Apply bottom friction only if needed, i.e. shelf
        IF (KB .LT. KB_tot) THEN 
        !       Logarithmic law, using bottom roughness length
        !        CBC=AMAX1(.0025,.16/ALOG((ZZ(KBM1)-Z(KB))*DH/.01)**2)
        !        CBC=CBC*SQRT(UB(KBM1)**2+(.25*(VB(KBM1) &
        !           +VB(KBM1)+VB(KBM1)+VB(KBM1)))**2)
        !
        !       Quadratic law
                CBC=0.0025 !constant drag value
                CBC=CBC*SQRT((.25*(UB(KBM1)+UB(KBM1) &
                   +UB(KBM1)+UB(KBM1)))**2+VB(KBM1)**2)
        ELSE ! open ocean point, no bottom friction
                CBC=0.
        ENDIF
        !
        UF(KBM1)=(C(KBM1)*VHP(KBM2)-UF(KBM1)) &
                /(CBC*DT2/(-DZ(KBM1)*DH)-1.-(VH(KBM2)-1.)*C(KBM1))

        DO K=2,KBM1
            KI=KB-K
            UF(KI)=VH(KI)*UF(KI+1)+VHP(KI)
        ENDDO
        WUBOT=-CBC*UF(KBM1)

        ! Reset reused intermediate variables
        DH   = 1./0.
        A    = 1./0.
        C    = 1./0.
        VH   = 1./0.
        VHP  = 1./0.
        CBC = 1./0.

    END SUBROUTINE PROFU

!-------------------------------------------------------------

    ! ======================================================
    ! PROFV() - Advance y-momentum equation.
    ! ======================================================
    !   Reused intermediate variables:
    !       DH, A, C, VH, VHP, CBC
    !   Returned variables:
    !       VF, WVBOT
    ! ======================================================

    SUBROUTINE PROFV()

        DH=H
        DO K=1,KB
            C(K)=0.5*(KM(K)+KM(K)+KMB(K)+KMB(K))
        ENDDO
        DO K=2,KBM1
            A(K-1)=-DT2*C(K)/(DZ(K-1)*DZZ(K-1)*DH*DH)
            C(K)=-DT2*C(K)/(DZ(K)*DZZ(K-1)*DH*DH)
        ENDDO

        ! Upper BC: Prescribed momentum flux
        VH(1)=A(1)/(A(1)-1.)
        VHP(1)=(-DT2*WVSURF/(-DZ(1)*DH)-VF(1))/(A(1)-1.)

        DO K=2,KBM2
            VHP(K)=1./(A(K)+C(K)*(1.-VH(K-1))-1.)
            VH(K)=A(K)*VHP(K)
            VHP(K)=(C(K)*VHP(K-1)-VF(K))*VHP(K)
        ENDDO

        ! Bottom BC: Calculate CBC to include bottom friction or
        ! set CBC = 0 for no bottom friction.
        ! Apply bottom friction only if needed, i.e. shelf
        IF (KB .LT. KB_tot) THEN 
        !       Logarithmic law, using bottom roughness length
        !        CBC=AMAX1(.0025,.16/ALOG((ZZ(KBM1)-Z(KB))*DH/.01)**2)
        !        CBC=CBC*SQRT((.25*(UB(KBM1)+UB(KBM1) &
        !           +UB(KBM1)+UB(KBM1)))**2+VB(KBM1)**2)
        !
        !       Quadratic law
                CBC=0.0025 ! constant drag value
                CBC=CBC*SQRT((.25*(UB(KBM1)+UB(KBM1) &
                   +UB(KBM1)+UB(KBM1)))**2+VB(KBM1)**2)
        ELSE ! open ocean point, no bottom friciton
                CBC=0.
        ENDIF
        VF(KBM1)=(C(KBM1)*VHP(KBM2)-VF(KBM1)) &
                /(CBC*DT2/(-DZ(KBM1)*DH)-1.-(VH(KBM2)-1.)*C(KBM1))

        DO K=2,KBM1
            KI=KB-K
            VF(KI)=VH(KI)*VF(KI+1)+VHP(KI)
        ENDDO
        WVBOT=-CBC*VF(KBM1)

        ! Reset reused intermediate variables
        DH   = 1./0.
        A    = 1./0.
        C    = 1./0.
        VH   = 1./0.
        VHP  = 1./0.
        CBC = 1./0.

    END SUBROUTINE PROFV

!-------------------------------------------------------------

    ! ======================================================
    ! PROFQ() - Advance Q2 and Q2L equations and calculate 
    !     new values for eddy diffusivities.
    ! ======================================================
    !   Reused intermediate variables:
    !       DH, A, C, VH, VHP
    !   Returned variables:
    !       Q2F, Q2LF, KM, KH, KQ, KMB, KHB
    ! ======================================================

    SUBROUTINE PROFQ()

        EPS=1.E-30    ! Different value used for EPS within PROFQ
        UMOLPR=UMOL/PR
        DH=H
        DO K=2,KBM1
            A(K)=-0.5*DT2*(KQ(K+1)+KQ(K)+KHB(K+1)+KHB(K))/ &
                 (DZZ(K-1)*DZ(K)*DH*DH)
            C(K)=-0.5*DT2*(KQ(K-1)+KQ(K)+KHB(K-1)+KHB(K))/ &
                 (DZZ(K-1)*DZ(K-1)*DH*DH)
        ENDDO

        ! 1. Solve the Q2 equation =========================
        !     DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B

        ! Upper BC
        USTAR=SQRT(SQRT(WUSURF**2+WVSURF**2))    ! u* on the water side
        L(1)=0.0
        VH(1)=0.0
        VHP(1)=CONST1*(USTAR**2)

        ! Lower BC
        Q2F(KB)=0.5*SQRT((WUBOT+WUBOT)*(WUBOT+WUBOT) &
               +(WVBOT+WVBOT)*(WVBOT+WVBOT))*CONST1
        Q2LF(KB)=0.
        L(KB)=0.

        DO K=2,KBM1
            ZIN=1.0/(DZZ(K-1)*DH)
            GDRDZ(K)=GEE*((RHO(K-1)-RHO(K))*ZIN)
            DUDZ=0.5*(U(K-1)+U(K-1)-U(K)-U(K))*ZIN
            DVDZ=0.5*(V(K-1)+V(K-1)-V(K)-V(K))*ZIN
            SPROD(K)=KM(K)*(DUDZ*DUDZ+DVDZ*DVDZ)*SEF
            BPROD(K)=KH(K)*GDRDZ(K)
        ENDDO

        DO K=2,KBM1
            Q2B(K)=ABS(Q2B(K))
            Q2LB(K)=ABS(Q2LB(K))
            DTEF(K)=Q2B(K)*SQRT(Q2B(K))/(B1*Q2LB(K)+EPS)
        ENDDO

        DO K=2,KBM1
            PROD(K)=BPROD(K)+SPROD(K)
            VHP(K)=1.0/(A(K)+C(K)*(1.0-VH(K-1))-(2.0*DT2*DTEF(K)+1.0))
            VH(K)=A(K)*VHP(K)
            VHP(K)=(-2.0*DT2*PROD(K)+C(K)*VHP(K-1)-Q2B(K))*VHP(K)
        ENDDO

        DO K=1,KBM1
            KI=KB-K
            Q2F(KI)=VH(KI)*Q2F(KI+1)+VHP(KI)
        ENDDO

        ! 2. Solve the Q2L equation =======================
        !     DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB

        VH(1)=0.
        VHP(1)=Q2F(1)*L(1)

        DO K=2,KBM1
            A(K)=-0.5*DT2*(E4*(KQ(K+1)+KQ(K))+KHB(K+1)+KHB(K))/ &
                 (DZZ(K-1)*DZ(K)*DH*DH)
            C(K)=-0.5*DT2*(E4*(KQ(K-1)+KQ(K))+KHB(K-1)+KHB(K))/ &
                 (DZZ(K-1)*DZ(K-1)*DH*DH)
        ENDDO

        DO K=2,KBM1
            DTEF(K)=DTEF(K)*(1.0+E2*((1./ABS(Z(K)-Z(1)) &
                   +1.0/ABS(Z(K)-Z(KB)))*L(K)/(DH*KAPPA))**2)
            VHP(K)=1.0/(A(K)+C(K)*(1.0-VH(K-1))-(DT2*DTEF(K)+1.0))
            VH(K)=A(K)*VHP(K)
            VHP(K)=(-DT2*(E1*L(K)*(SPROD(K)+E3*BPROD(K))) &
                  +C(K)*VHP(K-1)-Q2LB(K))*VHP(K)
        ENDDO

        DO K=1,KBM1
            KI=KB-K
            Q2LF(KI)=VH(KI)*Q2LF(KI+1)+VHP(KI)
        ENDDO

        DO K=2,KBM1
            IF(Q2F(K).LT.SMALL .OR. Q2LF(K).LT.SMALL) THEN
                Q2F(K)=SMALL
                Q2LF(K)=SMALL
            ENDIF
        ENDDO

        ! 3. Solve for new KM, KH, KQ, KMB, KHB ===========

        DO K=2,KBM1
            L(K)=Q2LF(K)/(Q2F(K)+EPS)
            QLI(K)=L(K)*L(K)/Q2F(K)
            ZIN=1.0/(DZZ(K-1)*DH)
            DUDZ=0.5*(U(K-1)+U(K-1)-U(K)-U(K))*ZIN
            DVDZ=0.5*(V(K-1)+V(K-1)-V(K)-V(K))*ZIN
            RII(K)=-GDRDZ(K)/((DUDZ**2+DVDZ**2+SMALL)*SEF)
            RII(K)=AMIN1(RII(K),0.7)

            ! Constant background mixing assumed.
            ! (Gargett background mixing option removed per CAC guidance.)
            KMB(K)=1.E-5
            KHB(K)=1.E-5

            IF(ABS(DUDZ).LT.EPS) DUDZ=EPS
            IF(ABS(DVDZ).LT.EPS) DVDZ=EPS
            GM(K)=QLI(K)*(DUDZ*DUDZ+DVDZ*DVDZ)*SEF
 
            GH(K)=QLI(K)*GDRDZ(K)
            GH(K)=AMAX1(GH(K),GHMIN)
            GH(K)=AMIN1(GH(K),GHMAX)
        ENDDO

        GH(1)=0.
        GH(KB)=0.
        GM(1)=GM(2)
        GM(KB)=GM(KBM1)
        RII(1)=RII(2)
        RII(KB)=RII(KBM1)
        KHB(1)=0.
        KMB(1)=0.
        KHB(KB)=0.
        KMB(KB-1)=0.

        DO K=1,KB
            SH(K)=CF1/(1.-CF2*GH(K))
            SM(K)=(CF3+CF4*SH(K)*GH(K))/(1.-CF5*GH(K))
        ENDDO

        DO K=1,KB
            POD(K)=B1*(SM(K)*GM(K)+SH(K)*GH(K))
            RIF(K)=-SH(K)*GH(K)/(SM(K)*GM(K)+EPS)

            ! SHEAR INSTABILITY RICHARDSON NUMBER DEPENDENT MIXING
            ! per LARGE et al.
            ! Peters et al. option removed per CAC guidance.
            IF(Q2(K).LT.1.E-6) THEN
                KMB(K)=1.E-4
                KHB(K)=5.E-5
                IF (RII(K).LT.0.0) THEN
                    KMB(K)=KMB(K)+5.E-3
                    KHB(K)=KHB(K)+5.E-3
                ENDIF
                IF (RII(K).GT.0.0.AND.RII(K).LT.0.7) THEN
                    KMB(K)=KMB(K)+5.E-3*(1.-(RII(K)/0.7)**2)**3
                    KHB(K)=KHB(K)+5.E-3*(1.-(RII(K)/0.7)**2)**3
                ENDIF
            ENDIF
        ENDDO

        DO K=1,KB
            KN(K)=L(K)*SQRT(Q2(K))
            KQ(K)=KN(K)*SQ1*(1.-FCT)*SH(K) + KQ(K)*FCT
            KM(K)=KN(K)*SM(K)*(1.-FCT) + KM(K)*FCT
            KH(K)=KN(K)*SH(K)*(1.-FCT) + KH(K)*FCT
        ENDDO

        EPS=1.E-10    ! Reset EPS to its value outside of PROFQ()

        ! Reset reused intermediate variables
        DH   = 1./0.
        A    = 1./0.
        C    = 1./0.
        VH   = 1./0.
        VHP  = 1./0.

    END SUBROUTINE PROFQ

!-------------------------------------------------------------

    ! ======================================================
    ! MLDPTH() - Calculate mixed-layer depths.  DP1 (temp-based),
    !     DP2 (salinity-based), and DP3 (turbulence-based)
    !     come from the KC 1994 code, and DP4 (density-based)
    !     comes from the KC 2004 code.
    ! ======================================================
    !   Reused intermediate variables:
    !       None
    !   Returned variables:
    !       DP1, DP2, DP3, DP4
    ! ======================================================

    SUBROUTINE MLDPTH()

        DP1=EPS    ! Temperature-based MLD
        DP2=EPS    ! Salinity-based MLD
        DP3=EPS    ! Turbulence-based MLD
        DP4=EPS    ! Density-based MLD

        DO K=2,KBM2
            IF(T(K).LT.(T(1)-DELT)) THEN
                DP1=-H*(ZZ(K-1)-(T(K-1)-T(1)+DELT)*DZZ(K-1)/(T(K-1)-T(K)+EPS))
                GO TO 120
            ENDIF
        ENDDO
    120 CONTINUE

        DO K=2,KBM2
            IF(S(K).GT.(S(1)+DELS)) THEN
                DP2=-H*(ZZ(K-1)-(S(K-1)-S(1)-DELS)*DZZ(K-1)/(S(K-1)-S(K)+EPS))
                GO TO 220
            ENDIF
        ENDDO
    220 CONTINUE

        DO K=2,KB
            IF(Q2(K).LT.10.*EPS) THEN
                DP3=-H*Z(K)
                GO TO 320
            ENDIF
        ENDDO
    320 CONTINUE

        DO K=2,KBM2
            IF(RHO(K).GT.(RHO(1)+DELR)) THEN
                DP4=-H*(ZZ(K-1)-(RHO(K-1)-RHO(1)-DELR)*DZZ(K-1)/(RHO(K-1)-RHO(K)+EPS))
                GO TO 420
            ENDIF
        ENDDO
    420 CONTINUE

    END SUBROUTINE MLDPTH

!-------------------------------------------------------------
!----------------- END NESTED SUBROUTINES --------------------
!-------------------------------------------------------------

    END SUBROUTINE KC1994OML

!-------------------------------------------------------------

END MODULE module_sf_kc1994oml


