module EDBtranMod

  !-------------------------------------------------------------------------------------
  ! Description:
  ! 
  ! ------------------------------------------------------------------------------------

  use EDPftvarcon       , only : EDPftvarcon_inst
  use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm 
  use FatesConstantsMod , only : itrue,ifalse,nearzero
  use FatesConstantsMod , only : nocomp_bareground
  use EDTypesMod        , only : ed_site_type
  use FatesPatchMod,      only : fates_patch_type
  use EDParamsMod,        only : maxpft
  use EDParamsMod,        only : soil_tfrz_thresh
  use FatesCohortMod,     only : fates_cohort_type
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use FatesInterfaceTypesMod , only : bc_in_type, &
       bc_out_type, &
       numpft
  use FatesInterfaceTypesMod , only : hlm_use_planthydro
  use FatesGlobals      , only : fates_log
  use FatesAllometryMod , only : set_root_fraction
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesAllometryMod, only : bsap_allom

  !
  implicit none
  private


  logical, parameter :: debug = .false.
  
  public :: btran_ed
  public :: get_active_suction_layers
  public :: check_layer_water

contains 

  ! ====================================================================================

  logical function check_layer_water(h2o_liq_vol, tempk)

    implicit none
    ! Arguments
    real(r8),intent(in) :: h2o_liq_vol
    real(r8),intent(in) :: tempk

    check_layer_water = .false.

    if ( h2o_liq_vol .gt. 0._r8 ) then
       if ( tempk .gt. soil_tfrz_thresh + tfrz) then
          check_layer_water = .true.
       end if
    end if
    return
  end function check_layer_water

  ! =====================================================================================

  subroutine get_active_suction_layers(nsites, sites, bc_in, bc_out)

    ! Arguments

    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    ! !LOCAL VARIABLES:
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    !------------------------------------------------------------------------------

    do s = 1,nsites
       if (bc_in(s)%filter_btran) then
          do j = 1,bc_in(s)%nlevsoil
             bc_out(s)%active_suction_sl(j) = check_layer_water( bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j) )
          end do
       else
          bc_out(s)%active_suction_sl(:) = .false.
       end if
    end do

  end subroutine get_active_suction_layers

  ! =====================================================================================

  subroutine btran_ed( nsites, sites, bc_in, bc_out)

    use FatesPlantHydraulicsMod, only : BTranForHLMDiagnosticsFromCohortHydr


    ! ---------------------------------------------------------------------------------
    ! Calculate the transpiration wetness function (BTRAN) and the root uptake
    ! distribution (ROOTR).
    ! Boundary conditions in: bc_in(s)%eff_porosity_sl(j)    unfrozen porosity
    !                         bc_in(s)%watsat_sl(j)          porosity
    !                         bc_in(s)%active_uptake_sl(j)   frozen/not frozen
    !                         bc_in(s)%smp_sl(j)             suction
    ! Boundary conditions out: bc_out(s)%rootr_pasl          root uptake distribution
    !                          bc_out(s)%btran_pa            wetness factor
    ! ---------------------------------------------------------------------------------

    ! Arguments

    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)

    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type),pointer             :: cpatch ! Current Patch Pointer
    type(fates_cohort_type),pointer            :: ccohort ! Current cohort pointer
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    integer  :: ifp               ! patch vector index for the site
    integer  :: ft                ! plant functional type index
    real(r8) :: smp_node          ! matrix potential
    real(r8) :: rresis            ! suction limitation to transpiration independent
    ! of root density
    real(r8) :: pftgs(maxpft)     ! pft weighted stomatal conductance m/s
    real(r8) :: temprootr
    real(r8) :: sum_pftgs         ! sum of weighted conductances (for normalization)
    real(r8), allocatable :: root_resis(:,:)  ! Root resistance in each pft x layer'
    real(r8), allocatable :: root_resis_cohort(:)  ! Root resistance in each cohort'

    !------------------------------------------------------------------------------
    ! added variables
  
    real(r8) :: mm_to_MPa, rho_w, g, f_wilt, f_star
    real(r8) :: Zm, RAI, dr
    real(r8) :: hc, LAI, Ks, Ps0
    real(r8) :: k_xl_max, k_max, Px50, Pg50
    real(r8) :: K_p_max, K_sr_max
    real(r8) :: pi_F, pi_R, pi_T, pi_S
    real(r8) :: x_ww, beta_ww
    real(r8) :: beta_c, y_c, xx_c, xxx_c, x_c, p_c, smpsc_t
    real(r8) :: beta_o, y_o, xx_o, xxx_o, x_o, p_o, smpso_t
    real(r8) :: btran_cohort
    real(r8) :: a_sapwood, c_sap_dummy
    real(r8) :: air_temp, netrad
    real(r8) :: Eo, lambda_v, delta, gamma_psy, tc, cp, p0, es

    !------------------------------------------------------------------------------

    associate(                                 &
    !     smpsc     => EDPftvarcon_inst%smpsc              , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
    !     smpso     => EDPftvarcon_inst%smpso                &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
          Pg50_pft      => EDPftvarcon_inst%hydr_p50_gs    , &
          Px50_pft      => EDPftvarcon_inst%hydr_p50_node  , &
          k_max_pft     => EDPftvarcon_inst%hydr_kmax_node  &
         )

    do s = 1,nsites ! grid cell

       allocate(root_resis(numpft,bc_in(s)%nlevsoil))
       allocate(root_resis_cohort(bc_in(s)%nlevsoil))

       bc_out(s)%rootr_pasl(:,:) = 0._r8

       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch)) ! patch ~ plot that has an age and can have multiple cohorts/pft

          ifp = cpatch%patchno
          
          if_bare: if(cpatch%nocomp_pft_label.ne.nocomp_bareground)then ! only for veg patches

            
             ! COHORT LOOP 
             cpatch%btran_ft(:) = 0.0_r8
             root_resis(:,:) = 0.0_r8
             ccohort => cpatch%tallest
             do while (associated(ccohort))
                  ft=ccohort%pft
                  call bsap_allom(ccohort%dbh,ft,ccohort%crowndamage, &
                       ccohort%canopy_trim, ccohort%efstem_coh, a_sapwood,c_sap_dummy)
                  call set_root_fraction(sites(s)%rootfrac_scr, ft, sites(s)%zi_soil, &
                       bc_in(s)%max_rooting_depth_index_col ) 

                !cpatch%btran_ft(ft) = 0.0_r8
                ! calculating by cohort verusus patch
                btran_cohort = 0.0_r8 
                do j = 1,bc_in(s)%nlevsoil

                   ! Calculations are only relevant where liquid water exists
                   ! see clm_fates%wrap_btran for calculation with CLM/ALM

                    ! US-Me2: Sandy loam / Ponderosa pine hydraulic thresholds
                    ! Computes smpsc_t (psi_wilt), smpso_t (psi_star), and rresis from smp_node
                    ! ---------------------------------------------------------------------------------
                    
                    ! set Parameters and Conversion factors
                    mm_to_MPa = 9.80665e-6_r8  ! Conversion factor: 1 mm H2O ≈ 9.80665e-6 MPa
                    rho_w    = 1000.0_r8      ! [kg m-3] Density of liquid water
                    g        = 9.81_r8        ! [m s-2] Gravitational acceleration
                    f_wilt   = 0.05_r8        ! [-] fraction of beta_ww at wilting point
                    f_star   = 0.95_r8        ! [-] fraction of beta_ww at incipient stomatal closure
                    cp   = 1004.0_r8          ! [J kg-1 K-1]
                    p0   = 101325.0_r8        ! [Pa]

                    ! not important to change but use clm parameter if they exist
                    RAI      = 10.0_r8        ! [m2 m-2] Root area index 
                    dr       = 0.0005_r8      ! [m] root diameter
                    Zm       = 0.5_r8         ! [m] Rooting depth

                    ! get FATES & CLM parameters 
                    k_max = k_max_pft(ft,2)    ! = 0.0001_r8 [kg/MPa/m/s] maximum xylem conductivity per unit conducting xylem area
                    Px50 = Px50_pft(ft,2)       ! = -2.25_r8 [MPa] leaf potential at 50% stem conductance loss (scale with sapwood area)  
                    Pg50 = Pg50_pft(ft)         ! = -1.5_r8  [MPa] leaf potential at 50% stomatal conductance loss
                    
                    ! get model variables
                    !LAI      = 2.0_r8        ! [m2 m-2] <- cohort specific, don't need because now kx scaled by hc
                    hc       = ccohort%height ! [m] <- cohort specific
                    air_temp = bc_in(s)%tgcm_pa(ifp) ! [K] air temperature to calculate Eo
                    netrad = bc_in(s)%netrad_net_pa(ifp) ! [M/m2] net radiation to calculate Eo 

                    ! calculate Eo with Priestley Taylor in [m /s]
                    ! to do ---> make variant with Gcmax*VPD
                    tc = air_temp - 273.15_r8                                   ! [C] air temp
                    lambda_v = (2.501_r8 - 2361.0_r8 * (air_temp - 273.15_r8))  ! [J/kg] latent heat of vaporization
                    gamma_psy = cp * p0 / (0.622_r8 * lambda_v)                 ! [Pa/K] psychrometric constant
                    es = 610.8_r8 * exp(17.27_r8 * tc/(tc + 237.3_r8))          ! [Pa]
                    delta = es * 17.27_r8 * 237.3_r8 / ((tc + 237.3_r8)**2)     ! [Pa/K] slope of saturation vapour pressure curve
                    if (netrad > 0.d0) then
                       Eo = 1.26_r8 * delta / (delta + gamma_psy) * netrad / (lambda_v * rho_w) ! [m /s ] 
                    else
                       Eo = 0.0_r8
                    end if
                    
                   if ( check_layer_water(bc_in(s)%h2o_liqvol_sl(j),bc_in(s)%tempk_sl(j)) )  then

                        ! get soil level parameters 
                        Ps0 = bc_in(s)%sucsat_sisl(j)  ! [mm in clm] Soil water potential near saturation 
                        Ks = bc_in(s)%hksat_sisl(j)    ! [mm/s] Saturated hydraulic conductivity ...  double check units?

                       ! Conductivities
                        K_p_max  = k_max * a_sapwood / hc / rho_w                                  ! [m/s/MPa]
                        K_sr_max = (Ks/1000.0_r8) * sqrt(RAI / (dr * Zm)) * 1.0e6_r8 / (rho_w * g) ! [m/s/]
                        
                        ! Pi-groups
                        pi_F = -Eo / (K_p_max * Pg50)
                        pi_R = Pg50 / Px50
                        pi_T = -(K_sr_max * Pg50) / Eo
                        pi_S = Pg50 / (Ps0 * mm_to_MPa) ! convert Ps0 from mm to MPa because Pg50 in MPa
                        
                        ! Well-watered beta - maybe not use
                        x_ww = (pi_F / 2.0_r8 + 1.0_r8)**2 - 2.0_r8 * pi_F * pi_R
                        beta_ww = 1.0_r8 - (1.0_r8 / (2.0_r8 * pi_R)) * (1.0_r8 + pi_F / 2.0_r8 - sqrt(x_ww))
                        
                        ! smpsc soil water potential close
                        beta_c = f_wilt !* beta_ww
                        y_c  = 1.0_r8 / (2.0_r8 * beta_c * pi_S / pi_T)
                        xx_c = 1.0_r8 - beta_c
                        xxx_c = pi_F * beta_c / (1.0_r8 - xx_c * pi_R)
                        x_c = 1.0_r8 + (2.0_r8 * xx_c - xxx_c) * 4.0_r8 * beta_c * pi_S * pi_S / pi_T
                        p_c = y_c * (sqrt(x_c) - 1.0_r8)
                        smpsc_t = p_c * Ps0 ! in mm because Ps0 and smp_sl in mm
                        
                        ! smpso soil water potential open 
                        beta_o = f_star !* beta_ww
                        y_o = 1.0_r8 / (2.0_r8 * beta_o * pi_S / pi_T)
                        xx_o = 1.0_r8 - beta_o
                        xxx_o = pi_F * beta_o / (1.0_r8 - xx_o * pi_R)
                        x_o = 1.0_r8 + (2.0_r8 * xx_o - xxx_o) * 4.0_r8 * beta_o * pi_S * pi_S / pi_T
                        p_o = y_o * (sqrt(x_o) - 1.0_r8)
                        smpso_t = p_o * Ps0 ! in mm because Ps0 and smp_sl in mm
 
                      smp_node = max(smpsc_t, bc_in(s)%smp_sl(j))

                      rresis = min( (bc_in(s)%eff_porosity_sl(j)/bc_in(s)%watsat_sl(j))*               &
                            (smp_node - smpsc_t) / (smpso_t - smpsc_t), 1._r8)

                      root_resis_cohort(j) = sites(s)%rootfrac_scr(j)*rresis
                     
                     
                      ! root water uptake is not linearly proportional to root density,
                      ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)

                      !cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + root_resis(ft,j)
                      ! calculating by cohort verusus patch
                      btran_cohort = btran_cohort + root_resis_cohort(j)

                   else
                      root_resis_cohort(j) = 0._r8
                   end if

                end do !j
               
                cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + btran_cohort
                root_resis(ft,:) = root_resis(ft,:) + root_resis_cohort(:)

                ccohort => ccohort%shorter
             end do ! cohort
             
             cpatch%btran_ft(ft) = cpatch%btran_ft(ft) / cpatch%num_cohorts
             root_resis = root_resis / cpatch%num_cohorts

             ! pft in the patch
             do ft = 1,numpft 
             
             ! Normalize root resistances to get layer contribution to ET
             do j = 1,bc_in(s)%nlevsoil  
                   if (cpatch%btran_ft(ft)  >  nearzero) then
                      root_resis(ft,j) = root_resis(ft,j)/cpatch%btran_ft(ft)
                   else
                      root_resis(ft,j) = 0._r8
                   end if
                end do
            end do


             ! PFT-averaged point level root fraction for extraction purposese.
             ! The cohort's conductance g_sb_laweighted, contains a weighting factor
             ! based on the cohort's leaf area. units: [m/s] * [m2]

             pftgs(1:maxpft) = 0._r8
             ccohort => cpatch%tallest
             do while(associated(ccohort))
                pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%g_sb_laweight
                ccohort => ccohort%shorter
             end do

             ! Process the boundary output, this is necessary for calculating the soil-moisture
             ! sink term across the different layers in driver/host.  Photosynthesis will
             ! pass the host a total transpiration for the patch.  This needs rootr to be
             ! distributed over the soil layers.
             sum_pftgs = sum(pftgs(1:numpft))

             do j = 1, bc_in(s)%nlevsoil
                bc_out(s)%rootr_pasl(ifp,j) = 0._r8
                do ft = 1,numpft
                   if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                      !bit-retart test as a result? FIX(RF,032414)  
                      bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                           root_resis(ft,j) * pftgs(ft)/sum_pftgs
                   else
                      bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j) + &
                           root_resis(ft,j) * 1._r8/real(numpft,r8)
                   end if
                end do
             end do

             ! Calculate the BTRAN that is passed back to the HLM
             ! used only for diagnostics. If plant hydraulics is turned off
             ! we are using the patchxpft level btran calculation

             if(hlm_use_planthydro.eq.ifalse) then
                !weight patch level output BTRAN for the
                bc_out(s)%btran_pa(ifp) = 0.0_r8
                do ft = 1,numpft
                   if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                      !bit-retart test as a result? FIX(RF,032414)   
                      bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft)  * pftgs(ft)/sum_pftgs
                   else
                      bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft) * 1./numpft
                   end if
                end do
             end if

             temprootr = sum(bc_out(s)%rootr_pasl(ifp,1:bc_in(s)%nlevsoil))

             if(abs(1.0_r8-temprootr) > 1.0e-10_r8 .and. temprootr > 1.0e-10_r8)then

                if(debug) write(fates_log(),*) 'error with rootr in canopy fluxes',temprootr,sum_pftgs
                
                do j = 1,bc_in(s)%nlevsoil
                   bc_out(s)%rootr_pasl(ifp,j) = bc_out(s)%rootr_pasl(ifp,j)/temprootr
                end do
                
             end if
          endif if_bare
          cpatch => cpatch%younger
       end do

       deallocate(root_resis)

    end do

    if(hlm_use_planthydro.eq.itrue) then
       call BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)
    end if

  end associate

end subroutine btran_ed


end module EDBtranMod
