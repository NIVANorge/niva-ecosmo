#include "fabm_driver.h"

! -----------------------------------------------------------------------------------------------------------
! This is a NIVA adaptation of the ECOSMO model, based on NERSC/ecosmo_operational.F90
! -----------------------------------------------------------------------------------------------------------
!
!
! Change Log ------------------------------------------------------------------------------------------------
! PJW 28/07/2025: Copied from nersc/ecosmo/ecosmo_operational.F90 (from FABM updated 27/07/2025).
!                 Removed temporal mean forcings and all code/parameters associated by cyanobacteria.
!                 Removed parameter vector BioC, replacing with named parameters.
!                 Tidied up formatting and deleted some commented-out experimental code.
!
! PJW 30/03/2026: Major revisions of phytoplankton growth parameterizations:
!                 1) Added simple Q10 temperature dependence and daylength dependence of maximum growth rate
!                    following Gilstad and Sakshaug (1990) (as in NERSEM, cf. PISCESv2, Aumont et al., 2015).
!                 2) Adopt modified Platt function (Platt et al., 1980) that accounts for photoinhibition
!                    within the Geider et al. (1997) C:Chl model, consistent with ERSEM (Blackford et al., 2004).
!                    This factor is combined multiplicatively with the Liebig minimum of nutrient limitation
!                    factors, rather than in one large minimum as previously --- this different treatment of
!                    light limitation seems to be more conventional (cf. ERSEM, PISCES, ECB models).
!                 3) Removed capping of light dependence within production Chl:C (max(0.01_rk,Ps_prod)).
!                    This was not consistent with Geider et al. (1997) and was not explained in
!                    Yumruktepe et al. (2022) (Y22).
!                 4) Generalized code to allow different nutrient half-saturation parameters for
!                    different phytoplankton functional groups.
!                 5) Corrected alkalinity dynamics to include negative contribution of phosphate changes.
!                 Various cleanup modifications:
!                 1) Removed secs_pr_day parameter duplicating sedy0.
!                 2) Converted d-1 to s-1 for alfa/beta parameters, avoiding use of sedy0 within subroutine do.
!                    Now we have *all* rate parameters converted to s-1 within subroutine initialize.
!                 3) Removed unused parameters: zpr
!
! PJW 13/04/2026: 1) Elaborated mortality fluxes for phytoplankton to include:
!                    a) quadratic mortality terms representing effects of e.g. aggregation and viral lysis,
!                       (new parameters m2Ps, m2Pl, default = 0 for consistency with legacy),
!                    b) enhancement factors due to nutrient limitation (similar to ERSEM formulation, new parameters
!                       amortPs, amortPl, default values 1 for consistency with legacy).
!                 2) Added missing maximum grazing rate parameters for feeding on detritus (GrZsD, GrZlD),
!                    assumed equal to (GrZsP, GrZlP) for now (as before). We also generalize to allow different
!                    maximum grazing rates on each prey type (GrZsP->GrZsPs,GrZsPl, etc.), and we include
!                    group-level cannibalism (GrZsZs, GrZlZl).
!                 3) Split growth efficiency on detritus gammaZD into two separate parameters (gammaZsD and gammaZlD,
!                    for now assumed equal).
!                 4) General format revision to calculate *fluxes* as intermediate quantities (mortPs etc.)
!                    that are subsequently re-used in the code, avoiding repeated calculations from rates.
!                    The rationale is that the small increase in memory demand should be outweighed by the
!                    increase in computation speed, and moreover the reduced risk of inconsistent repeat calculations
!                    which could lead to violations of mass conservation.
!                    NOTE! We now use grazing variables ZsonPs etc. to store *fluxes* rather than rates.
!                 5) Numbered switch functions are replaced with explicit unitary functions (f_*)
!                    combined with explicit stoichiometric ratio parameters (eta_*), rather than hard-coded values.
!                    For the denitrification stoichiometry parameter (eta_denit), the default value is set to 5.3.
!                    This was previously hard-coded to the value 5.0, with no explanation for the deviation from
!                    the original value of 5.3 in Daewel and Schrum (2013). The latter is the standard value used in
!                    (OXYDEP, ECB) models and is consistent with the rationale of Paulmier et al. (2009), also
!                    followed by (ERSEM, PISCES) and the original Baltic model of Neumann (2000).
!                    NOTE! This new default means that fabm.yaml scripts will need to set eta_denit=5.0 in order to
!                          reproduce legacy output.
!                    Further stoichiometric ratio parameters are added for the photosynthetic quotient
!                    (eta_PQ, previously assumed = 1.0), the oxic respiration quotient (eta_O2resp, previously assumed = 1.0),
!                    the anoxic sulphate respiration quotient (eta_SO4resp, previously assumed = 1.0),
!                    and the nitrification oxygen demand (eta_nitrif, previously assumed = 2.0).
!                    NOTE! Current literature suggests values larger than the default unitary values for:
!                          eta_PQ (range 1.0-1.8, Sanz-Martin et al., 2019; Laws, 1991; Burris, 1981),
!                          eta_O2resp (range 0.8-1.6, Stephens et al., 2025; Moreno et al., 2022; Tanioka and Matsumoto, 2020).
!                    NOTE! The default value of eta_SO4resp=1.0 assumes that negative oxygen concentration is
!                          interpreted as -0.5 x the H2S concentration in moles (see Neumann, 2000).
!                          This is a standard assumption in marine biogeochemical models.
!                    NOTE! The original Baltic model (Neumann, 2000) actually assumed a different values for
!                          eta_nitrif = 1.5, but the default value of 2.0 is more consistent with current
!                          standard sources (e.g. Bianchi et al., 2012; Middelburg, 2019).
!                    For nitrification we add parameters (nitrifmax, aTnitrif, KO2nitrif) replacing hard-coded
!                    values, but keeping the latter as default values (from Stigebrandt and Wulff, 1987).
!                 6) Zooplankton losses are elaborated following the sloppy-feeding/egestion/excretion paradigm
!                    of Steinberg and Landry (2017).
!                 6a) Zooplankton sloppy feeding is implemented only for large zooplankton assuming that such losses are only
!                    significant for metazoan predators (Steinberg and Landry, 2017, although MEDUSA (Yool et al., 2013)
!                    applies then uniformly to all zooplankton groups). The slopped fraction is assumed to be
!                    dependent on prey type due primarily to relative size dependence, with larger prey types
!                    resulting in larger slopped fractions (Møller, 2007). This is implemented with four new
!                    parameters (fGslpZlPs, fGslpZlPl, fGslpZlZs, fGslpZlD). Slopped material is routed entirely
!                    to labile DOM (Steinberg and Landry, 2017).
!                 6b) Zooplankton egestion is assumed proportional to ingestion via the 'gamma' assimilation parameters,
!                    so fractions (1-gammaZlP) etc. are egested. Egested material is routed to (det, dom) state
!                    variables with dissolved fractions set by new parameters (fregesZs, fregesZl) (the former
!                    uniform dom fraction parameter frr is deprecated).
!                 6c) Zooplankton excretion/respiration is parameterized as a basal maintenance term plus
!                    'activity' excretion/respiration terms proportional to food absorption (=ingestion x assim. efficiency).
!                    This is implemented using new fixed fraction parameters (fAexcZs, fAexcZl) describing the
!                    fraction of absorbed food that is excreted/respired. The total excreted/respired material
!                    is fractionated into DOM and mineral components using new parameters (fexcdomZs, fexcdomZl).
!                    The organic matter component fluxes are routed to labile dom, while the mineral components
!                    are routed to ammonium, phosphate, and DIC state variables assuming Redfield rations defined
!                    by redf. Oxygen consumption due to zooplankton respiration is assumed to be proportional to
!                    the total mineral excretion flux.
!                 7) Zooplankton mortality is elaborated to include additional temperature dependences to allow for
!                    different feeding activity dependence of higher trophic levels (new parameters q10mZs, q10mZl).
!                    Also, mesozooplankton mortality is given an additional light-dependent term to represent
!                    visual predation, following Aksnes et al. (1997) (new parameters mvZl, KE0mvZl).
!                 8) The source-independent dissolved fraction of organic waste parameter frr was split into a set of
!                    parameters allowing different fractions for different mortality/egestion sources (frmort1Ps etc.).
!                    This was motivated by an expected size-dependence: fractions for larger sources (e.g. fregesZl)
!                    and processes with contributions from aggregation mortality (e.g. frmort2Pl) may be expected
!                    to yield smaller dissolved fractions in the detrital products. Note that the dissolved fraction
!                    of mesozooplankton mortality should be reduced to account for the expected higher particulate
!                    fraction of egesta and mortality products of consumers from higher trophic levels.
!                 9) The sediment module was extensively recoded, removing the numbered switch variables and
!                    applying the new stoichiometry parameters in a manner consistent with the revised pelagic code.
!                 9a) Settling rate (Rds) is split into two parameters for detritus/opal (Rds/RdsO) and set equal
!                     to the respective sinking speeds (sinkD, sinkOPAL) when bottom stress is sub-critical, rather
!                     than applying an additional settling velocity parameter (sedimRt, now removed).
!                 9b) Consistent with the pelagic code, remineralization is set as a temperature-dependent flux
!                     independent of the mode (oxic, denitrification, sulphate reduction), but the mode determines
!                     remineralization impacts on (nitrate, oxygen, alkalinity) using the unitary switch variables
!                     (f_O2resp etc.). This is consistent with the old code, but the old code was rather unclear
!                     in the way it was coded e.g. using Rsdenit in the oxygen impacts, when denitrification does
!                     not impact oxygen.
!                     NOTE! We remove the unexplained factor of 2 in the application of the oxic remineralization
!                           rate to sed1. Instead we multiply the default value of reminSED by a factor of 2.
!                           This means that any fabm.yaml files where reminSED is specified should used double
!                           the value in application to the present code.
!                 9c) Removed term for bottom water oxygen consumption due to nitrification (2.0_rk*BioOM1*Rsa*sed1).
!                     We believe this was an erroneous term in Daewel and Schrum (2013), because ammonium in the
!                     sediments is not represented --- rather the remineralization flux is routed to ammonium
!                     in the bottom water, which means that this term will double-count the oxygen consumption
!                     due to nitrification of sediment-derived ammonium. Also it is not clear why this process
!                     should scale with (half of) the sedimentary remineralization flux (Rsa*sed1). 
!                     NOTE! Absence of this term in the new code means that the new code cannot exactly
!                           reproduce legacy simulation results. 
!                 9d) Burial fluxes are elaborated to include quadratic terms and distinguish organic matter vs. opal.
!                     This requires 3 new parameters (burialRt2, burialRtO, burialRt2O) but allows us to
!                     accommodate the hard-coded alterations in the operational code.
!                     NOTE! In order to reproduce legacy results of ecosmo_operational.F90, we need in fabm.yaml:
!                           burialRt   = 0
!                           burialRt2  = 2e-3 * 5e-5 = 1e-7
!                           burialRtO  = 0
!                           burialRt2O = 2e-3 * 5e-5 = 1e-7
!                 9e) We reverted the factor 2.0 inflation of the settling flux of opal. This was inconsistent with DS13,
!                     and appears to have already been reverted in the latest NERSC code.
!                 9f) Contribution to bottom water carbonate dynamics is reformulated for consistency with the pelagic
!                     code, and account for the contribution of phosphate release (previously neglected).
!                10) Nutrient limitation diagnostics were split between groups and overall nutrient limitation
!                    diagnostics (NutlimPs, NutlimPl) were added.
!                11) We decided not to correct phytoplankton production for excretion / exudation because:
!                    a) These fluxes are thought to have highly non-Redfield stoichiometry (e.g. in ERSEM, only the carbon
!                       productivity is corrected, so only DOC is produced). Included them within the Redfield stoichiometry
!                       assumed elsewhere in the ECOSMO model (e.g. in the remineralization of labile dom)
!                       seems difficult and may result in excessive nutrient recycling.
!                    b) These fluxes are thought to produce mainly semi-labile DOM (see e.g. ERSEM, ECB models).
!                       Routing them to labile DOM may again result in excessive nutrient recycling.
!                    c) It seems likely that the effect of these processes to reduce particulate primary
!                       production and phytoplankton net growth rates can be accounted for by adjusting the
!                       nutrient half saturation constants, which are anyway uncertain and generally in
!                       need of tuning to regional observations (e.g. Aumont et al., 2015).
!                    For similar reasons we also decided not to try to account for 'excess production', i.e.
!                    primary production that is routed directly to semi-labile DOC (see ECB model, Feng et al., 2015).
!                12) Renamed parameters (regenSi, reminSEDSi) to (dissO, dissSEDO) for opal dissolution.
!                    We note in passing that a more complex, temperature-dependent parameterization for dissolution
!                    of pelagic detrital opal is used in the NERSEM code (following Ridgwell et al., 2002), while
!                    a temperature-dependent opal dissolution rate is applied in the sediments in (N)ERSEM.
!                    For now we stick with the simple constant opal dissolution rates (dissO, dissSEDO).
!
! PJW 07/05/2026: Added detrital calcite state variable 'cal' with calcification / dissolution dynamics
!                 following the 'virtual calcite' approach of ERSEM, PISCES (Butenschon et al., 2016; Aumont et al., 2015).
!                 Note however that we follow the simpler PISCES approach to parameterizing dissolution
!                 losses inside zooplankton grazers (A15, Eqn. 76), not including the inefficiency factors that
!                 are presently included in the ERSEM code (B16, Eqn. 97).
!                 The code prior to this update was renamed niva_ecosmo->niva_ecosmo_nocalc.
!                 This also involved adding sedimentary calcite (sed4) to the sediment module, with
!                 dissolution dynamics following ERSEM/benthic_calcite.
!                 A consequence of including detrital calcite is that carbonate chemistry and Omega(calcite)
!                 now needs to be computed for the whole water column, not just the surface for air-sea CO2 flux.
!                 We therefore added back the do subroutine to the module niva_ecosmo_carbonate that
!                 should be used with this module, and renamed the old code to niva_ecosmo_carbonate_nocalc.
!
!
!References (short form)
!
! Aksnes et al. (1997), doi:10.1080/00364827.1997.10413647
! Aumont et al. (2015) (A15), doi:10.5194/gmd-8-2465-2015
! Bianchi et al. (2012), doi:10.1029/2011GB004209
! Blackford et al. (2004), doi:10.1016/j.jmarsys.2004.02.004
! Burris (1981), doi.org/10.1007/bf00397114
! Butenschon et al. (2016) (B16), doi:10.5194/gmd-9-1293-2016
! Daewel and Schrum (2013) (D13), doi:10.1016/j.jmarsys.2013.03.008
! Feng et al. (2015) (F15), doi:10.1002/2015JG002931
! Geider et al. (1997), doi:10.3354/meps148187
! Gilstad and Sakshaug (1990), Marine Ecology Progress Series 64: 169-173.
! Laws (1991), doi:10.1016/0198-0149(91)90059-O
! Middelburg (2019), doi:10.1007/978-3-030-10822-9
! Moreno et al. (2022), doi:10.1029/2022AV000679
! Møller (2007), doi:10.4319/lo.2007.52.1.0079
! Neumann (2000), Journal of Marine Systems 25 (2000) 405–419
! Paulmier et al. (2009), doi:10.5194/bg-6-923-2009
! Platt et al. (1980), Journal of Marine Research 38, 687–701.
! Ridgwell et al. (2002), doi:10.1029/2002GB001877
! Sanz-Martin et al. (2019), doi:10.3389/fmars.2019.00468
! Steinberg and Landry (2017), doi:10.1146/annurev-marine-010814-015924
! Stephens et al. (2025), doi.org/10.1038/s42003-025-07574-2
! Stigebrandt and Wulff (1987), doi:10.1357/002224087788326812
! Tanioka and Matsumoto (2020), doi:10.1029/2019GL085564
! Yool et al. (2013), doi:10.5194/gmd-6-1767-2013
! Yumruktepe et al. (2022) (Y22), doi:10.5194/gmd-15-3901-2022
! -----------------------------------------------------------------------------------------------------------

   module niva_ecosmo
!
! !DESCRIPTION:
!
! The ECOSMO model is based on Daewel & Schrum (JMS,2013)
!
! !USES:
   use fabm_types
   use fabm_expressions
   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_niva_ecosmo
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: sedy0 = 86400.0_rk
   real(rk), parameter :: mmolm3_in_mll = 44.6608009_rk
   real(rk)            :: redf(20)=0.0_rk
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_niva_ecosmo
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_opa, id_cal, id_det, id_dia, id_fla
      type (type_state_variable_id)         :: id_diachl, id_flachl
      type (type_state_variable_id)         :: id_mesozoo, id_microzoo, id_dom, id_oxy
      type (type_state_variable_id)         :: id_dic, id_alk
      type (type_bottom_state_variable_id)  :: id_sed1, id_sed2, id_sed3, id_sed4
      type (type_dependency_id)             :: id_temp, id_salt, id_E0, id_om_cal
      type (type_horizontal_dependency_id)  :: id_tbs, id_wnd, id_aice, id_day_length
      type (type_diagnostic_variable_id)    :: id_denit, id_primprod, id_secprod
      type (type_diagnostic_variable_id)    :: id_c2chl_fla, id_c2chl_dia
      type (type_diagnostic_variable_id)    :: id_NlimPs, id_NlimPl, id_PlimPs, id_PlimPl, id_SilimPl
      type (type_diagnostic_variable_id)    :: id_NutlimPs, id_NutlimPl
      type (type_horizontal_diagnostic_variable_id)    :: id_tbsout

!     Model parameters
      real(rk) :: Exphy, rtsom_cnp, rtsim_s, rtsim_c_mesozoo
      real(rk) :: rNH4Ps, rNO3Ps, psiPs, rPO4Ps
      real(rk) :: muPs, q10Ps, gammaDPs, alfaPs, betaPs, MINchl2cPs, MAXchl2cPs
      real(rk) :: amortPs, m1Ps, m2Ps, frmort1Ps, frmort2Ps
      real(rk) :: rNH4Pl, rNO3Pl, psiPl, rPO4Pl, rSiPl, minSiPl
      real(rk) :: muPl, q10Pl, gammaDPl, alfaPl, betaPl, MINchl2cPl, MAXchl2cPl
      real(rk) :: amortPl, m1Pl, m2Pl, frmort1Pl, frmort2Pl
      real(rk) :: GrZsPs, GrZsPl, GrZsZs, GrZsD
      real(rk) :: prefZsPs, prefZsPl, prefZsZs, prefZsD, RgZs, q10Zs
      real(rk) :: gammaZsP, gammaZsD, fregesZs, excZs, fAexcZs, fexcdomZs
      real(rk) :: mZs, q10mZs, frmortZs
      real(rk) :: GrZlPs, GrZlPl, GrZlZs, GrZlZl, GrZlD
      real(rk) :: prefZlPs, prefZlPl, prefZlZs, prefZlZl, prefZlD, RgZl, q10Zl
      real(rk) :: fGslpZlPs, fGslpZlPl, fGslpZlZs, fGslpZlZl, fGslpZlD
      real(rk) :: gammaZlP, gammaZlD, fregesZl, excZl, fAexcZl, fexcdomZl
      real(rk) :: mZl, mvZl, KE0mvZl, q10mZl, frmortZl
      real(rk) :: eta_PQ, eta_O2resp, eta_denit, eta_SO4resp, eta_nitrif
      real(rk) :: nitrifmax, aTnitrif, KO2nitrif, reminD, dissO
      real(rk) :: Rain0, Kcalom, fdissCZs, fdissCZl, dissCmax, ndissC
      real(rk) :: sinkD, sinkOPAL, sinkCAL
      real(rk) :: crBotStr, resuspRt, burialRt, burialRt2, reminSED, aTreminSED
      real(rk) :: dissSEDO, burialRtO, burialRt2O
      real(rk) :: dissSEDCmax, dissSEDCmin, ndissSEDC, burialRtC, burialRt2C
      real(rk) :: releaseP, aTreleaseP, RelSEDp1, RelSEDp2

!     ECOSMO modules
      logical  :: couple_co2 ! activates PML (Blackford, 2004) carbon module
      logical  :: use_aice ! activates use of ice area to limit surface fluxes

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom

   end type type_niva_ecosmo
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ECOSMO model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_niva_ecosmo),intent(inout),target  :: self
   integer,                 intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!  Caglar Yumruktepe:
!     Added fabm.yaml support: parameters from yaml file are copied to
!                           BioC array. Eventually, BioC array will be removed
!                           from the model where parameter names from the yaml
!                           file will be used.
!     Added dynamic chlorophyll-a from Geider etal., 1997
!     Added community dependent particle sinking rates
!     Added chlorophyll-a dependent light-limitation   
!
! !LOCAL VARIABLES:
! Everything else taken from yaml file
!
   integer :: i
   ! set Redfield ratios:
   redf(1) = 6.625_rk      !C_N
   redf(2) = 106.0_rk      !C_P
   redf(3) = 6.625_rk      !C_SiO
   redf(4) = 16.0_rk       !N_P
   redf(5) = 1.0_rk        !N_SiO
   redf(6) = 12.01_rk      !C_Cmg
   redf(7) = 44.6608009_rk !O2mm_ml
   redf(8) = 14.007_rk     !N_Nmg
   redf(9) = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk     !Si_Simg
   do i=1,10
     redf(i+10) = 1._rk/redf(i)
   end do


!EOP
!-----------------------------------------------------------------------
!BOC

   !! Input parameters, changing units 1/day to 1/sec and mmolN,P,Si to mgC

   ! aggregation scale parameters
   call self%get_parameter( self%Exphy,    'Exphy',      'm2/mgCHL',   'phyto extinction',                default=0.04_rk )
   call self%get_parameter( self%rtsom_cnp,'rtsom_cnp',  '-',          'ratio of total suspended organic mass (dry) to total CNP mass', default=1.5547_rk)
   call self%get_parameter( self%rtsim_s,  'rtsim_s',    '-',          'ratio of total suspended inorganic mass (dry) to total silicon mass', default=2.4596_rk)
   call self%get_parameter( self%rtsim_c_mesozoo,'rtsim_c_mesozoo','-','ratio of total suspended inorganic mass (dry) to total carbon mass in mesozooplankton', default=0.1583_rk)

   ! small phytoplankton parameters
   call self%get_parameter( self%rNH4Ps,  'rNH4Ps',      'mmolN/m3',   'NH4 half saturation for Ps',      default=0.20_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%rNO3Ps,  'rNO3Ps',      'mmolN/m3',   'NO3 half saturation for Ps',      default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%psiPs,   'psiPs',       'm3/mmolN',   'NH4 inhibition for Ps',           default=3.0_rk,   scale_factor=1.0_rk/(redf(1)*redf(6)))
   call self%get_parameter( self%rPO4Ps,  'rPO4Ps',      'mmolP/m3',   'PO4 half saturation for Ps',      default=0.05_rk,  scale_factor=redf(2)*redf(6))
   call self%get_parameter( self%muPs,    'muPs',        '1/day',      'max growth rate for Ps',          default=1.10_rk,  scale_factor=1.0_rk/sedy0) !Note: realized maximum growth rate at (T=10degC, D=0.5) mumax = muP * (alfaP/(alfaP+betaP)) * (betaP/(alfaP+betaP))^(betaP/alfaP)
   call self%get_parameter( self%q10Ps,   'q10Ps',       '-',          'Q_10 temperature coefficient for Ps', default=1.0_rk)
   call self%get_parameter( self%gammaDPs,'gammaDPs',    '-',          'day length adaptation coefficient for Ps', default=1.e8_rk)
   call self%get_parameter( self%alfaPs,  'alfaPs', 'mgC/(mgChl day W/m2)', 'initial slope P-I curve for Ps', default=3.127_rk, scale_factor=1.0_rk/sedy0) !From original default 0.0393*6.625*12.01
   call self%get_parameter( self%betaPs,  'betaPs', 'mgC/(mgChl day W/m2)', 'photoinhibition parameter for Ps', default=0.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%MINchl2cPs,'MINchl2cPs','mgChl/mgC',  'minimum Chl to C ratio Ps',       default=0.0063_rk) !From original default 0.5/(6.625*12.01)
   call self%get_parameter( self%MAXchl2cPs,'MAXchl2cPs','mgChl/mgC',  'maximum Chl to C ratio Ps',       default=0.0481_rk) !From original default 3.83/(6.625*12.01)
   call self%get_parameter( self%amortPs, 'amortPs',     '-',          'ratio max:min Ps mortality under nutrient limitation', default=1.0_rk)
   call self%get_parameter( self%m1Ps,    'm1Ps',        '1/day',      'linear Ps mortality rate',        default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%m2Ps,    'm2Ps',       'm3/(mgC day)','quadratic Ps mortality rate',     default=0.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%frmort1Ps,'frmort1Ps',  '-',          'dissolved fraction of linear mortality',    default=0.4_rk)
   call self%get_parameter( self%frmort2Ps,'frmort2Ps',  '-',          'dissolved fraction of quadratic mortality', default=0.4_rk)

   ! large phytoplankton parameters
   call self%get_parameter( self%rNH4Pl,  'rNH4Pl',      'mmolN/m3',   'NH4 half saturation for Pl',      default=0.20_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%rNO3Pl,  'rNO3Pl',      'mmolN/m3',   'NO3 half saturation for Pl',      default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%psiPl,   'psiPl',       'm3/mmolN',   'NH4 inhibition for Pl',           default=3.0_rk,   scale_factor=1.0_rk/(redf(1)*redf(6)))
   call self%get_parameter( self%rPO4Pl,  'rPO4Pl',      'mmolP/m3',   'PO4 half saturation for Pl',      default=0.05_rk,  scale_factor=redf(2)*redf(6))
   call self%get_parameter( self%rSiPl,   'rSiPl',       'mmolSi/m3',  'SiO2 half saturation for Pl',     default=0.50_rk,  scale_factor=redf(3)*redf(6))
   call self%get_parameter( self%minSiPl, 'minSiPl',     'mmolSi/m3',  'stop Si uptake below this concentration', default=1.0_rk, scale_factor=redf(3)*redf(6))
   call self%get_parameter( self%muPl,    'muPl',        '1/day',      'max growth rate for Pl',          default=1.30_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%q10Pl,   'q10Pl',       '-',          'Q_10 temperature coefficient for Pl', default=1.0_rk)
   call self%get_parameter( self%gammaDPl,'gammaDPl',    '-',          'day length adaptation coefficient for Pl', default=0.5_rk)
   call self%get_parameter( self%alfaPl,  'alfaPl', 'mgC/(mgChl day W/m2)', 'initial slope P-I curve for Pl', default=4.225_rk, scale_factor=1.0_rk/sedy0) !From original default 0.0531*6.625*12.01
   call self%get_parameter( self%betaPl,  'betaPl', 'mgC/(mgChl day W/m2)', 'photoinhibition parameter for Pl', default=0.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%MINchl2cPl,'MINchl2cPl','mgChl/mgC',  'minimum Chl to C ratio Pl',       default=0.0063_rk) !From original default 0.5/(6.625*12.01)
   call self%get_parameter( self%MAXchl2cPl,'MAXchl2cPl','mgChl/mgC',  'maximum Chl to C ratio Pl',       default=0.0370_rk) !From original default 2.94/(6.625*12.01)
   call self%get_parameter( self%amortPl, 'amortPl',     '-',          'ratio max:min Pl mortality under nutrient limitation', default=1.0_rk)
   call self%get_parameter( self%m1Pl,    'm1Pl',        '1/day',      'linear Pl mortality rate',        default=0.04_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%m2Pl,    'm2Pl',       'm3/(mgC day)','quadratic Pl mortality rate',     default=0.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%frmort1Pl,'frmort1Pl',  '-',          'dissolved fraction of linear mortality',    default=0.4_rk)
   call self%get_parameter( self%frmort2Pl,'frmort2Pl',  '-',          'dissolved fraction of quadratic mortality', default=0.0_rk)

   ! small zooplankton parameters
   call self%get_parameter( self%GrZsPs,  'GrZsPs',      '1/day',      'Grazing rate Zs on Ps',           default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZsPl,  'GrZsPl',      '1/day',      'Grazing rate Zs on Pl',           default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZsZs,  'GrZsZs',      '1/day',      'Grazing rate Zs on Zs',           default=0.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZsD,   'GrZsD',       '1/day',      'Grazing rate Zs on Det',          default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%prefZsPs,'prefZsPs',    '-',          'Grazing preference Zs on Ps',     default=0.70_rk)
   call self%get_parameter( self%prefZsPl,'prefZsPl',    '-',          'Grazing preference Zs on Pl',     default=0.25_rk)
   call self%get_parameter( self%prefZsZs,'prefZsZs',    '-',          'Grazing preference Zs on Zs',     default=0.00_rk)
   call self%get_parameter( self%prefZsD, 'prefZsD',     '-',          'Grazing preference Zs on Det.',   default=0.00_rk)
   call self%get_parameter( self%RgZs,    'RgZs',        'mgC/m3',     'Zs grazing half saturation',      default=40.0_rk) !From original default 0.5*6.625*12.01
   call self%get_parameter( self%q10Zs,   'q10Zs',       '-',          'Q_10 temperature coefficient for Zs grazing', default=1.0_rk)
   call self%get_parameter( self%gammaZsP,'gammaZsP',    '-',          'Zs assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%gammaZsD,'gammaZsD',    '-',          'Zs assim. eff. on det',           default=0.75_rk)
   call self%get_parameter( self%fregesZs,'fregesZs',    '-',          'dissolved fraction of Zs egestion', default=0.4_rk)
   call self%get_parameter( self%excZs,   'excZs',       '1/day',      'Zs basal excretion/respiration rate', default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%fAexcZs, 'fAexcZs',     '-',          'Zs fraction of assimilated food that is excreted/respired', default=0.0_rk)
   call self%get_parameter( self%fexcdomZs,'fexcdomZs',  '-',          'Zs fraction of excreted/respired matter that is dom', default=0.0_rk)
   call self%get_parameter( self%mZs,     'mZs',         '1/day',      'Zs mortality rate',               default=0.20_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%q10mZs,  'q10mZs',      '-',          'Q_10 temperature coefficient for Zs mortality', default=1.0_rk)
   call self%get_parameter( self%frmortZs,'frmortZs',    '-',          'dissolved fraction of Zs mortality', default=0.4_rk)

   ! large zooplankton parameters
   call self%get_parameter( self%GrZlPs,  'GrZlPs',      '1/day',      'Grazing rate Zl on Ps',           default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZlPl,  'GrZlPl',      '1/day',      'Grazing rate Zl on Pl',           default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZlZs,  'GrZlZs',      '1/day',      'Grazing rate Zl on Zs',           default=0.50_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZlZl,  'GrZlZl',      '1/day',      'Grazing rate Zl on Zl',           default=0.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZlD,   'GrZlD',       '1/day',      'Grazing rate Zl on Det',          default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%prefZlPs,'prefZlPs',    '-',          'Grazing preference Zl on Ps',     default=0.10_rk)
   call self%get_parameter( self%prefZlPl,'prefZlPl',    '-',          'Grazing preference Zl on Pl',     default=0.85_rk)
   call self%get_parameter( self%prefZlZs,'prefZlZs',    '-',          'Grazing preference Zl on Zs',     default=0.15_rk)
   call self%get_parameter( self%prefZlZl,'prefZlZl',    '-',          'Grazing preference Zl on Zl',     default=0.00_rk)
   call self%get_parameter( self%prefZlD, 'prefZlD',     '-',          'Grazing preference Zl on Det.',   default=0.00_rk)
   call self%get_parameter( self%RgZl,    'RgZl',        'mgC/m3',     'Zl grazing half saturation',      default=40.0_rk) !From original default 0.5*6.625*12.01
   call self%get_parameter( self%q10Zl,   'q10Zl',       '-',          'Q_10 temperature coefficient for Zl grazing', default=1.0_rk)
   call self%get_parameter( self%fGslpZlPs,'fGslpZlPs',  '-',          'Zl fraction of Ps-grazing that is slopped (not ingested)', default=0.0_rk)
   call self%get_parameter( self%fGslpZlPl,'fGslpZlPl',  '-',          'Zl fraction of Pl-grazing that is slopped (not ingested)', default=0.0_rk)
   call self%get_parameter( self%fGslpZlZs,'fGslpZlZs',  '-',          'Zl fraction of Zs-feeding that is slopped (not ingested)', default=0.0_rk)
   call self%get_parameter( self%fGslpZlZl,'fGslpZlZl',  '-',          'Zl fraction of Zl-feeding that is slopped (not ingested)', default=0.0_rk)
   call self%get_parameter( self%fGslpZlD, 'fGslpZlD',   '-',          'Zl fraction of Det-feeding that is slopped (not ingested)', default=0.0_rk)
   call self%get_parameter( self%gammaZlP,'gammaZlP',    '-',          'Zl assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%gammaZlD,'gammaZlD',    '-',          'Zl assim. eff. on det',           default=0.75_rk)
   call self%get_parameter( self%fregesZl,'fregesZl',    '-',          'dissolved fraction of Zl egestion', default=0.4_rk)
   call self%get_parameter( self%excZl,   'excZl',       '1/day',      'Zl basal excretion/respiration rate', default=0.06_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%fAexcZl, 'fAexcZl',     '-',          'Zl fraction of assimilated food that is excreted/respired', default=0.0_rk)
   call self%get_parameter( self%fexcdomZl,'fexcdomZl',  '-',          'Zl fraction of excreted/respired matter that is dom', default=0.0_rk)
   call self%get_parameter( self%mZl,     'mZl',         '1/day',      'Zl basal mortality rate',         default=0.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%mvZl,    'mvZl',        '1/day',      'Zl mortality due to visual predation', default=0.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%KE0mvZl, 'KE0mvZl',     'W/m2',       'Half-saturation PAR for Zl visual predation mortality', default=2.2_rk)
   call self%get_parameter( self%q10mZl,  'q10mZl',      '-',          'Q_10 temperature coefficient for Zl mortality', default=1.0_rk)
   call self%get_parameter( self%frmortZl,'frmortZl',    '-',          'dissolved fraction of Zl mortality', default=0.4_rk)

   ! stoichiometry parameters
   call self%get_parameter( self%eta_PQ,    'eta_PQ',    '-',          'photosynthetic quotient: moles O2 produced per mole C fixed', default=1.0_rk)
   call self%get_parameter( self%eta_O2resp,'eta_O2resp','-',          'O2 respiration quotient: moles O2 consumed per mole C respired', default=1.0_rk)
   call self%get_parameter( self%eta_denit, 'eta_denit', '-',          'denitrification quotient: moles NO3 used per mole N in organic matter', default=5.3_rk) !Default from 0.8*6.625, see Paulmier et al. (2009).
   call self%get_parameter( self%eta_SO4resp,'eta_SO4resp','-',        'SO4 respiration quotient: moles negative O2 produced per mole C respired', default=1.0_rk) !Default follows Neumann (2000) and standard assumption [H2S] = -0.5*[O2].
   call self%get_parameter( self%eta_nitrif,'eta_nitrif','-',          'nitrification quotient: moles O2 consumed per mole NH4 nitrified', default=2.0_rk) !Default from Middelburg (2019): NH4 + 2O2 -> NO3 + 2H + H2O

   ! remineralization/regeneration parameters
   call self%get_parameter( self%nitrifmax,'nitrifmax',  '1/day',      'max nitrification rate',          default=0.1_rk,   scale_factor=1.0_rk/sedy0) !Default follows Stigebrandt and Wulff (1987).
   call self%get_parameter( self%aTnitrif,'aTnitrif',    '1/degC',     'temp. control nitrification',     default=0.11_rk) !Default follows Stigebrandt and Wulff (1987).
   call self%get_parameter( self%KO2nitrif,'KO2nitrif',  'mmolO2/m3',  'O2 half saturation for nitrification', default=0.45_rk) !Default from 0.01*44.6608009, with 0.01 from Stigebrandt and Wulff (1987).
   call self%get_parameter( self%reminD,  'reminD',      '1/day',      'Detritus remin. rate',            default=0.003_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%dissO,   'dissO',       '1/day',      'Opal dissolution rate',           default=0.015_rk, scale_factor=1.0_rk/sedy0)

   ! calcite formation/dissolution parameters
   call self%get_parameter( self%Rain0,   'Rain0',       '-',          'maximum Rain Ratio (PIC:POC) within calcifiers', default=1.0_rk)
   call self%get_parameter( self%Kcalom,  'Kcalom',      '-',          'half-saturation constant for calcifier Rain Ratio dependence on calcite saturation state', default=1.0_rk)
   call self%get_parameter( self%fdissCZs,'fdissCZs',    '-',          'fraction of small zooplankton prey calcite that dissolves after ingestion', default=0.50_rk)
   call self%get_parameter( self%fdissCZl,'fdissCZl',    '-',          'fraction of large zooplankton prey calcite that dissolves after ingestion', default=0.25_rk)
   call self%get_parameter( self%dissCmax,'dissCmax',    '1/day',      'maximum specific dissolution rate', default=0.03_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%ndissC,  'ndissC',      '-',          'power of the dissolution law (Keir 1980)', default=2.22_rk)

   ! sinking parameters
   call self%get_parameter( self%sinkD,   'sinkD',       'm/day',      'Detritus sinking rate',           default=5.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkOPAL,'sinkOPAL',    'm/day',      'OPAL sinking rate',               default=5.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%sinkCAL, 'sinkCAL',     'm/day',      'Calcite sinking rate',            default=5.0_rk,   scale_factor=1.0_rk/sedy0)

   ! sediment module parameters
   call self%get_parameter( self%crBotStr,'crBotStr',    'N/m2',       'critic. bot. stress for resusp.', default=0.007_rk)
   call self%get_parameter( self%resuspRt,'resuspRt',    '1/day',      'resuspension rate',               default=25.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%burialRt,'burialRt',    '1/day',      'linear detritus burial rate',     default=1e-5_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%burialRt2,'burialRt2',  'm2/(mgC day)','quadratic detritus burial parameter', default=0.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%reminSED,'reminSED',    '1/day',      'remineralization rate',           default=0.002_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%aTreminSED,'aTreminSED','1/degC',     'temp. control remineralization',  default=0.15_rk)
   call self%get_parameter( self%dissSEDO,'dissSEDO',    '1/day',      'sed. opal dissolution rate',      default=0.0002_rk,scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%burialRtO,'burialRtO',  '1/day',      'linear opal burial rate',         default=1e-5_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%burialRt2O,'burialRt2O','m2/(mgC day)','quadratic opal burial parameter', default=0.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%dissSEDCmax,'dissSEDCmax','1/day',    'sed. max. dissolution rate calcite', default=30.0_rk,scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%dissSEDCmin,'dissSEDCmin','1/day',    'sed. min. dissolution rate calcite', default=0.05_rk,scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%ndissSEDC,'ndissSEDC',  '-',          'sed. power of calcite dissolution law (Keir 1980)', default=2.22_rk)
   call self%get_parameter( self%burialRtC,'burialRtC',  '1/day',      'linear calcite burial rate',      default=1e-5_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%burialRt2C,'burialRt2C','m2/(mgC day)','quadratic calcite burial parameter', default=0.0_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%releaseP,'releaseP',    '1/day',      'P sed. release rate',             default=0.002_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%aTreleaseP,'aTreleaseP','1/degC',     'temp. control P sed. release',    default=0.15_rk)
   call self%get_parameter( self%RelSEDp1,'RelSEDp1',    '-',          'P sed. release O2-dependence max. fractional reduction', default=0.15_rk)
   call self%get_parameter( self%RelSEDp2,'RelSEDp2',    '-',          'P sed. release O2-dependence half-saturation coefficient', default=0.10_rk)

   ! switches
   call self%get_parameter( self%couple_co2,'couple_co2', '',          'switch coupling to carbonate module', default=.false.)
   call self%get_parameter( self%use_aice,'use_aice',     '',          'use ice area to limit air-sea flux',default=.false.)

   ! Register state variables
   call self%register_state_variable( self%id_no3,      'no3',     'mgC/m3',    'nitrate',                   minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_nh4,      'nh4',     'mgC/m3',    'ammonium',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.1_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_pho,      'pho',     'mgC/m3',    'phosphate',                 minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.3_rk*redf(2)*redf(6)  )
   call self%register_state_variable( self%id_sil,      'sil',     'mgC/m3',    'silicate',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(3)*redf(6)  )
   call self%register_state_variable( self%id_oxy,      'oxy',     'mmolO2/m3', 'oxygen',                    minimum=0.0_rk,   vertical_movement=0.0_rk,  &
                                      initial_value=85.0_rk  )
   call self%register_state_variable( self%id_fla,      'fla',     'mgC/m3',    'small phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6))
   call self%register_state_variable( self%id_dia,      'dia',     'mgC/m3',    'large phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_flachl,   'flachl',  'mgChl/m3',  'small phytoplankton chl-a', minimum=1.0e-7_rk/20., vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/20.)
   call self%register_state_variable( self%id_diachl,   'diachl',  'mgChl/m3',  'large phytoplankton chl-a', minimum=1.0e-7_rk/27., vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/27.)
   call self%register_state_variable( self%id_microzoo, 'microzoo','mgC/m3',    'microzooplankton',          minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_mesozoo,  'mesozoo', 'mgC/m3',    'mesozooplankton',           minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_det,      'det',     'mgC/m3',    'detritus',                  minimum=0.0_rk, vertical_movement=-self%sinkD,   &
                                      initial_value=2.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_opa,      'opa',     'mgC/m3',    'opal',                      minimum=0.0_rk, vertical_movement=-self%sinkOPAL,&
                                      initial_value=2.0_rk*redf(3)*redf(6) )
   call self%register_state_variable( self%id_cal,      'cal',     'mgC/m3',    'calcite',                   minimum=0.0_rk, vertical_movement=-self%sinkCAL, &
                                      initial_value=0.05_rk ) !Initial value from ERSEM
   call self%register_state_variable( self%id_dom,      'dom',     'mgC/m3',    'labile dissolved om',       minimum=0.0_rk , &
                                      initial_value=3.0_rk*redf(1)*redf(6)   )
   call self%register_state_variable( self%id_sed1,     'sed1',    'mgC/m2',    'sediment detritus',         minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(1)*redf(6)*redf(18) )
   call self%register_state_variable( self%id_sed2,     'sed2',    'mgC/m2',    'sediment opal',             minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(3)*redf(6)*redf(20) )
   call self%register_state_variable( self%id_sed3,     'sed3',    'mgC/m2',    'sediment adsorbed pho.',    minimum=0.0_rk , &
                                      initial_value=2.0_rk*redf(2)*redf(6)*redf(19) )
   call self%register_state_variable( self%id_sed4,     'sed4',    'mgC/m2',    'sediment calcite',          minimum=0.0_rk , &
                                      initial_value=0.05_rk ) !Initial value from ERSEM

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m**3/s', &
         'primary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_secprod,'secprod','mgC/m**3/s', &
         'secondary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NlimPs,'NlimPs','-', &
         'N-limitation of small phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NlimPl,'NlimPl','-', &
         'N-limitation of large phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_PlimPs,'PlimPs','-', &
         'P-limitation of small phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_PlimPl,'PlimPl','-', &
         'P-limitation of large phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_SilimPl,'SilimPl','-', &
         'Si-limitation of large phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_nutlimPs,'NutlimPs','-', &
         'nutrient limitation of small phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_nutlimPl,'NutlimPl','-', &
         'nutrient limitation of large phytoplankton production', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmolN/m**3/s', &
         'denitrification rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tbsout,'botstrss','fill_later', &
         'total bottom stress', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c2chl_fla,'c2chl_fla','mgC/mgCHL', &
         'daily-mean C to CHL ratio for flagellates', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_c2chl_dia,'c2chl_dia','mgC/mgCHL', &
         'daily-mean C to CHL ratio for diatoms', output=output_time_step_averaged)

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_E0,type_bulk_standard_variable(name='E0_s')) !scalar PAR energy flux [W/m2]
   call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
   call self%register_dependency(self%id_wnd,standard_variables%wind_speed)
   if (self%use_aice) call self%register_dependency(self%id_aice,standard_variables%ice_area_fraction)
   call self%register_horizontal_dependency(self%id_day_length,type_horizontal_standard_variable(name='day_length_s'))
   if (self%couple_co2) then
     call self%register_state_dependency(self%id_dic, 'dic_target','mmol m-3','dic budget')
     call self%register_state_dependency(self%id_alk, 'alk_target','mmol m-3','alkalinity budget')
     call self%register_dependency(self%id_om_cal,'om_cal','-','calcite saturation')
   end if

   ! Register contributions to aggregate variables
   ! light attenuation due to chlorophyll:
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_flachl,scale_factor=self%Exphy,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_diachl,scale_factor=self%Exphy,include_background=.true.)
   ! total chlorophyll:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='Chl_s',units='mg/m^3',aggregate_variable=.true.), &
         self%id_flachl,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='Chl_s',units='mg/m^3',aggregate_variable=.true.), &
         self%id_diachl,include_background=.true.)
   ! total suspended organic matter (dry weight):
   !fla:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_fla,scale_factor=1.e-3_rk*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_fla,scale_factor=redf(11)*redf(16)*1.e-3_rk*redf(8)*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_fla,scale_factor=redf(12)*redf(16)*1.e-3_rk*redf(9)*self%rtsom_cnp,include_background=.true.)
   !dia:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_dia,scale_factor=1.e-3_rk*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_dia,scale_factor=redf(11)*redf(16)*1.e-3_rk*redf(8)*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_dia,scale_factor=redf(12)*redf(16)*1.e-3_rk*redf(9)*self%rtsom_cnp,include_background=.true.)
   !microzoo:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_microzoo,scale_factor=1.e-3_rk*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_microzoo,scale_factor=redf(11)*redf(16)*1.e-3_rk*redf(8)*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_microzoo,scale_factor=redf(12)*redf(16)*1.e-3_rk*redf(9)*self%rtsom_cnp,include_background=.true.)
   !mesozoo:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_mesozoo,scale_factor=1.e-3_rk*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_mesozoo,scale_factor=redf(11)*redf(16)*1.e-3_rk*redf(8)*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_mesozoo,scale_factor=redf(12)*redf(16)*1.e-3_rk*redf(9)*self%rtsom_cnp,include_background=.true.)
   !det:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_det,scale_factor=1.e-3_rk*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_det,scale_factor=redf(11)*redf(16)*1.e-3_rk*redf(8)*self%rtsom_cnp,include_background=.true.)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_det,scale_factor=redf(12)*redf(16)*1.e-3_rk*redf(9)*self%rtsom_cnp,include_background=.true.)
   ! total suspended inorganic matter (dry weight excluding sediments):
   !diatom frustules:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSIM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_dia,scale_factor=redf(13)*redf(16)*1.e-3_rk*redf(10)*self%rtsim_s,include_background=.true.)
   !opal:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSIM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_opa,scale_factor=redf(13)*redf(16)*1.e-3_rk*redf(10)*self%rtsim_s,include_background=.true.)
   !calcite:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSIM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_cal,scale_factor=redf(16)*1.e-3_rk*100.09_rk,include_background=.true.) !Here we assume a molecular weight or 100.09 g/mol for CaCO3
   !mesozooplankton mineral content:
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='TSIM_s',units='mg/L',aggregate_variable=.true.), &
         self%id_mesozoo,scale_factor=1.e-3_rk*self%rtsim_c_mesozoo,include_background=.true.)

   return

end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ECOSMO model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_niva_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: fla_loss=1.0_rk
   real(rk) :: dia_loss=1.0_rk
   real(rk) :: mic_loss=1.0_rk
   real(rk) :: mes_loss=1.0_rk
   real(rk) :: no3,nh4,pho,sil,t_silPl,oxy,fla,dia
   real(rk) :: flachl,diachl,chl2c_fla,chl2c_dia
   real(rk) :: microzoo,mesozoo,opa,cal,det,dom
   real(rk) :: temp,salt,D,E0
   real(rk) :: frem,remdet,fremDOM,remdom,remtot,denit,nitrif
   real(rk) :: up_no3Ps,up_nh4Ps,up_nPs,up_phoPs,NutlimPs
   real(rk) :: up_no3Pl,up_nh4Pl,up_nPl,up_phoPl,up_silPl,NutlimPl
   real(rk) :: muPsc,muPlc,ChlCPs,ChlCPl,ProdPs,ProdPl,Prod
   real(rk) :: NutfacmPs,NutfacmPl
   real(rk) :: mort1Ps,mort1Pl,mort2Ps,mort2Pl,mortPs,mortPl
   real(rk) :: Fs,etZs,ZsonPs,ZsonPl,ZsonZs,ZsonD
   real(rk) :: Fl,etZl,ZlonPs,ZlonPl,ZlonZs,ZlonZl,ZlonD
   real(rk) :: f_O2pos,f_O2resp,f_denit,f_SO4resp
   real(rk) :: egesZs,excrZs,etmZs,mortZs,prodZs
   real(rk) :: inZlPs,inZlPl,inZlZs,inZlZl,inZlPZ,inZlD
   real(rk) :: egesZl,excrZl,etmZl,mZlt,mortZl,prodZl
   real(rk) :: dxxt,dxxdom,slpdom,excrt,excrdom,excrmin,dissO
   real(rk) :: RainR,om_cal,t,calc,diss
   real(rk) :: rhs,rhs_nit,rhs_amm,rhs_pho,rhs_cal,rhs_oxy
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_HORIZONTAL_(self%id_day_length,D)
   _GET_(self%id_E0,E0) !scalar PAR energy flux [W/m2]
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_pho,pho)
   _GET_(self%id_sil,sil)
   _GET_(self%id_dia,dia)
   _GET_(self%id_fla,fla)
   _GET_(self%id_diachl,diachl)
   _GET_(self%id_flachl,flachl)
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_det,det)
   _GET_(self%id_dom,dom)
   _GET_(self%id_opa,opa)
   _GET_(self%id_cal,cal)
   _GET_(self%id_oxy,oxy)

   ! Retrieve current (local) diagnostic variable values.
   _GET_(self%id_om_cal,om_cal)

   ! CAGLAR
   ! checks - whether the biomass of plankton is below a predefined threshold,
   !          where below the threshold, loss terms are removed from the RHS of
   !          the equations. The idea is to keep plankton safe from extinction.
   ! loss terms are multiplied by the constants below, which can only be set
   ! by the model to 0 or 1.
   fla_loss = max(sign(-1.0_rk,fla-0.1_rk),0.0_rk)       ! flagellates
   dia_loss = max(sign(-1.0_rk,dia-0.1_rk),0.0_rk)       ! diatoms
   mic_loss = max(sign(-1.0_rk,microzoo-0.01_rk),0.0_rk) ! microzooplankton
   mes_loss = max(sign(-1.0_rk,mesozoo-0.01_rk),0.0_rk)  ! mesozooplankton

   ! nutrient limitation factors
   up_nh4Ps = nh4/(self%rNH4Ps+nh4)
   up_nh4Pl = nh4/(self%rNH4Pl+nh4)
   up_no3Ps = no3/(self%rNO3Ps+no3)*exp(-self%psiPs*nh4)
   up_no3Pl = no3/(self%rNO3Pl+no3)*exp(-self%psiPl*nh4)
   up_nPs   = up_nh4Ps+up_no3Ps
   up_nPl   = up_nh4Pl+up_no3Pl
   up_phoPs = pho/(self%rPO4Ps+pho)
   up_phoPl = pho/(self%rPO4Pl+pho)
   t_silPl  = max(sil-self%minSiPl,0.0_rk)
   up_silPl = t_silPl/(self%rSiPl+t_silPl)
   NutlimPs = min(up_nPs, up_phoPs)
   NutlimPl = min(up_nPl, up_phoPl, up_silPl)

   ! maximum phytoplankton growth rates, corrected for temperature and fractional daylength
   muPsc = self%q10Ps**((temp-10._rk)/10._rk) * (self%gammaDPs+0.5_rk)/(self%gammaDPs+D) * self%muPs
   muPlc = self%q10Pl**((temp-10._rk)/10._rk) * (self%gammaDPl+0.5_rk)/(self%gammaDPl+D) * self%muPl

   ! further correction for light limitation
   ChlCPs = flachl/max(fla,1e-8_rk) !to avoid 0/0
   ChlCPs = max(self%MINchl2cPs,ChlCPs)
   ChlCPs = min(self%MAXchl2cPs,ChlCPs)
   ChlCPl = diachl/max(dia,1e-8_rk) !to avoid 0/0
   ChlCPl = max(self%MINchl2cPl,ChlCPl)
   ChlCPl = min(self%MAXchl2cPl,ChlCPl)
   if (E0>1e-8_rk) then
     muPsc = muPsc * (1-exp(-self%alfaPs*E0*ChlCPs/muPsc)) * exp(-self%betaPs*E0*ChlCPs/muPsc)
     muPlc = muPlc * (1-exp(-self%alfaPl*E0*ChlCPl/muPlc)) * exp(-self%betaPl*E0*ChlCPl/muPlc)

     ! ratio of chlorophyll-a synthesis to carbon fixation rates
     chl2c_fla = self%MAXchl2cPs * (muPsc*fla)/(self%alfaPs*E0*flachl)
     chl2c_dia = self%MAXchl2cPl * (muPlc*dia)/(self%alfaPl*E0*diachl)
     chl2c_fla = max(self%MINchl2cPs,chl2c_fla)
     chl2c_fla = min(self%MAXchl2cPs,chl2c_fla)
     chl2c_dia = max(self%MINchl2cPl,chl2c_dia)
     chl2c_dia = min(self%MAXchl2cPl,chl2c_dia)
   else
     muPsc = 0.0_rk
     muPlc = 0.0_rk
     chl2c_fla = self%MAXchl2cPs
     chl2c_dia = self%MAXchl2cPl
     ! This result for chl2c follows from applying L'Hopital's Rule:
     ! Lim(E0->0) {(1-exp(alfa*E0*ChlC/mu))/E0} = alfa*ChlC/mu
   end if

   ! further correction for nutrient limitation
   muPsc = muPsc * NutlimPs
   muPlc = muPlc * NutlimPl

   ! primary productivity
   ProdPs = muPsc*fla
   ProdPl = muPlc*dia
   Prod = ProdPs + ProdPl

   ! phytoplankton mortality
   NutfacmPs = self%amortPs / (1._rk + NutlimPs*(self%amortPs-1._rk)) !enhancement factor due to nutrient limitation
   NutfacmPl = self%amortPl / (1._rk + NutlimPl*(self%amortPl-1._rk))
   mort1Ps = NutfacmPs * self%m1Ps*fla_loss * fla
   mort1Pl = NutfacmPl * self%m1Pl*dia_loss * dia
   mort2Ps = NutfacmPs * self%m2Ps*fla*fla_loss * fla
   mort2Pl = NutfacmPl * self%m2Pl*dia*dia_loss * dia
   mortPs  = mort1Ps + mort2Ps
   mortPl  = mort1Pl + mort2Pl

   ! grazing/feeding fluxes
   Fs = self%prefZsPs*fla + self%prefZsPl*dia + self%prefZsZs*microzoo + self%prefZsD*det
   Fl = self%prefZlPs*fla + self%prefZlPl*dia + self%prefZlZs*microzoo + self%prefZlZl*mesozoo &
       + self%prefZlD*det

   etZs = self%q10Zs**((temp-10._rk)/10._rk)
   etZl = self%q10Zl**((temp-10._rk)/10._rk)

   ZsonPs = fla_loss * self%GrZsPs * etZs * self%prefZsPs * fla/(self%RgZs + Fs) * microzoo
   ZsonPl = dia_loss * self%GrZsPl * etZs * self%prefZsPl * dia/(self%RgZs + Fs) * microzoo
   ZsonZs = mic_loss * self%GrZsZs * etZs * self%prefZsZs * microzoo/(self%RgZs + Fs) * microzoo
   ZsonD  =            self%GrZsD  * etZs * self%prefZsD  * det/(self%RgZs + Fs)  * microzoo

   ZlonPs = fla_loss * self%GrZlPs * etZl * self%prefZlPs * fla/(self%RgZl + Fl) * mesozoo
   ZlonPl = dia_loss * self%GrZlPl * etZl * self%prefZlPl * dia/(self%RgZl + Fl) * mesozoo
   ZlonZs = mic_loss * self%GrZlZs * etZl * self%prefZlZs * microzoo/(self%RgZl + Fl) * mesozoo
   ZlonZl = mes_loss * self%GrZlZl * etZl * self%prefZlZl * mesozoo/(self%RgZl + Fl) * mesozoo
   ZlonD =             self%GrZlD  * etZl * self%prefZlD  * det/(self%RgZl + Fl)  * mesozoo

   ! unitary switching functions
   f_O2pos   = 1.0_rk
   f_O2resp  = 1.0_rk
   f_denit   = 0.0_rk
   f_SO4resp = 0.0_rk
   if (oxy.le.0.0_rk) then
     f_O2pos   = 0.0_rk
     f_O2resp  = 0.0_rk
     if (no3.gt.0.0_rk) then
       f_denit   = 1.0_rk
     else
       f_SO4resp = 1.0_rk
     end if
   end if

! reaction rates

   ! phytoplankton
   _ADD_SOURCE_(self%id_fla, ProdPs - mortPs - ZsonPs - ZlonPs)
   _ADD_SOURCE_(self%id_dia, ProdPl - mortPl - ZsonPl - ZlonPl)

   ! phytoplankton chlorophyll-a
   _ADD_SOURCE_(self%id_flachl, ProdPs*chl2c_fla - (mortPs + ZsonPs + ZlonPs)*flachl/fla)
   _ADD_SOURCE_(self%id_diachl, ProdPl*chl2c_dia - (mortPl + ZsonPl + ZlonPl)*diachl/dia) 

   ! microzooplankton
   egesZs   =  (1._rk-self%gammaZsP)*(ZsonPs+ZsonPl+ZsonZs) &       !egestion -> DOM,POM
             + (1._rk-self%gammaZsD)*ZsonD
   excrZs   = self%excZs*etZs*mic_loss * microzoo &                 !basal excretion/respiration
             + self%fAexcZs*(self%gammaZsP*(ZsonPs+ZsonPl+ZsonZs) & !activity excretion/respiration as fraction of absorbed
             +               self%gammaZsD*ZsonD)                   ! -> DOM,DIC/NUT
   etmZs    = self%q10mZs**((temp-10._rk)/10._rk)
   mortZs   = self%mZs*etmZs*mic_loss * microzoo
   prodZs   = ZsonPs + ZsonPl + ZsonZs + ZsonD - egesZs             !contribution to 'secondary production'
   _ADD_SOURCE_(self%id_microzoo, prodZs - excrZs - mortZs - ZsonZs - ZlonZs)

   ! mesozooplankton
   inZlPs   = (1._rk-self%fGslpZlPs)*ZlonPs                         !ingestion, corrected for sloppy feeding -> DOM
   inZlPl   = (1._rk-self%fGslpZlPl)*ZlonPl
   inZlZs   = (1._rk-self%fGslpZlZs)*ZlonZs
   inZlZl   = (1._rk-self%fGslpZlZl)*ZlonZl
   inZlPZ   = inZlPs + inZlPl + inZlZs + inZlZl
   inZlD    = (1._rk-self%fGslpZlD) *ZlonD
   egesZl   =  (1._rk-self%gammaZlP)*inZlPZ &                       !egestion -> DOM,POM
             + (1._rk-self%gammaZlD)*inZlD
   excrZl   = self%excZl*etZl*mes_loss * mesozoo &                  !basal excretion/respiration
             + self%fAexcZl*(self%gammaZlP*inZlPZ &                 !activity excretion/respiration as fraction of absorbed
             +               self%gammaZlD*inZlD)                   ! -> DOM,DIC/NUT
   etmZl    = self%q10mZl**((temp-10._rk)/10._rk)
   mZlt     = self%mZl + self%mvZl*E0/(self%KE0mvZl+E0)             !total mortality rate inc. visual predation
   mortZl   = mZlt*etmZl*mes_loss * mesozoo
   prodZl   = inZlPZ + inZlD - egesZl                               !contribution to 'secondary production'
   _ADD_SOURCE_(self%id_mesozoo, prodZl - excrZl - mortZl - ZlonZl)

   ! detritus
   dxxt     =  mortPs + mortPl + mortZs + mortZl + egesZs + egesZl
   dxxdom   =  self%frmort1Ps*mort1Ps + self%frmort1Pl*mort1Pl &
             + self%frmort2Ps*mort2Ps + self%frmort2Pl*mort2Pl &
             + self%frmortZs*mortZs + self%frmortZl*mortZl &
             + self%fregesZs*egesZs + self%fregesZl*egesZl
   frem     = self%reminD * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))
   remdet   = frem * det
   _ADD_SOURCE_(self%id_det, (dxxt-dxxdom) - ZsonD - ZlonD - remdet)

   ! labile dissolved organic matter
   slpdom   = (self%fGslpZlPs*ZlonPs + self%fGslpZlPl*ZlonPl &
             + self%fGslpZlZs*ZlonZs + self%fGslpZlZl*ZlonZl + self%fGslpZlD*ZlonD)
   excrt    = excrZs + excrZl
   excrdom  = self%fexcdomZs*excrZs + self%fexcdomZl*excrZl
   excrmin  = excrt - excrdom
   fremDOM  = 10.0_rk * frem
   remdom   = fremDOM * dom
   _ADD_SOURCE_(self%id_dom, dxxdom + slpdom + excrdom - remdom)

   ! nitrate
   nitrif   = f_O2pos * self%nitrifmax * exp(self%aTnitrif*temp) * oxy/(self%KO2nitrif+oxy) * nh4 !nitrification
   remtot   = remdet + remdom
   denit    = f_denit * self%eta_denit * remtot
   rhs_nit  = -(up_no3Ps+0.5d-10)/(up_nPs+1.0d-10)*ProdPs &
              -(up_no3Pl+0.5d-10)/(up_nPl+1.0d-10)*ProdPl &
             + nitrif &
             - denit
   _ADD_SOURCE_(self%id_no3, rhs_nit)

   ! ammonium
   rhs_amm  = -(up_nh4Ps+0.5d-10)/(up_nPs+1.0d-10)*ProdPs &
              -(up_nh4Pl+0.5d-10)/(up_nPl+1.0d-10)*ProdPl &
             + excrmin & 
             + remtot &
             - nitrif
   _ADD_SOURCE_(self%id_nh4, rhs_amm)

   ! phosphate
   rhs_pho  = -Prod &
             + excrmin &
             + remtot
   _ADD_SOURCE_(self%id_pho, rhs_pho)

   ! silicate
   dissO    = self%dissO * opa
   _ADD_SOURCE_(self%id_sil, -ProdPl + dissO)

   ! opal
   _ADD_SOURCE_(self%id_opa, mortPl + ZsonPl + ZlonPl - dissO)

   ! calcite
   !First we calculate the calcifier 'rain ratio', i.e. the ratio PIC:POC
   !within coccolithophores as a function of Omega(calcite), based on experimental data
   !(Gehlen et al., 2007; Zondervan et al., 2002).
   RainR    = self%Rain0 * max(0._rk, (om_cal-1._rk)/(om_cal-1._rk+self%Kcalom))
   !Next we correct the rain ratio to represent a ratio of PIC:POC within the
   !'small phytoplankton' functional group as a whole. This makes use of semi-empirical
   !factors that reflect the dependency of the calcifying fraction of small phytoplankton 
   !on the environmental conditions. The approach here follows ERSEM (Butenschon et al., 2016).
   t        = max(0._rk,temp)  ! this is to avoid funny values of rain ratio when temp ~ -2 degrees
   RainR    = RainR * min((1._rk-up_phoPs),up_nPs) * (t/(2._rk+t))
   RainR    = max(RainR,0.005_rk)
   !Next we use the rain ratio to calculate fluxes to the detrital calcite pool
   !arising from particulate fractions of small phytoplankton mortality and grazer ingestion,
   !corrected for dissolution within zooplankton guts or vacuoles. Here we deviate from the
   !ERSEM B16 formulation in favour of the PISCES formulation (Aumont et al., 2015) which corrects
   !for internal dissolution but does not involve grazer inefficiency factors.
   calc     = ((1._rk-self%frmort1Ps)*mort1Ps + (1._rk-self%frmort2Ps)*mort2Ps) * RainR &
             + ZsonPs * RainR * (1._rk-self%fdissCZs) &
             + ZlonPs * RainR * (1._rk-self%fdissCZl)
   !Next we calculate dissolution flux following ERSEM B16, Kier (1980).
   diss     = self%dissCmax * max(1._rk-om_cal,0._rk)**self%ndissC * cal
   rhs_cal  = calc - diss
   _ADD_SOURCE_(self%id_cal, rhs_cal)

   ! oxygen
   rhs_oxy  = ((up_nh4Ps + 1.226_rk*up_no3Ps+1.d-10)/(up_nPs+1.d-10) * self%eta_PQ * ProdPs &
             + (up_nh4Pl + 1.226_rk*up_no3Pl+1.d-10)/(up_nPl+1.d-10) * self%eta_PQ * ProdPl &
             - f_O2resp * self%eta_O2resp * excrmin &
             - (f_O2resp*self%eta_O2resp + f_SO4resp*self%eta_SO4resp) * remtot &
             - redf(11)*self%eta_nitrif * nitrif) * redf(16)
   _ADD_SOURCE_(self%id_oxy, rhs_oxy)

   ! carbonate dynamics
   if (self%couple_co2) then
     _ADD_SOURCE_(self%id_dic, redf(16)*(excrmin + remtot - Prod - rhs_cal))

     rhs    = redf(16)*((rhs_amm-rhs_nit)*redf(11) &               !\Delta[TA] = + \Delta[NH4] - \Delta[NO3]
             - rhs_pho*redf(12) &                                  !             - \Delta[PO4]
             - 2.0_rk * rhs_cal) &                                 !             - 2*\Delta[CaCO3]
             - 0.5_rk * rhs_oxy * (1._rk-f_O2pos)                  !             + \Delta[H2S]   (O2<=0)
     ! The last term here accounts for alkalinity changes due to release/consumption
     ! of H2S, which correspond to -1/2 * oxygen changes on molar basis when oxygen is negative.
     ! \Delta[TA] = \Delta[H2S] = -1/2 * \Delta[O2]         (O2<0)
     _ADD_SOURCE_(self%id_alk, rhs)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_primprod, Prod)
   _SET_DIAGNOSTIC_(self%id_secprod, prodZs + prodZl)
   _SET_DIAGNOSTIC_(self%id_denit,denit*redf(11)*redf(16))
   _SET_DIAGNOSTIC_(self%id_NlimPs, up_nPs)
   _SET_DIAGNOSTIC_(self%id_NlimPl, up_nPl)
   _SET_DIAGNOSTIC_(self%id_PlimPs, up_phoPs)
   _SET_DIAGNOSTIC_(self%id_PlimPl, up_phoPl)
   _SET_DIAGNOSTIC_(self%id_SilimPl, up_silPl)
   _SET_DIAGNOSTIC_(self%id_NutlimPs, NutlimPs)
   _SET_DIAGNOSTIC_(self%id_NutlimPl, NutlimPl)
   _SET_DIAGNOSTIC_(self%id_c2chl_fla, 1.0_rk/chl2c_fla)
   _SET_DIAGNOSTIC_(self%id_c2chl_dia, 1.0_rk/chl2c_dia)


   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_niva_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   ! Code for surface O2 flux was copied from nersem_oxygen.F90 (PJW 02/08/2025).

   real(rk) :: O2o,ETW,T,X1X,wnd,aice
   real(rk) :: OSAT,sc,ko2o,FAIRO2

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_oxy,O2o) ! NOTE: Both NERSEM and ECOSMO codes use O2 in [mmolO2/m3]
   _GET_(self%id_temp,ETW)
   _GET_(self%id_salt,X1X)
   _GET_HORIZONTAL_(self%id_wnd,wnd)
   if (self%use_aice) then
     _GET_HORIZONTAL_(self%id_aice,aice)
   end if

   X1X = max(X1X, 0.0_rk)
   wnd = max(wnd, 0.0_rk)

   OSAT = oxygen_saturation_concentration(self,ETW,X1X)

   ! New formulation for the Schmidt number for O2 following Wanninkhof 2014
   T = max(min(ETW,40.0_rk), -2.0_rk)
   sc = 1920.4_rk - 135.6_rk*T + 5.2122_rk*T**2 - 0.10939_rk*T**3 + 0.00093777_rk*T**4

   ko2o = 0.251_rk * (wnd**2) * (sc/660._rk)**(-0.5_rk) ! Wanninkhof 2014

   ! units of ko2 converted from cm/hr to m/s
   ko2o = ko2o/360000._rk
   if (self%use_aice) then
     ko2o = max(0._rk, (1._rk-aice))*ko2o !Limit flux to area fraction not covered by ice
   end if

   FAIRO2 = ko2o*(OSAT-O2o)

   _ADD_SURFACE_FLUX_(self%id_oxy,FAIRO2)

   ! Leave spatial loops over the horizontal domain (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bottom fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_niva_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk) :: temp, oxy, no3, det, opa, cal, sed1, sed2, sed3, sed4
   real(rk) :: tbs, om_cal
   real(rk) :: Rds, RdsO, RdsC, Rsd
   real(rk) :: settle, resusp, remSED, burSED, denitSED
   real(rk) :: yt1, yt2, releaseP
   real(rk) :: rhs, flux, flux_oxy
   real(rk) :: f_O2pos, f_O2resp, f_denit, f_SO4resp
   real(rk) :: settleO, resuspO, dissSEDO, burSEDO
   real(rk) :: settleC, resuspC, rdissC, dissSEDC, burSEDC
   ! add community sinking local variables
   real(rk) :: dsnk
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_det,det)
   _GET_(self%id_opa,opa)
   _GET_(self%id_cal,cal)
   _GET_(self%id_no3,no3)
   _GET_HORIZONTAL_(self%id_sed1,sed1)
   _GET_HORIZONTAL_(self%id_sed2,sed2)
   _GET_HORIZONTAL_(self%id_sed3,sed3)
   _GET_HORIZONTAL_(self%id_sed4,sed4)

   _GET_(self%id_om_cal, om_cal)
   om_cal    = max(om_cal, 0._rk)
   _GET_HORIZONTAL_(self%id_tbs,tbs)


   ! unitary redox switching functions
   f_O2pos   = 1.0_rk
   f_O2resp  = 1.0_rk
   f_denit   = 0.0_rk
   f_SO4resp = 0.0_rk
   if (oxy.le.0.0_rk) then
     f_O2pos   = 1.0_rk
     f_O2resp  = 0.0_rk
     if (no3.gt.0.0_rk) then
       f_denit   = 1.0_rk
     else
       f_SO4resp = 1.0_rk
     end if
   end if

   !----suspension, settling, remineralization rates
   if (tbs.ge.self%crBotStr) then
     Rsd    = min(self%resuspRt, self%resuspRt * tbs**2 * 100._rk) !sets to max=self%resuspRt when tbs~=0.070, else rapid increase from 0.
     Rds    = 0.0_rk
     RdsO   = 0.0_rk
     RdsC   = 0.0_rk
   else
     Rsd    = 0.0_rk
     Rds    = self%sinkD
     RdsO   = self%sinkOPAL
     RdsC   = self%sinkCAL
   end if

   ! sediment organic matter (sed1)
   settle   = Rds * det
   resusp   = Rsd * sed1
   remSED   = self%reminSED * exp(self%aTreminSED*temp) * sed1
   burSED   = (self%burialRt + self%burialRt2*sed1) * sed1
   denitSED = f_denit * self%eta_denit * remSED
   _ADD_BOTTOM_FLUX_(self%id_det, resusp - settle)
   _ADD_BOTTOM_FLUX_(self%id_no3, -denitSED)
   _ADD_BOTTOM_FLUX_(self%id_nh4, remSED)
   _ADD_BOTTOM_SOURCE_(self%id_sed1, settle - resusp - remSED - burSED)

   ! sediment opal (sed2)
   settleO  = RdsO * opa
   resuspO  = Rsd  * sed2
   dissSEDO = self%dissSEDO * sed2
   burSEDO  = (self%burialRtO + self%burialRt2O*sed2) * sed2
   _ADD_BOTTOM_FLUX_(self%id_opa, resuspO - settleO)
   _ADD_BOTTOM_FLUX_(self%id_sil, dissSEDO)
   _ADD_BOTTOM_SOURCE_(self%id_sed2, settleO - resuspO - dissSEDO - burSEDO)

   ! sediment phosphate (sed3)
   yt2      = max(oxy/375.0_rk, 0.0_rk)   !normieren des wertes wie in Neumann et al 2002
   yt1      = yt2**2.0_rk/(self%RelSEDp2**2.0_rk+yt2**2.0_rk)
   releaseP = self%releaseP * exp(self%aTreleaseP*temp) * (1._rk-self%RelSEDp1*yt1) * sed3 !release of phosphate from sediments
   _ADD_BOTTOM_FLUX_(self%id_pho, releaseP)
   _ADD_BOTTOM_SOURCE_(self%id_sed3, remSED - releaseP)

   ! sediment calcite (sed4)
   settleC  = RdsC * cal
   resuspC  = Rsd  * sed4
   rdissC   = self%dissSEDCmax * (max(1._rk-om_cal,0._rk))**self%ndissSEDC
   rdissC   = max(rdissC, self%dissSEDCmin)
   dissSEDC = rdissC * sed4
   burSEDC  = (self%burialRtC + self%burialRt2C*sed4) * sed4
   _ADD_BOTTOM_FLUX_(self%id_cal, resuspC - settleC)
   _ADD_BOTTOM_SOURCE_(self%id_sed4, settleC - resuspC - dissSEDC - burSEDC)

   ! bottom oxygen
   flux_oxy = -(f_O2resp*self%eta_O2resp  + f_SO4resp*self%eta_SO4resp) * remSED * redf(16)
   _ADD_BOTTOM_FLUX_(self%id_oxy, flux_oxy)

   ! bottom carbonate dynamics
   if (self%couple_co2) then
     _ADD_BOTTOM_FLUX_(self%id_dic, redf(16)*(remSED + dissSEDC))

     flux   = redf(16)*((remSED+denitSED)*redf(11) &      !\Delta[TA] = + \Delta[NH4] - \Delta[NO3]
             - releaseP*redf(12) &                        !             - \Delta[PO4]
             + 2.0_rk * dissSEDC) &                       !             - 2*\Delta[CaCO3]
             - 0.5_rk * flux_oxy * (1._rk-f_O2pos)        !             + \Delta[H2S]   (O2<=0)
     ! The last term here accounts for alkalinity changes due to release/consumption
     ! of H2S, which correspond to -1/2 * oxygen changes on molar basis when oxygen is negative.
     ! \Delta[TA] = \Delta[H2S] = -1/2 * \Delta[O2]         (O2<0)
     _ADD_BOTTOM_FLUX_(self%id_alk, flux)
   end if

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tbsout, tbs)

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

   ! Code for OSAT function was copied from nersem_oxygen.F90 (PJW 02/08/2025).
   function oxygen_saturation_concentration(self,ETW,X1X) result(OSAT)
      class (type_niva_ecosmo), intent(in) :: self
      real(rk),                 intent(in) :: ETW,X1X
      real(rk)                             :: OSAT

      real(rk),parameter :: A1 = -173.4292_rk
      real(rk),parameter :: A2 = 249.6339_rk
      real(rk),parameter :: A3 = 143.3483_rk
      real(rk),parameter :: A4 = -21.8492_rk
      real(rk),parameter :: B1 = -0.033096_rk
      real(rk),parameter :: B2 = 0.014259_rk
      real(rk),parameter :: B3 = -0.0017_rk
      real(rk),parameter :: R = 8.3145_rk
      real(rk),parameter :: P = 101325_rk
      real(rk),parameter :: T = 273.15_rk

      ! volume of an ideal gas at standard temp (0C) and pressure (1 atm)
      real(rk),parameter :: VIDEAL = (R * 273.15_rk / P) *1000._rk

      real(rk)           :: ABT

      ! calc absolute temperature
      ABT = ETW + T

      ! calc theoretical oxygen saturation for temp + salinity
      ! From WEISS 1970 DEEP SEA RES 17, 721-735.
      ! units of ln(ml(STP)/l)
      OSAT = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
               + A4 * (ABT/100._rk) &
               + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

      ! convert units to ml(STP)/l then to mMol/m3
      OSAT = exp( OSAT )
      OSAT = OSAT * 1000._rk / VIDEAL
   end function

   end module niva_ecosmo