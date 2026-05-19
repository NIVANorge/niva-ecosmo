#include "fabm_driver.h"

! -----------------------------------------------------------------------------------------------------------
! This is a NIVA adaptation of ERSEM carbonate module for use with the ECOSMO model
! in applications where calcite formation/dissolution is not modelled.
! -----------------------------------------------------------------------------------------------------------
!
!
! Change Log ------------------------------------------------------------------------------------------------
! PJW 05/08/2025: Copied from niva/nersem/nersem_carbonate.F90, renamed to niva_ecosmo_carbonate.F90.
!                 Removed subroutine do (not currently needed since no pH etc. sensitivity yet in ECOSMO).
!                 Removed subroutine CaCO3_SATURATION (only called from do).
!                 Harmonized capitalization of CO2DYN.
!                 Updated _SET_SURFACE_EXCHANGE_ to _ADD_SURFACE_FLUX_ and removed diagnostic fair.
!                 (redundant since the latter automatically makes a diagnostic dic_sfl).
!                 Renamed (O3c,TA) -> (dic,alk) to match ECOSMO conventions.
!                 Renamed instances of ctot to Ctot.
!                 Removed unused: (bioalk,Carb_in,BiCarb_in,CarbA_in,pH*,Carb,BiCarb,CarbA,Om_*,alk_diag).
!
! PJW 07/05/2025: Renamed niva_ecosmo_carbonate->niva_ecosmo_carbonate_nocalc.
!                 Renamed PWA->PJW.
!                 This module is now reserved for use with niva_ecosmo_nocalc in cases where modelling
!                 calcite formation/dissolution is not considered necessary. It is significantly
!                 faster than niva_ecosmo_carbonate because the carbonate system only needs to
!                 be calculated for the surface layer (to calculate air-sea CO2 exchange).
! -----------------------------------------------------------------------------------------------------------


module niva_ecosmo_carbonate_nocalc

   use fabm_types
   use fabm_builtin_models

   implicit none

   private

   type,extends(type_base_model),public :: type_niva_ecosmo_carbonate_nocalc
!     Variable identifiers
      type (type_state_variable_id)     :: id_dic,id_alk
      type (type_dependency_id)         :: id_ETW, id_X1X, id_dens, id_pres
      type (type_dependency_id)         :: id_pco2_in
      type (type_horizontal_dependency_id) :: id_wnd,id_PCO2A,id_aice

      type (type_diagnostic_variable_id) :: id_pco2
      type (type_horizontal_diagnostic_variable_id) :: id_wnd_diag

      integer :: phscale
      logical :: use_aice
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

   public :: CO2DYN

contains

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_niva_ecosmo_carbonate_nocalc), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%phscale,'pHscale','','pH scale (1: total, 0: SWS, -1: SWS backward compatible)',default=1,minimum=-1,maximum=1)
      call self%get_parameter(self%use_aice,'use_aice','','use ice area to limit air-sea flux',default=.false.)

      call self%register_state_variable(self%id_dic,'dic','mmol C/m^3','total dissolved inorganic carbon', 2200._rk,minimum=0._rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_dic)

      ! Total alkalinity is a state variable.
      call self%register_state_variable(self%id_alk,'alk','mmol/m^3','total alkalinity',2300._rk,minimum=1.e-4_rk, &
         standard_variable=standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_diagnostic_variable(self%id_pco2,  'pCO2',  '1e-6',    'partial pressure of CO2',missing_value=0._rk)
      call self%register_diagnostic_variable(self%id_wnd_diag,'wind','m/s','surface wind speed',source=source_do_surface)

      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_X1X, standard_variables%practical_salinity)
      call self%register_dependency(self%id_dens,standard_variables%density)
      call self%register_dependency(self%id_pres,standard_variables%pressure)
      call self%register_dependency(self%id_pco2_in,'pCO2','1e-6','previous pCO2')

      call self%register_dependency(self%id_wnd,  standard_variables%wind_speed)
      call self%register_dependency(self%id_PCO2A,standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      if (self%use_aice) call self%register_dependency(self%id_aice,standard_variables%ice_area_fraction)

   end subroutine


   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_niva_ecosmo_carbonate_nocalc), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: dic,T,S,PRSS,density
      real(rk) :: wnd,PCO2A,aice
      real(rk) :: sc,fwind,UPTAKE,FAIRCO2

      real(rk) :: Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2
      logical  :: success


      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_dic,dic)
         _GET_(self%id_ETW,T)
         _GET_(self%id_X1X,S)
         _GET_(self%id_pres,PRSS)
         _GET_(self%id_dens,density)
         _GET_HORIZONTAL_(self%id_wnd,wnd)
         _GET_HORIZONTAL_(self%id_PCO2A,PCO2A)
         if (self%use_aice) then
            _GET_HORIZONTAL_(self%id_aice,aice)
         end if

         S = max(S, 0.0_rk) ! Inserted PJW 06/01/2020
         wnd = max(wnd, 0.0_rk)

         ! Alkalinity is a state variable.
         _GET_(self%id_alk,TA)

         TA = TA / 1.e3_rk / density     ! from mmol m-3 to mol kg-1
         Ctot  = dic / 1.e3_rk / density ! from mmol m-3 to mol kg-1

!  for surface box only calculate air-sea flux
!..Only call after 2 days, because the derivation of instability in the
!..
         CALL CO2DYN(T, S, PRSS*0.1_rk,Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success,self%phscale)
         if (.not.success) then
            _GET_(self%id_pco2_in,PCO2)
            PCO2 = PCO2*1.e-6_rk
         end if

         ! New formulation for the Schmidt number for CO2 following Wanninkhof 2014
         T = max(min(T,40.0_rk), -2.0_rk)
         sc = 2116.8_rk - 136.25_rk*T + 4.7353_rk*T**2._rk - 0.092307_rk*T**3 + 0.0007555_rk*T**4._rk

         fwind = 0.251_rk*wnd**2.0_rk*(sc/660._rk)**(-0.5_rk) ! Wanninkhof 2014

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wnd_diag,wnd) ! diagnostic in m/s

         fwind = fwind/360000._rk   ! convert from cm/hr to m/s
         UPTAKE = fwind * k0co2 * ( PCO2A/1.e6_rk - PCO2 )
         if (self%use_aice) then
            UPTAKE = max(0._rk, (1._rk-aice))*UPTAKE !Limit flux to area fraction not covered by ice
         end if

         FAIRCO2 = UPTAKE * 1.e3_rk * density

         _ADD_SURFACE_FLUX_(self%id_dic,FAIRCO2)

      _HORIZONTAL_LOOP_END_
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2DYN \label{sec:CO2DYN}
!
! !DESCRIPTION:
!  TODO - description
!
!     This subroutine calculates the partial pressure of CO2 (pCO2) at
!     the ambient salinity, temperature, alkalinity and total CO2 and
!     hence the CO2 exchange across the air sea interface.
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2DYN ( T, S, PRSS,Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,k0co2,success,hscale)
!
! !LOCAL VARIABLES:
!     ! TODO - SORT THESE!
      real(rk),intent(in) :: T, S, PRSS   ! NB PRSS is pressure in bar
      real(rk),intent(inout) :: Ctot,TA
      real(rk),intent(out) :: pH,PCO2,H2CO3,HCO3,CO3,k0co2
      logical, intent(out) :: success
      integer, intent(in)  :: hscale

      real(rk) :: k1co2,k2co2,kb
      real(rk) :: Tmax, Btot
      INTEGER  :: ICALC
      LOGICAL  :: BORON
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC

      ICALC = 1               ! use DIC and TA to calculate CO2sys
      Tmax = max(T,0._rk)

      BORON=.True.
      IF(BORON) THEN
        BTOT=0.0004128_rk*S/35._rk
      ENDIF

      ! Initialize CO2CLC inputs that are marked intent(inout) to avoid compiler warnings
      ! These are not used in practice because we set ICALC=1
      pH = 0._rk
      PCO2 = 0._rk
      H2CO3 = 0._rk
      HCO3 = 0._rk
      CO3 = 0._rk

      CALL CO2SET(PRSS,Tmax,S,k0co2,k1co2,k2co2,kb,hscale)
      CALL CO2CLC(k0co2,k1co2,k2co2,kb,ICALC,BORON,BTOT,Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,success)

      END SUBROUTINE CO2DYN
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2SET \label{sec:CO2SET}
!
! !DESCRIPTION:
!  TODO - CHECK THIS!
!
! Routine to calculate CO2 system constants under the conditions set by P,S,T
! This subroutine has been rewritten to calculate constants in total scale 
! using Millero, Marine and Freshwater Research, 2010 for k1 and k2 to cover low salinity areas
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2SET(P,T,S,k0co2,k1co2,k2co2,kb,hscale)
         real(rk),intent(in)  :: P,T,S
         integer, intent (in) :: hscale
         real(rk),intent(out) :: k0co2,k1co2,k2co2,kb
!
! !USES:
!
! !LOCAL VARIABLES
!
      real(rk)              :: TK, delta, kappa
      real(rk)              :: dlogTK, S2, S15, sqrtS,TK100
      real(rk),parameter    :: Rgas = 83.131_rk
      real(rk)              :: ST, FT, kS,kF, Cl
      real(rk)              :: is,sqrtis,invtk,total2free_surface,total2free_depth,free2sws_surface,total2sws_surface
      real (rk)             :: free2sws_depth,total2sws_depth,pk1co2,pk2co2

!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!

!  Derive simple terms used more than once
       TK=T+273.15_rk
       invtk=1._rk/TK
       dlogTK = log(TK)
       TK100=TK/100._rk
       S2 = S*S
       sqrtS = sqrt(S)
       S15 = S**1.5_rk


! is : ionic strength, needed to calculate ks 
        is = 19.924_rk*S/(1000._rk-1.005_rk*S)
        sqrtis=sqrt(is)

! cl : Chloride concentration, used to calculate total sulphate and total fluoride
       cl = S / 1.80655_rk

! st : total sulfate using Morris & Riley, Deep Sea Research, 1966
! 0.14 is the S:Cl ratio observed, 96.065 is the molecular weight of SO4--
        st= 0.14_rk * cl / 96.065_rk

! ks = [H][SO4]/[HSO4] in free scale from Dickson, J. chem. Thermodynamics, 1990 and Perez & Frega, Mar Chem, 1987
        ks=exp(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk &
     &      + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk) * sqrtis &
     &      + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is &
     &      - 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is**2._rk &
     &      + log(1.0d0 - 0.001005d0*s))

! ft : total fluoride using Riley, Analytica chimica Acta, 1966
! 0.000067 is the F:Cl ratio observed, 19.9984 is the molecular weight of F-
        ft= 0.000067_rk * cl / 19.9984_rk

! kf = [H][F]/[HF] in total scale from Perez and Fraga (1987)
!          Formulation as given in Dickson et al. (2007)
         kf = exp(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)

! this is the conversion factor from total scale to free scale at surface
        total2free_surface = 1._rk/(1._rk + st/ks)
! this is the conversion factor from free to SWS at surface
        free2sws_surface= 1._rk + st/ks + ft/kf
! this is the conversion factor from total to SWS at surface
        total2sws_surface= total2free_surface*free2sws_surface

! Correction for high pressure (from Mocsy)
        delta=-18.03_rk+0.0466_rk*T+0.000316_rk*T**2._rk
        kappa=-4.53_rk+0.00009_rk*T
        ks=ks*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))
! this is the conversion factor from total scale to free scale at depth
        total2free_depth = 1._rk/(1._rk + st/ks)

! Correction for high pressure (from Mocsy) - this requires kf being in free scale, final value still in total scale
        delta=-9.78_rk-0.009_rk*T-0.000942_rk*T**2._rk
        kappa=-3.91_rk+0.000054_rk*T
        kf=kf*total2free_surface*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))/total2free_depth
! this is the conversion factor from free to SWS at depth
        free2sws_depth= 1._rk + st/ks + ft/(kf*total2free_depth)
! this is the conversion factor from total to SWS at surface
        total2sws_depth= total2free_depth*free2sws_depth

!  Calculation of constants as used in the OCMIP process for ICONST = 3 or 6
!  see http://www.ipsl.jussieu.fr/OCMIP/
! k0co2 = CO2/fCO2 from Weiss 1974
        k0co2 = exp(93.4517_rk/tk100 - 60.2409_rk + 23.3585_rk * log(tk100) + &
        &       s * (.023517_rk - 0.023656_rk * tk100 + 0.0047036_rk * tk100 ** 2._rk))
! correction for high pressure from Weiss 1974
!        vbarCO2 = 32.3      partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
!        P is in bar, hence the reference is 1.01325 instead of 1 as in Weiss 1974
        k0co2 = k0co2 * exp( ((1.01325_rk-P)*32.3_rk)/(Rgas*tk) )
! kb = [H][BO2]/[HBO2]
! Millero p.669 (1995) using data from Dickson (1990)
        kb=exp((-8966.9_rk - 2890.53_rk*sqrtS - 77.942_rk*S + &
     &      1.728_rk*S15 - 0.0996_rk*S2)/TK + &
     &      (148.0248_rk + 137.1942_rk*sqrtS + 1.62142_rk*S) + &
     &      (-24.4344_rk - 25.085_rk*sqrtS - 0.2474_rk*S) * &
     &      dlogTK + 0.053105_rk*sqrtS*TK)

! k1co2 = [H][HCO3]/[H2CO3]
! k2co2 = [H][CO3]/[HCO3]
     if (hscale==-1) then
           ! if phscale = -1 then we use old formulation for backward compatibility
           ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
           ! kb is left in total scale because this was in the old formulation
        k1co2=10._rk**(-1._rk*(3670.7_rk/TK - 62.008_rk + 9.7944_rk*dlogTK - &
     &     0.0118_rk * S + 0.000116_rk*S2))
        k2co2=10._rk**(-1._rk*(1394.7_rk/TK + 4.777_rk - &
     &     0.0184_rk*S + 0.000118_rk*S2))
     else if (hscale==0) then
           ! if phscale = 0 then we use  Millero 2010 on seawater scale
           ! kb is converted in seawater scale
           pk1co2 = -126.34048_rk+6320.813_rk*invtk+19.568224_rk*dlogtk + &
                &  (13.4038_rk*sqrtS + 0.03206_rk*S - 0.00005242_rk*S**2._rk) + &
                &  (-530.659_rk *sqrtS - 5.8210_rk *S) * invTK + &
                &  (-2.0664_rk *sqrts)*dlogTK
           k1co2 =10._rk**(-pk1co2)
           pk2co2 = -90.18333_rk+5143.692_rk*invtk+14.613358_rk*dlogtk +&
                &  (21.3728_rk*sqrtS + 0.1218_rk*S - 0.0003688_rk*S**2._rk) + &
                &  (-788.289_rk *sqrtS - 19.189_rk *S) * invTK + &
                &  (-3.374_rk *sqrts)*dlogTK
           k2co2=10._rk**(-pk2co2)
           kb=kb*total2sws_depth
      else if (hscale==1) then
           ! if phscale = 1 then we use  Millero 2010 on total scale
           pk1co2 = -126.34048_rk+6320.813_rk*invtk+19.568224_rk*dlogtk + &
                &  (13.4051_rk*sqrtS + 0.03185_rk*S - 0.00005218_rk*S**2._rk) + &
                &  (-531.095_rk *sqrtS - 5.7789 *S) * invTK + &
                &  (-2.0663_rk *sqrts)*dlogTK
           k1co2 = 10._rk**(-pk1co2)
           pk2co2 = -90.18333_rk+5143.692_rk*invtk+14.613358_rk*dlogtk +&
                &  (21.5724_rk*sqrtS + 0.1212_rk*S - 0.0003714_rk*S**2._rk) + &
                &  (-798.292_rk *sqrtS - 18.951_rk *S) * invTK + &
                &  (-3.403_rk *sqrts)*dlogTK
           k2co2=10._rk**(-pk2co2)
      endif

! here k1, k2, kb are corrected for high pressure using Millero 1995
! please not that MOCSY assume that this correction is valid for SWSscale, so it first convert everything to SWS then back to total
! differences are minimal
! correction of k1co2
        delta=-25.5_rk+0.1271_rk*T
        kappa=(-3.08_rk+0.0877_rk*T)/1000._rk
        k1co2=k1co2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

! Correction for k2co2
        delta=-15.82_rk-0.0219_rk*T
        kappa=(1.13_rk-0.1475_rk*T)/1000._rk
        k2co2=k2co2*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

! Correction for kb
        delta=-29.48_rk+0.1622_rk*T-0.002608_rk*T**2._rk
        kappa=-2.84_rk/1000._rk
        kb=kb*exp((-delta+0.5_rk*kappa*P)*P/(Rgas*TK))

      END SUBROUTINE CO2SET
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CO2CLC \label{sec:COCCLC}
!
! !DESCRIPTION:
!  TODO - check this.
!
! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.
!\\
!\\
! !INTERFACE:
      SUBROUTINE CO2CLC(k0co2,k1co2,k2co2,kb,ICALC,BORON,BTOT,Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3,success)
!
! !USES:
!
! !LOCAL VARIABLES:
      real(rk),intent(in)    :: k0co2,k1co2,k2co2,kb
      real(rk),intent(inout) :: Ctot,TA,pH,PCO2,H2CO3,HCO3,CO3
      logical,intent(out)    :: success
      INTEGER ICALC, KARL, LQ
!put counter in to check duration in convergence loop
      INTEGER                :: COUNTER,C_CHECK,C_SW
      real(rk)              :: ALKC, ALKB,BTOT
      real(rk)              :: AKR,AHPLUS
      real(rk)              :: PROD,tol1,tol2,tol3,tol4,steg,fak
      real(rk)              :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
      LOGICAL               :: BORON,DONE
!
! !REVISION HISTORY:
!  Original author(s) TODO
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
!  DERIVING PH REQUIRES FOLLOWING LOOP TO CONVERGE.
!  THIS SUBROUTINE RELIES ON CONVERGENCE.  IF THE ENVIRONMENTAL
!  CONDITIONS DO NOT ALLOW FOR CONVERGENCE (IN 3D MODEL THIS IS
!  LIKELY TO OCCUR NEAR LOW SALINITY REGIONS) THE MODEL WILL
!  BE STUCK IN THE LOOP.  TO AVOID THIS A CONVERGENCE CONDITION
!  IS PUT IN PLACE TO SET A FLAGG OF -99 IN THE PH VAR FOR NON CONVEGENCE.
!  THE MODEL IS THEN ALLOWED TO CONTINUE. 'COUNTER, C_SW,C_CHECK' ARE
!  THE LOCAL VARS USED.
! C_SW = condition of convergence 0=yes, 1= no
! COUNTER = number of iterations
! C_CHECK = maximum number of iterations
      success = .true.
      DONE = .false.
! SET COUNTER AND SWITCH TO ZERO AND OFF
      COUNTER=0
      C_SW=0
! FROM EXPERIENCE IF THE ITERATIONS IN THE FOLLOWING DO LOOP
! EXCEEDS 15 CONVERGENCE WILL NOT OCCUR.  THE OVERHEAD OF 25 ITERATIONS
! IS OK FOR SMALL DOMAINS WITH 1/10 AND 1/15 DEG RESOLUTION.
! I RECOMMEND A LOWER VALUE OF 15 FOR HIGHER RESOLUTION OR LARGER DOMAINS.
      C_CHECK=25

!      DO II=1,NKVAL
 !       AKVAL2(II)=AKVAL(II)
  !    END DO

      AKR = k1co2/k2co2
      AHPLUS=10._rk**(-PH)
      PROD=AKR*k0co2*PCO2

      IF(BORON) THEN

        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** TA, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE

          H2CO3=PCO2*k0co2
          CO3=TA/10._rk
          AHPLUS=1.e-8_rk
          ALKB=BTOT
          TOL1=TA/1.e5_rk
          TOL2=H2CO3/1.e5_rk
          TOL3=CTOT/1.e5_rk
          TOL4=BTOT/1.e5_rk
!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
          KARL=1
          STEG=2._rk
          FAK=1._rk
          STEGBY=0.4_rk

          DO WHILE (.not.DONE)
            DONE=.TRUE.

!
! SET COUNTER UPDATE.
            COUNTER=COUNTER+1

! CHECK IF CONVERGENCE HAS OCCURED IN THE NUMBER OF
! ACCEPTABLE ITTERATIONS.
            if(counter.ge.c_check)then
!!        IF(MASTER)THEN
!!! LOG FILE TO SHOW WHEN AND WHERE NON CONVERGENCE OCCURS.
!!           PPWRITELOG 'ERRORLOG ',III,' ',(CONCS2(II),II=1,NCONC)
!!        ENDIF
! IF NON CONVERGENCE, THE MODEL REQUIRES CONCS TO CONTAIN USABLE VALUES.
! BEST OFFER BEING THE OLD CONCS VALUES WHEN CONVERGENCE HAS BEEN
! ACHIEVED
              success = .false.

!RESET SWITCH FOR NEXT CALL TO THIS SUBROUTINE
              C_SW=0
              RETURN
!
            endif

            IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***
              Y=AHPLUS*AHPLUS*CO3/(k1co2*k2co2)
              IF(ABS(Y-H2CO3).GT.TOL2) THEN
                CO3=CO3*H2CO3/Y
                DONE=.FALSE.
              ENDIF
            ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED ***
              Y=CO3*(1._rk+AHPLUS/k2co2+AHPLUS*AHPLUS/(k1co2*k2co2))
              IF(ABS(Y-CTOT).GT.TOL3) THEN
                CO3=CO3*CTOT/Y
                DONE=.FALSE.
              ENDIF
            ENDIF
            Y=ALKB*(1._rk+AHPLUS/kb)
            IF(ABS(Y-BTOT).GT.TOL4) THEN
              ALKB=ALKB*BTOT/Y
              DONE=.FALSE.
            ENDIF

! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
! to that normally used

            Y=CO3*(2._rk+AHPLUS/k2co2)+ALKB
            IF(ABS(Y-TA).GT.TOL1) THEN
              DONE=.FALSE.
              X=LOG(AHPLUS)
              W=SIGN(1._rk,Y-TA)
              IF(W.GE.0._rk) THEN
                X1=X
                Y1=Y
              ELSE
                X2=X
                Y2=Y
              ENDIF
              LQ=KARL
              IF(LQ.EQ.1) THEN
                KARL=2*NINT(W)
              ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0._rk) THEN
                FAK=0.5_rk
                KARL=3
              ENDIF
              IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN
                W=(X2-X1)/(Y2-Y1)
                X=X1+W*(TA-Y1)
              ELSE
                STEG=STEG*FAK
                X=X+STEG*W
              ENDIF
              AHPLUS=EXP(X)
            ENDIF
! LOOP BACK UNTIL CONVERGENCE HAS BEEN ACHIEVED
! OR MAX NUMBER OF ITERATIONS (C_CHECK) HAS BEEN REACHED.
          ENDDO

          HCO3=CO3*AHPLUS/k2co2
          IF(ICALC.EQ.4) THEN
            CTOT=H2CO3+HCO3+CO3
          ELSEIF(ICALC.EQ.1) THEN
            H2CO3=HCO3*AHPLUS/k1co2
            PCO2=H2CO3/k0co2
          ENDIF
          PH=-LOG10(AHPLUS)
          ALKC=TA-ALKB
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
          Y=SQRT(PROD*(PROD-4._rk*k0co2*PCO2+4._rk*CTOT))
          H2CO3=PCO2*k0co2
          HCO3=(Y-PROD)/2._rk
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2._rk*CO3
          AHPLUS=k1co2*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          TA=ALKC+ALKB
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+k1co2*AHPLUS+k1co2*k2co2)
          CO3=FACTOR*k1co2*k2co2
          HCO3=FACTOR*k1co2*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/k0co2
          ALKC=HCO3+2._rk*CO3
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          TA=ALKC+ALKB
        ELSEIF(ICALC.EQ.5) THEN
!         *** TA, PH AND BTOT FIXED ***
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          ALKC=TA-ALKB
          HCO3=ALKC/(1._rk+2._rk*k2co2/AHPLUS)
          CO3=HCO3*k2co2/AHPLUS
          H2CO3=HCO3*AHPLUS/k1co2
          PCO2=H2CO3/k0co2
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
          ALKB=BTOT/(1._rk+AHPLUS/kb)
          H2CO3=PCO2*k0co2
          HCO3=H2CO3*k1co2/AHPLUS
          CO3=HCO3*k2co2/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALKC=HCO3+2._rk*CO3
          TA=ALKC+ALKB
        ENDIF
      ELSE
        IF(ICALC.EQ.1) THEN
!         *** CTOT AND TA FIXED ***
          TERM=4._rk*TA+CTOT*AKR-TA*AKR
          Z=SQRT(TERM*TERM+4._rk*(AKR-4._rk)*TA*TA)
          CO3=(TA*AKR-CTOT*AKR-4._rk*TA+Z)/(2._rk*(AKR-4._rk))
          HCO3=(CTOT*AKR-Z)/(AKR-4._rk)
          H2CO3=CTOT-TA+CO3
          PCO2=H2CO3/k0co2
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4._rk*k0co2*PCO2+4._rk*CTOT))
          H2CO3=PCO2*k0co2
          HCO3=(Y-PROD)/2._rk
          CO3=CTOT-H2CO3-HCO3
          TA=HCO3+2._rk*CO3
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+k1co2*AHPLUS+k1co2*k2co2)
          CO3=FACTOR*k1co2*k2co2
          HCO3=FACTOR*k1co2*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/k0co2
          TA=HCO3+2._rk*CO3
        ELSEIF(ICALC.EQ.4) THEN
!         *** TA AND PCO2 FIXED ***
          TERM=SQRT((8._rk*TA+PROD)*PROD)
          CO3=TA/2._rk+PROD/8._rk-TERM/8._rk
          HCO3=-PROD/4._rk+TERM/4._rk
          H2CO3=PCO2*k0co2
          CTOT=CO3+HCO3+H2CO3
          PH=-LOG10(k1co2*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.5) THEN
!         *** TA AND PH FIXED ***
          HCO3=TA/(1._rk+2._rk*k2co2/AHPLUS)
          CO3=HCO3*k2co2/AHPLUS
          H2CO3=HCO3*AHPLUS/k1co2
          PCO2=H2CO3/k0co2
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2 AND PH FIXED ***
          H2CO3=PCO2*k0co2
          HCO3=H2CO3*k1co2/AHPLUS
          CO3=HCO3*k2co2/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          TA=HCO3+2._rk*CO3
        ENDIF
      ENDIF

      END SUBROUTINE CO2CLC

   end module
