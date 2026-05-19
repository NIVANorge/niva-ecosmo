#include "fabm_driver.h"

! -----------------------------------------------------------------------------
! This is a two-band model based on the Aksnes (2015) empirical model for
! absorption at 440 nm, the Lee et al. (2005) model to relate Kd to IOPs,
! a power law to relate Kd(440 nm) to Kd(400-550 nm) (Lee et al., 2014),
! and a chlorophyll-based model for Kd(550-700 nm) (Foujols et al., 2000;
! Manizza et al. 2005).
!
! It calculates profiles of blue PAR (400-560 nm), red PAR (560-700 nm), and
! total PAR (quantum and energy units).
!
! Change Log (last updated: 25/09/2018)----------------------------------------
! By: Phil Wallhead, branched from ersem_edge/src/light_iop.F90 (last updated 21/04/2017)
!
! 09/09/2019: Carried across update to base code ersem-edge/src/light_iop.F90 (git pulled 09/09/2019)
!             (application of max() instead of if statement for I_0).
!
! 06/01/2020: Cap input salinity at 0 psu.
!
! 09/04/2026: Refined model for blue band attenuation, and update default parameter
!             values using optimizations based on ECOLIGHT data (see fit_light_iop_aksnes_to_ecolight.m).
!             Simplified model for bbp, using TSM-bp scaling from Sørensen et al. (2007)
!             and bbp:bp ratio from Petzold (1972), Bowers and Binding (2006).
! -----------------------------------------------------------------------------

module light_iop_aksnes

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_light_iop_aksnes
      ! Identifiers for state variables of other models
      type (type_state_variable_id) :: id_O2o

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)   :: id_Chl,id_TSOM,id_TSIM,id_TSM
      type (type_diagnostic_variable_id)   :: id_a440,id_bb440,id_Kd440,id_K0b,id_K0r
      type (type_diagnostic_variable_id)   :: id_E0,id_PAR0,id_PAR0b,id_PAR0r

      ! Environmental dependencies
      type (type_dependency_id)            :: id_dz,id_X1X,id_TSS,id_Chl_s,id_TSOM_s,id_TSIM_s
      type (type_horizontal_dependency_id) :: id_I_0,id_zenithA

      ! Parameters
      real(rk) :: facTSS,rE00_I_0,fE00b,k0,k1,k2,k3,a440min
      real(rk) :: bbp440_per_TSM,bb440w
      real(rk) :: kb0,kb1,kb2,kb3,kr0,kr1,kr2,Qb,Qr
      logical  :: use_TSS
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_light_iop_aksnes

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_light_iop_aksnes),intent(inout),target :: self
      integer,                      intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
      !Parameters for input of total suspended sediments data
      call self%get_parameter(self%use_TSS,'use_TSS','','use input Total Suspended Sediments to supplement TSM',default=.false.)
      call self%get_parameter(self%facTSS,'facTSS','','scale factor for input Total Suspended Sediments',default=1._rk)

      !Surface transmission parameters
      call self%get_parameter(self%rE00_I_0,'rE00_I_0','-','ratio of scalar PAR energy flux just below the sea surface to input energy flux',default=0.54_rk) !Default assumes input downwelling above sea surface, PAR = 0.43*shortwave (Aumont et al., 2015), and surface transmission factors from Mobley and Boss (2012)
      call self%get_parameter(self%fE00b,'fE00b','-','blue fraction of scalar PAR energy flux just below sea surface',default=0.5_rk) !Manizza et al. (2005)

      !Absorption model parameters (Aksnes, 2015)
      call self%get_parameter(self%k0,'k0','1/m','absorption (440 nm) at zero (S,O2,Chl)',default=1.7_rk)
      call self%get_parameter(self%k1,'k1','1/m/PSU','absorption (440 nm) increase per unit salinity increase',default=-0.041_rk)
      call self%get_parameter(self%k2,'k2','1/m/(mmol O2/m3)','absorption (440 nm) increase per unit oxygen increase',default=-7.165e-4_rk)
      call self%get_parameter(self%k3,'k3','1/m/(mg Chla/m3)','absorption (440 nm) increase per unit chlorophyll a increase',default=0.073_rk)
      call self%get_parameter(self%a440min,'a440min','1/m','minimum absorption at 440 nm (clear water)',default=0.006_rk) !Pope and Fry (1997)

      !Backscatter model parameters (Sørensen et al., 2007)
      call self%get_parameter(self%bbp440_per_TSM,'bbp440_per_TSM','m^2/g','linear scale parameter relating bbp(440nm) to TSM',default=0.0145_rk) !Assuming TSM=1.73*bp (Sørensen et al., 2007) and bbp:bp = 0.025 (Petzold, 1972; Bowers and Binding, 2006).
      call self%get_parameter(self%bb440w,'bb440w','1/m','backscatter of clear water at 440 nm', default=.0022_rk) !Buiteveld and Donze (1994), default assumes T = 10 degC, S = 35 psu, and depolarization ratio rho = 0.039 (Reynolds et al., 2016)

      !Broadband attenuation model parameters
      call self%get_parameter(self%kb0,'kb0','1/m (1/m)^(-kb1)','power-law scale parameter relating Ko(blue) to Kd(440) at zero (Chla, depth)',default=0.52_rk) !see fit_light_iop_aknes_to_ecolight.m
      call self%get_parameter(self%kb1,'kb1','-','power-law exponent relating Ko(blue) to Kd(440)',default=0.67_rk) !see fit_light_iop_aknes_to_ecolight.m
      call self%get_parameter(self%kb1,'kb2','-','power-law exponent relating Ko(blue) to (1+Chla)',default=0.14_rk) !see fit_light_iop_aknes_to_ecolight.m
      call self%get_parameter(self%kb3,'kb3','-','power-law exponent relating Ko(blue) to (1+depth)',default=-0.083_rk) !see fit_light_iop_aknes_to_ecolight.m
      call self%get_parameter(self%kr0,'kr0','1/m','red PAR attenuation due to water',default=0.225_rk) !Foujols et al. (2000), Manizza et al. (2005)
      call self%get_parameter(self%kr1,'kr1','1/m (mg Chla/m3)^(-kr2)','power-law scale parameter relating Ko(red) to Chla',default=0.037_rk) !Foujols et al. (2000), Manizza et al. (2005)
      call self%get_parameter(self%kr2,'kr2','-','power-law exponent relating Kd(red) to Chla',default=0.629_rk) !Foujols et al. (2000), Manizza et al. (2005)

      !Unit conversion parameters (energy to quantum units)
      call self%get_parameter(self%Qb,'Qb','umol quanta/s/W','factor to convert blue PAR energy flux to quantum units',default=4.01_rk) !Estimate as Q = 1e6*mean(lambda)/(N_A*h*c)
      call self%get_parameter(self%Qr,'Qr','umol quanta/s/W','factor to convert red PAR energy flux to quantum units',default=5.26_rk)


      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_Chl,'Chl','mg/m^3','total chlorophyll a', source=source_do_column)
      call self%register_diagnostic_variable(self%id_TSOM,'TSOM','mg/L','total suspended organic matter', source=source_do_column)
      call self%register_diagnostic_variable(self%id_TSIM,'TSIM','mg/L','total suspended inorganic matter', source=source_do_column)
      call self%register_diagnostic_variable(self%id_TSM,'TSM','mg/L','total suspended matter', source=source_do_column)
      !Note: The aggregate results will also be automatically stored in diagnostics e.g. Chl_calculator_result
      !      HOWEVER: These diagnostics may not be calculated if only one variable is contributing (this may be addressed in future FABM updates)
      !      Therefore it is safest to define a separate diagnostic, for now.
      call self%register_diagnostic_variable(self%id_a440,'a440','1/m','total absorption at 440 nm', source=source_do_column)
      call self%register_diagnostic_variable(self%id_bb440,'bb440','1/m','total backscatter at 440 nm', source=source_do_column)
      call self%register_diagnostic_variable(self%id_Kd440,'Kd440','1/m','diffuse attenuation of downwelling irradiance at 440 nm', source=source_do_column)
      call self%register_diagnostic_variable(self%id_K0b,'K0b','1/m','diffuse attenuation of scalar blue PAR', source=source_do_column)
      call self%register_diagnostic_variable(self%id_K0r,'K0r','1/m','diffuse attenuation of scalar red PAR', source=source_do_column)

      ! Register diagnostics to be used elsewhere as standard variables (simplifies the fabm.yaml)
      call self%register_diagnostic_variable(self%id_E0,'E0','W/m2', &
              'scalar photosynthetically active radiation energy flux', &
              standard_variable=type_bulk_standard_variable(name='E0_s'), source=source_do_column)
      call self%register_diagnostic_variable(self%id_PAR0,'PAR0','umol quanta/m2/s', &
              'quantum scalar photosynthetically active radiation', &
              standard_variable=type_bulk_standard_variable(name='PAR0_s'), source=source_do_column)
      call self%register_diagnostic_variable(self%id_PAR0b,'PAR0b','umol quanta/m2/s', &
              'quantum scalar irradiance in blue PAR band', &
              standard_variable=type_bulk_standard_variable(name='PAR0b_s'), source=source_do_column)
      call self%register_diagnostic_variable(self%id_PAR0r,'PAR0r','umol quanta/m2/s', &
              'quantum scalar irradiance in red PAR band', &
              standard_variable=type_bulk_standard_variable(name='PAR0r_s'), source=source_do_column)

      ! State variable dependencies
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen')

      ! Environmental dependencies
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_X1X, standard_variables%practical_salinity)
      call self%register_horizontal_dependency(self%id_zenithA, type_horizontal_standard_variable(name='zenith_angle'))
      if (self%use_TSS) then
        call self%register_dependency(self%id_TSS,type_bulk_standard_variable(name='tss',units='kg/m^3'))
      end if

      ! Aggregate variable dependencies
      call self%register_dependency(self%id_Chl_s, type_bulk_standard_variable(name='Chl_s',units='mg/m^3',aggregate_variable=.true.))
      call self%register_dependency(self%id_TSOM_s,type_bulk_standard_variable(name='TSOM_s',units='mg/L',aggregate_variable=.true.))
      call self%register_dependency(self%id_TSIM_s,type_bulk_standard_variable(name='TSIM_s',units='mg/L',aggregate_variable=.true.))
   end subroutine

   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_light_iop_aksnes),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: I_0,zenithA,Edair,E00,E0b,E0r,z,dz,X1X,O2o,Chl,TSOM,TSS,TSIM,TSM
      real(rk) :: a440,bbp440,bb440
      real(rk) :: Kd440,K0b,K0r,xtncb,xtncr,E0bav,E0rav,E0av,PAR0bav,PAR0rav,PAR0av

      _GET_HORIZONTAL_(self%id_I_0,I_0)
      _GET_HORIZONTAL_(self%id_zenithA,zenithA)   ! Zenith angle

      ! Propagate light through sea surface
      I_0 = max(I_0, 0.0_rk)
      E00 = self%rE00_I_0*I_0         ! Scalar PAR energy flux just below the sea surface [W/m2]

      E0b = self%fE00b*E00            ! Scalar blue PAR energy flux, initialized just below sea surface [W/m2]
      E0r = (1._rk - self%fE00b)*E00  ! Scalar red PAR energy flux, initialized just below sea surface [W/m2]

      z = 0.0_rk
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)                  ! Layer height [m]

         _GET_(self%id_X1X,X1X)                ! Practical salinity [psu]
         _GET_(self%id_O2o,O2o)                ! Dissolved oxygen [mmol O_2/m^3]
         _GET_(self%id_Chl_s,Chl)               ! Total chlorophyll a [mg/m3]
         _GET_(self%id_TSOM_s,TSOM)            ! Total Suspended Organic Matter (mg/L)
         _GET_(self%id_TSIM_s,TSIM)            ! Total Suspended Inorganic Matter (mg/L)

         z = z + dz * 0.5_rk                   ! Depth of layer centre
         X1X = max(X1X, 0.0_rk)
         TSOM = max(TSOM, 0.0_rk)
         TSIM = max(TSIM, 0.0_rk)
         if (self%use_TSS) then
            _GET_(self%id_TSS,TSS)                     ! Total Suspended Sediments (kg/m3) from sediment module
            TSS = max(self%facTSS*1.e3_rk*TSS, 0.0_rk) ! Convert TSS from (kg/m3=g/L) to (mg/L) and scale by factor facTSS
            TSIM = TSIM + TSS                          ! Assume suspended sediments are entirely inorganic
         end if
         TSM = TSOM + TSIM

         _SET_DIAGNOSTIC_(self%id_Chl,Chl)     ! Total chlorophyll a [mg/m3]
         _SET_DIAGNOSTIC_(self%id_TSOM,TSOM)   ! Total Suspended Organic Matter (mg/L)
         _SET_DIAGNOSTIC_(self%id_TSIM,TSIM)   ! Total Suspended Inorganic Matter (mg/L)
         _SET_DIAGNOSTIC_(self%id_TSM,TSM)     ! Total Suspended Matter (mg/L)


         ! Total absorption at 440 nm [1/m] using empirical formula from (Aksnes, 2015)
         a440 = max(self%a440min, self%k0 + self%k1*X1X + self%k2*O2o + self%k3*Chl)

         ! Particulate backscatter at 440 nm by simple scaling with TSM (Sørensen et al., 2007)
         bbp440 = self%bbp440_per_TSM * TSM

         ! Total backscatter including water component
         bb440 = self%bb440w + bbp440

         ! Downwelling attenuation at 440 nm using formula from Lee et al. (2005)
         Kd440 = (1._rk+.005_rk*zenithA)*a440 + 4.18_rk*(1._rk-.52_rk*exp(-10.8_rk*a440))*bb440

         ! Power law to derive scalar blue-band attenuation from Kd(440 nm), Chl, depth
         K0b = self%kb0 * Kd440**self%kb1 * (1._rk+Chl)**self%kb2 * (1._rk+z)**self%kb3

         ! Simple biooptical model for scalar red band attenuation (Foujols et al., 2000)
         K0r = self%kr0 + self%kr1*(Chl**self%kr2)

         xtncb = K0b*dz
         xtncr = K0r*dz

         E0bav = E0b/xtncb*(1.0_rk-exp(-xtncb))  ! Note: this computes the vertical average, not the value at the layer centre.
         E0rav = E0r/xtncr*(1.0_rk-exp(-xtncr))
         E0av  = E0bav + E0rav
         PAR0bav = self%Qb * E0bav
         PAR0rav = self%Qr * E0rav
         PAR0av  = PAR0bav + PAR0rav

         E0b = E0b*exp(-xtncb)
         E0r = E0r*exp(-xtncr)

         _SET_DIAGNOSTIC_(self%id_a440,a440)       ! Local total absorption at 440 nm [1/m]
         _SET_DIAGNOSTIC_(self%id_bb440,bb440)     ! Local total backscatter at 440 nm [1/m]
         _SET_DIAGNOSTIC_(self%id_Kd440,Kd440)     ! Local attenuation of downwelling irradiance at 440 nm [1/m]
         _SET_DIAGNOSTIC_(self%id_K0b,K0b)         ! Local attenuation of scalar energy flux over blue PAR band [1/m]
         _SET_DIAGNOSTIC_(self%id_K0r,K0r)         ! Local attenuation of scalar energy flux over red PAR band [1/m]
         _SET_DIAGNOSTIC_(self%id_E0,E0av)         ! Local scalar PAR energy flux [W/m2]
         _SET_DIAGNOSTIC_(self%id_PAR0,PAR0av)     ! Local quantum scalar PAR [umol quanta/m2/s]
         _SET_DIAGNOSTIC_(self%id_PAR0b,PAR0bav)   ! Local quantum scalar irradiance in blue PAR band [umol quanta/m2/s]
         _SET_DIAGNOSTIC_(self%id_PAR0r,PAR0rav)   ! Local quantum scalar irradiance in red PAR band [umol quanta/m2/s]
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module
