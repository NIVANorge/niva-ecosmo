#include "fabm_driver.h"

! This module calculates the day length as a function of longitude, latitude, and day of year.
! The module uses the NOAA formula for declination (following ersem/zenith_angle.F90) and
! the corresponding NOAA formula for day length (see sunrisesunset.m).
!
! Change Log --------------------------------------------------------------------------------------
! By: Phil Wallhead, NIVA, branched from ersem-edge/src/zenith_angle.F90
!
! 27/10/2020: Renamed modules and classes to format niva_ersem_*.
!
! 25/11/2021: Renamed daylength_niva.F90 -> nersem_daylength.F90
!                     niva_ersem_... -> nersem_ in modules and classes
!
! 16/06/2022: Replaced ersem_shared->nersem_shared.
!
! 31/03/2026: Renamed to nersem->niva_ecosmo.
!
! 09/04/2026: Added pi as parameter (since not taken from nersem_shared).
! -------------------------------------------------------------------------------------------------

module niva_ecosmo_daylength

   use fabm_types

   implicit none

   private

   real(rk), parameter :: pi=3.141592653589793_rk

   type,extends(type_base_model),public :: type_niva_ecosmo_daylength
      type (type_horizontal_dependency_id)          :: id_lat
      type (type_horizontal_diagnostic_variable_id) :: id_day_length
      type (type_global_dependency_id)              :: id_yday
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self,configunit)
    class (type_niva_ecosmo_daylength), intent(inout), target :: self
    integer,                            intent(in)            :: configunit

      call self%register_global_dependency(self%id_yday,standard_variables%number_of_days_since_start_of_the_year)
      call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
      call self%register_horizontal_diagnostic_variable(self%id_day_length,'day_length','-','day length as fraction', &
              standard_variable=type_horizontal_standard_variable(name='day_length_s'),source=source_do_surface)

   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_niva_ecosmo_daylength), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: D,lat,latr,th0,th02,th03,sundec,yday
      integer :: iday

      ! Retrieve time since beginning of the year
      _GET_GLOBAL_(self%id_yday,yday)

      !Get day integer and day fraction in hours:
      iday=int(yday)+1

      ! Leap year not considered:
      th0 = pi*iday/182.5_rk
      th02 = 2._rk*th0
      th03 = 3._rk*th0

      ! Sun declination :
      sundec = 0.006918_rk - 0.399912_rk*cos(th0) + 0.070257_rk*sin(th0) &
            - 0.006758_rk*cos(th02) + 0.000907_rk*sin(th02)              &
            - 0.002697_rk*cos(th03) + 0.001480_rk*sin(th03)

      ! Enter spatial loops (if any)
      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve latitude
         _GET_HORIZONTAL_(self%id_lat,lat)
         latr = lat*pi/180._rk

         ! Calculate day length
         D = acos(min(1._rk,max(-1._rk,                                  &
           cos(90.833_rk*pi/180._rk)/(cos(latr)*cos(sundec))             &
           - tan(latr)*tan(sundec))))/pi

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_day_length,D)

      ! Leave spatial loops (if any)
      _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module