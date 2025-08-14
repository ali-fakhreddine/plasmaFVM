!=======================================================================
!  File       : m_constants.f90
!  Author     : Ali Fakhreddine
!  Created    : 14-07-2025
!  Description: Defines fundamental physical constants for use in
!               simulations.
!=======================================================================
module m_constants
  use m_precision
  implicit none

  real(dp), parameter :: qe = -1.602176634e-19   ! electron charge [C]
  real(dp), parameter :: me   = 9.10938356e-31   ! electron mass   [kg]
  real(dp), parameter :: qi =  1.602176634e-19   ! ion charge [C]
  real(dp), parameter :: mi   = 1.67262192e-27   ! ion mass   [kg]
  real(dp), parameter :: eps0 = 8.854187817e-9  ! permittivity    [F/m]
  real(dp), parameter :: kB   = 1.380649e-23     ! Boltzmann const [J/K]
  real(dp), parameter :: mu_e  = -0.01_dp          ! Mobility m2/V.s
  real(dp), parameter :: mu_i  = 0.0001_dp        ! Mobility m2/V.s 
  real(dp), parameter :: charge = 1.602176634e-19

end module m_constants
