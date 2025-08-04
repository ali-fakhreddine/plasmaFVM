!=======================================================================
!  File       : m_varspace.f90
!  Author     : Ali Fakhreddine
!  Created    : 14-07-2025
!  Description: Provides variable-space definitions and management
!               routines for grid/field arrays.
!=======================================================================
module m_varspace
  ! Use associations
  use m_precision
  implicit none

  public
  ! Grid variables
  integer(ii) :: nx, ny, nz
  real(dp) :: Lx, Ly, Lz, dx, dy, dz
  real(dp), allocatable :: xc(:), yc(:), zc(:)  ! Cell-center grid coordinates
  real(dp), allocatable :: xf(:), yf(:), zf(:)  ! Cell-face grid coordinates

  ! Field variables
  real(dp) :: input_ne, input_ni, Ts
  real(dp), allocatable :: phi(:,:), rhs(:,:)
  real(dp), allocatable :: ue(:,:), ve(:,:), ne(:,:), pe(:,:) 
  real(dp), allocatable :: rho_e(:,:), rhou_e(:,:), rhov_e(:,:)
  real(dp), allocatable :: rhouu_e(:,:), rhouv_e(:,:), rhovu_e(:,:), rhovv_e(:,:) 
  real(dp), allocatable :: ui(:,:), vi(:,:), ni(:,:), pi(:,:)
  real(dp), allocatable :: rho_i(:,:), rhou_i(:,:), rhov_i(:,:)
  real(dp), allocatable :: Ex(:,:), Ey(:,:), dpedx(:,:), dpedy(:,:) 

contains

  subroutine allocate_variables
    implicit none

    ! Allocate grid metrics
    allocate(xc(1:nx),yc(1:ny),zc(1:nz))
    allocate(xf(0:nx),yf(0:ny),zf(0:nz))

    ! Allocate field variables
    allocate(phi(0:nx+1,0:ny+1), rhs(0:nx+1,0:ny+1))
    allocate(pe(0:nx+1,0:ny+1),pi(0:nx+1,0:ny+1))
    allocate(ne(0:nx+1,0:ny+1),ni(0:nx+1,0:ny+1))
    allocate(rho_e(0:nx+1,0:ny+1),rho_i(0:nx+1,0:ny+1))
    allocate(ue(0:nx+1,0:ny+1), ve(0:nx+1,0:ny+1))
    allocate(ui(0:nx+1,0:ny+1), vi(0:nx+1,0:ny+1))
    allocate(rhou_e(0:nx+1,0:ny+1), rhov_e(0:nx+1,0:ny+1))
    allocate(rhouu_e(0:nx+1,0:ny+1), rhouv_e(0:nx+1,0:ny+1))
    allocate(rhovu_e(0:nx+1,0:ny+1), rhovv_e(0:nx+1,0:ny+1))
    allocate(rhou_i(0:nx+1,0:ny+1), rhov_i(0:nx+1,0:ny+1))
    allocate(Ex(0:nx,0:ny+1), Ey(0:nx+1,0:ny))
    allocate(dpedx(0:nx,0:ny+1), dpedy(0:nx+1,0:ny))
  end subroutine allocate_variables

  subroutine initialize_variables
    implicit none

    xc = 0.0_dp; yc = 0.0_dp; zc = 0.0_dp
    xf = 0.0_dp; yf = 0.0_dp; zf = 0.0_dp

    Ts = 0.0_dp
    input_ne = 0.0_dp
    input_ni = 0.0_dp

    phi = 0.0_dp; rhs = 0.0_dp
    pe  = 0.0_dp; pi = 0.0_dp
    rho_e = 0.0_dp; rhou_e = 0.0_dp; rhov_e = 0.0_dp
    rhouu_e = 0.0_dp; rhouv_e = 0.0_dp
    rhovu_e = 0.0_dp; rhovv_e = 0.0_dp
    Ex = 0.0_dp; Ey = 0.0_dp

  end subroutine initialize_variables

end module m_varspace
