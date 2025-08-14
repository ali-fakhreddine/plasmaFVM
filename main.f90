!=======================================================================
!  File       : main.f90
!  Author     : Ali Fakhreddine
!  Created    : 14-07-2025
!  Description: Main driver program that sets up parameters, allocates
!               and initializes data structures, and advances the
!               simulation time loop.
!=======================================================================
program main
  ! Use associations
  use m_varspace
  use m_precision
  use m_constants
  use m_grid
  use m_dataio
  use m_init
  use m_solver
  implicit none
  integer(ii) :: iter, save_iter, bc_option
  real(dp) :: t, dt, t_final, n0
  real(dp) :: phi_L, phi_R, phi_T, phi_B
  real(dp) :: sigma, delta_n0, urf

  ! -------------------------
  ! Grid input parameters
  ! -------------------------
  Lx = 1d-3; Ly = 1d-3
  nx = 64; ny = 64; nz = 1
  ! -------------------------
  call allocate_variables
  call initialize_variables
  ! -------------------------
  ! Sim. input parameters
  ! -------------------------
  save_iter = 10000
  input_ne = 1.0d17
  input_ni = 1.0d17
  n0 = 1.0d17
  t  = 0.0_dp
  dt = 1d-9
  t_final = 5d-4

  Ts = 5000_dp
  diff_e = mu_e * kB * Ts/qe
  diff_i = mu_i * kB * Ts/qi
 
  phi_L = 0.0_dp
  phi_R = 0.0_dp
  phi_T = 0.0_dp
  phi_B = 0.0_dp
  urf = 0.5_dp

  sigma = 0.0001_dp
  delta_n0 = 0.02_dp
  ! -------------------------
  call create_grid
  call initialize_species_Gaussian(ne, ni, n0, sigma, delta_n0)
  rho_e = me*ne
  rho_i = mi*ni

  !bc_option = 1 --> all Dirichlet
  !bc_option = 2 --> left/right Dirichlet, top/bottom Neumann
  !bc_option = 3 --> left/right Neumann, top/bottom Dirichlet
  !bc_option = 4 --> all Neumann

  bc_option = 1
  phi = 0.0_dp
  call initialize_field(bc_option, phi, phi_L, phi_R, phi_B, phi_T)
  call compute_rhs(ne, ni, rhs)
  call solve_poisson_SOR(phi, rhs, 1d-8, 5000, urf)
  call compute_gradient(Ex, Ey, -phi)
  call compute_velocity()
  call write_solution_vtk('result_'//trim(filenostr(0))//'.vtk')

  iter = 1
  do while (t < t_final)

     call advance_conserved_quantity_MUSCL(dt, rho_e, ue, ve, rhou_e, rhov_e)
     call advance_conserved_quantity_MUSCL(dt, rho_i, ui, vi, rhou_i, rhov_i)
     ne =  rho_e/me 
     ni = rho_i/mi

     call compute_rhs(ne,ni,rhs)
     call solve_poisson_SOR(phi, rhs, 1d-8, 5000, urf)
     call compute_gradient(Ex, Ey, -phi)
     call compute_velocity()

     if(mod(iter,save_iter) .eq. 0) then
        call write_solution_vtk('result_'//trim(filenostr(iter))//'.vtk')
     endif

     t = t + dt
     iter = iter + 1
     print*, 'Simulation at t =', t, ', Time iterations:', iter
  end do

  print*, '--------------------------------------------------------'
  print*, 'Simulation complete at t =', t
  print*, '--------------------------------------------------------'

contains
  function filenostr(i) result(str)
    implicit none
    integer(ii), intent(in) :: i
    character(len=10) :: str
    write(str, '(i10.10)') i
  end function filenostr

end program main
