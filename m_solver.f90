!=======================================================================
!  File       : m_solver.f90
!  Author     : Ali Fakhreddine
!  Created    : 16-07-2025
!  Description: Defines general numerical solver routines for field and
!               variable operations within the simulation framework.
!=======================================================================
module m_solver
  ! Use associations
  use m_precision
  use m_varspace
  use m_constants
  implicit none

contains

  subroutine compute_rhs(s1, s2, rhsVal)
    implicit none
    integer(ii) :: i, j
    real(dp), intent(in) :: s1(0:nx+1,0:ny+1), s2(0:nx+1,0:ny+1)
    real(dp), intent(inout) :: rhsVal(0:nx+1,0:ny+1)
    do j = 0, ny+1
       do i = 0, nx+1
          rhsVal(i,j) = -charge*(s2(i,j) - s1(i,j))
          rhsVal(i,j) = rhsVal(i,j)/eps0
       enddo
    enddo

  end subroutine compute_rhs

  subroutine update_primitive_variables(Tconst, p1, n1, p2, n2)
    implicit none
    integer(ii) :: i, j
    real(dp), intent(in) :: Tconst
    real(dp), intent(in), optional :: n1(0:nx+1,0:ny+1)
    real(dp), intent(in), optional :: n2(0:nx+1,0:ny+1)
    real(dp), intent(inout), optional :: p1(0:nx+1,0:ny+1)
    real(dp), intent(inout), optional :: p2(0:nx+1,0:ny+1)


    do j = 0, ny+1
       do i = 0, nx+1
          if (present(p1) .and. present(n1)) p1(i,j)  = n1(i,j) * kB * Tconst
          if (present(p2) .and. present(n2)) p2(i,j)  = n2(i,j) * kB * Tconst
       enddo
    enddo

  end subroutine update_primitive_variables

  subroutine advance_momentum(dt)
    implicit none
    integer(ii) :: i, j
    real(dp), intent(in) :: dt
    real(dp):: force_x, force_y

    force_x = 0.0_dp
    force_y = 0.0_dp
    ! Advance x-momentum
    call advance_conserved_quantity_MUSCL(dt, rhou_e, ue, ve, rhouu_e, rhouv_e)
    ! Advance y-momentum
    call advance_conserved_quantity_MUSCL(dt, rhov_e, ue, ve, rhovu_e, rhovv_e)

    ! Add isothermal source: Lorentz force and pressure gradient
    call compute_gradient(dpedx,dpedy,pe)
    do j = 1, ny
       do i = 1, nx
          force_x = (qe/me) * rho_e(i,j) * Ex(i,j)
          force_y = (qe/me) * rho_e(i,j) * Ey(i,j)

          ! Update momentum
          rhou_e(i,j) = rhou_e(i,j) + dt*(force_x - dpedx(i,j))
          rhov_e(i,j) = rhov_e(i,j) + dt*(force_y - dpedy(i,j))
       enddo
    enddo

    ! Update boundaries (Neumann)
    do i = 0, nx+1
       rhou_e(i,0    ) = rhou_e(i,1)
       rhou_e(i,ny+1) = rhou_e(i,ny)
       rhov_e(i,0    ) = rhov_e(i,1)
       rhov_e(i,ny+1) = rhov_e(i,ny)
    end do
    do j = 0, ny+1
       rhou_e(0    ,j) = rhou_e(1,j)
       rhou_e(nx+1,j) = rhou_e(nx,j)
       rhov_e(0    ,j) = rhov_e(1,j)
       rhov_e(nx+1,j) = rhov_e(nx,j)
    end do

    ! Compute velocities from updated momentum
    ue = rhou_e/(rho_e + 1d-9)
    ve = rhov_e/(rho_e + 1d-9)

  end subroutine advance_momentum

  subroutine advance_conserved_quantity_upwind(dt, cq, velx, vely, mdotx, mdoty)
    implicit none
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: cq(0:nx+1,0:ny+1)
    real(dp), intent(inout) :: mdotx(0:nx+1,0:ny+1), mdoty(0:nx+1,0:ny+1)
    real(dp), intent(in) :: velx(0:nx+1,0:ny+1), vely(0:nx+1,0:ny+1)

    real(dp), allocatable :: cq_old(:,:), flux_u(:,:), flux_v(:,:)
    real(dp), allocatable :: u_face(:,:), v_face(:,:)
    integer(ii) :: i, j

    ! Save old cq and allocate fluxes
    allocate(cq_old(0:nx+1,0:ny+1)); cq_old = cq
    allocate(flux_u(0:nx,0:ny+1))
    allocate(flux_v(0:nx+1,0:ny))
    allocate(u_face(0:nx,0:ny+1))
    allocate(v_face(0:nx+1,0:ny))

    ! Compute face-centered velocities via central interpolation
    do j = 0, ny+1
       do i = 0, nx
          u_face(i,j) = 0.5_dp * (velx(i,j) + velx(i+1,j))
       end do
    end do

    do j = 0, ny
       do i = 0, nx+1
          v_face(i,j) = 0.5_dp * (vely(i,j) + vely(i,j+1))
       end do
    end do

    ! Compute face-centered mdots
    mdotx = cq_old * velx
    mdoty = cq_old * vely

    ! Compute fluxes using upwind scheme with face velocities
    do j = 0, ny+1
       do i = 0, nx
          if (u_face(i,j) >= 0.0_dp) then
             flux_u(i,j) = u_face(i,j) * cq_old(i,j)
          else
             flux_u(i,j) = u_face(i,j) * cq_old(i+1,j)
          end if
       end do
    end do

    do j = 0, ny
       do i = 0, nx+1
          if (v_face(i,j) >= 0.0_dp) then
             flux_v(i,j) = v_face(i,j) * cq_old(i,j)
          else
             flux_v(i,j) = v_face(i,j) * cq_old(i,j+1)
          end if
       end do
    end do

    ! Update cell-centered cq via divergence of fluxes
    do j = 1, ny
       do i = 1, nx
          cq(i,j) = cq_old(i,j) &
               - (dt/dx) * (flux_u(i,j) - flux_u(i-1,j)) &
               - (dt/dy) * (flux_v(i,j) - flux_v(i,j-1))
       end do
    end do

    ! Apply zero-gradient (Neumann) BCs
    do i = 0, nx+1
       cq(i,0    ) = cq(i,1)
       cq(i,ny+1) = cq(i,ny)
    end do
    do j = 0, ny+1
       cq(0    ,j) = cq(1,j)
       cq(nx+1,j) = cq(nx,j)
    end do

    ! Clean up
    deallocate(cq_old, flux_u, flux_v, u_face, v_face)
  end subroutine advance_conserved_quantity_upwind

  subroutine advance_conserved_quantity_MUSCL(dt, cq, velx, vely, mdotx, mdoty)
    implicit none
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: cq(0:nx+1,0:ny+1)
    real(dp), intent(inout) :: mdotx(0:nx+1,0:ny+1), mdoty(0:nx+1,0:ny+1)
    real(dp), intent(in) :: velx(0:nx+1,0:ny+1), vely(0:nx+1,0:ny+1)

    real(dp), allocatable :: cq_old(:,:), flux_u(:,:), flux_v(:,:)
    real(dp), allocatable :: slope_x(:,:), slope_y(:,:)
    real(dp), allocatable :: u_face(:,:), v_face(:,:)
    real(dp) :: rL, rR
    integer(ii) :: i, j

    allocate(cq_old(0:nx+1,0:ny+1)); cq_old = cq
    allocate(slope_x(0:nx+1,0:ny+1), slope_y(0:nx+1,0:ny+1))
    allocate(flux_u(0:nx,  0:ny+1))
    allocate(flux_v(0:nx+1,0:ny  ))
    allocate(u_face(0:nx,0:ny+1), v_face(0:nx+1,0:ny))

    ! Interpolate face velocities
    do j = 0, ny+1
       do i = 0, nx
          u_face(i,j) = 0.5_dp * (velx(i,j) + velx(i+1,j))
       end do
    end do
    do j = 0, ny
       do i = 0, nx+1
          v_face(i,j) = 0.5_dp * (vely(i,j) + vely(i,j+1))
       end do
    end do

    ! Compute limited slopes using minmod limiter
    do j = 0, ny+1
       do i = 1, nx
          rL = (cq_old(i,j) - cq_old(i-1,j)) / dx
          rR = (cq_old(i+1,j) - cq_old(i,j)) / dx
          slope_x(i,j) = 0.5_dp * (sign(1.0_dp, rL) + sign(1.0_dp, rR)) * min(abs(rL), abs(rR))
       end do
    end do

    do j = 1, ny
       do i = 0, nx+1
          rL = (cq_old(i,j) - cq_old(i,j-1)) / dy
          rR = (cq_old(i,j+1) - cq_old(i,j)) / dy
          slope_y(i,j) = 0.5_dp * (sign(1.0_dp, rL) + sign(1.0_dp, rR)) * min(abs(rL), abs(rR))
       end do
    end do

    ! Compute fluxes using MUSCL reconstruction
    do j = 0, ny+1
       do i = 0, nx
          if (i == 0 .or. i == nx) then
             rL = cq_old(i,j)
             rR = cq_old(i+1,j)
          else
             rL = cq_old(i,j) + 0.5_dp*dx*slope_x(i,j)
             rR = cq_old(i+1,j) - 0.5_dp*dx*slope_x(i+1,j)
          end if

          if (u_face(i,j) >= 0.0_dp) then
             flux_u(i,j) = u_face(i,j) * rL
          else
             flux_u(i,j) = u_face(i,j) * rR
          end if
       end do
    end do

    do j = 0, ny
       do i = 0, nx+1
          if (j == 0 .or. j == ny) then
             rL = cq_old(i,j)
             rR = cq_old(i,j+1)
          else
             rL = cq_old(i,j) + 0.5_dp*dy*slope_y(i,j)
             rR = cq_old(i,j+1) - 0.5_dp*dy*slope_y(i,j+1)
          end if

          if (v_face(i,j) >= 0.0_dp) then
             flux_v(i,j) = v_face(i,j) * rL
          else
             flux_v(i,j) = v_face(i,j) * rR
          end if
       end do
    end do

    ! Update cell-centered value
    do j = 1, ny
       do i = 1, nx
          cq(i,j) = cq_old(i,j) &
               - (dt/dx)*(flux_u(i,j) - flux_u(i-1,j)) &
               - (dt/dy)*(flux_v(i,j) - flux_v(i,j-1))
       end do
    end do

    ! Apply zero-gradient (Neumann) BCs
    do i = 0, nx+1
       cq(i,0    ) = cq(i,1)
       cq(i,ny+1) = cq(i,ny)
    end do
    do j = 0, ny+1
       cq(0    ,j) = cq(1,j)
       cq(nx+1,j) = cq(nx,j)
    end do

    ! Clean up
    deallocate(cq_old, slope_x, slope_y, flux_u, flux_v, u_face, v_face)
  end subroutine advance_conserved_quantity_MUSCL

  subroutine solve_poisson_SOR(f, rhsVal, tol, max_iter, omega)
    implicit none
    real(dp), intent(inout) :: f(0:nx+1,0:ny+1)
    real(dp), intent(in)    :: rhsVal(0:nx+1,0:ny+1)
    real(dp), intent(in)    :: tol, omega
    integer(ii), intent(in) :: max_iter
    integer(ii) :: i, j, iter
    real(dp) :: dx2, dy2, denom, f_old, frac
    real(dp) :: diff, max_diff

    dx2   = dx*dx
    dy2   = dy*dy
    denom = 2.0*(dx2 + dy2)
    frac    = omega/denom

    do iter = 1, max_iter

       max_diff = 0.0_dp

!!$       ! Zero-gradient BC on top/bottom
!!$       do i = 0, nx+1
!!$          f(i,   0) = f(i,   1)    ! bottom ghost row
!!$          f(i,ny+1) = f(i,  ny)    ! top ghost row
!!$       end do

       ! Single lexicographic sweep
       do j = 1, ny
          do i = 1, nx
             f_old    = f(i,j)
             f(i,j) = (1.0-omega)*f_old + frac * ( &
                  (f(i+1,j) + f(i-1,j))*dy2 + &
                  (f(i,j+1) + f(i,j-1))*dx2 - &
                  dx2*dy2*rhsVal(i,j) )
             diff = abs(f(i,j) - f_old)
             if (diff > max_diff) max_diff = diff
          end do
       end do

       if (max_diff < tol) exit
    end do
    print *, 'Poisson SOR: iterations =', iter, 'residual =', max_diff
  end subroutine solve_poisson_SOR

  subroutine solve_poisson_RBSOR(f, rhsVal, tol, max_iter, omega)
    implicit none
    real(dp), intent(inout) :: f(0:nx+1,0:ny+1)
    real(dp), intent(in)    :: rhsVal(0:nx+1,0:ny+1)
    real(dp), intent(in)    :: tol, omega
    integer(ii), intent(in) :: max_iter
    integer(ii) :: i, j, iter
    real(dp) :: dx2, dy2, denom, f_old, diff, max_diff

    dx2 = dx*dx
    dy2 = dy*dy
    denom = 2.0*(dx2 + dy2)

    do iter = 1, max_iter
       max_diff = 0.0


       ! Zero-gradient BC on top/bottom
       do i = 0, nx+1
          f(i,   0) = f(i,   1)    ! bottom ghost row
          f(i,ny+1) = f(i,  ny)    ! top ghost row
       end do

       ! Red-black SOR
       do j = 1, ny
          do i = 1 + mod(j,2), nx, 2
             f_old = f(i,j)
             f(i,j) = (1.0-omega)*f_old + (omega/denom)*( &
                  (f(i+1,j)+f(i-1,j))*dy2 + (f(i,j+1)+f(i,j-1))*dx2 - dx2*dy2*rhsVal(i,j) )
             diff = abs(f(i,j) - f_old)
             if (diff > max_diff) max_diff = diff
          end do
       end do
       do j = 1, ny
          do i = 1 + 1-mod(j,2), nx, 2
             f_old = f(i,j)
             f(i,j) = (1.0-omega)*f_old + (omega/denom)*( &
                  (f(i+1,j)+f(i-1,j))*dy2 + (f(i,j+1)+f(i,j-1))*dx2 - dx2*dy2*rhsVal(i,j) )
             diff = abs(f(i,j) - f_old)
             if (diff > max_diff) max_diff = diff
          end do
       end do
       if (max_diff < tol) exit
       print *, 'Poisson SOR: iterations =', iter, 'residual =', max_diff
    end do

  end subroutine solve_poisson_RBSOR

  subroutine compute_gradient(fx, fy, f)
    implicit none
    integer(ii) :: i, j
    real(dp), intent(inout) :: fx(0:nx,0:ny+1), fy(0:nx+1,0:ny)
    real(dp), intent(in) :: f(0:nx+1,0:ny+1)

    do j = 0, ny+1
       do i = 0, nx
          fx(i,j) = (f(i+1,j) - f(i,j))/dx
       end do
    end do

    do j = 0, ny
       do i = 0, nx+1
          fy(i,j) = (f(i,j+1) - f(i,j))/dy
       end do
    end do

  end subroutine compute_gradient

  subroutine compute_velocity()
    implicit none
    integer(ii) :: i, j

    call compute_gradient(dnedx,dnedy,ne)
    call compute_gradient(dnidx,dnidy,ni)
    !===============================
    ! Interior: compute center velocities from E-field
    !===============================
    do j = 1, ny
       do i = 1, nx
          ue(i,j) = mu_e * 0.5_dp * (Ex(i,j) + Ex(i-1,j))
          ue(i,j) = ue(i,j) - diff_e/ne(i,j)*dnedx(i,j)

          ve(i,j) = mu_e * 0.5_dp * (Ey(i,j) + Ey(i,j-1))
          ve(i,j) = ve(i,j) - diff_e/ne(i,j)*dnedy(i,j)

          ui(i,j) = mu_i * 0.5_dp * (Ex(i,j) + Ex(i-1,j))
          ui(i,j) = ui(i,j) - diff_i/ni(i,j)*dnidx(i,j)

          vi(i,j) = mu_i * 0.5_dp * (Ey(i,j) + Ey(i,j-1))
          vi(i,j) = vi(i,j) - diff_i/ni(i,j)*dnidy(i,j)
       end do
    end do

    call apply_symmetry_bc_velocity(ue, ve)
    call apply_symmetry_bc_velocity(ui, vi)

  end subroutine compute_velocity

  subroutine apply_symmetry_bc_velocity(u, v)
    implicit none
    real(dp), intent(inout) :: u(0:nx+1, 0:ny+1)
    real(dp), intent(inout) :: v(0:nx+1, 0:ny+1)
    integer :: i, j

    ! Left and right boundaries
    do j = 1, ny
       u(0,j)     = 0.0_dp
       u(nx+1,j)  = 0.0_dp
       v(0,j)     = v(1,j)
       v(nx+1,j)  = v(nx,j)
    end do

    ! Bottom and top boundaries
    do i = 1, nx
       v(i,0)     = 0.0_dp                ! normal
       v(i,ny+1)  = 0.0_dp
       u(i,0)     = u(i,1)                ! tangential
       u(i,ny+1)  = u(i,ny)
    end do

    ! Corner extrapolation
    u(0,0)       = 0.0_dp
    u(0,ny+1)    = 0.0_dp
    u(nx+1,0)    = 0.0_dp
    u(nx+1,ny+1) = 0.0_dp

    v(0,0)       = v(1,1)
    v(0,ny+1)    = v(1,ny)
    v(nx+1,0)    = v(nx,1)
    v(nx+1,ny+1) = v(nx,ny)
  end subroutine apply_symmetry_bc_velocity

end module m_solver
