!=======================================================================
!  File       : m_grid.f90
!  Author     : Ali Fakhreddine
!  Created    : 14-07-2025
!  Description: Provides routines to construct the computational grid—
!               computing cell-face and cell-center coordinates—
!               and to perform basic grid extent checks.
!=======================================================================
module m_grid
  ! Use associations
  use m_precision
  use m_varspace
  implicit none

contains

  subroutine create_grid
    implicit none
    integer(ii) :: i, j, k
    character(len=60) :: sep
    sep = repeat('-', len(sep))

    dx = Lx/nx; dy = Ly/ny;  dz = Lz/nz

    ! Face positions
    do i = 0, nx
       xf(i) = i * dx
    enddo
    do j = 0, ny
       yf(j) = j * dy
    enddo
    do k = 0, nz
       zf(k) = k * dz
    enddo

    ! Cell-centered positions
    do i = 1, nx
       xc(i) = 0.5_dp * (xf(i-1) + xf(i))
    enddo
    do j = 1, ny
       yc(j) = 0.5_dp * (yf(j-1) + yf(j))
    enddo
    do k = 1, nz
       zc(k) = 0.5_dp * (zf(k-1) + zf(k))
    enddo

    ! Grid check
    write(*,'(A)') sep
    write(*,'(A)') '      G R I D   E X T E N T   C H E C K      '
    write(*,'(A)') sep

    ! Number of cells
    write(*,'(A,I0,A,I0,A,I0)') &
         ' Grid: ', nx, ' x ', ny, ' x ', nz

    ! Domain extents
    write(*,'(A,F10.4,A,F10.4,A,F10.4,A)') &
         ' Domain: [0,', Lx, '] x [0,', Ly, '] x [0,', Lz, ']'

    ! Spacing
    write(*,'(A,F10.4,A,F10.4,A,F10.4)') &
         ' Spacing: dx =', dx, ',  dy =', dy, ',  dz =', dz
    write(*,'(A)') sep

  end subroutine create_grid
end module m_grid
