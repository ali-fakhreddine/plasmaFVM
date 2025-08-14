!=======================================================================
!  File       : m_dataio.f90
!  Author     : Ali Fakhreddine
!  Created    : 14-07-2025
!  Description: I/O routines for writing simulation data and grid
!               geometry in VTK format + directory management.
!=======================================================================
module m_dataio
  ! Use associations
  use m_precision
  use m_grid
  implicit none

contains

  subroutine write_point_data(unit)
    implicit none
    integer(ii), intent(in) :: unit
    integer(ii) :: npts

    !npts = (nx+1)*(ny+1)!*(nz+1)
    npts = (nx)*(ny)!*(nz)
    write(unit,'(A,I0)') 'POINT_DATA ',npts
    ! Write registered fields below
    call write_scalar(unit, 'ne', ne)
    call write_scalar(unit, 'ni', ni)
    call write_scalar(unit, 'phi', phi)
    call write_vector_face_field(unit, 'E', Ex, Ey)
    call write_vector_cv_field(unit,'electron_velocity', ue, ve)
    call write_vector_cv_field(unit,'ion_velocity', ui, vi)
  end subroutine write_point_data
  
  subroutine ensure_output_dir()
    implicit none
    logical(ll) :: data_dir_exists
    logical(ll) :: grid_dir_exists
    logical(ll) :: restart_dir_exists

    data_dir_exists = .FALSE.
    inquire(file='data_files', exist=data_dir_exists)
    if (.not. data_dir_exists) then
       call execute_command_line('mkdir -p data_files')
    end if

    grid_dir_exists = .FALSE.
    inquire(file='grid_files', exist=grid_dir_exists)
    if (.not. grid_dir_exists) then
       call execute_command_line('mkdir -p grid_files')
    end if

    restart_dir_exists = .FALSE.
    inquire(file='restart_files', exist=restart_dir_exists)
    if (.not. restart_dir_exists) then
       call execute_command_line('mkdir -p restart_files')
    end if

  end subroutine ensure_output_dir

  subroutine write_vtk_header(unit, title)
    implicit none
    integer(ii), intent(in) :: unit
    character(len=*), intent(in) :: title
    integer(ii) :: npts

    !npts = (nx+1)*(ny+1)!*(nz+1)
    npts = (nx)*(ny)!*(nz)

    write(unit,'(A)') '# vtk DataFile Version 3.0'
    write(unit,'(A)') trim(title)
    write(unit,'(A)') 'ASCII'
    write(unit,'(A)') 'DATASET STRUCTURED_GRID'
    write(unit,'(A,I0,1X,I0,1X,I0)') 'DIMENSIONS ',nx,ny,1!nz+1
    write(unit,'(A,I0,1X,A)') 'POINTS ',npts,' double'
  end subroutine write_vtk_header

  subroutine write_grid_points(unit)
    implicit none
    integer(ii), intent(in) :: unit
    integer(ii) :: i, j, k


    !do k = 0,nz
       do j = 1,ny
          do i = 1,nx
             write(unit,'(3F14.6)') xc(i), yc(j), zc(0)!zf(k)
          enddo
       enddo
    !enddo
  end subroutine write_grid_points

  subroutine write_scalar(unit, name, field)
    implicit none
    integer(ii), intent(in) :: unit
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: field(0:nx+1,0:ny+1)
    integer(ii) :: i, j, k

    write(unit,'(A,1X,A,1X,A)') 'SCALARS', trim(name), 'double'
    write(unit,'(A)') 'LOOKUP_TABLE default'
    !do k = 0,nz
       do j = 1,ny
          do i = 1,nx
             write(unit,'(E14.6E3)') field(i,j)
          enddo
       enddo
    !enddo
  end subroutine write_scalar

  subroutine write_vector_cv_field(unit, name, fx, fy)
    implicit none
    integer(ii), intent(in) :: unit
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: fx(0:nx+1,0:ny+1), fy(0:nx+1,0:ny+1)
    integer(ii) :: i, j, k

    write(unit,'(A,1X,A,1X,A)') 'VECTORS', trim(name), 'double'
    !do k = 0, nz
       do j = 1, ny
          do i = 1, nx
             write(unit,'(3E14.6E3)') fx(i,j), fy(i,j), 0.0
          enddo
       enddo
    !enddo
  end subroutine write_vector_cv_field

  subroutine write_vector_face_field(unit, name, fx, fy)
    implicit none
    integer(ii), intent(in) :: unit
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: fx(0:nx,0:ny+1), fy(0:nx+1,0:ny)
    integer(ii) :: i, j, k

    write(unit,'(A,1X,A,1X,A)') 'VECTORS', trim(name), 'double'
    !do k = 0, nz
       do j = 1, ny
          do i = 1, nx
             write(unit,'(3E14.6E3)') 0.5_dp*(fx(i,j)+fx(i-1,j)), 0.5_dp*(fy(i,j)+fy(i,j-1)), 0.0d0
             !write(unit,'(3E14.6E3)') fx(i,j), fy(i,j), 0.0
          enddo
       enddo
    !enddo
  end subroutine write_vector_face_field

  subroutine write_solution_vtk(name)
    implicit none
    character(len=*), intent(in) :: name
    integer(ii) :: unit
    character(len=256) :: path

    call ensure_output_dir()
    path = trim('data_files/'//name)
    open(newunit=unit, file=path, status='replace', action='write', form='formatted')
    call write_vtk_header(unit, 'Solution data')
    call write_grid_points(unit)
    call write_point_data(unit)
    close(unit)
  end subroutine write_solution_vtk

  subroutine write_solution_vtk_basic(name)
    implicit none
    character(len=*), intent(in) :: name
    integer(ii) :: i,j,k, unit, npts
    character(len=256) :: path

    call ensure_output_dir()
    path = trim('data_files/'//name)
    open(newunit=unit, file=path, status='replace', action='write', form='formatted')
    npts = (nx+1)*(ny+1)*(nz+1)
    write(unit,'(A)') '# vtk DataFile Version 3.0'
    write(unit,'(A)') 'Solution data'
    write(unit,'(A)') 'ASCII'
    write(unit,'(A)') 'DATASET STRUCTURED_GRID'
    write(unit,'(A,I0,1X,I0,1X,I0)') 'DIMENSIONS ', nx+1, ny+1, nz+1
    write(unit,'(A,I0,1X,A)') 'POINTS ', npts, ' double'
    do k = 0, nz
       do j = 0, ny
          do i = 0, nx
             write(unit,'(3F14.6)') xf(i), yf(j), zf(k)
          end do
       end do
    end do
    write(unit,'(A,I0)') 'POINT_DATA ', npts
    write(unit,'(A)') 'SCALARS ne double'
    write(unit,'(A)') 'LOOKUP_TABLE default'
    do k = 0, nz
       do j = 0, ny
          do i = 0, nx
             write(unit,'(F14.6)') ne(i,j)
          end do
       end do
    end do
    close(unit)
  end subroutine write_solution_vtk_basic

  subroutine write_grid_vtk(name)
    implicit none
    character(len=*), intent(in) :: name
    integer(ii) :: i,j,k, unit, npts
    character(len=256) :: path

    call ensure_output_dir()
    path = 'grid_files/'//trim(name)
    open(newunit=unit, file=path, status='replace', action='write')
    npts = (nx+1)*(ny+1)*(nz+1)
    write(unit,'(A)') '# vtk DataFile Version 3.0'
    write(unit,'(A)') 'Grid geometry'
    write(unit,'(A)') 'ASCII'
    write(unit,'(A)') 'DATASET STRUCTURED_GRID'
    write(unit,'(A,I0,1X,I0,1X,I0)') 'DIMENSIONS ', nx+1, ny+1, nz+1
    write(unit,'(A,I0,A)') 'POINTS ', npts, ' double'
    do k = 0, nz
       do j = 0, ny
          do i = 0, nx
             write(unit,'(3F12.5)') xf(i), yf(j), zf(k)
          end do
       end do
    end do
    close(unit)
  end subroutine write_grid_vtk

end module m_dataio
