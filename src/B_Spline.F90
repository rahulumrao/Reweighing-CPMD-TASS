MODULE B_Spline
CONTAINS
SUBROUTINE bspline(nr,gridmin,gridmax,pcons,fes)
!*****************************************************************************************
    USE bspline_module
    USE bspline_kinds_module, only: wp, ip
!    USE pyplot_module

    IMPLICIT NONE

    integer(ip) :: nx    !! number of points in x
    integer(ip) :: nxv   !! number of points to evaluate interpolant
    integer(ip),parameter :: kx    = 4    !! order in x
    integer(ip),parameter :: iknot = 0    !! automatically select the knots

    real(wp),ALLOCATABLE :: x(:),f1(:),fval(:)
    real(wp),ALLOCATABLE :: xval(:)
    real(wp),ALLOCATABLE :: tx(:)
    real(wp) :: val,tru,err,errmax,gridmin(*),gridmax(*)
    integer(ip) :: i,j,idx,iflag,inbvx,iloy
    logical :: extrap
!    type(pyplot) :: plt
    integer :: istat  !! pyplot-fortran status flag
    real(wp),dimension(3*kx) :: w1_1d !! work array
    INTEGER :: n,nr,ios
    REAL*8  :: pcons(*),fes(*)
    REAL*8 :: dummy,dummy1,nb

    idx = 0
    nx = nr ; nxv = 15*(nr)
    ALLOCATE(x(nx))
    ALLOCATE(xval(nxv))
    ALLOCATE(tx(nx+kx))
    ALLOCATE(f1(nx))
    ALLOCATE(fval(nxv))
    
    DO i = 1,nr
      x(i)  = pcons(i)
      f1(i) = fes(i)
!      WRITE(*,*)pcons(i),fes(i)
    ENDDO

    !nb = (gridmax(1) - gridmin(1))/nxv
    nb = (pcons(nr) - pcons(2))/nxv
    DO i = 1,nxv
      xval(i) = pcons(2)+nb*i
    ENDDO
    !have to set these before the first evaluate call:
    inbvx = 1
    iloy  = 1
    ! initialize
    call db1ink(x,nx,f1,kx,iknot,tx,f1,iflag)

    if (iflag/=0) then
        write(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
    end if

    !initialize the plot:
!    call plt%initialize(grid=.true.,xlabel='x (deg)',ylabel='f(x)',&
!                        title='Extrapolation Test',legend=.true.)
!    call plt%add_plot(x,f1,label='free_energy$',&
!                        linestyle='ko',markersize=5,linewidth=2,istat=istat)

        OPEN(21,FILE='interp_free_energy.out')
        errmax = 0.0_wp
        do i=1,nxv
            call db1val(xval(i),idx,tx,nx,kx,f1,val,iflag,inbvx,w1_1d,extrap=extrap)
!            write(*,*) xval(i), val
            IF (xval(i) .lt. pcons(1)) val = 0.0
            fval(i) = val  ! save it for plot
            WRITE(21,*)xval(i), val
        end do

        write(*,*) ''
        write(*,*) 'interpolated free energy written in : interp_free_energy.dat'
        write(*,*) ''

!        if (extrap) then
!            call plt%add_plot(xval,fval,&
!                    label='Interpolated',&
!                    linestyle='g.-',linewidth=1,istat=istat)
!            call plt%savefig('bspline_extrap_test.png',istat=istat)
!        end if


DEALLOCATE(x)
DEALLOCATE(xval)
DEALLOCATE(f1)
CLOSE(21)
END
END MODULE B_Spline
