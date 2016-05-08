! Test for "rkf45.f" solver

! Solve a dimension-3 problem on 5 time-points.
!
!   x'_i(t) = i x_i(t)
!   x_i(0) = i
!
! Use contained derivative function.
! Use no module.

program main

  implicit none

  interface
    subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr,abserr,work(1)
      external f
    end subroutine rkf45
  end interface

  integer &
    :: ind

  integer , parameter &
    :: dim_ode = 3
  integer , parameter &
    :: dim_time = 5
  double precision , dimension(dim_time) &
    :: time_point = (/ ( ind , ind=0,(dim_time-1) ) /)
  double precision &
    :: tol = 1e-8

  double precision , dimension(dim_ode,dim_time) &
    :: ode_res

  double precision , dimension(dim_ode) &
    :: state
  double precision &
    :: time
  integer &
    :: iflag = 1
  double precision , dimension(3+6*dim_ode) &
    :: work
  integer , dimension(5) &
    :: iwork

  state = (/ ( ind , ind=1,dim_ode ) /)
  time = time_point(1)

  do ind = 1 , dim_time!{{{

    call rkf45 ( &
      derivative &
      , dim_ode &
      , state &
      , time &
      , time_point(ind) &
      , tol &
      , tol &
      , iflag &
      , work &
      , iwork &
    )

    ode_res(:,ind) = state
    time = time_point(ind)

  end do!}}}

  write(*,*) ode_res

  stop

  contains

! Derivative: y'_i = i * y_i

  subroutine derivative ( t , y , d )!{{{

    implicit none

    double precision , intent(in) &
      :: t
    double precision , intent(in) , dimension(:) &
      :: y
    double precision , intent(out) , dimension(:) &
      :: d

    d = (/ ( ind*y(ind) , ind=1,size(y) ) /)

    return

  end subroutine derivative!}}}

end program main
