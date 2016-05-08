program main

  implicit none

  integer , parameter &
    :: dim_fun = 5

  double precision , dimension(dim_fun) &
    :: x

  integer , parameter &
    :: num_point = 2*dim_fun+1
  double precision &
    :: rhobeg = 2
  double precision &
    :: rhoend = 1e-8
  integer &
    :: iprint = 0
  integer &
    :: maxfun = 1000
  double precision , &
    dimension((num_point+13)*(num_point+dim_fun)+3*dim_fun*(dim_fun+3)/2) &
    :: w

  interface
    SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),W(*)
      external CALFUN
    end SUBROUTINE NEWUOA
  end interface

  x = 0

  call newuoa &
    ( &
      dim_fun &
      , num_point &
      , x &
      , rhobeg &
      , rhoend &
      , iprint &
      , maxfun &
      , w &
      , calfun &
    )

  write(*,*) x

  stop

  contains

! Testing function:
! F = \sum (x_i - i*x_{i-1})^2 + (x_1-1)^2
! Optimal solution: (1,2,6,24,120,...)

  subroutine calfun &
    ( &
      n , x , f &
    )

    implicit none

    integer , intent(in) &
      :: n
    double precision , intent(in) , dimension(n) &
      :: x
    double precision , intent(out) &
      :: f

    integer &
      :: ind

    f = ( x(1)-1 ) ** 2

    do ind = 2 , n
      f = f + (x(ind)-ind*x(ind-1))**2
    end do

  end subroutine calfun

end program main
