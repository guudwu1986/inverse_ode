! Test program for "ResidualSumOfSquares" of module
! "Linear_Ode_Inverse_Eigen_mod"

! Start from true eigenvalues to reconstruct the system,
! the calculate the discrepancy between the system curve and
! true observation,
! which should be close to 0.

! Or use 0-eigenvalues, with 2 time-points,
! then uncomment "zero_rss" block in "ResidualSumOfSquares"
! should give two close values.

program main

  use Linear_Ode_Inverse_Eigen_mod

  implicit none

  integer &
    :: NUM_UNIT

  integer &
    :: dim_ode
  integer &
    :: dim_time
  double precision , dimension(:) , allocatable &
    :: timepoint
  double precision , dimension(:) , allocatable &
    :: linear
  double precision , dimension(:) , allocatable &
    :: eigen
  double precision , dimension(:) , allocatable &
    :: initial
  double precision , dimension(:) , allocatable &
    :: observation
  double precision &
    :: ridge_parameter = 1e-5
  double precision &
    :: rss

! Read!{{{

  open ( newunit = NUM_UNIT , file = 'data10' ,&
    action = 'read' , status = 'old' )

  read ( NUM_UNIT , * ) dim_ode
  read ( NUM_UNIT , * ) dim_time

  allocate ( timepoint(dim_time) )
  read ( NUM_UNIT , * ) timepoint

  allocate ( observation(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) observation

  allocate ( initial(dim_ode) )
  read ( NUM_UNIT , * ) initial

  allocate ( linear(dim_ode*dim_ode) )
  read ( NUM_UNIT , * ) linear

  allocate ( eigen(dim_ode) )
  read ( NUM_UNIT , * ) eigen ( 1 : ((dim_ode+1)/2) )
  read ( NUM_UNIT , * ) eigen ( ((dim_ode+1)/2+1) : dim_ode )

  close ( unit = NUM_UNIT )
!}}}

!  eigen = 0

  call ResidualSumOfSquares &
    ( &
      dim_ode &
      , dim_time &
      , eigen &
      , timepoint &
      , observation &
      , ridge_parameter &
      , rss &
    )

  write(*,*) rss

  stop

end program main
