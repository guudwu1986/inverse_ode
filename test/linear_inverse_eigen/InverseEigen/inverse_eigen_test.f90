! Test program for "InverseEigen" in module
! "Linear_Ode_Inverse_Eigen_mod"

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
    :: scaling
  double precision , dimension(:) , allocatable &
    :: initial
  double precision , dimension(:) , allocatable &
    :: observation
  double precision &
    :: ridge_parameter = 1e-3
  double precision &
    :: rhobeg = 2e-1
  double precision &
    :: rhoend = 1e-1
  integer &
    :: maxfun = 100000

  double precision , dimension(:) , allocatable &
    :: linear_e
  double precision , dimension(:) , allocatable &
    :: initial_e
  double precision , dimension(:) , allocatable &
    :: curve_e

  integer &
    :: t1
  integer &
    :: t2

! Read!{{{

  open ( newunit = NUM_UNIT , file = 'data' ,&
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
  eigen = 1
  allocate ( scaling(dim_ode) )
  read ( NUM_UNIT , * ) scaling ( 1 : ((dim_ode+1)/2) )
  read ( NUM_UNIT , * ) scaling ( ((dim_ode+1)/2+1) : dim_ode )

  close ( unit = NUM_UNIT )
!}}}

  allocate ( linear_e ( dim_ode*dim_ode ) )
  allocate ( initial_e ( dim_ode ) )
  allocate ( curve_e ( dim_ode*dim_time ) )

  call system_clock(t1)

  call InverseEigen &
    ( &
      dim_ode &
      , dim_time &
      , eigen &
      , scaling &
      , timepoint &
      , observation &
      , ridge_parameter &
      , rhobeg &
      , rhoend &
      , maxfun &
      , linear_e &
      , initial_e &
      , curve_e &
    )

  call system_clock(t2)

!  write(*,*) maxval ( abs ( linear_e - linear ) )
!  write(*,*) maxval ( abs ( initial_e - initial ) )
!  write(*,*) maxval ( abs ( curve_e - observation ) )

  write(*,*) sum((linear_e-linear)**2) / sum(linear**2)
  write(*,*) sum((curve_e-observation)**2) / sum(observation**2)

  write(*,*) t2-t1

!  open ( newunit = NUM_UNIT , file = 'regenerate_c' ,&
!    action = 'write' )
!  write(NUM_UNIT,*) linear_e
!  close ( unit = NUM_UNIT )

!  open ( newunit = NUM_UNIT , file = 'regenerate_i' ,&
!    action = 'write' )
!  write(NUM_UNIT,*) curve_e(1:dim_ode)
!  close ( unit = NUM_UNIT )

  stop

end program main
