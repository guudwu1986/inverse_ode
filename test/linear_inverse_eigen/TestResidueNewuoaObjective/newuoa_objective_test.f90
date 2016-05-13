! Test program for "ResidueNewuoaObjective" in module
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
    :: ridge_parameter = 1e-5

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
  read ( NUM_UNIT , * ) scaling ( 1 : (dim_ode/2) )
  read ( NUM_UNIT , * ) scaling ( (dim_ode/2+1) : dim_ode )

  close ( unit = NUM_UNIT )
!}}}

  call TestResidueNewuoaObjective &
    ( &
      dim_ode &
      , dim_time &
      , eigen &
      , scaling &
      , timepoint &
      , observation &
      , ridge_parameter &
    )

  stop

end program main
