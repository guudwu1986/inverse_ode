! Test program for "ConstructODE" in module
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
  :: initial
double precision , dimension(:) , allocatable &
  :: observation
double precision &
  :: ridge_parameter = 1e-3
double precision , dimension(:) , allocatable &
  :: linear_e
double precision , dimension(:) , allocatable &
  :: initial_e
double precision , dimension(:) , allocatable &
  :: curve_e
integer &
  :: info

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
read ( NUM_UNIT , * ) eigen ( 1 : ((dim_ode+1)/2) )
read ( NUM_UNIT , * ) eigen ( ((dim_ode+1)/2+1) : dim_ode )

close ( unit = NUM_UNIT )
!}}}

allocate ( linear_e(dim_ode*dim_ode) )
allocate ( initial_e(dim_ode) )
allocate ( curve_e(dim_ode*dim_time) )

call ConstructODE &
( &
  dim_ode &
  , dim_time &
  , eigen &
  , timepoint &
  , observation &
  , ridge_parameter &
  , linear_e &
  , initial_e &
  , curve_e &
  , info &
)

write(*,*) info
write(*,*) maxval ( abs ( linear_e - linear ) )
write(*,*) maxval ( abs ( initial_e - initial ) )
write(*,*) maxval ( abs ( curve_e - observation ) )

stop

end program main
