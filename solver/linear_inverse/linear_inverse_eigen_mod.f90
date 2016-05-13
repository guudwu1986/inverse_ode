! A module to solve a linear-ODE inverse problem
! via eigenvalue-updating scheme.

! Used by module Linear_Ode_Inverse_mod

! ResidualSumOfSquares:
!   A private subroutine which computes the residue between observation
!   and system curves when eigenvalues of the system is given.
! ResidueNewuoaObjective:
!   A private subroutine which is a wrapper for "ResidualSumOfSquares",
!   to be called by NEWUOA-solver.
! TestResidueNewuoaObjective:
!   A wrapper for "ResidueNewuoaObjective" testing.
!   Should be commented.
! ConstructODE:
!   A private subroutine to construct an ODE system
!   from given eigenvalues.
! InverseEigen:
!   A public subroutine to compute an ODE inverse problem
!   via eigenvalue-updating.
! InverseEigenLasso:
! ResidualSumOfSquaresLasso:
! ResidueNewuoaObjectiveLasso:
!   With LASSO penalty term.

module Linear_Ode_Inverse_Eigen_mod

  implicit none
  private

!  public &
!    :: ResidualSumOfSquares
!  public &
!    :: ResidueNewuoaObjective
!  public &
!    :: ConstructODE
  public &
    :: InverseEigen
  public &
    :: InverseEigenLasso

!  public &
!    :: TestResidueNewuoaObjective

  integer &
    :: m_dim_time
  double precision &
    :: m_ridge_parameter

  double precision , dimension(:) , pointer &
    :: mp_timepoint
  double precision , dimension(:) , pointer &
    :: mp_observation
  double precision , dimension(:) , pointer &
    :: mp_scaling
  double precision , dimension(:) , pointer &
    :: mp_lasso

  contains

subroutine ResidualSumOfSquares &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Rss &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(in) &
    :: Eigen
  double precision , dimension(Dim_Time) , intent(in) &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , intent(out) &
    :: Rss

  integer &
    :: ind
  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: basis
  double precision , &
    dimension ( (Dim_Ode+Dim_Time) , Dim_Ode ) &
    :: q
!  double precision &
!    :: zero_rss

  double precision , dimension(:) , allocatable &
    :: work
  integer &
    :: lwork
  double precision , dimension(Dim_Ode) &
    :: tau
  integer &
    :: info

  integer &
    :: itemp
  double precision &
    :: temp

  double precision , dimension(:,:) , pointer &
    :: p_observation

  p_observation ( 1:Dim_Ode , 1:Dim_Time ) => Observation

! Generate basis:
!   Each column is a basis.
!   A diagonal square matrix is binded to the bottom,
!   used for ridge regression.

  basis = 0

  itemp = (Dim_Ode+1)/2
  allocate ( work(Dim_Time) )
  do ind = 1 , Dim_Ode/2
    work = dexp ( Eigen(ind) * Timepoint )
    basis ( 1:Dim_Time , 2*ind-1 ) &
      = work * dsin ( Eigen(ind+itemp) * Timepoint )
    basis ( 1:Dim_Time , 2*ind ) &
      = work * dcos ( Eigen(ind+itemp) * Timepoint )
  end do
  if ( mod(Dim_Ode,2) .eq. 1 ) then
    basis ( 1:Dim_Time , Dim_Ode ) &
      = dexp ( Eigen(itemp) * Timepoint )
  end if
  deallocate ( work )

  temp = sum ( basis**2 )
  temp = sqrt(temp)
  do ind = 1 , Dim_Ode
    basis ( ind+Dim_Time , ind ) = Ridge_Parameter * temp
  end do

! QR-decompose basis

  allocate ( work(1) )
  lwork = -1
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info &
  )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info &
  )
  deallocate ( work )

! Reconstruct Q

  q = basis

  allocate ( work(1) )
  lwork = -1
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info &
    )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info &
    )
  deallocate ( work )

! Compute RSS:
!   Residue is ||Y||^2 - ||Q^T Y||^2.

  Rss = sum ( p_observation**2 ) &
    - sum ( matmul ( p_observation , q(1:Dim_Time,:) ) **2 )

! Uncomment following block to check
!   cases with 0-eigenvalue with 2 time-points.

!  zero_rss = 0
!  do ind = 1 , Dim_Ode
!    zero_rss = zero_rss &
!      + sum ( &
!          ( p_observation(ind,:) &
!            - sum(p_observation(ind,:))/Dim_Time &
!          ) **2 &
!        )
!  end do
!  write(*,*) zero_rss

  return

end subroutine ResidualSumOfSquares!}}}

subroutine ResidueNewuoaObjective &!{{{
  ( &
    N , X , F &
  )

  implicit none

  integer , intent(in) :: N
  double precision , intent(in) , dimension(N) :: X
  double precision , intent(out) :: F

  call ResidualSumOfSquares &
    ( &
      N &
      , m_dim_time &
      , X * mp_scaling &
      , mp_timepoint &
      , mp_observation &
      , m_ridge_parameter &
      , F &
    )

  return

end subroutine ResidueNewuoaObjective!}}}

!  subroutine TestResidueNewuoaObjective &!{{{
!    ( &
!      Dim_Ode &
!      , Dim_Time &
!      , Eigen &
!      , Scaling &
!      , Timepoint &
!      , Observation &
!      , Ridge_Parameter &
!    )
!
!    implicit none
!
!    integer , intent(in) &
!      :: Dim_Ode
!    integer , intent(in) &
!      :: Dim_Time
!    double precision , dimension(Dim_Ode) , intent(inout) &
!      :: Eigen
!    double precision , dimension(Dim_Ode) , intent(in) , target &
!      :: Scaling
!    double precision , dimension(Dim_Time) , intent(in) , target &
!      :: Timepoint
!    double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
!      , target &
!      :: Observation
!    double precision , intent(in) &
!      :: Ridge_Parameter
!
!    double precision &
!      :: f
!
!    m_dim_time = Dim_Time
!    m_ridge_parameter = Ridge_Parameter
!
!    mp_scaling => Scaling
!    mp_timepoint => Timepoint
!    mp_observation => Observation
!
!    call ResidueNewuoaObjective ( Dim_Ode , Eigen , f )
!
!    write(*,*) f
!
!    return
!
!  end subroutine TestResidueNewuoaObjective!}}}

subroutine ConstructODE &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Linear &
    , Initial &
    , Curve &
    , Info &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(in) &
    :: Eigen
  double precision , dimension(Dim_Time) , intent(in) , target &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(out) &
    , target &
    :: Linear
  double precision , dimension(Dim_Ode) , intent(out) &
    :: Initial
  double precision , dimension(Dim_Ode*Dim_Time) , intent(out) &
    , target &
    :: Curve
  integer , intent(out) &
    :: Info
! 0: Normal
! 1: LU-factorize Q^T Y
! 2: Inverse Q^T Y
! 3: Inverse R^T
! 255: QR-factorize basis
! 254: Reconstruct Q
! 253: LU-factorize Q^T Y
! 252: Inverse Q^T Y
! 251: Inverse R^T

  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: basis
  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: q
  double precision , dimension(Dim_Ode,Dim_Ode) &
    :: r

  double precision , dimension(:) , allocatable &
    :: work
  double precision , dimension(:) , allocatable &
    :: work2
  integer , dimension(:) , allocatable &
    :: iwork
  integer &
    :: ind
  double precision , dimension(Dim_Ode) &
    :: tau
  integer &
    :: info_i

  double precision , dimension(:,:) , pointer &
    :: p_linear
  double precision , dimension(:,:) , pointer &
    :: p_observation
  double precision , dimension(:,:) , pointer &
    :: p_curve

  integer &
    :: itemp
  double precision &
    :: temp

! Main target:
!   Compute similarity transformation operator S and inverse
!   by S = Y Q R^{-T}.
!   Compute Linear output by: S Sigma S^{-1}

  p_linear ( 1:Dim_Ode , 1:Dim_Ode ) => Linear
  p_observation ( 1:Dim_Ode , 1:Dim_Time ) => Observation
  p_curve ( 1:Dim_Ode , 1:Dim_Time ) => Curve
  Info = 0

! Generate basis:
!   Each column is a basis.
!   A diagonal square matrix is binded to the bottom,
!   used for ridge regression.
!   Generate Sigma, stored in Linear.
!   Generate Initial.
!   Generate Curve.

  basis = 0
  p_linear = 0
  Initial = 0
  p_curve = 0

  itemp = (Dim_Ode+1)/2
  allocate ( work(Dim_Time) )
  allocate ( work2(Dim_Time) )
  do ind = 1 , Dim_Ode/2
    work = dexp ( Eigen(ind) * Timepoint )

    work2 = work * dsin ( Eigen(ind+itemp) * Timepoint )
    basis ( 1:Dim_Time , 2*ind-1 ) = work2
    p_curve ( 2*ind-1 , 1:Dim_Time ) = work2

    work2 = work * dcos ( Eigen(ind+itemp) * Timepoint )
    basis ( 1:Dim_Time , 2*ind ) = work2
    p_curve ( 2*ind , 1:Dim_Time ) = work2

    p_linear ( 2*ind-1 , 2*ind-1 ) = Eigen(ind)
    p_linear ( 2*ind , 2*ind ) = Eigen(ind)
    p_linear ( 2*ind-1 , 2*ind ) = Eigen(ind+itemp)
    p_linear ( 2*ind , 2*ind-1 ) = -Eigen(ind+itemp)

    Initial ( 2*ind ) = 1
  end do
  if ( mod(Dim_Ode,2) .eq. 1 ) then
    work = dexp ( Eigen(itemp) * Timepoint )
    basis ( 1:Dim_Time , Dim_Ode ) = work
    p_curve ( Dim_Ode , 1:Dim_Time ) = work
    p_linear ( Dim_Ode , Dim_Ode ) = Eigen(itemp)
    Initial ( Dim_Ode ) = 1
  end if
  deallocate ( work )
  deallocate ( work2 )

  temp = sum ( basis**2 )
  temp = sqrt(temp)
  do ind = 1 , Dim_Ode
    basis ( ind+Dim_Time , ind ) = Ridge_Parameter * temp
  end do

! QR-decompose basis

  allocate ( work(1) )
  ind = -1
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , ind &
    , info_i &
  )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , ind &
    , info_i &
  )
  deallocate ( work )

  if ( info_i .lt. 0 ) then
    Info = 255
    stop
  end if

! Construct R^T stored in "r", then compute
! R^{-T} Sigma R^T
! R^{-T} Initial
! R^{-T} Curve

  r = 0

  do ind = 1 , Dim_Ode
    r ( ind , 1 : ind ) = basis ( 1 : ind , ind )
  end do

  p_linear = matmul ( p_linear , r )

  call dtrtri &
    ( &
      'L' &
      , 'N' &
      , Dim_Ode &
      , r &
      , Dim_Ode &
      , info_i &
    )

  if ( info_i .lt. 0 ) then
    Info = 251
    stop
  end if

  if ( info_i .gt. 0 ) then
    Info = 3
  end if

  p_linear = matmul ( r , p_linear )

  Initial = matmul ( r , Initial )

  p_curve = matmul ( r , p_curve )

! Reconstruct Q

  q = basis

  allocate ( work(1) )
  ind = -1
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , ind &
      , info_i &
    )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , ind &
      , info_i &
    )
  deallocate ( work )

  if ( info_i .lt. 0 ) then
    Info = 254
    stop
  end if

! Compute YQ, saved in "r"

  r = matmul ( p_observation , q(1:Dim_Time,:) )

! Compute YQ Sigma (YQ)^{-1}, YQ Initial

  Initial = matmul ( r , Initial )

  p_linear = matmul ( r , p_linear )

  p_curve = matmul ( r , p_curve )

  allocate ( iwork(Dim_Ode) )
  call dgetrf &
    ( &
      Dim_Ode &
      , Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , info_i &
    )

  if ( info_i .lt. 0 ) then
    Info = 253
    stop
  end if

  if ( info_i .gt. 0 ) then
    Info = 1
  end if

  allocate ( work(1) )
  ind = -1
  call dgetri &
    ( &
      Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dgetri &
    ( &
      Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  deallocate ( work )
  deallocate ( iwork )

  if ( info_i .lt. 0 ) then
    Info = 252
    stop
  end if

  if ( info_i .gt. 0 ) then
    Info = 2
  end if

  p_linear = matmul ( p_linear , r )

  return

end subroutine ConstructODE!}}}

subroutine InverseEigen &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Scaling &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Rhobeg &
    , Rhoend &
    , Maxfun &
    , Linear &
    , Initial &
    , Curve &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(inout) &
    :: Eigen
  double precision , dimension(Dim_Ode) , intent(in) , target &
    :: Scaling
  double precision , dimension(Dim_Time) , intent(in) , target &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , intent(in) &
    :: Rhobeg
  double precision , intent(in) &
    :: Rhoend
  integer , intent(in) &
    :: Maxfun
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(out) &
    :: Linear
  double precision , dimension(Dim_Ode) , intent(out) &
    :: Initial
  double precision , dimension(Dim_Ode*Dim_Time) , intent(out) &
    :: Curve

  integer &
    :: iprint = 3
  integer &
    :: info

  double precision , &
    dimension ( (15*Dim_Ode+97)*Dim_Ode/2+14 ) &
    :: work

  interface
    SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),W(*)
      external CALFUN
    end SUBROUTINE NEWUOA
  end interface

  m_dim_time = Dim_Time
  m_ridge_parameter = Ridge_Parameter

  mp_scaling => Scaling
  mp_timepoint => Timepoint
  mp_observation => Observation

  call newuoa &
    ( &
      Dim_Ode &
      , 2*Dim_Ode+1 &
      , Eigen &
      , Rhobeg &
      , Rhoend &
      , iprint &
      , Maxfun &
      , work &
      , ResidueNewuoaObjective &
    )

  call ConstructODE &
    ( &
      Dim_Ode &
      , Dim_Time &
      , Eigen*Scaling &
      , Timepoint &
      , Observation &
      , Ridge_Parameter &
      , Linear &
      , Initial &
      , Curve &
      , info &
    )

  return

end subroutine InverseEigen!}}}

subroutine ResidualSumOfSquaresLasso &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Lasso_Parameter &
    , Rss &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(in) &
    :: Eigen
  double precision , dimension(Dim_Time) , intent(in) &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(in) &
    , target &
    :: Lasso_Parameter
  double precision , intent(out) &
    :: Rss

  integer &
    :: ind
  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: basis
  double precision , dimension( Dim_Ode , Dim_Ode ) &
    :: linear
  double precision , &
    dimension ( (Dim_Ode+Dim_Time) , Dim_Ode ) &
    :: q
  double precision , dimension( Dim_Ode , Dim_Ode ) &
    :: r

  double precision , dimension(:) , allocatable &
    :: work
  integer , dimension(:) , allocatable &
    :: iwork
  integer &
    :: lwork
  double precision , dimension(Dim_Ode) &
    :: tau
  integer &
    :: info_i

  double precision , dimension(:,:) , pointer &
    :: p_observation
  double precision , dimension(:,:) , pointer &
    :: p_lasso_parameter

  p_observation ( 1:Dim_Ode , 1:Dim_Time ) => Observation
  p_lasso_parameter ( 1:Dim_Ode , 1:Dim_Ode ) => Lasso_Parameter

! Generate basis:
!   Each column is a basis.
!   A diagonal square matrix is binded to the bottom,
!   used for ridge regression.

  basis = 0
  linear = 0

  allocate ( work(Dim_Time) )
  do ind = 1 , Dim_Ode/2
    work = dexp ( Eigen(ind) * Timepoint )
    basis ( : , 2*ind-1 ) &
      = work * dsin ( Eigen(ind+Dim_Ode/2) * Timepoint )
    basis ( : , 2*ind ) &
      = work * dcos ( Eigen(ind+Dim_Ode/2) * Timepoint )
    linear ( 2*ind-1 , 2*ind-1 ) = Eigen(ind)
    linear ( 2*ind , 2*ind ) = Eigen(ind)
    linear ( 2*ind-1 , 2*ind ) = Eigen(ind+Dim_Ode/2)
    linear ( 2*ind , 2*ind-1 ) = -Eigen(ind+Dim_Ode/2)
  end do
  deallocate ( work )

  do ind = 1 , Dim_Ode
    basis ( ind+Dim_Time , ind ) = Ridge_Parameter
  end do

! QR-decompose basis

  allocate ( work(1) )
  lwork = -1
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info_i &
  )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info_i &
  )
  deallocate ( work )

! Reconstruct Q

  q = basis

  allocate ( work(1) )
  lwork = -1
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info_i &
    )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info_i &
    )
  deallocate ( work )

! Construct linear term
! Compute YQ, saved in "r"
! Compute RSS:
!   Residue is ||Y||^2 - ||Q^T Y||^2.
! Compute final return value.

  r = 0

  do ind = 1 , Dim_Ode
    r ( ind , 1 : ind ) = basis ( 1 : ind , ind )
  end do

  linear = matmul ( linear , r )

  call dtrtri &
    ( &
      'L' &
      , 'N' &
      , Dim_Ode &
      , r &
      , Dim_Ode &
      , info_i &
    )

  linear = matmul ( r , linear )

  r = matmul ( p_observation , q(1:Dim_Time,:) )

  Rss = sum ( p_observation**2 ) - sum ( r**2 )

  linear = matmul ( r , linear )

  allocate ( iwork(Dim_Ode) )
  call dgetrf &
    ( &
      Dim_Ode &
      , Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , info_i &
    )

  allocate ( work(1) )
  ind = -1
  call dgetri &
    ( &
      Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dgetri &
    ( &
      Dim_Ode &
      , r &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  deallocate ( work )
  deallocate ( iwork )

  linear = matmul ( linear , r )

  Rss = Rss + sum ( abs ( linear * p_lasso_parameter ) )
!  write(*,*) sum ( abs ( linear * p_lasso_parameter ) )
!write(*,*) p_lasso_parameter * linear

  return

end subroutine ResidualSumOfSquaresLasso!}}}

subroutine ResidueNewuoaObjectiveLasso &!{{{
  ( &
    N , X , F &
  )

  implicit none

  integer , intent(in) :: N
  double precision , intent(in) , dimension(N) :: X
  double precision , intent(out) :: F

  call ResidualSumOfSquaresLasso &
    ( &
      N &
      , m_dim_time &
      , X * mp_scaling &
      , mp_timepoint &
      , mp_observation &
      , m_ridge_parameter &
      , mp_lasso &
      , F &
    )

  return

end subroutine ResidueNewuoaObjectiveLasso!}}}

subroutine InverseEigenLasso &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Scaling &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Lasso_Parameter &
    , Rhobeg &
    , Rhoend &
    , Maxfun &
    , Linear &
    , Initial &
    , Curve &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(inout) &
    :: Eigen
  double precision , dimension(Dim_Ode) , intent(in) , target &
    :: Scaling
  double precision , dimension(Dim_Time) , intent(in) , target &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(in) &
    , target &
    :: Lasso_Parameter
  double precision , intent(in) &
    :: Rhobeg
  double precision , intent(in) &
    :: Rhoend
  integer , intent(in) &
    :: Maxfun
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(out) &
    :: Linear
  double precision , dimension(Dim_Ode) , intent(out) &
    :: Initial
  double precision , dimension(Dim_Ode*Dim_Time) , intent(out) &
    :: Curve

  integer &
    :: iprint = 3
  integer &
    :: info

  double precision , &
    dimension ( (15*Dim_Ode+97)*Dim_Ode/2+14 ) &
    :: work

  interface
    SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),W(*)
      external CALFUN
    end SUBROUTINE NEWUOA
  end interface

  m_dim_time = Dim_Time
  m_ridge_parameter = Ridge_Parameter

  mp_scaling => Scaling
  mp_timepoint => Timepoint
  mp_observation => Observation
  mp_lasso => Lasso_Parameter

  call newuoa &
    ( &
      Dim_Ode &
      , 2*Dim_Ode+1 &
      , Eigen &
      , Rhobeg &
      , Rhoend &
      , iprint &
      , Maxfun &
      , work &
      , ResidueNewuoaObjectiveLasso &
    )

  call ConstructODE &
    ( &
      Dim_Ode &
      , Dim_Time &
      , Eigen*Scaling &
      , Timepoint &
      , Observation &
      , Ridge_Parameter &
      , Linear &
      , Initial &
      , Curve &
      , info &
    )

  return

end subroutine InverseEigenLasso!}}}

end module Linear_Ode_Inverse_Eigen_mod
