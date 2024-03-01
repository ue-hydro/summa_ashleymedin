! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module summaSolve4cvode_module


!======= Inclusions ===========
USE, intrinsic :: iso_c_binding
USE nrtype
USE type4cvode

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! access matrix information
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! constants
USE multiconst,only: Tfreeze        ! temperature at freezing              (K)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookDERIV     ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
  qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                      & ! a big bucket (lumped aquifer model)
  noExplicit                         ! no explicit groundwater parameterization
 
! look-up values for method used to compute derivative
USE mDecisions_module,only:       &
  numerical,                      & ! numerical solution
  analytical                        ! analytical solution

! privacy
 implicit none
 private::setInitialCondition
 private::setSolverParams
 private::find_rootdir
 public::layerDisCont4cvode
 private::getErrMessage
 public::summaSolve4cvode

contains


! ************************************************************************************
! * public subroutine summaSolve4cvode: solve F(y,y') = 0 by CVODE (y is the state vector)
! ************************************************************************************
subroutine summaSolve4cvode(&
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    data time step
                      atol,                    & ! intent(in):    absolute tolerance
                      rtol,                    & ! intent(in):    relative tolerance
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nStat,                   & ! intent(in):    total number of state variables
                      ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      computMassBalance,       & ! intent(in):    flag to compute mass balance
                      computNrgBalance,        & ! intent(in):    flag to compute energy balance
                      ! input: state vectors
                      stateVecInit,            & ! intent(in):    initial state vector
                      sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                      dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                      ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt_cur
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                      ! output
                      ixSaturation,            & ! intent(inout)  index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      cvodeSucceeds,           & ! intent(out):   flag to indicate if CVODE successfully solved the problem in current data step
                      tooMuchMelt,             & ! intent(inout): lag to denote that there was too much melt
                      nSteps,                  & ! intent(out):   number of time steps taken in solver
                      stateVec,                & ! intent(out):   model state vector
                      stateVecPrime,           & ! intent(out):   derivative of model state vector
                      balance,                 & ! intent(inout): balance per state
                      err,message)               ! intent(out):   error control

  !======= Inclusions ===========
  USE fcvode_mod                                              ! Fortran interface to CVODE
  USE fsundials_core_mod                                      ! Fortran interface to SUNContext
  USE fnvector_serial_mod                                     ! Fortran interface to serial N_Vector
  USE fsunmatrix_dense_mod                                    ! Fortran interface to dense SUNMatrix
  USE fsunmatrix_band_mod                                     ! Fortran interface to banded SUNMatrix
  USE fsunlinsol_dense_mod                                    ! Fortran interface to dense SUNLinearSolver
  USE fsunlinsol_band_mod                                     ! Fortran interface to banded SUNLinearSolver
  USE fsunnonlinsol_newton_mod                                ! Fortran interface to Newton SUNNonlinearSolver
  USE allocspace_module,only:allocLocal                       ! allocate local data structures
  USE getVectorz_module, only:checkFeas                       ! check feasibility of state vector
  USE eval8summaWithPrime_module,only:eval8summa4ida        ! DAE/ODE functions
  USE eval8summaWithPrime_module,only:eval8summaWithPrime     ! residual of DAE
  USE computJacobWithPrime_module,only:computJacob4ida      ! system Jacobian
  USE tol4ida_module,only:computWeight4ida                ! weight required for tolerances
  USE var_lookup,only:maxvarDecisions                         ! maximum number of decisions
  !======= Declarations =========
  implicit none

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! calling variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  real(qp),intent(in)             :: dt                     ! data time step
  real(qp),intent(inout)          :: atol(:)                ! vector of absolute tolerances
  real(qp),intent(inout)          :: rtol(:)                ! vector of relative tolerances
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nStat                  ! total number of state variables
  integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)         :: computMassBalance      ! flag to compute mass balance
  logical(lgt),intent(in)         :: computNrgBalance       ! flag to compute energy balance
  ! input: state vectors
  real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
  real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
  real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
  type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  type(var_d),        intent(in)  :: forc_data              ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
   ! input-output: data structures
  type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt_cur
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
  ! output: state vectors
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  integer(i4b),intent(out)        :: nSteps                 ! number of time steps taken in solver
  real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
  real(rkind),intent(inout)       :: stateVecPrime(:)       ! model state vector (y')
  logical(lgt),intent(out)        :: cvodeSucceeds          ! flag to indicate if CVODE is successful
  logical(lgt),intent(inout)      :: tooMuchMelt            ! flag to denote that there was too much melt
  ! output: residual terms and balances
  real(rkind),intent(inout)       :: balance(:)             ! balance per state
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  type(N_Vector),           pointer :: sunvec_y                               ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_yp                              ! sundials derivative vector
  type(N_Vector),           pointer :: sunvec_av                              ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A                               ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS                           ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS                          ! sundials nonlinear solver
  type(c_ptr)                       :: cvode_mem                              ! CVODE memory
  type(c_ptr)                       :: sunctx                                 ! SUNDIALS simulation context
  type(data4ida),           target  :: eqns_data                            ! CVODE type
  integer(i4b)                      :: retval, retvalr                        ! return value
  logical(lgt)                      :: feasible                               ! feasibility flag
  real(qp)                          :: t0                                     ! starting time
  real(qp)                          :: dt_last(1)                             ! last time step
  real(qp)                          :: dt_diff                                ! difference from previous timestep
  integer(c_long)                   :: mu, lu                                 ! in banded matrix mode in SUNDIALS type
  integer(c_long)                   :: nState                                 ! total number of state variables in SUNDIALS type
  integer(i4b)                      :: iVar, i                                ! indices
  integer(i4b)                      :: nRoot                                  ! total number of roots (events) to find
  real(qp)                          :: tret(1)                                ! time in data window
  real(qp)                          :: tretPrev                               ! previous time in data window
  integer(i4b),allocatable          :: rootsfound(:)                          ! crossing direction of discontinuities
  integer(i4b),allocatable          :: rootdir(:)                             ! forced crossing direction of discontinuities
  logical(lgt)                      :: tinystep                               ! if step goes below small size
  type(var_dlength)                 :: flux_prev                              ! previous model fluxes for a local HRU
  character(LEN=256)                :: cmessage                               ! error message of downwind routine
  real(rkind),allocatable           :: mLayerMatricHeadPrimePrev(:)           ! previous derivative value for total water matric potential (m s-1)
  real(rkind),allocatable           :: resVecPrev(:)                          ! previous value for residuals
  real(rkind),allocatable           :: dCompress_dPsiPrev(:)                  ! previous derivative value soil compression
  ! flags
  logical(lgt)                      :: use_fdJac                              ! flag to use finite difference Jacobian, controlled by decision fDerivMeth
  logical(lgt),parameter            :: offErrWarnMessage = .true.             ! flag to turn CVODE warnings off, default true
  logical(lgt),parameter            :: detect_events = .true.                 ! flag to do event detection and restarting, default true
  ! -----------------------------------------------------------------------------------------------------
  ! link to the necessary variables
  associate(&
    ! number of state variables of a specific type
    nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1) ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1) ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    ! model indices
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)      ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)      ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)      ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)       ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat    ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat    ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat    ,& ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
    ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat    ,& ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
    ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat    ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
    layerType               => indx_data%var(iLookINDEX%layerType)%dat         & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
    ) ! association to necessary variables for the residual computations

    ! initialize error control
    err=0; message="summaSolve4cvode/"
    
    ! choose Jacobian type
    select case(model_decisions(iLookDECISIONS%fDerivMeth)%iDecision) 
      case(numerical);  use_fdJac =.true.
      case(analytical); use_fdJac =.false.
      case default; err=20; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian'; return
    end select
    
    nState = nStat ! total number of state variables in SUNDIALS type
    cvodeSucceeds = .true.
    
    ! fill eqns_data which will be required later to call eval8summaWithPrime
    eqns_data%dt             = dt
    eqns_data%nSnow          = nSnow
    eqns_data%nSoil          = nSoil
    eqns_data%nLayers        = nLayers
    eqns_data%nState         = nState
    eqns_data%ixMatrix       = ixMatrix
    eqns_data%firstSubStep   = firstSubStep
    eqns_data%computeVegFlux = computeVegFlux
    eqns_data%scalarSolution = scalarSolution
    eqns_data%type_data      = type_data
    eqns_data%attr_data      = attr_data
    eqns_data%mpar_data      = mpar_data
    eqns_data%forc_data      = forc_data
    eqns_data%bvar_data      = bvar_data
    eqns_data%prog_data      = prog_data
    eqns_data%indx_data      = indx_data
    eqns_data%diag_data      = diag_data
    eqns_data%flux_data      = flux_data
    eqns_data%deriv_data     = deriv_data
    eqns_data%ixSaturation   = ixSaturation
    
    ! allocate space and fill
    allocate( eqns_data%model_decisions(maxvarDecisions) ); eqns_data%model_decisions = model_decisions
    allocate( eqns_data%atol(nState) ); eqns_data%atol = atol
    allocate( eqns_data%rtol(nState) ); eqns_data%rtol = rtol
    allocate( eqns_data%sMul(nState) ); eqns_data%sMul = sMul
    allocate( eqns_data%dMat(nState) ); eqns_data%dMat = dMat
    
    ! allocate space for the to save previous fluxes
    call allocLocal(flux_meta(:),flux_prev,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
    flux_prev                         = eqns_data%flux_data
    
    ! allocate space for other variables
    if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel)then
      allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
    else
      allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
    end if
    allocate( eqns_data%mLayerTempPrev(nLayers) )
    allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
    allocate( eqns_data%mLayerTempTrial(nLayers) )
    allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
    allocate( eqns_data%mLayerMatricHeadLiqPrime(nSoil) )
    allocate( eqns_data%mLayerVolFracWatTrial(nLayers) )
    allocate( eqns_data%mLayerTempPrime(nLayers) )       
    allocate( eqns_data%mLayerMatricHeadPrime(nSoil) )
    allocate( eqns_data%mLayerVolFracWatPrime(nLayers) ) 
    allocate( eqns_data%mLayerVolFracIcePrime(nLayers) )
    allocate( mLayerMatricHeadPrimePrev(nSoil) )
    allocate( dCompress_dPsiPrev(nSoil) )
    allocate( eqns_data%fluxVec(nState) )
    allocate( eqns_data%resVec(nState) )
    allocate( eqns_data%resSink(nState) )
    allocate( resVecPrev(nState) )

    ! need the following values for the first substep
    eqns_data%scalarCanopyTempPrev    = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)
    eqns_data%mLayerTempPrev(:)       = prog_data%var(iLookPROG%mLayerTemp)%dat(:)
    eqns_data%mLayerMatricHeadPrev(:) = prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
    dCompress_dPsiPrev(:)             = 0._rkind
    resVecPrev(:)                     = 0._rkind
    
    retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
    
    ! create serial vectors
    sunvec_y => FN_VMake_Serial(nState, stateVec, sunctx)
    if (.not. associated(sunvec_y)) then; err=20; message=trim(message)//'sunvec = NULL'; return; endif
    sunvec_yp => FN_VMake_Serial(nState, stateVecPrime, sunctx)
    if (.not. associated(sunvec_yp)) then; err=20; message=trim(message)//'sunvec = NULL'; return; endif
    
    ! initialize solution vectors
    call setInitialCondition(nState, stateVecInit, sunvec_y, sunvec_yp)
    
    ! create memory
    cvode_mem = FCVODECreate(sunctx)
    if (.not. c_associated(cvode_mem)) then; err=20; message=trim(message)//'cvode_mem = NULL'; return; endif
    
    ! Attach user data to memory
    retval = FCVODESetUserData(cvode_mem, c_loc(eqns_data))
    if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODESetUserData'; return; endif
    
    ! Set the function CVODE will use to advance the state
    t0 = 0._rkind
    retval = FCVODEInit(cvode_mem, c_funloc(eval8summa4ida), t0, sunvec_y, sunvec_yp)
    if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODEInit'; return; endif
    
    ! set tolerances
    retval = FCVODEWFtolerances(cvode_mem, c_funloc(computWeight4ida))
    if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODEWFtolerances'; return; endif
    
    ! initialize rootfinding problem and allocate space, counting roots
    if(detect_events)then
      nRoot = 0
      if(ixVegNrg/=integerMissing) nRoot = nRoot+1
      if(nSnow>0)then
        do i = 1,nSnow
          if(ixSnowOnlyNrg(i)/=integerMissing) nRoot = nRoot+1
        enddo
      endif
      if(nSoil>0)then
        do i = 1,nSoil
          if(ixSoilOnlyHyd(i)/=integerMissing) nRoot = nRoot+1
          if(ixSoilOnlyNrg(i)/=integerMissing) nRoot = nRoot+1
        enddo
      endif
      allocate( rootsfound(nRoot) )
      allocate( rootdir(nRoot) )
      rootdir = 0
      retval = FCVODERootInit(cvode_mem, nRoot, c_funloc(layerDisCont4cvode))
      if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODERootInit'; return; endif
    else ! will not use, allocate at something
      nRoot = 1
      allocate( rootsfound(nRoot) )
      allocate( rootdir(nRoot) )
    endif
    
    ! define the form of the matrix
    select case(ixMatrix)
      case(ixBandMatrix)
        mu = ku; lu = kl;
        ! Create banded SUNMatrix for use in linear solves
        sunmat_A => FSUNBandMatrix(nState, mu, lu, sunctx)
        if (.not. associated(sunmat_A)) then; err=20; message=trim(message)//'sunmat = NULL'; return; endif
    
        ! Create banded SUNLinearSolver object
        sunlinsol_LS => FSUNLinSol_Band(sunvec_y, sunmat_A, sunctx)
        if (.not. associated(sunlinsol_LS)) then; err=20; message=trim(message)//'sunlinsol = NULL'; return; endif
    
      case(ixFullMatrix)
        ! Create dense SUNMatrix for use in linear solves
        sunmat_A => FSUNDenseMatrix(nState, nState, sunctx)
        if (.not. associated(sunmat_A)) then; err=20; message=trim(message)//'sunmat = NULL'; return; endif
    
        ! Create dense SUNLinearSolver object
        sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
        if (.not. associated(sunlinsol_LS)) then; err=20; message=trim(message)//'sunlinsol = NULL'; return; endif
    
        ! check
      case default;  err=20; message=trim(message)//'error in type of matrix'; return
    
    end select  ! form of matrix
    
    ! Attach the matrix and linear solver
    ! For the nonlinear solver, CVODE uses a Newton SUNNonlinearSolver-- it is not necessary to create and attach it
    retval = FCVODESetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A);
    if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODESetLinearSolver'; return; endif
    
    ! Set the user-supplied Jacobian routine
    if(.not.use_fdJac)then
      retval = FCVODESetJacFn(cvode_mem, c_funloc(computJacob4ida))
      if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODESetJacFn'; return; endif
    endif
    
    ! Enforce the solver to stop at end of the time step
    retval = FCVODESetStopTime(cvode_mem, dt_cur)
    if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODESetStopTime'; return; endif
    
    ! Set solver parameters at end of setup
    call setSolverParams(dt_cur, cvode_mem, retval)
    if (retval /= 0) then; err=20; message=trim(message)//'error in setSolverParams'; return; endif
    
    ! Disable error messages and warnings
    if(offErrWarnMessage) then
      retval = FSUNLogger_SetErrorFilename(cvode_mem, c_null_char)
      retval = FSUNLogger_SetWarningFilename(cvode_mem, c_null_char)
      retval = FCVODESetNoInactiveRootWarn(cvode_mem)
    endif
    
    !*********************** Main Solver * loop on one_step mode *****************************
    tinystep = .false.
    tret(1) = t0 ! initial time
    tretPrev = tret(1)
    nSteps = 0 ! initialize number of time steps taken in solver
    
    do while(tret(1) < dt_cur)
    
      ! call this at beginning of step to reduce root bouncing (only looking in one direction)
      if(detect_events .and. .not.tinystep)then
        call find_rootdir(eqns_data, rootdir)
        retval = FCVODESetRootDirection(cvode_mem, rootdir)
        if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODESetRootDirection'; return; endif
      endif
    
      eqns_data%firstFluxCall = .false. ! already called for initial
      eqns_data%firstSplitOper = .true. ! always true at start of dt_cur since no splitting
    
      ! call CVODESolve, advance solver just one internal step
      retvalr = FCVODESolve(cvode_mem, dt_cur, tret, sunvec_y, sunvec_yp, CVODE_ONE_STEP)
      ! early return if CVODESolve failed
      if( retvalr < 0 )then
        cvodeSucceeds = .false.
        call getErrMessage(retvalr,cmessage)
        message=trim(message)//trim(cmessage)
        !if(retvalr==-1) err = -20 ! max iterations failure, exit and reduce the data window time in varSubStep
        exit
      end if
    
      tooMuchMelt = .false.
      ! loop through non-missing energy state variables in the snow domain to see if need to merge
      do concurrent (i=1:nSnow,ixSnowOnlyNrg(i)/=integerMissing)
        if (stateVec(ixSnowOnlyNrg(i)) > Tfreeze) tooMuchMelt = .true. !need to merge
      enddo
      if(tooMuchMelt)exit
    
      ! get the last stepsize and difference from previous end time, not necessarily the same
      retval = FCVODEGetLastStep(cvode_mem, dt_last)
      dt_diff = tret(1) - tretPrev
      nSteps = nSteps + 1 ! number of time steps taken in solver
    
      ! check the feasibility of the solution
      feasible=.true.
      call checkFeas(&
                      ! input
                      stateVec,            & ! intent(in):    model state vector (mixed units)
                      eqns_data%mpar_data, & ! intent(in):    model parameters
                      eqns_data%prog_data, & ! intent(in):    model prognostic variables for a local HRU
                      eqns_data%indx_data, & ! intent(in):    indices defining model states and layers
                      ! output: feasibility
                      feasible,            & ! intent(inout):   flag to denote the feasibility of the solution
                    ! output: error control
                      err,cmessage)           ! intent(out):   error control
    
      ! early return for non-feasible solutions, right now will just fail if goes infeasible
      if(.not.feasible)then
        cvodeSucceeds = .false.
        message=trim(message)//trim(cmessage)//'non-feasible' ! err=0 is already set, could make this a warning and reduce the data window time in varSubStep
        exit
      end if
    
      ! sum of fluxes smoothed over the time step, average from instantaneous values
      do iVar=1,size(flux_meta)
        flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + ( eqns_data%flux_data%var(iVar)%dat(:) &
                                                                  + flux_prev%var(iVar)%dat(:) ) *dt_diff/2._rkind
      end do
      mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + ( eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:) * eqns_data%mLayerMatricHeadPrime(:) &
                                                      + dCompress_dPsiPrev(:)  * mLayerMatricHeadPrimePrev(:) ) * dt_diff/2._rkind
    
      ! ----
      ! * compute energy balance, from residuals
      !  formulation with prime variables would cancel to closedForm version, so does not matter which formulation is used
      !------------------------
      if(computNrgBalance)then    
    
        ! compute energy balance mean, resVec is the instanteous residual vector from the solver
        if(ixCasNrg/=integerMissing) balance(ixCasNrg) = balance(ixCasNrg) + ( eqns_data%resVec(ixCasNrg) + resVecPrev(ixCasNrg) )*dt_diff/2._rkind/dt
        if(ixVegNrg/=integerMissing) balance(ixVegNrg) = balance(ixVegNrg) + ( eqns_data%resVec(ixVegNrg) + resVecPrev(ixVegNrg) )*dt_diff/2._rkind/dt
        if(nSnowSoilNrg>0)then
          do concurrent (i=1:nLayers,ixSnowSoilNrg(i)/=integerMissing) 
            balance(ixSnowSoilNrg(i)) = balance(ixSnowSoilNrg(i)) + ( eqns_data%resVec(ixSnowSoilNrg(i)) + resVecPrev(ixSnowSoilNrg(i)) )*dt_diff/2._rkind/dt
          enddo
        endif
      endif
    
      ! ----
      ! * compute mass balance, from residuals
      !------------------------
      if(computMassBalance)then
    
        ! compute mass balance mean, resVec is the instanteous residual vector from the solver
        if(ixVegHyd/=integerMissing) balance(ixVegHyd) = balance(ixVegHyd) + ( eqns_data%resVec(ixVegHyd) + resVecPrev(ixVegHyd) )*dt_diff/2._rkind/dt
        if(nSnowSoilHyd>0)then
          do concurrent (i=1:nLayers,ixSnowSoilHyd(i)/=integerMissing) 
            balance(ixSnowSoilHyd(i)) = balance(ixSnowSoilHyd(i)) + ( eqns_data%resVec(ixSnowSoilHyd(i)) + resVecPrev(ixSnowSoilHyd(i)) )*dt_diff/2._rkind/dt
          enddo
        endif
        if(ixAqWat/=integerMissing) balance(ixAqWat) = balance(ixAqWat) + ( eqns_data%resVec(ixAqWat) + resVecPrev(ixAqWat) )*dt_diff/2._rkind/dt
      endif
    
      ! save required quantities for next step
      eqns_data%scalarCanopyTempPrev         = eqns_data%scalarCanopyTempTrial
      eqns_data%mLayerTempPrev(:)            = eqns_data%mLayerTempTrial(:)
      eqns_data%mLayerMatricHeadPrev(:)      = eqns_data%mLayerMatricHeadTrial(:)
      mLayerMatricHeadPrimePrev(:)           = eqns_data%mLayerMatricHeadPrime(:) 
      dCompress_dPsiPrev(:)                  = eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:)
      tretPrev                               = tret(1)
      resVecPrev(:)                          = eqns_data%resVec(:)
    
      ! Restart for where vegetation and layers cross freezing point
      if(detect_events)then
        if (retvalr .eq. CVODE_ROOT_RETURN) then !CVODESolve succeeded and found one or more roots at tret(1)
          ! rootsfound[i]= +1 indicates that gi is increasing, -1 g[i] decreasing, 0 no root
          !retval = FCVODEGetRootInfo(cvode_mem, rootsfound)
          !if (retval < 0) then; err=20; message=trim(message)//'error in FCVODEGetRootInfo'; return; endif
          !print '(a,f15.7,2x,17(i2,2x))', "time, rootsfound[] = ", tret(1), rootsfound
          ! Reininitialize solver for running after discontinuity and restart
          retval = FCVODEReInit(cvode_mem, tret(1), sunvec_y, sunvec_yp)
          if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODEReInit'; return; endif
          if(dt_last(1) < 0.1_rkind)then ! don't keep calling if step is small (more accurate with this tiny but getting hung up)
            retval = FCVODERootInit(cvode_mem, 0, c_funloc(layerDisCont4cvode))
            tinystep = .true.
          else
            retval = FCVODERootInit(cvode_mem, nRoot, c_funloc(layerDisCont4cvode))
            tinystep = .false.
          endif
          if (retval /= 0) then; err=20; message=trim(message)//'error in FCVODERootInit'; return; endif
        endif
      endif
    
    enddo ! while loop on one_step mode until time dt_cur
    !****************************** End of Main Solver ***************************************
    
    if(cvodeSucceeds)then
      ! copy to output data
      diag_data     = eqns_data%diag_data
      flux_data     = eqns_data%flux_data
      deriv_data    = eqns_data%deriv_data
      ixSaturation  = eqns_data%ixSaturation
      indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = eqns_data%indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) !only number of flux calculations changes in indx_data
      err           = eqns_data%err
      message       = eqns_data%message
    endif
    
    ! free memory
    deallocate( eqns_data%model_decisions)
    deallocate( eqns_data%sMul )
    deallocate( eqns_data%dMat )
    deallocate( eqns_data%dBaseflow_dMatric )
    deallocate( eqns_data%mLayerTempPrev )
    deallocate( eqns_data%mLayerMatricHeadPrev )
    deallocate( eqns_data%mLayerTempTrial )
    deallocate( eqns_data%mLayerMatricHeadTrial )
    deallocate( eqns_data%mLayerVolFracWatTrial )
    deallocate( eqns_data%mLayerTempPrime )       
    deallocate( eqns_data%mLayerMatricHeadPrime )
    deallocate( eqns_data%mLayerMatricHeadLiqPrime)
    deallocate( eqns_data%mLayerVolFracWatPrime ) 
    deallocate( eqns_data%mLayerVolFracIcePrime )
    deallocate( mLayerMatricHeadPrimePrev )
    deallocate( dCompress_dPsiPrev )
    deallocate( eqns_data%resVec )
    deallocate( eqns_data%resSink )
    deallocate( rootsfound )
    deallocate( rootdir )
    
    call FCVODEFree(cvode_mem)
    retval = FSUNLinSolFree(sunlinsol_LS)
    if(retval /= 0)then; err=20; message=trim(message)//'unable to free the linear solver'; return; endif
    call FSUNMatDestroy(sunmat_A)
    call FN_VDestroy(sunvec_y)
    call FN_VDestroy(sunvec_yp)
    retval = FSUNContext_Free(sunctx)
    if(retval /= 0)then; err=20; message=trim(message)//'unable to free the SUNDIALS context'; return; endif

  end associate

end subroutine summaSolve4cvode

! ----------------------------------------------------------------
! SetInitialCondition: routine to initialize u and up vectors.
! ----------------------------------------------------------------
subroutine setInitialCondition(neq, y, sunvec_u, sunvec_up)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fsundials_core_mod
  USE fnvector_serial_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)  :: sunvec_u  ! solution N_Vector
  type(N_Vector)  :: sunvec_up ! derivative N_Vector
  integer(c_long) :: neq
  real(rkind)     :: y(neq)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  real(c_double), pointer :: up(:)

  ! get data arrays from SUNDIALS vectors
  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
  up(1:neq) => FN_VGetArrayPointer(sunvec_up)

  uu = y
  up = 0._rkind

end subroutine setInitialCondition

! ----------------------------------------------------------------
! setSolverParams: private routine to set parameters in CVODE solver
! ----------------------------------------------------------------
subroutine setSolverParams(dt_cur,cvode_mem,retval)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fcvode_mod   ! Fortran interface to CVODE

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind),intent(in)      :: dt_cur             ! current whole time step
  type(c_ptr),intent(inout)   :: cvode_mem            ! CVODE memory
  integer(i4b),intent(out)    :: retval             ! return value

  !======= Internals ============
  integer,parameter           :: nonlin_iter = 4    ! maximum number of nonlinear iterations before reducing step size, default = 4
  integer,parameter           :: max_order = 5      ! maximum BDF order,  default and max = 5
  real(qp),parameter          :: coef_nonlin = 0.33 ! coefficient in the nonlinear convergence test, default = 0.33
  integer(c_long),parameter   :: max_step = 999999  ! maximum number of steps,  default = 500
  integer,parameter           :: fail_iter = 50     ! maximum number of error test and convergence test failures, default 10
  real(qp)                    :: h_max              ! maximum stepsize,  default = infinity
  real(qp),parameter          :: h_init = 0         ! initial stepsize
 
  ! Set the maximum BDF order
  retval = FCVODESetMaxOrd(cvode_mem, max_order)
  if (retval /= 0) return

  ! Set coefficient in the nonlinear convergence test
  retval = FCVODESetNonlinConvCoef(cvode_mem, coef_nonlin)
  if (retval /= 0) return

  ! Set maximun number of nonliear iterations, maybe should just make 4 (instead of SUMMA parameter)
  retval = FCVODESetMaxNonlinIters(cvode_mem, nonlin_iter)
  if (retval /= 0) return

  !  Set maximum number of convergence test failures
  retval = FCVODESetMaxConvFails(cvode_mem, fail_iter)
  if (retval /= 0) return

  !  Set maximum number of error test failures
  retval = FCVODESetMaxErrTestFails(cvode_mem, fail_iter)
  if (retval /= 0) return

  ! Set maximum number of steps
  retval = FCVODESetMaxNumSteps(cvode_mem, max_step)
  if (retval /= 0) return

  ! Set maximum stepsize
  h_max = dt_cur
  retval = FCVODESetMaxStep(cvode_mem, h_max)
  if (retval /= 0) return

  ! Set initial stepsize
  retval = FCVODESetInitStep(cvode_mem, h_init)
  if (retval /= 0) return

end subroutine setSolverParams

! ----------------------------------------------------------------------------------------
! find_rootdir: private routine to determine which direction to look for the root, by
!  determining if the variable is greater or less than the root. Need to do this to prevent
!  bouncing around solution
! ----------------------------------------------------------------------------------------
subroutine find_rootdir(eqns_data,rootdir)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
  use globalData,only:integerMissing     ! missing integer
  use var_lookup,only:iLookINDEX         ! named variables for structure elements
  use multiconst,only:Tfreeze            ! freezing point of pure water (K)

  !======= Declarations =========
  implicit none

  ! calling variables
  type(data4ida),intent(in)  :: eqns_data  ! equations data
  integer(i4b),intent(inout) :: rootdir(:)   ! root function directions to search

  ! local variables
  integer(i4b)               :: i,ind        ! indices
  integer(i4b)               :: nState       ! number of states
  integer(i4b)               :: nSnow        ! number of snow layers
  integer(i4b)               :: nSoil        ! number of soil layers
  real(rkind)                :: xPsi         ! matric head at layer (m)
  real(rkind)                :: TcSoil       ! critical point when soil begins to freeze (K)

  ! get equations data variables
  nState = eqns_data%nState
  nSnow = eqns_data%nSnow
  nSoil = eqns_data%nSoil

  ! initialize
  ind = 0

  ! identify the critical point when vegetation begins to freeze
  if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
    ind = ind+1
    rootdir(ind) = 1
    if(eqns_data%scalarCanopyTempPrev > Tfreeze) rootdir(ind) = -1
  endif

  if(nSnow>0)then
    do i = 1,nSnow
      ! identify the critical point when the snow layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        rootdir(ind) = 1
        if(eqns_data%mLayerTempPrev(i) > Tfreeze) rootdir(ind) = -1
      endif
    end do
  endif

  if(nSoil>0)then
    do i = 1,nSoil
      xPsi = eqns_data%mLayerMatricHeadPrev(i)
      ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
      if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
        ind = ind+1
        rootdir(ind) = 1
        if(xPsi > 0._rkind ) rootdir(ind) = -1
      endif
      ! identify the critical point when the soil layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        TcSoil = crit_soilT(xPsi)
        rootdir(ind) = 1
        if(eqns_data%mLayerTempPrev(i+nSnow) > TcSoil) rootdir(ind) = -1
      endif
    end do
  endif

end subroutine find_rootdir

! ----------------------------------------------------------------------------------------
! layerDisCont4cvode: The root function routine to find soil matrix potential = 0,
!  soil temp = critical frozen point, and snow and veg temp = Tfreeze
! ----------------------------------------------------------------------------------------
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------------------------------
integer(c_int) function layerDisCont4cvode(t, sunvec_u, sunvec_up, gout, user_data) &
      result(ierr) bind(C,name='layerDisCont4cvode')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
  use globalData,only:integerMissing     ! missing integer
  use var_lookup,only:iLookINDEX         ! named variables for structure elements
  use multiconst,only:Tfreeze            ! freezing point of pure water (K)

  !======= Declarations =========
  implicit none

  ! calling variables
  real(c_double), value      :: t         ! current time
  type(N_Vector)             :: sunvec_u  ! solution N_Vector
  type(N_Vector)             :: sunvec_up ! derivative N_Vector
  real(c_double)             :: gout(999) ! root function values, if (nVeg + nSnow + 2*nSoil)>999, problem
  type(c_ptr),    value      :: user_data ! user-defined data

  ! local variables
  integer(i4b)               :: i,ind     ! indices
  integer(i4b)               :: nState    ! number of states
  integer(i4b)               :: nSnow     ! number of snow layers
  integer(i4b)               :: nSoil     ! number of soil layers
  real(rkind)                :: xPsi      ! matric head at layer (m)
  real(rkind)                :: TcSoil    ! critical point when soil begins to freeze (K)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  type(data4ida), pointer :: eqns_data  ! equations data

  !======= Internals ============
  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)
  nState = eqns_data%nState
  nSnow = eqns_data%nSnow
  nSoil = eqns_data%nSoil

  ! get data array from SUNDIALS vector
  uu(1:nState) => FN_VGetArrayPointer(sunvec_u)

  ! initialize
  ind = 0

  ! identify the critical point when vegetation begins to freeze
  if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
    ind = ind+1
    gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)) - Tfreeze
  endif

  if(nSnow>0)then
    do i = 1,nSnow
      ! identify the critical point when the snow layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)) - Tfreeze
      endif
    end do
  endif

  if(nSoil>0)then
    do i = 1,nSoil
      ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
      if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
        ind = ind+1
        xPsi = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
      else
        xPsi = eqns_data%prog_data%var(iLookPROG%mLayerMatricHead)%dat(i)
      endif
      ! identify the critical point when the soil layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        TcSoil = crit_soilT(xPsi)
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)) - TcSoil
      endif
    end do
  endif

  ! return success
  ierr = 0
  return

end function layerDisCont4cvode

! ----------------------------------------------------------------
! getErrMessage: private routine to get error message for CVODE solver
! ----------------------------------------------------------------
subroutine getErrMessage(retval,message)

  !======= Declarations =========
  implicit none

  ! calling variables
  integer(i4b),intent(in)    :: retval              ! return value from CVODE
  character(*),intent(out)   :: message             ! error message

  ! get message
   if( retval==-1 ) message = 'CV_TOO_MUCH_WORK'    ! The solver took mxstep internal steps but could not reach tout.
   if( retval==-2 ) message = 'CV_TOO_MUCH_ACC'     ! The solver could not satisfy the accuracy demanded by the user for some internal step.
   if( retval==-3 ) message = 'CV_ERR_FAIL'         ! Error test failures occurred too many times during one internal timestep or minimum step size was reached.
   if( retval==-4 ) message = 'CV_CONV_FAIL'        ! Convergence test failures occurred too many times during one internal time step or minimum step size was reached.
   if( retval==-5 ) message = 'CV_LINIT_FAIL'       ! The linear solver’s initialization function failed.
   if( retval==-6 ) message = 'CV_LSETUP_FAIL'      ! The linear solver’s setup function failed in an unrecoverable manner.
   if( retval==-7 ) message = 'CV_LSOLVE_FAIL'      ! The linear solver’s solve function failed in an unrecoverable manner.
   if( retval==-8 ) message = 'CV_RES_FAIL'         ! The user-provided residual function failed in an unrecoverable manner.
   if( retval==-9 ) message = 'CV_REP_RES_FAIL'     ! The user-provided residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.
   if( retval==-10) message = 'CV_RTFUNC_FAIL'      ! The rootfinding function failed in an unrecoverable manner.
   if( retval==-11) message = 'CV_CONSTR_FAIL'      ! The inequality constraints were violated and the solver was unable to recover.
   if( retval==-12) message = 'CV_FIRST_RES_FAIL'   ! The user-provided residual function failed recoverably on the first call.
   if( retval==-13) message = 'CV_LINESEARCH_FAIL'  ! The line search failed.
   if( retval==-14) message = 'CV_NO_RECOVERY'      ! The residual function, linear solver setup function, or linear solver solve function had a recoverable failure, but CVCalcIC could not recover.
   if( retval==-15) message = 'CV_NLS_INIT_FAIL'    ! The nonlinear solver’s init routine failed.
   if( retval==-16) message = 'CV_NLS_SETUP_FAIL'   ! The nonlinear solver’s setup routine failed.
   if( retval==-20) message = 'CV_MEM_NULL'         ! The cvode_mem argument was NULL.
   if( retval==-21) message = 'CV_MEM_FAIL'         ! A memory allocation failed.
   if( retval==-22) message = 'CV_ILL_INPUT'        ! One of the function inputs is illegal.
   if( retval==-23) message = 'CV_NO_MALLOC'        ! The CVODE memory was not allocated by a call to CVODEInit.
   if( retval==-24) message = 'CV_BAD_EWT'          ! Zero value of some error weight component.
   if( retval==-25) message = 'CV_BAD_K'            ! The k-th derivative is not available.
   if( retval==-26) message = 'CV_BAD_T'            ! The time t is outside the last step taken.
   if( retval==-27) message = 'CV_BAD_DKY'          ! The vector argument where derivative should be stored is NULL.

end subroutine getErrMessage


end module summaSolve4cvode_module