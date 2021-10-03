!----------------------------------------------------------------------------- best with 100 columns

!> environment for fepic
module modFEPIC
  use modParticle
  use modPICGrid
  use modCondition
  use modSparse
  use modSource
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  double precision,parameter::AMU=1.6605390666050d-27 !< atomic mass unit (Dalton) [kg]
  double precision,parameter::QE=1.602176634d-19 !< elementary charge [C]
  double precision,parameter::EPS0=8.85418782d-12 !< vacuum permittivity [F/m]
  double precision,parameter::KB=1.38065d-23 !< Boltzmann constant [J/K]
  double precision,parameter::EV2K=11604.52d0 !< 1eV in Kelvin [QE/K]
  
  integer,parameter::BC_PHI_DEFAULT=0 !< default phi BC: Neumann with zero normal phi gradient
  integer,parameter::BC_PHI_DIRICHLET=1 !< Dirichlet BC with given potential
  
  integer,parameter::BC_PTCL_DEFAULT=0 !< default particle BC: removed upon impact
  
  integer::iProc,nProc,ierr !< mpi variables
  
  type(ptcls),allocatable::p(:) !< particles of each species
  type(PICGrid)::grid !< the grid
  type(condTab)::phibc !< electric field boundary conditions
  type(condTab)::ptclbc !< particle boundary conditions
  type(multiFront)::PoissonPhi !< FEM Poisson equation for phi
  type(pSrc),allocatable::ptclSrc(:) !< particle source
  
  integer::nSp !< number of particle species
  
  double precision,allocatable::phi(:) !< electric potential [V]
  double precision,allocatable::ef(:,:) !< electric field strength [V/m]
  double precision,allocatable::rhsPhi(:) !< RHS of the phi equation [C], and [V] (Dirichlet nodes)
  double precision,allocatable::rhsPhiDi(:) !< RHS of the phi equation, only Dirichlet nodes [V]
  double precision,allocatable::rhsPhiLocal(:) !< RHS of the phi equation, local version [C]
  double precision,allocatable::rho(:) !< plasma density [QE/m^3]
  double precision,allocatable::rhoLocal(:) !< plasma density, local version [QE/m^3]
  double precision,allocatable::nVol(:) !< FEM nodal volume
  integer,allocatable::iPhiBC(:) !< indexes of electric field boundary conditions
  integer,allocatable::iPtclBC(:) !< indexes of particle boundary conditions
  logical,allocatable::isDirichlet(:) !< whether a point is a Dirichlet point
  
  double precision::t !< current time [s]
  double precision::dt !< time step size [s]
  double precision::tFinal !< total simulation time [s]
  
  double precision::tOutNext !< time to write the next output [s]
  double precision::dtOut !< interval to write output [s]
  integer::iOut !< output index
  
  ! Boltzmann relationship electron model parameters
  double precision::ne0BR !< reference electron density [m^-3]
  double precision::phi0BR !< reference electric potential [V]
  double precision::kbTeBR !< electron temperature [eV]
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    use modCondition
    use modSparse
    use modBasicFEM
    use modSource
    integer,parameter::FID=10
    double precision::mass,charge
    
    ! read grid
    open(FID,file='grid.msh',action='read')
      call readGMSH(FID,grid)
    close(FID)
    call grid%up()
    
    ! work space and initial state
    allocate(phi(grid%nN))
    allocate(rhsPhi(grid%nN))
    allocate(rhsPhiDi(grid%nN))
    allocate(rhsPhiLocal(grid%nN))
    allocate(rho(grid%nN))
    allocate(rhoLocal(grid%nN))
    allocate(nVol(grid%nN))
    allocate(isDirichlet(grid%nN))
    allocate(ef(3,grid%nN))
    
    ! read particle species
    open(FID,file='ptclSp',action='read')
      read(FID,*)nSp
      allocate(p(nSp))
      do i=1,nSp
        read(FID,*)
        read(FID,*)mass
        read(FID,*)charge
        call p(i)%init(mass*AMU,charge*QE)
      end do
    close(FID)
    
    ! read particle sources
    open(FID,file='ptclSrc',action='read')
      read(FID,*)n
      allocate(ptclSrc(n))
      do i=1,n
        call ptclSrc(i)%load(FID)
      end do
    close(FID)
    
    ! read electric field BCs
    open(FID,file='phiBC',action='read')
      call readCondTab(FID,phibc)
    close(FID)
    call mapCondTab(grid,phibc,iPhiBC)
    rhsPhiDi(:)=0d0
    isDirichlet(:)=.false.
    do i=grid%nC+1,grid%nE
      if(iPhiBC(i)>0)then
        select case(phibc%t(iPhiBC(i)))
        case(BC_PHI_DIRICHLET) !< Dirichlet with given phi
          do j=1,grid%nNE(i)
            isDirichlet(grid%iNE(j,i))=.true.
            rhsPhiDi(grid%iNE(j,i))=phibc%p(1,iPhiBC(i))
          end do
        case default
        end select
      end if
    end do
    
    ! read particle BCs
    open(FID,file='ptclBC',action='read')
      call readCondTab(FID,ptclbc)
    close(FID)
    call mapCondTab(grid,ptclbc,iPtclBC)
    
    ! TODO: read electron model parameters
    ne0BR=1d11
    phi0BR=0d0
    kbTeBR=2d0
    
    ! TODO: read simulation parameters
    tFinal=1d-4
    dt=2d-7
    dtOut=4d-7
    
    t=0d0
    iOut=0
    tOutNext=dtOut
    
  end subroutine
  
  !> step particles
  subroutine stepPtcls()
    use modPush
    integer::f,info
    double precision::h
    integer::n0(size(p))
    logical,allocatable::toBeRemoved(:)
    
    n0(:)=p(:)%n
    ! emission from particle sources
    do i=1,size(ptclSrc)
      call ptclSrc(i)%emit(p,grid,dt)
    end do
    ! push particles
    do j=1,size(p)
      allocate(toBeRemoved(p(j)%n))
      toBeRemoved(:)=.false.
      do i=1,p(j)%n
        h=dt
        if(i<=n0(j))then ! old particles
          call push(grid,p(j),i,h,phi,f,LF_NORMAL,info)
        else ! new particles
          call push(grid,p(j),i,h,phi,f,LF_REWIND,info)
        end if
        if(info==PUSH_IMPACT)then ! TODO add other particle BCs
          toBeRemoved(i)=.true.
        end if
        if(info==PUSH_LOST)then
          write(*,*)'[W] PARTICLE LOST AT <x,y,z>:'
          write(*,*)p(j)%x(:,i)
          toBeRemoved(i)=.true.
        end if
      end do
      do i=p(j)%n,1,-1 ! reversed loop to ensure correct removal
        if(toBeRemoved(i))then
          call p(j)%rm(i)
        end if
      end do
      deallocate(toBeRemoved)
    end do
    
  end subroutine
  
  !> prepare data before writing the state (must run on all processes)
  subroutine preOut()
    use modMGS
    use mpi
    
    rhoLocal(:)=0d0
    do j=1,size(p)
      do i=1,p(j)%n
        call scatter(grid,p(j)%w(i)*p(j)%q/QE,p(j)%iC(i),p(j)%xx(:,i),rhoLocal)
      end do
    end do
    rhoLocal(:)=rhoLocal(:)/nVol(:)
    call mpi_reduce(rhoLocal,rho,size(rho),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  end subroutine
  
  !> write the state to post-processing file
  subroutine writeState(fName)
    use modFileIO
    character(*),intent(in)::fName
    integer,parameter::FID=10
    
    open(FID,file=trim(fName),action='write')
    call writeVTK(FID,grid)
    call writeVTK(FID,grid,N_DATA)
    call writeVTK(FID,'phi',phi)
    call writeVTK(FID,'rho',rho)
    call writeVTK(FID,'ef',ef)
    close(FID)
  end subroutine
  
end module
