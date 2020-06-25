!----------------------------------------------------------------------------- best with 100 columns

!> environment for fepic
module modFEPIC
  use modParticle
  use modPICGrid
  use modCondition
  use modSparse
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  integer,parameter::BC_PHI_DEFAULT=0 !< default BC: Neumann with zero normal phi gradient
  integer,parameter::BC_PHI_DIRICHLET=1 !< Dirichlet BC with given potential
  
  integer::iProc,nProc,ierr !< mpi variables
  
  type(ptcls),allocatable::p(:) !< particles of each species
  type(PICGrid)::grid !< the grid
  type(condTab)::ebc !< electric field boundary conditions
  type(multiFront)::PoissonPhi !< FEM Poisson equation for phi
  
  double precision,allocatable::phi(:) !< electric potential [V]
  double precision,allocatable::ef(:,:) !< electric field strength [V/m]
  double precision,allocatable::rhsPhi(:) !< RHS of the phi equation [C], and [V] (Dirichlet nodes)
  double precision,allocatable::rhsPhiDi(:) !< RHS of the phi equation, only Dirichlet nodes [V]
  double precision,allocatable::rhsPhiLocal(:) !< RHS of the phi equation, local version [C]
  double precision,allocatable::nVol(:) !< FEM nodal volume
  integer,allocatable::iEBC(:) !< indexes electric field boundary conditions
  logical,allocatable::isDirichlet(:) !< whether a point is a Dirichlet point
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    use modCondition
    use modSparse
    use modBasicFEM
    integer,parameter::FID=10
    
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
    allocate(nVol(grid%nN))
    allocate(isDirichlet(grid%nN))
    allocate(ef(3,grid%nN))
    ! read BCs
    open(FID,file='ebc',action='read')
      call readCondTab(FID,ebc)
    close(FID)
    call mapCondTab(grid,ebc,iEBC)
    rhsPhiDi(:)=0d0
    isDirichlet(:)=.false.
    do i=grid%nC+1,grid%nE
      if(iEBC(i)>0)then
        select case(ebc%t(iEBC(i)))
        case(BC_PHI_DIRICHLET) !< Dirichlet with given phi
          do j=1,grid%nNE(i)
            isDirichlet(grid%iNE(j,i))=.true.
            rhsPhiDi(grid%iNE(j,i))=ebc%p(1,iEBC(i))
          end do
        case default
        end select
      end if
    end do
    
  end subroutine
  
end module
