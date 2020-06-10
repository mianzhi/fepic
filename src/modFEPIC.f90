!----------------------------------------------------------------------------- best with 100 columns

!> environment for fepic
module modFEPIC
  use modPICGrid
  use modSparse
  
  public
  
  type(PICGrid)::grid !< the grid
  
  double precision,allocatable::phi(:) !< electric potential [V]
  double precision,allocatable::ef(:,:) !< electric field strength [V/m]
  
  type(multiFront)::PoissonPhi !< FEM Poisson equation for phi
  
  double precision,allocatable::rhsPhi(:) !< RHS of the phi equation [C], and [V] (Dirichlet nodes)
  double precision,allocatable::nVol(:) !< FEM nodal volume
  logical,allocatable::isDirichlet(:) !< whether a point is a Dirichlet point
  
contains
  
  !> initialize the simulation
  subroutine init()
    use modFileIO
    integer,parameter::FID=10
    
    ! read inputs
    open(FID,file='grid.msh',action='read')
      call readGMSH(FID,grid)
    close(FID)
    call grid%up()
    ! initialize algebraic solver
    call PoissonPhi%init(grid%nN,size(grid%iNE,1)**2*grid%nE)
    ! work space and initial state
    allocate(phi(grid%nN))
    allocate(rhsPhi(grid%nN))
    allocate(nVol(grid%nN))
    allocate(isDirichlet(grid%nN))
    allocate(ef(3,grid%nN))
  end subroutine
  
end module
