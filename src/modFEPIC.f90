!----------------------------------------------------------------------------- best with 100 columns

!> environment for fepic
module modFEPIC
  use modParticle
  use modPICGrid
  use modCondition
  use modSUNDIALS
  use modSparse
  use modSource
  
  public
  
  integer,parameter::DIMS=3 !< three dimensions
  
  double precision,parameter::AMU=1.6605390666050d-27 !< atomic mass unit (Dalton) [kg]
  double precision,parameter::QE=1.602176634d-19 !< elementary charge [C]
  double precision,parameter::EPS0=8.85418782d-12 !< vacuum permittivity [F/m]
  
  integer,parameter::BC_PHI_DEFAULT=0 !< default phi BC: Neumann with zero normal phi gradient
  integer,parameter::BC_PHI_DIRICHLET=1 !< Dirichlet BC with given potential
  
  integer,parameter::BC_PTCL_DEFAULT=0 !< default particle BC: removed upon impact
  
  integer,parameter::PHI_PREC_LFILL=50 !< phi equation ILU preconditioning fill limit
  
  integer::iProc,nProc,ierr !< mpi variables
  
  type(ptcls),allocatable::p(:) !< particles of each species
  type(PICGrid)::grid !< the grid
  type(condTab)::phibc !< electric field boundary conditions
  type(condTab)::ptclbc !< particle boundary conditions
  type(NewtonKrylov)::phiEq !< non-linear FEM phi equation
  type(ILU)::phiLinEq !< linearized FEM phi equation
  type(ILU)::negLaPhi !< negative Laplacian of phi (phiLinEq but without electron model part)
  type(pSrc),allocatable::ptclSrc(:) !< particle source
  
  integer::nSp !< number of particle species
  
  double precision,allocatable::phi(:) !< electric potential [V]
  double precision,allocatable::ef(:,:) !< electric field strength [V/m]
  double precision,allocatable::rhsPhi(:) !< RHS of the phi equation [C], and [V] (Dirichlet nodes)
  double precision,allocatable::rhsPhiDi(:) !< RHS of the phi equation, only Dirichlet nodes [V]
  double precision,allocatable::qDepoLocal(:) !< discretized ionDensity/EPS0, local version [C]
  double precision,allocatable::qDepo(:) !< discretized ionDensity/EPS0, reduced version [C]
  double precision,allocatable::den(:,:) !< species density [m^-3]
  double precision,allocatable::denLocal(:,:) !< species density, local version [m^-3]
  double precision,allocatable::wallIDen(:) !< wall hit current density [A/m^2]
  double precision,allocatable::wallHitDenLocal(:) !< wall hit charge density, local version [C/m^2]
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
  integer::nAvgOut !< number of time steps for averaged output
  
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
    
    ! work space
    allocate(phi(grid%nN))
    allocate(rhsPhi(grid%nN))
    allocate(rhsPhiDi(grid%nN))
    allocate(qDepoLocal(grid%nN))
    allocate(qDepo(grid%nN))
    allocate(den(grid%nN,size(p)))
    allocate(denLocal(grid%nN,size(p)))
    allocate(wallIDen(grid%nE))
    allocate(wallHitDenLocal(grid%nE))
    allocate(nVol(grid%nN))
    allocate(isDirichlet(grid%nN))
    allocate(ef(3,grid%nN))
    
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
    
    ! read electron model parameters
    open(FID,file='eMod',action='read')
      read(FID,*)ne0BR
      read(FID,*)phi0BR
      read(FID,*)kbTeBR
    close(FID)
    
    ! read simulation parameters
    open(FID,file='sim',action='read')
      read(FID,*)tFinal
      read(FID,*)dt
      read(FID,*)dtOut
    close(FID)
    
    t=0d0
    iOut=0
    tOutNext=dtOut
    nAvgOut=0
    
  end subroutine
  
  !> step particles
  subroutine stepPtcls()
    use modPush
    use modPICGrid
    use modPolyGrid
    use modGeometry
    integer::f,info
    double precision::h,a
    integer::n0(size(p))
    logical,allocatable::toBeRemoved(:)
    
    n0(:)=p(:)%n
    ! emission from particle sources
    do i=1,size(ptclSrc)
      call ptclSrc(i)%emit(p,grid,dt,nProc)
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
          select case(grid%sE(f))
          case(TRI)
            a=a3p(grid%pN(:,grid%iNE(1:3,f)))
          case default
            a=0d0
          end select
          if(a>0d0)then
            wallHitDenLocal(f)=wallHitDenLocal(f)+p(j)%w(i)*p(j)%q/a
          end if
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
    use modUglyFEM
    use mpi
    
    ! species density
    denLocal(:,:)=0d0
    do j=1,size(p)
      do i=1,p(j)%n
        call scatter(grid,p(j)%w(i),p(j)%iC(i),p(j)%xx(:,i),denLocal(:,j))
      end do
    end do
    forall(j=1:size(p))
      denLocal(:,j)=denLocal(:,j)/nVol(:)
    end forall
    call mpi_reduce(denLocal,den,size(den),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    ! electric field
    if(iProc==0)then
      call findNodalGrad(grid,phi,nVol,ef)
      ef(:,:)=-ef(:,:)
    end if
    
    ! wall hit current density
    call mpi_reduce(wallHitDenLocal,wallIDen,size(wallIDen),MPI_DOUBLE_PRECISION,MPI_SUM,0,&
    &               MPI_COMM_WORLD,ierr)
    if(nAvgOut>0)then
      wallIDen(:)=wallIDen(:)/dt/nAvgOut
    else
      wallIDen(:)=0d0
    end if
    wallHitDenLocal(:)=0d0
  end subroutine
  
  !> write the state to post-processing file
  subroutine writeState(fName)
    use modFileIO
    character(*),intent(in)::fName
    integer,parameter::FID=10
    character(20)::tmpStr
    double precision,allocatable::tmp1(:)
    
    if(.not.allocated(tmp1)) allocate(tmp1(grid%nE))
    
    open(FID,file=trim(fName),action='write')
    call writeVTK(FID,grid)
    call writeVTK(FID,grid,N_DATA)
    call writeVTK(FID,'phi',phi)
    call writeVTK(FID,'ef',ef)
    do i=1,size(p)
      write(tmpStr,*)i
      call writeVTK(FID,'den'//trim(adjustl(tmpStr)),den(:,i))
    end do
    call writeVTK(FID,grid,E_DATA)
    call writeVTK(FID,'geoID',dble(grid%gid))
    forall(i=1:grid%nE)
      tmp1(i)=merge(0d0,1d0,i<=grid%nC) ! 0 for cell, 1 for facet
    end forall
    call writeVTK(FID,'geoType',tmp1)
    call writeVTK(FID,'wallIDen',wallIDen)
    close(FID)
  end subroutine
  
  !> residual function of the FEM phi equation
  function phiRes(oldVector,newVector,dat)
    use iso_c_binding
    use ieee_arithmetic
    type(C_PTR),value::oldVector !< old N_Vector
    type(C_PTR),value::newVector !< new N_Vector
    type(C_PTR),value::dat !< optional user data object
    integer(C_INT)::phiRes !< error code
    double precision,pointer::x(:) !< fortran pointer associated with oldVector
    double precision,pointer::y(:) !< fortran pointer associated with newVector
    double precision,allocatable,save::tmp(:)
    
    call associateVector(oldVector,x)
    call associateVector(newVector,y)
    
    if(.not.allocated(tmp)) allocate(tmp(grid%nN))
    
    ! update rhsPhi to include electron model
    do i=1,grid%nN
      rhsPhi(i)=qDepo(i)-QE/EPS0*ne0BR*exp((x(i)-phi0BR)/kbTeBR)*nVol(i)
    end do
    rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
    
    ! construct residual vector
    call negLaPhi%mulVec(x,y)
    y(1:grid%nN)=y(1:grid%nN)-rhsPhi(1:grid%nN)
    
    phiRes=merge(1,0,any(ieee_is_nan(y).or.(.not.ieee_is_finite(y))))
    if(c_associated(dat))then
    end if
  end function
  
  !> phi preconditioning matrix setup/factor
  function phiPSet(c_x,c_xScale,c_f,c_fScale,dat)
    use iso_c_binding
    type(C_PTR),value::c_x,c_xScale,c_f,c_fScale,dat
    integer(C_INT)::phiPSet
    double precision,pointer::x(:),xScale(:),f(:),fScale(:)
    type(C_PTR)::foo
    double precision,parameter::DROPTOL=1e-6
    double precision,allocatable::dEm_dPhi(:) !< d (electron model term) / d phi
    
    call associateVector(c_x,x)
    call associateVector(c_xScale,xScale)
    call associateVector(c_f,f)
    call associateVector(c_fScale,fScale)
    
    if(.not.allocated(dEm_dPhi)) allocate(dEm_dPhi(grid%nN))
    
    ! construct preconditioning matrix
    call phiLinEq%setCSR(negLaPhi%iA,negLaPhi%jA,negLaPhi%A)
    forall(i=1:grid%nN)
      dEm_dPhi(i)=merge(0d0,QE/EPS0*ne0BR/kbTeBR*exp((x(i)-phi0BR)/kbTeBR)*nVol(i),isDirichlet(i))
    end forall
    call phiLinEq%addDiag(dEm_dPhi)
    ! factorize preconditioning matrix
    call phiLinEq%fact(DROPTOL)
    
    phiPSet=0
    foo=dat
  end function
  
  !> phi preconditioning problem solving
  function phiPSol(c_x,c_xScale,c_f,c_fScale,c_v,dat)
    use iso_c_binding
    type(C_PTR),value::c_x,c_xScale,c_f,c_fScale,c_v,dat
    integer(C_INT)::phiPSol
    double precision,pointer::x(:),xScale(:),f(:),v(:),fScale(:)
    type(C_PTR)::foo
    double precision,allocatable,save::tmp1(:),tmp2(:)
    
    call associateVector(c_x,x)
    call associateVector(c_xScale,xScale)
    call associateVector(c_f,f)
    call associateVector(c_fScale,fScale)
    call associateVector(c_v,v)
    
    if(.not.allocated(tmp1)) allocate(tmp1(grid%nN))
    if(.not.allocated(tmp2)) allocate(tmp2(grid%nN))
    
    ! find x increment from f increment v and save back to v
    tmp1(1:grid%nN)=v(1:grid%nN)
    call phiLinEq%solve(tmp1,tmp2)
    v(1:grid%nN)=tmp2(1:grid%nN)
    
    phiPsol=0
    foo=dat
  end function
  
end module
