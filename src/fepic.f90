program fepic
  use mpi
  use modFEPIC
  use modFileIO
  use modSparse
  use modUglyFEM
  use modBasicFEM
  
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,iProc,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nProc,ierr)
  
  ! initialization at all processes
  call init()
  
  ! initial FEM setup and factorization at process 0
  if(iProc==0)then
    call PoissonPhi%init(grid%nN,size(grid%iNE,1)**2*grid%nE)
    call findLaplacian(grid,PoissonPhi,isDirichlet)
    call PoissonPhi%fact()
    call findVolSrc(grid,nVol)
  end if
  
  ! initial particle loading and deposition
  rhsPhiLocal(:)=0d0
  call mpi_reduce(rhsPhiLocal,rhsPhi,size(rhsPhi),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,&
  &               ierr)
  
  ! solve field at process 0
  if(iProc==0)then
    rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
    call PoissonPhi%solve(rhsPhi,phi)
    call findNodalGrad(grid,phi,nVol,ef)
    ef(:,:)=-ef(:,:)
  end if
  call mpi_bcast(ef,size(ef),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  ! time loop
  !do l=1,1
  !  ! push particles and deposit again
  !  call mpi_reduce(rhsPhiLocal,rhsPhi,size(rhsPhi),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,&
  !  &               ierr)
  !  ! solve field at process 0
  !  if(iProc==0)then
  !    rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
  !    call PoissonPhi%solve(rhsPhi,phi)
  !    call findNodalGrad(grid,phi,nVol,ef)
  !    ef(:,:)=-ef(:,:)
  !  end if
  !  call mpi_bcast(ef,size(ef),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !end do
  
  if(iProc==0)then
    open(10,file='rst.vtk',action='write')
    call writeVTK(10,grid)
    call writeVTK(10,grid,N_DATA)
    call writeVTK(10,'phi',phi)
    call writeVTK(10,'ef',ef)
    close(10)
  end if
  
  call mpi_finalize(ierr)
  
end program
