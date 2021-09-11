program fepic
  use mpi
  use modFEPIC
  use modFileIO
  use modSparse
  use modUglyFEM
  use modBasicFEM
  
  use modParticle
  use modPush
  type(ptcls)::pp
  integer::info,f
   
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
  do while(t<tFinal)
    ! push particles and deposit again
    call stepPtcls()
  !  call mpi_reduce(rhsPhiLocal,rhsPhi,size(rhsPhi),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,&
  !  &               ierr)
    ! solve field at process 0
  !  if(iProc==0)then
  !    rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
  !    call PoissonPhi%solve(rhsPhi,phi)
  !    call findNodalGrad(grid,phi,nVol,ef)
  !    ef(:,:)=-ef(:,:)
  !  end if
  !  call mpi_bcast(ef,size(ef),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    t=t+dt
  end do
  
  if(iProc==0)then
    call pp%init(16*1.660538921d-27,1.602d-19)
    call pp%add([0.08d0,0.08d0,0d0],[0d0,0d0,7000d0],1d0)
    t=1d-7
    call push(grid,pp,1,t,phi,f,LF_REWIND,info)
    do while(info==PUSH_DONE)
      t=1d-7
      call push(grid,pp,1,t,phi,f,LF_NORMAL,info)
      write(*,*)pp%x(:,1)
    end do
    write(*,*)info,f,t,iPtclBC(f)==BC_PTCL_DEFAULT
    
    open(10,file='rst.vtk',action='write')
    call writeVTK(10,grid)
    call writeVTK(10,grid,N_DATA)
    call writeVTK(10,'phi',phi)
    call writeVTK(10,'ef',ef)
    close(10)
  end if
  
  call mpi_finalize(ierr)
  
end program
