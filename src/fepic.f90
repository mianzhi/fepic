program fepic
  use mpi
  use modFEPIC
  use modFileIO
  use modSparse
  use modBasicFEM
  use modMGS
  
  character(20)::tmpStr
  integer::info
  
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,iProc,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nProc,ierr)
  
  ! initialization at all processes
  call init()
  
  ! initial FEM setup at process 0
  if(iProc==0)then
    call phiEq%init(grid%nN,phiRes,pSet=phiPSet,pSol=phiPSol)
    call phiLinEq%init(grid%nN,size(grid%iNE,1)**2*grid%nE,PHI_PREC_LFILL)
    call negLaPhi%init(grid%nN,size(grid%iNE,1)**2*grid%nE,PHI_PREC_LFILL)
    call findLaplacian(grid,negLaPhi,isDirichlet)
    call findVolSrc(grid,nVol)
  end if
  call mpi_bcast(nVol,size(nVol),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  ! initial particle loading and deposition
  qDepoLocal(:)=0d0
  call mpi_reduce(qDepoLocal,qDepo,size(qDepo),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  ! solve field at process 0
  if(iProc==0)then
    phi(:)=rhsPhiDi(:) ! initial guess of phi
    call phiEq%solve(phi,info=info)
  end if
  call mpi_bcast(phi,size(phi),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  ! write initial state
  call preOut()
  if(iProc==0)then
    write(tmpStr,*)iOut
    write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
    call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
  end if
  
  ! time loop
  do while(t<tFinal)
  
    ! push particles
    call stepPtcls()
    call mpi_reduce(sum(p(:)%n),n,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(iProc==0) write(*,'(a,g12.6,a,i9)')"[i] t = ",t,": particle count ",n
    
    ! deposit particle charge
    qDepoLocal(:)=0d0
    do j=1,size(p)
      do i=1,p(j)%n
        call scatter(grid,p(j)%w(i)*p(j)%q/EPS0,p(j)%iC(i),p(j)%xx(:,i),qDepoLocal)
      end do
    end do
    call mpi_reduce(qDepoLocal,qDepo,size(qDepo),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    ! solve field
    if(iProc==0)then
      call phiEq%solve(phi,info=info)
    end if
    call mpi_bcast(phi,size(phi),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    t=t+dt
    
    ! write visualization
    call preOut()
    if(iProc==0.and.t+tiny(1d0)>=tOutNext)then
      iOut=iOut+1
      write(tmpStr,*)iOut
      write(*,'(a)')'[i] writing: rst_'//trim(adjustl(tmpStr))//'.vtk'
      call writeState('rst_'//trim(adjustl(tmpStr))//'.vtk')
      tOutNext=tOutNext+dtOut
    end if
    
  end do
  
  call mpi_finalize(ierr)
  
end program
