program fepic
  use mpi
  use modFEPIC
  use modFileIO
  use modSparse
  use modUglyFEM
  use modBasicFEM
  use modMGS
  
  character(20)::tmpStr
  
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,iProc,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nProc,ierr)
  
  ! initialization at all processes
  call init()
  
  ! initial FEM setup and factorization at process 0
  if(iProc==0)then
    call phiLinEq%init(grid%nN,size(grid%iNE,1)**2*grid%nE)
    call findLaplacian(grid,phiLinEq,isDirichlet)
    call phiLinEq%fact()
    call findVolSrc(grid,nVol)
  end if
  call mpi_bcast(nVol,size(nVol),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  ! initial particle loading and deposition
  qDepoLocal(:)=0d0
  call mpi_reduce(qDepoLocal,qDepo,size(qDepo),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(iProc==0)then
    do i=1,grid%nN
      rhsPhi(i)=qDepo(i) ! no electron model initially
    end do
  end if
  
  ! solve field at process 0
  if(iProc==0)then
    rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
    call phiLinEq%solve(rhsPhi,phi)
    call findNodalGrad(grid,phi,nVol,ef)
    ef(:,:)=-ef(:,:)
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
    
    ! deposit ions, reduce RHS, apply electron model
    qDepoLocal(:)=0d0
    do j=1,size(p)
      do i=1,p(j)%n
        call scatter(grid,p(j)%w(i)*p(j)%q/EPS0,p(j)%iC(i),p(j)%xx(:,i),qDepoLocal)
      end do
    end do
    call mpi_reduce(qDepoLocal,qDepo,size(qDepo),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(iProc==0)then
      do i=1,grid%nN
        rhsPhi(i)=qDepo(i)-QE/EPS0*ne0BR*exp((phi(i)-phi0BR)/kbTeBR)*nVol(i)
      end do
    end if
    
    ! solve field
    if(iProc==0)then
      rhsPhi=merge(rhsPhiDi,rhsPhi,isDirichlet)
      call phiLinEq%solve(rhsPhi,phi)
      call findNodalGrad(grid,phi,nVol,ef)
      ef(:,:)=-ef(:,:)
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
