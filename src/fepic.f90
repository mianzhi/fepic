program fepic
  use modFEPIC
  use modPICGrid
  use modFileIO
  use modSparse
  use modBasicFEM
  use modUglyFEM
  
  call init()
  
  rhsPhi(:)=0d0
  isDirichlet(:)=.false.
  do i=grid%nC+1,grid%nE
    if(grid%gid(i)==4)then
      isDirichlet(grid%iNE(1:grid%nNE(i),i))=.true.
      rhsPhi((grid%iNE(1:grid%nNE(i),i)))=0d0
    end if
    if(grid%gid(i)==7)then
      isDirichlet(grid%iNE(1:grid%nNE(i),i))=.true.
      rhsPhi((grid%iNE(1:grid%nNE(i),i)))=-100d0
    end if
  end do
  call findLaplacian(grid,PoissonPhi,isDirichlet)
  call findVolSrc(grid,nVol)
  
  call PoissonPhi%fact()
  call PoissonPhi%solve(rhsPhi,phi)
  call findNodalGrad(grid,phi,nVol,ef)
  ef(:,:)=-ef(:,:)
  
  open(10,file='rst.vtk',action='write')
  call writeVTK(10,grid)
  call writeVTK(10,grid,N_DATA)
  call writeVTK(10,'phi',phi)
  call writeVTK(10,'ef',ef)
  close(10)
  
end program
