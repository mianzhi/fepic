program fepic
  use modFEPIC
  use modFileIO
  use modSparse
  use modUglyFEM
  
  call init()
  
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
