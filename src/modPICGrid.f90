!----------------------------------------------------------------------------- best with 100 columns

!> FEM PIC grid
module modPICGrid
  use modPolyEdgeFeGrid
  private
  
  !> edge-based (and nodal) FEM grid with auxiliary features for PIC
  type,extends(polyEdgeFeGrid),public::PICGrid
  end type
  
end module
