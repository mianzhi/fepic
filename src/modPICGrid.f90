!----------------------------------------------------------------------------- best with 100 columns

!> FEM PIC grid
module modPICGrid
  use modPolyFeGrid
  private
  
  !> FEM grid with auxiliary features for PIC
  type,extends(polyFeGrid),public::PICGrid
  end type
  
end module
