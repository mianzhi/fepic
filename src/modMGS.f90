!----------------------------------------------------------------------------- best with 100 columns

!> matching, gathering, and scattering between particle and cell
module modMGS
  public
  
  integer,parameter,private::DIMS=3 !< three dimensions
  
  !> generic match
  interface match
    module procedure::matchDefault
  end interface
  
  !> generic gather
  interface gather
    module procedure::gatherScal
  end interface
  
  !> generic scatter
  interface scatter
    module procedure::scatterScal
  end interface
  
contains
  
  !> match (particle) location x with cell iC in grid
  !>   iC=0 at input: test all cells, otherwise test only neighbors of iC
  !>   iC=0 at output: no matching cell found
  pure subroutine matchDefault(grid,x,iC,xx)
    use modPICGrid
    use modUglyFEM
    class(PICGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::x(DIMS) !< global location [m]
    integer,intent(inout)::iC !< initial cell at input, matched cell at output
    double precision,intent(inout),optional::xx(DIMS) !< location at reference cell of iC
    double precision::xxx(DIMS)
    logical::neibOnly,isInside
    
    neibOnly=iC/=0
    if(neibOnly)then
      m=iC
      iC=0
      call x2xx(grid,m,x,xxx,isInside)
      if(isInside)then
        iC=m
      else
        do i=1,grid%nNeibC(m)
          n=grid%iNeibC(i,m)
          call x2xx(grid,n,x,xxx,isInside)
          if(isInside)then
            iC=n
            exit
          end if
        end do
      end if
    else
      iC=0
      do i=1,grid%nC
        call x2xx(grid,i,x,xxx,isInside)
        if(isInside)then
          iC=i
          exit
        end if
      end do
    end if
    if(iC/=0.and.present(xx))then
      xx=xxx
    end if
  end subroutine
  
  !> gather grid nodal scalar f at reference coordinate xx of cell iC
  pure subroutine gatherScal(grid,f,iC,xx,rst)
    use modPICGrid
    class(PICGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::f(grid%nN) !< the nodal scalar field
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::xx(DIMS) !< location at reference cell of iC
    double precision,intent(inout)::rst !< the local value of f
    
    
  end subroutine
  
  !> scatter extensive quantity ap at reference coordinate xx of cell iC onto nodal grid value a
  pure subroutine scatterScal(grid,ap,iC,xx,a)
    use modPICGrid
    class(PICGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::ap !< the point quantity
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::xx(DIMS) !< location at reference cell of iC
    double precision,intent(inout)::a(grid%nN) !< the nodal grid value
    
    
  end subroutine
  
end module
