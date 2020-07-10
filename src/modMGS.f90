!----------------------------------------------------------------------------- best with 100 columns

!> matching, gathering, and scattering between particle and cell
module modMGS
  public
  
  integer,parameter,private::DIMS=3 !< three dimensions
  
  !> generic match
  interface match
    module procedure::matchSimple
  end interface
  
  !> generic gather
  interface gather
    module procedure::gatherScal
    module procedure::gatherGradPhi
  end interface
  
  !> generic scatter
  interface scatter
    module procedure::scatterScal
  end interface
  
contains
  
  !> match (particle) location x with cell iC in grid
  !>   iC=0 at input: test all cells, otherwise test only neighbors of iC
  !>   iC=0 at output: no matching cell found
  pure subroutine matchSimple(grid,x,iC,xx)
    use modPICGrid
    use modUglyFEM
    class(PICGrid),intent(in)::grid !< the grid
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
    use modPolyGrid
    use modSMQ
    class(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::f(grid%nN) !< the nodal scalar field
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::xx(DIMS) !< location at reference cell of iC
    double precision,intent(inout)::rst !< the local value of f
    integer,parameter::MAX_N=30 !< maximum number of node per element
    double precision::s(MAX_N)
    
    s(:)=0d0
    rst=0d0
    select case(grid%sE(iC))
    case(TET)
      s(1:grid%nNE(iC))=shapeTet(xx)
      rst=dot_product(f(grid%iNE(1:grid%nNE(iC),iC)),s(1:grid%nNE(iC)))
    case(TET10)
      s(1:grid%nNE(iC))=shapeTet10(xx)
      rst=dot_product(f(grid%iNE(1:grid%nNE(iC),iC)),s(1:grid%nNE(iC)))
    case default
    end select
  end subroutine
  
  !> gather gradient of a nodal potential phi at reference coordinate xx of cell iC
  pure subroutine gatherGradPhi(grid,phi,iC,xx,rst)
    use modPICGrid
    use modPolyGrid
    use modSMQ
    class(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::phi(grid%nN) !< the nodal potential field
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::xx(DIMS) !< location at reference cell of iC
    double precision,intent(inout)::rst(DIMS) !< the local grad(phi)
    
    rst(:)=0d0
    select case(grid%sE(iC))
    case(TET)
      rst=matmul(grid%invJ(:,:,1,iC),matmul(gradShapeTet(xx),phi(grid%iNE(1:grid%nNE(iC),iC))))
      ! making use of the fact that invJ is constant in the TET
    case(TET10)
      ! TODO
    case default
    end select
  end subroutine
  
  !> scatter extensive quantity ap at reference coordinate xx of cell iC onto nodal grid value a
  pure subroutine scatterScal(grid,ap,iC,xx,a)
    use modPICGrid
    use modPolyGrid
    use modSMQ
    class(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::ap !< the point quantity
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::xx(DIMS) !< location at reference cell of iC
    double precision,intent(inout)::a(grid%nN) !< the nodal grid value
    
    select case(grid%sE(iC))
    case(TET)
      a(grid%iNE(1:grid%nNE(iC),iC))=a(grid%iNE(1:grid%nNE(iC),iC))+ap*shapeTet(xx)
    case(TET10)
      a(grid%iNE(1:grid%nNE(iC),iC))=a(grid%iNE(1:grid%nNE(iC),iC))+ap*shapeTet10(xx)
    case default
    end select
  end subroutine
  
end module
