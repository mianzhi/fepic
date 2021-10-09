!----------------------------------------------------------------------------- best with 100 columns

!> FEM operations that you would avoid be needed in this PIC code
module modUglyFEM
  public
  
  integer,private,parameter::DIMS=3
  
contains
  
  !> find nodal gradient of phi
  subroutine findNodalGrad(grid,phi,nVol,grad)
    use modPICGrid
    use modPolyFeGrid
    use modPolyGrid
    use modSMQ
    class(PICGrid),intent(inout)::grid !< the grid
    double precision,intent(in)::phi(:) !< the potential
    double precision,intent(in)::nVol(:) !< nodal volume, which can be found by findVolSrc()
    double precision,allocatable,intent(inout)::grad(:,:) !< the gradient of phi
    integer,parameter::MAX_QP=12 !< max number of quadrature points per element
    double precision::gradQP(DIMS,MAX_QP)
    
    call grid%up()
    if(.not.allocated(grad))then
      allocate(grad(DIMS,grid%nN))
    end if
    grad(:,:)=0d0
    do i=1,grid%nC
      select case(grid%sE(i))
      case(TET)
        gradQP(:,:)=0d0
        do j=1,TET_N
          do l=1,size(TET_QW)
            gradQP(:,l)=gradQP(:,l)+matmul(grid%invJ(:,:,l,i),TET_GRAD_QP(:,j,l))&
            &                       *phi(grid%iNE(j,i))
          end do
        end do
        do j=1,TET_N
          do l=1,size(TET_QW)
            grad(:,grid%iNE(j,i))=grad(:,grid%iNE(j,i))+gradQP(:,l)*TET_SHAPE_QP(j,l)&
            &                                           *grid%detJ(l,i)*TET_QW(l)
          end do
        end do
      case(TET10)
        gradQP(:,:)=0d0
        do j=1,TET10_N
          do l=1,size(TET10_QW)
            gradQP(:,l)=gradQP(:,l)+matmul(grid%invJ(:,:,l,i),TET10_GRAD_QP(:,j,l))&
            &                       *phi(grid%iNE(j,i))
          end do
        end do
        do j=1,TET10_N
          do l=1,size(TET10_QW)
            grad(:,grid%iNE(j,i))=grad(:,grid%iNE(j,i))+gradQP(:,l)*TET10_SHAPE_QP(j,l)&
            &                                           *grid%detJ(l,i)*TET10_QW(l)
          end do
        end do
      case default
      end select
    end do
    forall(i=1:grid%nN)
      grad(:,i)=grad(:,i)/nVol(i)
    end forall
  end subroutine
  
  !> inverse iso-parametric map
  !>   NOTE: assuming grid is updated, i.e. called grid%up(), before calling this subroutine
  pure subroutine x2xx(grid,iC,x,xx,isInside)
    use modPICGrid
    use modPolyGrid
    use modSMQ
    class(PICGrid),intent(in)::grid !< the grid
    integer,intent(in)::iC !< cell index
    double precision,intent(in)::x(DIMS) !< global location to be mapped [m]
    double precision,intent(inout)::xx(DIMS) !< location in the reference element of cell iC
    logical,optional,intent(inout)::isInside !< whether x is inside cell iC
    double precision,parameter::TOL=0d0 !< tolerance in the reference coordinate
    
    select case(grid%sE(iC))
    case(TET)
      xx=matmul(transpose(grid%invJ(:,:,1,iC)),x-grid%pN(:,grid%iNE(1,iC)))
      ! matrix identity: inv(A^T)=inv(A)^T
      if(present(isInside))then
        isInside=xx(1)>=-TOL.and.xx(2)>=-TOL.and.xx(3)>=-TOL.and.xx(1)+xx(2)+xx(3)<=1d0+TOL
      end if
    case(TET10)
      ! TODO Newton iteration: matmul(grid%pN(:,grid%iNE(1:TET10_N,iC)),shapeTet10(xx))-x==[0,0,0]
    case default
    end select
  end subroutine
  
end module
