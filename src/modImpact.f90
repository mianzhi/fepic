!----------------------------------------------------------------------------- best with 100 columns

!> impact of particle on facets
module modImpact
  private
  
  integer,parameter::DIMS=3 !< three dimensions
  
  public::testImpact
  
contains
  
  !> see if particle k impacts with any neighbor facets if moved to location newX
  pure subroutine testImpact(grid,p,k,newX,f,h)
    use modPICGrid
    use modPolyGrid
    use modParticle
    class(PICGrid),intent(in)::grid !< the grid
    class(ptcls),intent(in)::p !< the particle collection
    integer,intent(in)::k !< the particle to be tested
    double precision,intent(inout)::newX(DIMS) !< new location of particle
    integer,intent(out)::f !< facet index
    double precision,intent(out)::h !< fraction of step at impact
    double precision::s(DIMS,2),t(DIMS,3),tmp,tmpX(DIMS)
    logical::isPenetrate
    
    f=0
    h=2d0 ! any value greater than 1
    s(:,1)=p%x(:,k)
    s(:,2)=newX(:)
    do i=1,grid%nNeibF(p%iC(k))
      n=grid%iNeibF(i,p%iC(k))
      select case(grid%sE(n))
      case(TRI)
        t(:,:)=grid%pN(:,grid%iNE(1:N_TRI,n))
        call testPenetration(s,t,isPenetrate,tmp,tmpX)
        if(isPenetrate)then
          if(tmp<h)then
            f=n
            h=tmp
            newX=tmpX
          end if
        end if
      case(TRI6)
      case default
      end select
    end do
  end subroutine
  
  !> see if line segment s penetrates triangle t, return the point of penetration
  pure subroutine testPenetration(s,t,isPenetrate,h,p)
    double precision,intent(in)::s(DIMS,2) !< two points defining the line segment
    double precision,intent(in)::t(DIMS,3) !< three points defining the line segment
    logical,intent(out)::isPenetrate !< whether s penetrates t
    double precision,intent(inout)::h !< location of the penetration point on the segment (0 to 1)
    double precision,intent(inout)::p(DIMS) !< location of the penetration point
    double precision::A(DIMS,DIMS),invA(DIMS,DIMS),rhs(DIMS),sol(DIMS)
    logical::isSingular
    
    isPenetrate=.false.
    A(:,1)=s(:,2)-s(:,1)
    A(:,2)=t(:,1)-t(:,2)
    A(:,3)=t(:,1)-t(:,3)
    rhs(:)=t(:,1)-s(:,1)
    call findInv3by3(A,invA,isSingular)
    if(isSingular)then
      isPenetrate=.false.
    else
      sol=matmul(invA,rhs)
      if(all(sol>=0d0).and.sol(1)<=1d0.and.sol(2)+sol(3)<=1d0)then
        isPenetrate=.true.
        h=sol(1)
        p=s(:,1)+h*(s(:,2)-s(:,1))
      else
        isPenetrate=.false.
      end if
    end if
  end subroutine
  
  !> find inverse of 3 by 3 matrix
  pure subroutine findInv3by3(A,invA,isSingular)
    double precision,intent(in)::A(3,3) !< 3 by 3 matrix
    double precision,intent(inout)::invA(3,3) !< inverse of A
    logical,intent(out)::isSingular !< whether detA is 0
    double precision::detA
    
    detA=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)&
    &    +A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
    if(abs(detA)>tiny(1d0))then
      invA(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      invA(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      invA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      invA(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      invA(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      invA(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      invA(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      invA(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      invA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
      invA(:,:)=invA(:,:)/detA
      isSingular=.false.
    else
      isSingular=.true.
    end if
  end subroutine
  
end module
