!----------------------------------------------------------------------------- best with 100 columns

!> particle source
module modSource
  private
  
  integer,parameter::DIMS=3 !< three dimensions
  
  integer,public,parameter::SRC_SURFACE_FLUX=1 !< particle flux from surface
  
  !> structure to hold information about a particle source
  type,public::pSrc
    integer,public::t !< type
    integer,public::iSp !< species index
    double precision,public::w !< target specific weight
    double precision,public::f !< flux [1/m^2/s] or flow rate [1/s]
    double precision,public::u(DIMS) !< drifting velocity [m/s]
    double precision,public::temp !< temperature [K] if Maxwellian
    integer,public::gid !< geometric id
  contains
    procedure,public::load=>readPSrc
    procedure,public::emit=>emitPSrc
  end type
  
contains
  
  !> read into this particle sources from fid
  subroutine readPSrc(this,fid)
    class(pSrc),intent(inout)::this !< this particle source
    integer,intent(in)::fid !< file id
    
    read(fid,*)
    read(fid,*)this%t
    select case(this%t)
    case(SRC_SURFACE_FLUX)
      read(fid,*)this%iSp
      read(fid,*)this%w
      read(fid,*)this%gid
      read(fid,*)this%u(:)
      read(fid,*)this%f
      read(fid,*)this%temp
    case default
    end select
  end subroutine
  
  !> emit particles from this particle source
  subroutine emitPSrc(this,p,grid,dt)
    use modParticle
    use modPICGrid
    use modPolyGrid
    use modGeometry
    class(pSrc),intent(in)::this !< this particle source
    type(ptcls),intent(inout)::p(:) !< the particles of all species
    type(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::dt !< duration [s]
    double precision,parameter::INTO_DOMAIN=1d-6 !< to keep new particles inside the domain
    double precision::x(DIMS),u(DIMS),w,total,a,norm(DIMS),aa,r1,r2,r3,r4
    
    select case(this%t)
    case(SRC_SURFACE_FLUX)
      do i=grid%nC+1,grid%nE
        if(grid%gid(i)==this%gid)then
          select case(grid%sE(i))
          case(TRI)
            a=a3p(grid%pN(:,grid%iNE(1:3,i)))
            aa=a
            norm(:)=n3p(grid%pN(:,grid%iNE(1:3,i)))
          case default
            a=0d0
            aa=a
            norm=0d0
          end select
          a=a*abs(dot_product(norm,this%u)/norm2(this%u))
          total=this%f*dt*a
          if(total<tiny(1d0)) cycle
          n=ceiling(total/this%w)
          w=total/n
          do j=1,n
            select case(grid%sE(i))
            case(TRI)
              call random_number(r1)
              call random_number(r2)
              r3=min(r1,r2)
              r4=max(r1,r2)
              r1=r4-r3
              r2=1d0-r4
              x=r1*grid%pN(:,grid%iNE(1,i))+r2*grid%pN(:,grid%iNE(2,i))+r3*grid%pN(:,grid%iNE(3,i))
            case default
              x(:)=0d0
            end select
            x=x-norm(:)*sqrt(aa)*INTO_DOMAIN
            u=this%u ! TODO: Maxiwellian u
            call p(this%iSp)%add(x,u,w)
          end do
        end if
      end do
    case default
    end select
  end subroutine
  
end module
