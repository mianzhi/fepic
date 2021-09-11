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
  pure subroutine emitPSrc(this,p,grid,dt)
    use modParticle
    use modPICGrid
    use modPolyGrid
    use modGeometry
    class(pSrc),intent(in)::this !< this particle source
    type(ptcls),intent(inout)::p(:) !< the particles of all species
    type(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::dt !< duration [s]
    double precision::x(DIMS),u(DIMS),w,total,a,norm(DIMS)
    
    select case(this%t)
    case(SRC_SURFACE_FLUX)
      do i=grid%nC+1,grid%nE
        if(grid%gid(i)==this%gid)then
          select case(grid%sE(i))
          case(TRI)
            a=a3p(grid%pN(:,grid%iNE(1:3,i)))
            norm(:)=n3p(grid%pN(:,grid%iNE(1:3,i)))
          case default
            a=0d0
            norm=0d0
          end select
          a=a*abs(dot_product(norm,this%u)/norm2(this%u))
          total=this%f*dt*a
          if(total<tiny(1d0)) cycle
          n=ceiling(total/this%w)
          w=total/n
          do j=1,n
            x=grid%p(:,i) ! FIXME: randomized x and Maxiwellian u
            u=this%u
            call p(this%iSp)%add(x,u,w)
          end do
        end if
      end do
    case default
    end select
  end subroutine
  
end module
