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
  pure subroutine emitPSrc(this,p)
    use modParticle
    class(pSrc),intent(in)::this !< this particle source
    type(ptcls),intent(inout)::p(:) !< the particles of all species
    
  end subroutine
  
end module
