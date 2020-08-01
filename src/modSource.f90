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
  end type
  
  ! public procedures
  public::readPSrc
  
contains
  
  !> read particle sources from fid
  subroutine readPSrc(fid,src)
    integer,intent(in)::fid !< file id
    class(pSrc),intent(inout)::src !< particle source
    
    read(fid,*)
    read(fid,*)src%t
    select case(src%t)
    case(SRC_SURFACE_FLUX)
      read(fid,*)src%iSp
      read(fid,*)src%w
      read(fid,*)src%gid
      read(fid,*)src%u(:)
      read(fid,*)src%f
      read(fid,*)src%temp
    case default
    end select
  end subroutine
  
end module
