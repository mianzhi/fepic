!----------------------------------------------------------------------------- best with 100 columns

!> particles
module modParticle
  private
  
  integer,parameter::DIMS=3 !< three dimensions
  
  !> collection of (macro) particles of a same species
  type,public::ptcls
    integer,public::n !< number of particles
    integer,public::sz !< storage size
    double precision,allocatable,public::x(:,:) !< particle locations [m]
    double precision,allocatable,public::v(:,:) !< particle velocities [m/s]
    double precision,allocatable,public::w(:) !< weight (# of real particles per macro particle)
    double precision,allocatable,public::xx(:,:) !< particle locations in reference element
    integer,allocatable,public::iC(:) !< index of the cell in which the particle is last found
    double precision,public::m !< particle mass [kg]
    double precision,public::q !< charge [Coulomb]
  contains
    procedure,public::init=>initPtcls
    procedure,public::clear=>clearPtcls
    final::purgePtcls
  end type
  
contains
  
  !> initialize this ptcls
  subroutine initPtcls(this,m,q,sz)
    class(ptcls),intent(inout)::this !< this ptcls
    double precision,intent(in)::m !< particle mass [kg]
    double precision,intent(in)::q !< charge [Coulomb]
    integer,intent(in),optional::sz !< storage size
    integer,parameter::DEFAULT_SZ=10000
    
    call this%clear()
    if(present(sz))then
      this%sz=sz
    else
      this%sz=DEFAULT_SZ
    end if
    allocate(this%x(DIMS,this%sz),source=0d0)
    allocate(this%v(DIMS,this%sz),source=0d0)
    allocate(this%w(this%sz),source=0d0)
    allocate(this%xx(DIMS,this%sz),source=0d0)
    allocate(this%iC(this%sz),source=0)
    this%n=0
    this%m=m
    this%q=q
  end subroutine
  
  !> clear this ptcls
  subroutine clearPtcls(this)
    class(ptcls),intent(inout)::this !< this ptcls
    
    if(allocated(this%x)) deallocate(this%x)
    if(allocated(this%v)) deallocate(this%v)
    if(allocated(this%w)) deallocate(this%w)
    if(allocated(this%xx)) deallocate(this%xx)
    if(allocated(this%iC)) deallocate(this%iC)
    this%n=0
    this%sz=0
  end subroutine
  
  !> purge this ptcls
  subroutine purgePtcls(this)
    type(ptcls),intent(inout)::this !< this ptcls
    
    call this%clear()
  end subroutine
  
end module