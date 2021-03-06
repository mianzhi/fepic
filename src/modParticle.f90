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
    integer,allocatable,public::iC(:) !< index of the cell in which the particle is found
    double precision,public::m !< particle mass [kg]
    double precision,public::q !< charge [Coulomb]
  contains
    procedure,public::init=>initPtcls
    procedure,public::clear=>clearPtcls
    procedure,public::expand=>expandPtcls
    procedure,public::add=>addPtcls
    procedure,public::rm=>rmPtcls
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
  
  !> expand the storage of this ptcls
  pure subroutine expandPtcls(this,sz)
    class(ptcls),intent(inout)::this !< this ptcls
    integer,intent(in)::sz !< the new storage size
    double precision,allocatable::tmpD1(:,:),tmpD2(:)
    integer,allocatable::tmpI1(:)
    
    allocate(tmpD1(DIMS,sz),source=0d0)
    allocate(tmpD2(sz),source=0d0)
    allocate(tmpI1(sz),source=0)
    tmpD1(:,1:this%sz)=this%x(:,:)
    deallocate(this%x)
    allocate(this%x,source=tmpD1)
    tmpD1(:,1:this%sz)=this%v(:,:)
    deallocate(this%v)
    allocate(this%v,source=tmpD1)
    tmpD1(:,1:this%sz)=this%xx(:,:)
    deallocate(this%xx)
    allocate(this%xx,source=tmpD1)
    tmpD2(1:this%sz)=this%w(:)
    deallocate(this%w)
    allocate(this%w,source=tmpD2)
    tmpI1(1:this%sz)=this%iC(:)
    deallocate(this%iC)
    allocate(this%iC,source=tmpI1)
    this%sz=sz
    deallocate(tmpD1,tmpD2,tmpI1)
  end subroutine
  
  !> add a particle to the back of this ptcls
  pure subroutine addPtcls(this,x,v,w,iC,xx)
    class(ptcls),intent(inout)::this !< this ptcls
    double precision,intent(in)::x(DIMS) !< location [m]
    double precision,intent(in)::v(DIMS) !< velocity [m/s]
    double precision,intent(in)::w !< weight
    integer,intent(in),optional::iC !< cell index
    double precision,intent(in),optional::xx(DIMS) !< location in reference element
    real,parameter::EXPAND_FACT=1.5 !< expansion factor
    
    if(this%n==this%sz) call this%expand(int(this%sz*EXPAND_FACT))
    this%n=this%n+1
    m=this%n
    this%x(:,m)=x(:)
    this%v(:,m)=v(:)
    this%w(m)=w
    if(present(iC))then
      this%iC(m)=iC
      if(present(xx))then
        this%xx(:,m)=xx(:)
      else
        this%xx(:,m)=[0d0,0d0,0d0]
      end if
    else
      this%iC(m)=0
    end if
  end subroutine
  
  !> remove a particle from this ptcls, and move the last particle to the empty space
  pure subroutine rmPtcls(this,m)
    class(ptcls),intent(inout)::this !< this ptcls
    integer,intent(in)::m !< particle to be removed
    
    n=this%n
    if(m/=n)then
      this%x(:,m)=this%x(:,n)
      this%v(:,m)=this%v(:,n)
      this%w(m)=this%w(n)
      this%iC(m)=this%iC(n)
      this%xx(:,m)=this%xx(:,n)
    end if
    this%n=n-1
  end subroutine
  
  !> purge this ptcls
  subroutine purgePtcls(this)
    type(ptcls),intent(inout)::this !< this ptcls
    
    call this%clear()
  end subroutine
  
end module
