!----------------------------------------------------------------------------- best with 100 columns

!> particle source
module modSource
  private
  
  integer,parameter::DIMS=3 !< three dimensions
  double precision,parameter::PI=4d0*atan(1d0) !< Pi
  double precision,parameter::EV2J=1.602176634d-19 !< 1 eV in J
  
  integer,public,parameter::SRC_SURFACE_FLUX=1 !< particle flux from surface
  
  !> structure to hold information about a particle source
  type,public::pSrc
    integer,public::t !< type
    integer,public::iSp !< species index
    double precision,public::w !< target specific weight
    double precision,public::f !< flux [1/m^2/s] or flow rate [1/s]
    double precision,public::u(DIMS) !< drifting velocity [m/s]
    double precision,public::temp !< temperature [eV] if Maxwellian
    integer,public::gid !< geometric id
    logical,allocatable,public::mask(:) !< cell mask for efficient initial cell match
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
  subroutine emitPSrc(this,p,grid,dt,nProc)
    use modParticle
    use modPICGrid
    use modPolyGrid
    use modGeometry
    use modMGS
    class(pSrc),intent(inout)::this !< this particle source
    type(ptcls),intent(inout)::p(:) !< the particles of all species
    type(PICGrid),intent(in)::grid !< the grid
    double precision,intent(in)::dt !< duration [s]
    integer,intent(in)::nProc !< number of processes (new particles equally split among processes)
    double precision::x(DIMS),u(DIMS),w,total,a,norm(DIMS),aa,r1,r2,r3,r4,xx(DIMS),vThermal
    integer::iC
    
    iC=0
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
          total=this%f*dt*a/nProc
          if(total<tiny(1d0)) cycle
          n=ceiling(total/this%w)
          w=total/n
          j=1
          do while(j<=n)
            vThermal=sqrt(2*this%temp*EV2J/p(this%iSp)%m)
            u=this%u+MaxwellianSpeed(vThermal)*randomDirection()
            if(dot_product(u,-norm)<0)then
              cycle
            end if
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
            call random_number(r1)
            x=x+r1*dt*u
            call match(grid,x,iC,xx,this%mask)
            if(iC>0)then
              call p(this%iSp)%add(x,u,w,iC=iC,xx=xx)
              j=j+1
            end if
          end do
        end if
      end do
    case default
    end select
  end subroutine
  
  !> return a Maxwellian random speed according to vThermal
  function MaxwellianSpeed(vThermal)
    double precision,intent(in)::vThermal !< thermal speed defined as sqrt(2kbT/m) [m/s]
    double precision::MaxwellianSpeed !< the result [m/s]
    integer,parameter::V_RANGE=6 !< considering speed up to 6*vThermal
    double precision::a,fMax,r,f
    
    a=4/sqrt(PI)*vThermal**(-3)
    fMax=a*vThermal**2*exp(-1d0)
    do while(.true.)
      call random_number(r)
      MaxwellianSpeed=r*vThermal*V_RANGE
      f=a*MaxwellianSpeed**2*exp(-MaxwellianSpeed**2/vThermal**2)
      call random_number(r)
      if(f/fMax>r) exit
    end do
  end function
  
  !> return a random direction unit vector in 3-D
  function randomDirection(vec)
    double precision::randomDirection(DIMS) !< the unit vector result
    double precision,optional,intent(in)::vec(DIMS) !< limit to the half dome selected by vec
    double precision::theta,r,a,d
    
    call random_number(theta)
    theta=theta*2*PI
    call random_number(r)
    r=-1d0+r*2
    a=sqrt(1d0-r**2)
    randomDirection(:)=[cos(theta)*a,sin(theta)*a,r]
    if(present(vec))then
      d=dot_product(vec,randomDirection)
      if(d<0)then
        randomDirection(:)=randomDirection(:)-2*d*vec(:)/dot_product(vec,vec)
      end if
    end if
  end function
  
  
end module
