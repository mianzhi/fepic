!----------------------------------------------------------------------------- best with 100 columns

!> particle pushing
module modPush
  public
  
  integer,parameter,private::DIMS=3 !< three dimensions
  
  integer,parameter::PUSH_DONE=0 !< pushed through the required time without impact
  integer,parameter::PUSH_IMPACT=1 !< impact with a facet
  
  integer,parameter::LF_NORMAL=0 !< no rewind or synchronize for leap-frog
  integer,parameter::LF_REWIND=1 !< do rewind for leap-frog
  integer,parameter::LF_SYNC=2 !< do synchronize for leap-frog
  integer,parameter::LF_REWIND_SYNC=3 !< do both rewind and synchronize for leap-frog
  
  !> generic push
  interface push
    module procedure::pushPhiLF
  end interface
  
contains
  
  !> push the particle k under the effect of electric potential phi with leap-frog scheme
  pure subroutine pushPhiLF(grid,p,k,t,phi,f,job,info)
    use modPICGrid
    use modMGS
    use modParticle
    class(PICGrid),intent(in)::grid !< the grid
    class(ptcls),intent(inout)::p !< the particle collection
    integer,intent(in)::k !< the particle to be pushed
    double precision,intent(inout)::t !< time left [s]
    double precision,intent(in)::phi(grid%nN) !< the nodal electric potential [V]
    integer,intent(out)::f !< the facet of impact
    integer,intent(in)::job !< rewind and synchronize
    integer,intent(out)::info !< returning status
    logical::doRewind,doSync
    double precision::x(DIMS),v(DIMS),a(DIMS),h
    
    f=0
    doRewind=(job==LF_REWIND.or.job==LF_REWIND_SYNC)
    doSync=(job==LF_SYNC.or.job==LF_REWIND_SYNC)
    if(p%iC(k)==0)then
      call match(grid,p%x(:,k),p%iC(k),p%xx(:,k))
    end if
    call gather(grid,phi,p%iC(k),p%xx(:,k),a)
    a=-a*p%q/p%m
    if(doRewind)then
      v=p%v(:,k)+0.5d0*a*t
    else
      v=p%v(:,k)+a*t
    end if
    x=p%x(:,k)+v*t
    call stride(grid,p,k,x,f,h,info)
    if(info==PUSH_DONE)then ! no impact
      p%v(:,k)=v
      if(doSync)then
        call gather(grid,phi,p%iC(k),p%xx(:,k),a)
        a=-a*p%q/p%m
        p%v(:,k)=p%v(:,k)+0.5d0*a*t
      end if
      t=0d0
    else ! impact
      p%v(:,k)=v*h+p%v(:,k)*(1d0-h)
      t=(1d0-h)*t
    end if
  end subroutine
  
  !> stride the particle k along a line segment while keeping track of iC and impact
  pure recursive subroutine stride(grid,p,k,targetX,f,h,info)
    use modPICGrid
    use modMGS
    use modImpact
    use modParticle
    class(PICGrid),intent(in)::grid !< the grid
    class(ptcls),intent(inout)::p !< the particle collection
    integer,intent(in)::k !< the particle to be tested
    double precision,intent(in)::targetX(DIMS) !< target location of particle
    integer,intent(out)::f !< facet index
    double precision,intent(out)::h !< fraction of step at impact
    integer,intent(out)::info !< returning status
    integer::iC
    double precision::midX(DIMS),x(DIMS),xx(DIMS)
    
    f=0
    h=0d0
    x(:)=targetX(:)
    iC=p%iC(k)
    call match(grid,x,iC,xx)
    if(iC>0)then ! target is within the neighborhood
      p%x(:,k)=x(:)
      p%iC(k)=iC
      p%xx(:,k)=xx(:)
      info=PUSH_DONE
    else
      call testImpact(grid,p,k,x,f,h)
      if(f>0)then ! impact detected within the neighborhood
        info=PUSH_IMPACT
        p%x(:,k)=x
      else
        midX(:)=0.5d0*(p%x(:,k)+targetX(:))
        call stride(grid,p,k,midX,f,h,info)
        if(info==PUSH_DONE)then
          call stride(grid,p,k,targetX,f,h,info)
        end if
      end if
    end if
  end subroutine
  
end module
