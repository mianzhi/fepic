!----------------------------------------------------------------------------- best with 100 columns

!> particle pushing
module modPush
  public
  
  integer,parameter,private::DIMS=3 !< three dimensions
  
  integer,parameter::PUSH_DONE=0 !< pushed through the required time without impact
  integer,parameter::PUSH_LOST=-1 !< particle lost, impact facet is not identified
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
    use modParticle
    use modMGS
    class(PICGrid),intent(in)::grid !< the grid
    class(ptcls),intent(inout)::p !< the particle collection
    integer,intent(in)::k !< the particle to be pushed
    double precision,intent(inout)::t !< time left [s]
    double precision,intent(in)::phi(grid%nN) !< the nodal electric potential [V]
    integer,intent(inout)::f !< the facet of impact
    integer,intent(in)::job !< rewind and synchronize
    integer,intent(out)::info !< returning status
    logical::doRewind,doSync
    double precision::x(DIMS),xx(DIMS),v(DIMS),a(DIMS)
    integer::iC
    
    doRewind=(job==LF_REWIND.or.job==LF_REWIND_SYNC)
    doSync=(job==LF_SYNC.or.job==LF_REWIND_SYNC)
    call gather(grid,phi,p%iC(k),p%xx(:,k),a)
    a=-a*p%q/p%m
    if(doRewind)then
      v=p%v(:,k)+0.5d0*a*t
    else
      v=p%v(:,k)+a*t
    end if
    x=p%x(:,k)+v*t
    iC=p%iC(k)
    call match(grid,x,iC,xx)
    if(iC>0)then ! push ended within the neighborhood (allows briefly flying out of neighborhood)
      p%x(:,k)=x
      p%v(:,k)=v
      p%iC(k)=iC
      p%xx(:,k)=xx
      if(doSync)then
        call gather(grid,phi,p%iC(k),p%xx(:,k),a)
        a=-a*p%q/p%m
        p%v(:,k)=p%v(:,k)+0.5d0*a*t
      end if
      info=PUSH_DONE
    else ! push ended out of the neighborhood
      info=PUSH_LOST
    end if
  end subroutine
  
end module
