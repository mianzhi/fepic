!----------------------------------------------------------------------------- best with 100 columns

!> FEM PIC grid
module modPICGrid
  use modPolyFeGrid
  private
  
  !> FEM grid with auxiliary features for PIC
  type,extends(polyFeGrid),public::PICGrid
    ! note: *Neib* (neighbor) means sharing 1 or more nodes
    integer,allocatable::nNeibC(:),nNeibF(:) !< number of neighbor cells and facets of each cell
    integer,allocatable::iNeibC(:,:),iNeibF(:,:) !< index of neighbor cells and facets of each cell
  contains
    procedure,public::clear=>clearPICGrid
    procedure,public::up=>upPICGrid
    final::purgePICGrid
  end type
  
contains
  
  !> update this PICGrid
  subroutine upPICGrid(this)
    class(PICGrid),intent(inout)::this !< this PICGrid
    integer,allocatable::iEN(:,:),nEN(:),iNeibC(:,:),nNeibC(:),iNeibF(:,:),nNeibF(:)
    
    if(.not.this%isUp)then
      call this%polyFeGrid%up()
      ! indices of elements containing each node
      allocate(nEN(this%nN),source=0)
      do i=1,this%nE
        do j=1,this%nNE(i)
          nEN(this%iNE(j,i))=nEN(this%iNE(j,i))+1
        end do
      end do
      allocate(iEN(maxval(nEN),this%nN),source=0)
      nEN(:)=0
      do i=1,this%nE
        do j=1,this%nNE(i)
          nEN(this%iNE(j,i))=nEN(this%iNE(j,i))+1
          iEN(nEN(this%iNE(j,i)),this%iNE(j,i))=i
        end do
      end do
      ! neighbor cells and facets
      allocate(nNeibC(this%nC),source=0)
      allocate(nNeibF(this%nC),source=0)
      allocate(iNeibC(size(iEN,1)*maxval(this%nNE),this%nC),source=0)
      allocate(iNeibF(size(iEN,1)*maxval(this%nNE),this%nC),source=0)
      do i=1,this%nN
        do j=1,nEN(i)
          m=iEN(j,i)
          if(m<=this%nC)then
            do k=1,nEN(i)
              n=iEN(k,i)
              if(n<=this%nC.and.n/=m.and.all(n/=iNeibC(1:nNeibC(m),m)))then
                nNeibC(m)=nNeibC(m)+1
                iNeibC(nNeibC(m),m)=n
              elseif(n>this%nC.and.all(n/=iNeibF(1:nNeibF(m),m)))then
                nNeibF(m)=nNeibF(m)+1
                iNeibF(nNeibF(m),m)=n
              end if
            end do
          end if
        end do
      end do
      if(allocated(this%nNeibC)) deallocate(this%nNeibC)
      if(allocated(this%nNeibF)) deallocate(this%nNeibF)
      if(allocated(this%iNeibC)) deallocate(this%iNeibC)
      if(allocated(this%iNeibF)) deallocate(this%iNeibF)
      allocate(this%nNeibC,source=nNeibC)
      allocate(this%nNeibF,source=nNeibF)
      allocate(this%iNeibC,source=iNeibC(1:maxval(nNeibC),:))
      allocate(this%iNeibF,source=iNeibF(1:maxval(nNeibF),:))
      deallocate(iEN,nEN,iNeibC,nNeibC,iNeibF,nNeibF)
    end if
  end subroutine
  
  !> clear this PICGrid
  elemental subroutine clearPICGrid(this)
    class(PICGrid),intent(inout)::this !< this PICGrid
    
    call this%polyFeGrid%clear()
    if(allocated(this%nNeibC)) deallocate(this%nNeibC)
    if(allocated(this%nNeibF)) deallocate(this%nNeibF)
    if(allocated(this%iNeibC)) deallocate(this%iNeibC)
    if(allocated(this%iNeibF)) deallocate(this%iNeibF)
  end subroutine
  
  !> destructor of PICGrid
  elemental subroutine purgePICGrid(this)
    type(PICGrid),intent(inout)::this !< this PICGrid
    
    call this%clear()
  end subroutine
  
end module
