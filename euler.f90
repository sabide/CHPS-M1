program mpi_df
  
  use mpi
  implicit none
  integer :: comm_1D, ierr
  integer :: rank,rank_left,rank_right

  integer :: i,is,ie
  
  integer :: N,NN
  
  character(len=1024) :: filename
  character(len=3) :: num
  integer,parameter :: nprocs=1

  
  real(kind=8) ,dimension(:),allocatable ::  X
  real(kind=8) ,dimension(:),allocatable ::  Uold,unew
  real(kind=8) :: DX,DT
  integer :: iteration

  !> init mpi
  call mpi_init(ierr)
  !> nb node per cpu
  N=10
  call create_comm()
  !> init communication hal
  call get_index(is,ie)

  
  PRINT*,"check istar iend",IS,IE

  !> allocation
  ALLOCATE( X   (IS-1:IE+1))
  ALLOCATE( UOLD(IS-1:IE+1))
  ALLOCATE( UNEW(IS-1:IE+1))
  
  NN = NPROCS*N
  DX  = 1./REAL(NN)
  do i=is,ie
     x(i)=(i-0.5)*dx
  end do
  call halo_update(x)
  !> 
  do i=is,ie
     uold(i)=0 
  end do
  if (is==1 ) uOLD(is-1)=1
  if (Ie==nn) uOLD(ie+1)=1
  call halo_update(uold)


  DT = 1.5*DX**2
  do iteration=1,4
     
     do i=is,ie
        unew(i) = uold(i) + dt/dx**2*(uold(i+1)-2*uold(i)+uold(i-1))
     end do
     call halo_update(unew)
     if (is==1 ) uNEW(is-1)=1
     if (Ie==nn) uNEW(ie+1)=1
     uold=unew
  end do
  
  num(1:3)='000'
  write (num, '(I3.3)')rank
  open(unit=10+rank,file="out"//num//".dat")
  do i=is,ie
     write(10+rank,'(4(e15.8,1x))')X(I),UOLD(I)
  end do

  close(10+rank)
  


  call mpi_finalize(ierr)


contains
  real(kind=8) function f(x)
    implicit none
    real(kind=8) :: x
    real(kind=8) :: pi=acos(-1d0)
    f=sin(2*pi*x)
  end function f
  real(kind=8) function df(x)
    implicit none
    real(kind=8) :: x
    real(kind=8) :: pi=acos(-1d0)
    df=2*pi*cos(2*pi*x)
  end function df
  subroutine get_index(is,ie)
    implicit none
    integer :: is,ie
    integer :: rank,ierr
    call mpi_comm_rank(comm_1d,rank,ierr)
    is = rank*N+1
    ie = (rank+1)*(N)

  end subroutine get_index

  subroutine create_comm()
    implicit none
    integer, parameter ::ndims = 1
    integer, dimension(ndims) ::dims
    logical, dimension(ndims) ::periods
    logical :: reorganisation

    
    dims(1) = nprocs
    periods(1) = .false.
    reorganisation = .false.
    call MPI_CART_CREATE( MPI_COMM_WORLD ,ndims,dims,periods,reorganisation,comm_1D,ierr)
    call mpi_comm_rank(comm_1d,rank,ierr)
    
    CALL MPI_CART_SHIFT(comm=comm_1d, direction=0,disp=1, rank_source=rank_left, rank_dest=rank_right, ierror=ierr)

    
  end subroutine create_comm

  
  subroutine halo_update(tab)
    implicit none
    real(kind=8),dimension(:),allocatable :: tab
    logical :: reorganisation
    integer, dimension( MPI_STATUS_SIZE) :: statut
    integer :: tag=100
    
    CALL MPI_SEND (TAB(IS)  ,1, MPI_DOUBLE_PRECISION ,RANK_LEFT ,TAG, COMM_1D, IERR )
    CALL MPI_RECV (TAB(IE+1),1, MPI_DOUBLE_PRECISION ,RANK_RIGHT,TAG, COMM_1D,STATUT,IERR)
    
    CALL MPI_SEND (TAB(IE)  ,1, MPI_DOUBLE_PRECISION ,RANK_RIGHT,TAG, COMM_1D, IERR )
    CALL MPI_RECV (TAB(IS-1),1, MPI_DOUBLE_PRECISION ,RANK_LEFT ,TAG, COMM_1D,STATUT,IERR)
    
    
  end subroutine halo_update


end program mpi_df
