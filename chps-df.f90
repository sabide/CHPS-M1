program mpi_df
  
  use mpi
  implicit none
  integer :: comm_1D, ierr
  integer :: rank,rank_left,rank_right

  integer :: i,is,ie
  
  integer :: N,NN
  
  character(len=1024) :: filename
  character(len=3) :: num
  integer,parameter :: nprocs=2

  
  real(kind=8) ,dimension(:),allocatable ::  X
  real(kind=8) ,dimension(:),allocatable ::  U 
  real(kind=8) ,dimension(:),allocatable :: DU
  real(kind=8) ,dimension(:),allocatable ::  E
  real(kind=8) :: H
  
  call mpi_init(ierr)
  N=200
  call create_comm()
  call get_index(is,ie)

  !> allocation
  ALLOCATE( X(IS-1:IE+1))
  ALLOCATE( U(IS-1:IE+1))
  ALLOCATE(DU(IS-1:IE+1))
  ALLOCATE( E(IS-1:IE+1))
  NN = NPROCS*N
  H  = 1./REAL(NN)
  do i=is,ie
     x(i)=(i-0.5)*h
  end do
  call halo_update(x)

  do i=is,ie
     u(i)=f(x(i))
  end do
  call halo_update(u)

  do i=is,ie
     du(i)=(u(i+1)-u(i-1))/(2*h)
  end do
  call halo_update(du)
  
  do i=is,ie
     E(I)=abs(du(i)-df(x(i)))
  end do
  call halo_update(E)
  
  
  num(1:3)='000'
  write (num, '(I3.3)')rank
  open(unit=10+rank,file="out"//num//".dat")
  do i=is-1,ie+1
     write(10+rank,'(4(e15.8,1x))')X(I),U(I),DU(I),E(I)
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
    periods(1) = .true.
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
