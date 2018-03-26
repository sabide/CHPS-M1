program mpi_df
  
  use mpi
  implicit none
  integer :: comm_1D, ierr
  integer :: rank,rank_left,rank_right

  integer :: i,is,ie
  
  integer :: N
  real(kind=8) ,dimension(:),allocatable :: U 
  character(len=1024) :: filename
  character(len=3) :: num
  
  call mpi_init(ierr)
  N=4
  call create_comm()
  call get_index(is,ie)
  allocate(U(IS-1:IE+1))
  call mpi_comm_rank(comm_1d,rank,ierr)
  U(IS:IE)=rank
  call halo_update(U)




  num(1:3)='000'
  write (num, '(I3.3)')rank
  open(unit=10+rank,file="out"//num//".dat")
  do i=is-1,ie+1
     write(10+rank,*)U(I)
  end do

  close(10+rank)
  


  call mpi_finalize(ierr)


contains
  
  
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
    
    dims(1) = 4
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
