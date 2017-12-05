 
program hw2
implicit none
include "mpif.h"
integer ninp/5/, nres/6/
integer(4) :: mpiErr, mpiRank

call mpi_init(mpiErr)

call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)


!процесс 0 хочет поделиться данными с процессом 1
if(mpiRank == 0) then
    call main()
else
    call helper()
endif
    
call mpi_finalize(mpiErr)


contains

subroutine main()
implicit none
include "mpif.h"
integer ninp/5/, nres/6/
real(8), dimension(:,:), allocatable :: a
real(8), dimension(5) :: res
integer(4) :: nrow, ncol, i, nextRowToSolve, rowSolved, activeHelpers, maxLeft, maxRight, maxTop, maxBottom
integer(4), dimension(2) :: nsize
integer(4) :: mpiErr, mpiSize, mpiRank
integer(4), dimension(MPI_STATUS_SIZE) :: status
real(8) :: currSum, maxSum

! read matrix
open(unit=ninp, file='main.input')
read(ninp,*) nrow, ncol
nsize(1) = nrow
nsize(2) = ncol
allocate(a(nrow,ncol))
read(ninp,*) a
close(ninp)

! send matrix to helpers
call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
activeHelpers = mpiSize - 1

!write(*,*) "master sends matrix"
do i=1,activeHelpers
    call mpi_send(nsize, 2, MPI_INTEGER4, i, 1, MPI_COMM_WORLD, mpiErr)
    call mpi_send(a, nrow * ncol, MPI_REAL8, i, 2, MPI_COMM_WORLD, mpiErr)
end do


! main loop

maxSum = a(1,1)
maxLeft = 1
maxRight = 1
maxTop = 1
maxBottom = 1

nextRowToSolve = 1
rowSolved = 0
do while (rowSolved < nrow .or. activeHelpers > 0)
    ! get response from helper
    !write(*,*) "master gets response"
    call mpi_recv(res, 5, MPI_REAL8, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpiErr)
    
    ! check if response is empty and calc max
    if (status(MPI_TAG) == 4) then
        !write(*,*) "master gets 4 from ", status(MPI_SOURCE)
        rowSolved = rowSolved + 1
        
        currSum = res(1)
        
        if (currSum > maxSum) then
            maxSum = currSum
            maxLeft = int(res(2))
            maxTop = int(res(3))
            maxRight = int(res(4))
            maxBottom = int(res(5))
        end if
            
    end if
    
    ! send next task
    if (nextRowToSolve <= nrow) then
        !write(*,*) "master sends 5 to ", status(MPI_SOURCE)
        call mpi_send(nextRowToSolve, 1, MPI_INTEGER4, status(MPI_SOURCE), 5, MPI_COMM_WORLD, mpiErr)
        
        nextRowToSolve = nextRowToSolve + 1
    else
        !write(*,*) "master sends 6 to ", status(MPI_SOURCE)
        call mpi_send(nextRowToSolve, 1, MPI_INTEGER4, status(MPI_SOURCE), 6, MPI_COMM_WORLD, mpiErr)
        
        activeHelpers = activeHelpers - 1
    end if

end do

open(unit=nres, file='main.output')
write(nres,*) maxTop, maxLeft, maxBottom, maxRight, maxSum
close(nres)

deallocate(a)

end subroutine



subroutine kadane1d(a, x1, x2, maxSum)
implicit none
real(8), intent(in), dimension(:) :: a
integer(4), intent(out) :: x1, x2
real(8), intent(out) :: maxSum
integer(4) :: i, leftIndex, n
real(8) :: MaxEndingHere, candidateSum1, candidateSum2

n = size(a)

! initialization
maxSum = a(1); x1 = 1; x2 = 1

MaxEndingHere=a(1); leftIndex=1
do i=2,n

  candidateSum1 = a(i)
  candidateSum2 = MaxEndingHere + a(i)

  if (candidateSum1 .gt. candidateSum2) then
    MaxEndingHere = candidateSum1
    leftIndex = i
  else
    MaxEndingHere = candidateSum2
  endif

  if (MaxEndingHere .gt. maxSum) then
    maxSum = MaxEndingHere
    x1 = leftIndex
    x2 = i
  endif

enddo


end subroutine





subroutine helper()
implicit none
include "mpif.h"
integer(4) :: mpiErr, mpiSize, mpiRank
integer(4), dimension(MPI_STATUS_SIZE) :: status
integer(4), dimension(2) :: nsize
real(8), dimension(:,:), allocatable :: a
real(8), dimension(:), allocatable :: p
real(8), dimension(5) :: res
integer(4) :: ncol, nrow, maxLeft, maxBottom, maxRight, top, bottom, left, right, j
real(8) :: currSum, maxSum


! get matrix
call mpi_recv(nsize, 2, MPI_INTEGER4, 0, 1, MPI_COMM_WORLD, status, mpiErr)
nrow = nsize(1)
ncol = nsize(2)
allocate(a(nrow, ncol))
allocate(p(ncol))
call mpi_recv(a, nrow * ncol, MPI_REAL8, 0, 2, MPI_COMM_WORLD, status, mpiErr)

! send empty response
!write(*,*) "helper sends 3 to master"
call mpi_send(currSum, 1, MPI_REAL8, 0, 3, MPI_COMM_WORLD, mpiErr)


do

    ! get next task
    !write(*,*) "helper gets task"
    call mpi_recv(top, 1, MPI_INTEGER4, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpiErr)
    
    ! check if there is something to do
    if (status(MPI_TAG) == 6) then
        !write(*,*) "no work"
        exit
    end if
    
    ! calculate task
    p = 0
    do bottom = top,nrow

        ! add another row to p
        do j = 1,ncol
            p(j) = p(j) + a(bottom, j)
        enddo

        ! call kadane for p
        call kadane1d(p, left, right, currSum)

        ! check if new candidate is better than current max
        if (currSum > maxSum .or. top == bottom) then
            maxSum = currSum;
            maxLeft = left
            maxRight = right
            maxBottom = bottom
        endif

    enddo

    
    ! send response
    res(1) = maxSum
    res(2) = maxLeft
    res(3) = top
    res(4) = maxRight
    res(5) = maxBottom
    !write(*,*) "helper sends result: ", top, maxLeft, maxBottom, maxRight, maxSum
    call mpi_send(res, 5, MPI_REAL8, 0, 4, MPI_COMM_WORLD, mpiErr)

end do

deallocate(p)
deallocate(a)

end subroutine

end program


