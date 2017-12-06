module homework

contains

 
subroutine FindMaxCoordinates(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
real(8), dimension(:,:), intent(out) :: a
integer(4), intent(out) :: x1, y1, x2, y2
integer(4) :: mpiErr, mpiRank

call mpi_init(mpiErr)

call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)


!процесс 0 хочет поделиться данными с процессом 1
if(mpiRank == 0) then
    call master(a, x1, y1, x2, y2)
else
    call helper(a, x1, y1, x2, y2)
endif
    
call mpi_finalize(mpiErr)

end subroutine

subroutine master(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
real(8), dimension(:,:) :: a
integer(4), intent(out) :: x1, y1, x2, y2
real(8), dimension(5) :: res
integer(4) :: nrow, ncol, i, nextRowToSolve, rowSolved, activeHelpers
integer(4) :: mpiErr, mpiSize, mpiRank
integer(4), dimension(MPI_STATUS_SIZE) :: status
real(8) :: currSum, maxSum
integer(4), dimension(4) :: answer

nrow = size(a, dim = 1)
ncol = size(a, dim = 2)

! send matrix to helpers
call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
activeHelpers = mpiSize - 1

! main loop


maxSum = a(1,1)
x1 = 1
y1 = 1
x2 = 1
y2 = 1

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
            x1 = int(res(2))
            y1 = int(res(3))
            x2 = int(res(4))
            y2 = int(res(5))
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

! send result back to helpers
answer(1) = x1
answer(2) = y1
answer(3) = x2
answer(4) = y2
do i=1,mpiSize-1
    call mpi_send(answer, 4, MPI_INTEGER4, i, 1, MPI_COMM_WORLD, mpiErr)
end do

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





subroutine helper(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
integer(4) :: mpiErr, mpiSize, mpiRank
integer(4), dimension(MPI_STATUS_SIZE) :: status
real(8), dimension(:,:) :: a
real(8), dimension(:), allocatable :: p
real(8), dimension(5) :: res
integer(4), intent(out) :: x1, y1, x2, y2
integer(4) :: ncol, nrow, top, bottom, left, right, j
integer(4), dimension(4) :: answer
real(8) :: currSum, maxSum

nrow = size(a, dim = 1)
ncol = size(a, dim = 2)
allocate(p(ncol))
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
    x1 = top
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
            y1 = left
            y2 = right
            x2 = bottom
        endif

    enddo

    
    ! send response
    res(1) = maxSum
    res(2) = x1
    res(3) = y1
    res(4) = x2
    res(5) = y2
    !write(*,*) "helper sends result: ", top, maxLeft, maxBottom, maxRight, maxSum
    call mpi_send(res, 5, MPI_REAL8, 0, 4, MPI_COMM_WORLD, mpiErr)

end do

deallocate(p)

call mpi_recv(answer, 4, MPI_INTEGER4, 0, 1, MPI_COMM_WORLD, status, mpiErr)
x1 = answer(1)
y1 = answer(2)
x2 = answer(3)
y2 = answer(4)

end subroutine

end module



