!===============================================================================
module stack
!===============================================================================
!
! Routines for stacking seismograms stored in SAC format.
!
! Andy Nowacki
! a.nowacki@leeds.ac.uk
!

   use, intrinsic :: iso_c_binding, only: C_FLOAT, C_DOUBLE
   use f90sac
   ! Use the OpenMP module if needed
!$ use omp_lib

   implicit none

   private

   ! Precision selectors--based on C types because we deal with SAC data
   integer, parameter :: r4 = C_FLOAT, r8 = C_DOUBLE
   integer, parameter :: rs = r4

   ! I/O
   integer, parameter :: lu_stdin = 5, lu_stdout = 6, lu_stderr = 0, &
                         lu_in = 10, lu_out = 11

   ! Useful constants
   real(rs), parameter :: pi = real(4.d0*atan2(1.d0, 1.d0), rs)
   real(rs), parameter :: deg2km = 111.19_rs

   interface stack_allocate
      module procedure :: stack_allocate_r4_1d
   end interface stack_allocate

   public :: &
      stack_sum

contains

!===============================================================================
subroutine stack_sum(s, out, t1, t2, pick)
!===============================================================================
! Take an array of SAC traces and return the summed trace between times t1 and t2
! (if given), or just the whole trace.
!  INPUT:
!     s(:)   : SACtrace array of input traces
!  INPUT (OPTIONAL):
!     t1, t2 : Start and stop times for output stack.  Relative to O time, unless
!              pick array is supplied
!     pick(:): Array of pick times relative to which stack is performed.  Must be
!              an array of size(s)
!  OUTPUT:
!     out(:) : SACtrace containing the stack
!
   type(SACtrace), intent(in) :: s(:)
   type(SACtrace), intent(inout) :: out
   real(rs), intent(in), optional :: t1, t2, pick(:)
   real(rs) :: w1, w2
   integer :: npts, i, j, iw1, iw2

   call stack_check(s)

   ! Check input
   w1 = maxval(s%b)
   w2 = minval(s%e)
   if (present(t1)) w1 = t1
   if (present(t2)) w2 = t2
   if (w2 <= w1) call stack_error('stack_sum: Start time must be before end time')
   if (present(pick)) then
      if (size(pick) /= size(s)) &
         call stack_error('stack_sum: Pick array must be same length as number of traces')
      if (any(pick > s%e .or. pick < s%b)) &
         call stack_error('stack_sum: Picks must not be outside trace')
      if (any(pick + t2 > s%e .or. pick - t1 < s%b)) &
         call stack_error('stack_sum: Window around picks is outside data range')
   endif

   ! Make space for stack
   npts = int((w2 - w1)/s(1)%delta) + 1
   call f90sac_newtrace(npts, s(1)%delta, out)

   out%b = w1
   out%e = w2

!$omp parallel do default(none) shared(s, npts, out, pick, w1, w2) private(i, j, iw1, iw2)
   do i = 1, size(s)
      ! If using pick times, then the times are relative to this
      if (present(pick)) then
         iw1 = nint((pick(i) + w1 - s(i)%b)/s(1)%delta) + 1
         iw2 = iw1 + npts - 1
      ! Otherwise use time relative to O marker
      else
         iw1 = nint((w1 - s(i)%b)/s(1)%delta) + 1
         iw2 = nint((w2 - s(i)%b)/s(1)%delta) + 1
      endif
      do j = iw1, iw2
         out%trace(j-iw1+1) = out%trace(j-iw1+1) + s(i)%trace(j)
      enddo
   enddo
!$omp end parallel do
end subroutine stack_sum
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_check(s)
!===============================================================================
! Check that an array of SAC traces conforms to the criteria necessary for use
! in routines in this module.
!
   type(SACtrace), intent(in) :: s(:)
   integer :: n

   n = size(s)
   if (n == 1) then
      call stack_warning('stack_check: Only one trace in stack')
      return
   endif
   if (any(s(2:n)%delta /= s(1)%delta)) &
      call stack_error('stack_check: All traces must have same delta')
end subroutine stack_check
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_error(string)
!===============================================================================
   character(len=*), intent(in) :: string
   write(lu_stderr,'(a)') 'Error: '// string
   error stop
end subroutine stack_error
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_warning(string)
!===============================================================================
   character(len=*), intent(in) :: string
   write(lu_stderr,'(a)') 'Warning: '// string
end subroutine stack_warning
!-------------------------------------------------------------------------------

!===============================================================================
subroutine stack_allocate_r4_1d(a, n)
!===============================================================================
   real(r4), intent(inout), allocatable, dimension(:) :: a
   integer, intent(in) :: n
   if (allocated(a)) then
      if (size(a) == n) return
      deallocate(a)
      allocate(a(n))
      return
   else
      allocate(a(n))
   endif
end subroutine stack_allocate_r4_1d
!-------------------------------------------------------------------------------

end module stack
!-------------------------------------------------------------------------------
