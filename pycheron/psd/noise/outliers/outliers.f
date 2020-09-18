!program test
!    implicit none
!    integer n
!    integer*8 tr(432040)
!    real*8 out(432040)
!
!
!    ! spike_example.mseed
!    open(1, file = "data.dat")
!    read(1, *) tr
!    close(1)
!
!    n = 432040
!
!    call hampelFilter(tr,n,42,1,out)
!!    print *, out
!
!end program test

      subroutine hampelFilter(data,n, nwin, increment, out)
      implicit none
      integer*8 data(n)
      real*8 out(n), hampel
      integer nwin, increment, k, i,n
Cf2py intent(out) out(n)


       k = nwin /2

C initializing to NaN
      do i =1, n
        out(i) = -99999
      end do

      do i = k, n - k, increment
        out(i) = hampel(data, nwin, i)
      end do

      end subroutine hampelFilter


      function hampel(data,nwin,ind) result(out)
      implicit none
      integer*8 data(*)
      integer :: nwin, ind, k, x0, n, medianAbsMinusMedian
      real*8 :: L, S0, out
      integer*8, allocatable :: tmp(:), absMinusMedian(:)
Cf2py intent(out) out

C Constants
      L = 1.4826

C Get half window width
      k = nwin/2

      n = size(data(ind-k: (ind-k+nwin) -2))
C Calcualting x0 = median(data[nwin])
      allocate(tmp(n))
      allocate(absMinusMedian(n))
      tmp = data(ind-k: (ind-k+nwin)-2)

      call quicksort(tmp, 1, n)


      if (mod(nwin,2) == 0) then
          x0 = (tmp((nwin/2) -1) + tmp(nwin/2))/2

          absMinusMedian = abs(tmp - x0)

          call quicksort(absMinusMedian, 1,n)

          medianAbsMinusMedian = (absMinusMedian((nwin/2)-1)
     &     + absMinusMedian(nwin/2))/2
      else
          x0 = tmp(nwin/2)
          absMinusMedian = abs(tmp - x0)

          call quicksort(absMinusMedian, 1,n)

          medianAbsMinusMedian = absMinusMedian(nwin/2)

      end if

      S0 = L * medianAbsMinusMedian

      out = abs(data(ind) - x0)/ S0

      deallocate(tmp)
      deallocate(absMinusMedian)

      end function hampel

      recursive subroutine quicksort(a, first, last)
      implicit none
      integer*8  a(*), x, t
      integer first, last
      integer i, j
C f2py intent(out) a

      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         i=i+1
         j=j-1
      end do
      if (first < i-1) call quicksort(a, first, i-1)
      if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort
