program test
    integer*8 :: tr(800000), out_e(800000), out_v(800000), out_s(800000)
    character out_t(800000)
    integer n, minRep, bound, trunc
    real fs
    integer*8 :: x(64)
!    data x/2,3,4,60,60,60,60,60,60,60,60,60,60,60,60,60,3,2,1,5,4,3,2,5,5,1,2,4,30,30,30,30,30,6,3,4,6,7,8,0,1,5,6, &
!    3,5,3,1,3,1,3,4,2, 20, 25, 22, 23, 27, 30, 31, 29, 31, 30, 30, 31/

    ! see repAmpsTest.py for how this was made
    open(1, file = "/Users/jbobeck/pycheron/pycheron/metrics/repAmps/jitter.dat")
    read(1, *) tr
    close(1)
    n= 800000
    minRep = 100
    bound = 45
    call repAmps(tr, n, minRep, bound, trunc, out_s, out_e, out_v, out_t)
!    print*, out_t(1:trunc+1)

end program test

subroutine repAmps(x, n, minRep, bound, trunc, out_s, out_e, out_v, out_t)

    implicit none
    integer*8, allocatable :: tmp1(:)
    character out_t(n)
    integer*8 x(n), out_s(n), out_v(n), out_e(n), tmp(minRep)
    integer :: n, minRep, i, j, bound, count, ind
    integer :: trunc, inc

    !f2py intent(out) trunc
    !f2py intent(out) out_s
    !f2py intent(out) out_e
    !f2py intent(out) out_v
    !f2py intent(out) out_t
    i = 1
    ind = 1
    do while (i< (size(x) - minRep))
        count = 1
        inc = 1
        !case 1 + 2
        if (x(i) == x(i+minRep)) then
            tmp = x(i+minRep)
            if (all([(tmp(j) == x(i),j = 1, minRep)])) then
                allocate(tmp1(1))
                tmp1 = x(i:i)
                do while (all([(tmp1(j) == x(i), j=1,count)]))
!                    print *, "count: ", count
!                    print *, x(i:i+(count-1))
                    count = count + 1
                    deallocate(tmp1)
                    allocate(tmp1(count))
                    ! need to subtract one from count because count is inc above before it fails
                    tmp1 = x(i:i+(count -1))
                end do
                ! need to subtract one from count because count is inc above before it fails
                if (count-1 >= minRep) then
                    ! start is -1 because of python indexes
                    out_s(ind) = i - 1
                    out_t(ind) = "f"
                    out_v(ind) = x(i)
                    print *, "--flat--"
                    print *, "start: ", i - 1
!                    print *, x(i:i+count-2)
                    i=i+(count-2)
                    print *, "end: ", i
                    out_e(ind) = i
                    ind = ind + 1
                    deallocate(tmp1)
                    ! repAmp found, skip jitter and inc i
                    goto 20
                else
                    deallocate(tmp1)
                    goto 10
                end if
             else
                goto 10
             end if
        else
            goto 10
        end if

        ! case 3+4
 10       if( (x(i) == x(i+minRep)) .OR. (abs(x(i) - x(i+minRep)) <= bound) ) then
            tmp = x(i+minRep)
            if (all([(abs(tmp(j) - x(i)) <=bound , j =1, minRep)])) then
                allocate(tmp1(1))
                tmp1 = x(i:i)
                do while ((all([(abs(tmp1(j) - x(i)) <= bound , j =1, count)])) .and. (i + count < size(x)))
                    count = count + 1
                    deallocate(tmp1)
                    allocate(tmp1(count))
                    tmp1 = x(i:i+count-1)
                end do

                if (count-1 >= minRep) then
                    ! start is -1 because of python indexes
                    out_s(ind) = i -1
                    out_t(ind) = "j"
                    out_v(ind) = -99999
                    print *, "--jitter--"
                    print *, "start: ", i - 1
                    i=i+count-2
                    print *, "end: ", i
                    out_e(ind) = i
                    ind = ind + 1
                end if
                deallocate(tmp1)
            end if
        end if
 20      i= i +1
    end do

    trunc = ind -1

end subroutine repAmps