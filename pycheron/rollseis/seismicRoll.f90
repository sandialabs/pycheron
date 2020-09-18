!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
! Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
! certain rights in this software.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTICE:
! For five (5) years from 10/21/2019 the United States Government is granted for
! itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
! license in this data to reproduce, prepare derivative works, and perform publicly and
! display publicly, by or on behalf of the Government. There is provision for the
! possible extension of the term of this license. Subsequent to that period or any
! extension granted, the United States Government is granted for itself and others
! acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
! data to reproduce, prepare derivative works, distribute copies to the public,
! perform publicly and display publicly, and to permit others to do so. The specific
! term of the license can be identified by inquiry made to National Technology and
! Engineering Solutions of Sandia, LLC or DOE. NEITHER THE UNITED STATES GOVERNMENT,
! NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING
! SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR
! IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
! USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
! THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. Any licensee of this software
! has the obligation and responsibility to abide by the applicable export control laws,
! regulations, and general prohibitions relating to the export of technical data.
! Failure to obtain an export control license or other authority from the Government
! may result in criminal liability under U.S. laws.
! (End of Notice)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     !----------------------------------------------------

subroutine stalta_filter(x, n, n_sta, n_lta, ind, out)
    implicit real (a-z)
    integer n_sta, n_lta, ind, n, i
    real out, x
    dimension x(n)
    !f2py   intent(out) out
    !       !sta
    total = 0.0
    do i = 1, n_sta, 1
        total = total + x(ind + i)
    end do
    sta = total / n_sta

    !       !lta
    total = 0.0
    ind = ind + 1
    do i = 0, n_lta - 1, 1
        if ((ind - i) .lt. 0) then
            total = total + x((ind - i) + n)
        else if((ind - i) .eq. 0) then
            total = total + x(n)
        else
            total = total + x(ind - i)
        end if
    end do
    lta = total / n_lta

    out = sta / lta
end subroutine stalta_filter
!     !----------------------------------------------------
subroutine roll_stalta(x, n, n_sta, n_lta, increment, arr)
    implicit real(a-z)
    integer n_sta, n_lta, ind, n, increment
    real arr, x
    !f2py   intent(out) arr
    dimension x(n), arr(n)

    if (increment .lt. 1) then
        increment = 1
    end if
    ind = n_lta

    do while (ind < (n - n_sta))
        call stalta_filter(x, n, n_sta, n_lta, ind, filter)
        arr(ind) = filter
    end do
end subroutine roll_stalta
!     !----------------------------------------------------