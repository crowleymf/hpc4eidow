module pm_subs
  use pmtypes
  implicit none

  ! Save Directory Name
  character(len=:),allocatable :: dir_name

contains
  !---------------------------------------
  !        get_pm_save_dir
  !---------------------------------------
  subroutine get_pm_save_dir
    use param, only:outu
    character(len=255):: arg

    !Command Line Arguments Needed ==> save directory
    if(command_argument_count() /= 1)then
       write(outu,'("pmsubs.f90:get_pm_save ",a,a)') &
            'ERROR, ONE COMMAND-LINE ARGUMENTS REQUIRED,', &
            'STOPPING; MISSING SAVE DIRECTORY NAME'
       stop
    endif

    call get_command_argument(1, arg)   !first, read in the value
   
    dir_name = trim(adjustl(arg))
    write (outu,'("PMSUBS> ",3a)') 'saving code to----', dir_name,'----'
    return
  end subroutine get_pm_save_dir

  !---------------------------------------
  !       RAN2
  !---------------------------------------
  !This function is the random number generator used for all simulation.
  !Function from numerical recipies      
  real(kind=pm_single) FUNCTION ran2(idum)
    use pmtypes
    implicit none
    integer :: idum
    !integer idum, im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, ndiv
    !real(kind=pm_single) am, eps, rnmx
    integer,parameter :: im1=2147483563,im2=2147483399,imm1=im1-1
    integer,parameter :: ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211
    integer,parameter :: ir2=3791,ntab=32,ndiv=1+imm1/ntab
    real(kind=pm_dbl),parameter :: am=1./im1     
    real(kind=pm_dbl),parameter :: eps=1.2e-7,rnmx=1.-eps
    integer idum2, j, k, iv(ntab), iy
    save iv, iy, idum2
    data idum2/123456789/, iv/ntab*0/, iy/0/

    if (idum<=0) then
       idum = max(-idum,1)
       idum2 = idum
       do j = ntab + 8, 1, -1

          k = idum/iq1
          idum = ia1*(idum-k*iq1) - k*ir1
          if (idum<0) idum = idum + im1
          if (j<=ntab) iv(j) = idum
       end do
       iy = iv(1)
    end if
    k = idum/iq1
    idum = ia1*(idum-k*iq1) - k*ir1
    if (idum<0) idum = idum + im1
    k = idum2/iq2
    idum2 = ia2*(idum2-k*iq2) - k*ir2
    if (idum2<0) idum2 = idum2 + im2
    j = 1 + iy/ndiv
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy<1) iy = iy + imm1
    ran2 = min(am*iy,rnmx)
    return
  end function ran2
  
  !---------------------------------------
  !       BIASD
  !---------------------------------------
  ! make the random direction moves in this monte carlo scheme.
  subroutine biasd(x,d,dir_name)

    use pmtypes
    use param
!    use pm_subs,only:ran2

    implicit none
    character(len=*) :: dir_name
    integer :: x,d
    !---real :: ran2
    double precision :: dtest
    double precision :: pxzx,ppyx,pmyx

    pxzx=pxz(x)
    ppyx=ppy(x)
    pmyx=pmy(x)

    iseed=10

    dtest=ran2(iseed)

    !move backward
    if (dtest<=pmyx) then
       d=int(ran2(iseed)*4)+5
       if (d==6) d=9
       if (d==8) d=12
       !move laterally
    else if ((dtest>pmyx).and.(dtest<=(pmyx+pxzx))) then
       d=int(ran2(iseed)*4)+1
       !move forward
    else if (dtest>(pmyx+pxzx)) then
       d=int(ran2(iseed)*4)+5
       if (d==5) d=10
       if (d==7) d=11
    end if

  end subroutine biasd

  subroutine pm_close_all_files
    use param, only:outu

    write(outu,'(a)')"PM> Closing all files"
    call close_chaindyn_files
    call close_chaincalcs_files
    call close_boxcalcs_files
    call close_autocorr_files
    return
  end subroutine pm_close_all_files

    subroutine close_chaincalcs_files
    !--- Formatted files --------
    close(unit = 40)
    close(unit = 41)
    close(unit = 42)
    close(unit = 43)
    close(unit = 44)
    close(unit = 45)
    close(unit = 46)
    close(unit = 47)

    !--- Unformatted files --------
    close(unit = 30)
    close(unit = 31)
    close(unit = 32)
    close(unit = 33)
    close(unit = 34)
    close(unit = 35)
    return
  end subroutine close_chaincalcs_files

subroutine close_boxcalcs_files
    close(unit=51)
    close(unit=52)
    close(unit=53)
    close(unit=54)

    close(unit=55)
    close(unit=56)
    close(unit=57)
    close(unit=58)
    return  
  end subroutine close_boxcalcs_files

  subroutine close_chaindyn_files
    close (unit=90)
    close (unit=91)
    close (unit=92)
    close (unit=93)
    close (unit=94)
    return
  end subroutine close_chaindyn_files

  subroutine close_autocorr_files
    close(unit=70)
    close(unit=71)
    return
  end subroutine close_autocorr_files


end module pm_subs
