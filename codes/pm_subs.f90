module pm_subs
  use pmtypes
  implicit none

  ! Save Directory Name
  character(len=:), allocatable:: dir_name

contains
  subroutine get_pm_save_dir
    character(len=255):: arg

    !FIX -- Need to define output units and format, do not use stars
    !Command Line Arguments Needed ==> save directory
    if(command_argument_count() /= 1)then
       write(*,*)'ERROR, ONE COMMAND-LINE ARGUMENTS REQUIRED,', &
            'STOPPING; MISSING SAVE DIRECTORY NAME'
       stop
    endif

    call get_command_argument(1, arg)   !first, read in the value
    dir_name = trim(arg)
    write (*,*) 'saving code to----', dir_name,'----'
    return
  end subroutine get_pm_save_dir

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

    IF (idum<=0) THEN
       idum = max(-idum,1)
       idum2 = idum
       DO j = ntab + 8, 1, -1

          k = idum/iq1
          idum = ia1*(idum-k*iq1) - k*ir1
          IF (idum<0) idum = idum + im1
          IF (j<=ntab) iv(j) = idum
       END DO
       iy = iv(1)
    END IF
    k = idum/iq1
    idum = ia1*(idum-k*iq1) - k*ir1
    IF (idum<0) idum = idum + im1
    k = idum2/iq2
    idum2 = ia2*(idum2-k*iq2) - k*ir2
    IF (idum2<0) idum2 = idum2 + im2
    j = 1 + iy/ndiv
    iy = iv(j) - idum2
    iv(j) = idum
    IF (iy<1) iy = iy + imm1
    ran2 = min(am*iy,rnmx)
    RETURN
  END FUNCTION ran2


end module pm_subs
