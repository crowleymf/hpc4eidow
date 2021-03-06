module chaindyn
  use pmtypes

  implicit none

contains

  SUBROUTINE chaindynamics(dir_name)
    !The chaindynamics subroutine presented here is used to calculate the mean squared
    !displacement variables called g1, g2, g3, g4, g5. These variables track the
    !motion of of various monomers. From these correlation functions we can
    !obtain the dynamic "finger print" of the system

    USE param
    use pmtypes
    USE correlation

    implicit none

    character(len=*) :: dir_name
    character(len=20) :: fmt200="(3(a20,x))"
    character(len=20) :: fmt201="(i20,x,2(f20.10,x))"

    integer x, y, z
    integer cm, xcm, ycm, zcm

    real(kind=pm_dbl) dnk
    real(kind=pm_dbl) dxcom, dycom, dzcom
    real(kind=pm_dbl) xcom, ycom, zcom
    real(kind=pm_dbl) xnwf, xnwl, ynwf
    real(kind=pm_dbl) ynwl, znwf, znwl
    real(kind=pm_dbl) g1, g1sum, g1avg
    real(kind=pm_dbl) g2, g2sum, g2avg
    real(kind=pm_dbl) g3, g3sum, g3avg
    real(kind=pm_dbl) g4f, g4fsum, g4favg
    real(kind=pm_dbl) g4l, g4lsum, g4lavg, g4avg
    real(kind=pm_dbl) g5f, g5fsum, g5favg
    real(kind=pm_dbl) g5l, g5lsum, g5lavg, g5avg

    !Initialization
    if (t==0) then
       call open_chaindyn_files(dir_name)
       !$OMP parallel sections
       !$OMP section
       write (90,fmt200) 'l','t','g1'
       !$OMP section
       write (91,fmt200) 'l','t','g2'
       !$OMP section
       write (92,fmt200) 'l','t','g3'
       !$OMP section
       write (93,fmt200) 'l','t','g4'
       !$OMP section
       write (94,fmt200) 'l','t','g5'
       !$OMP section
       !allocation of initial posistion variables
       allocate (x0cm(nkt), y0cm(nkt), z0cm(nkt))
       !$OMP section
       allocate (x0nwf(nkt), y0nwf(nkt), z0nwf(nkt))
       !$OMP section
       allocate (x0nwl(nkt), y0nwl(nkt), z0nwl(nkt))
       !$OMP section
       allocate (x0com(nkt), y0com(nkt), z0com(nkt))
       !$OMP end parallel sections
    end if

    !Initialize variables
    call initalize

    do k=1,nkt
       call zero_var
       call chain_step_down

       do while (atemp/=code)
          call calc_central_monomer
       end do

       call increment_position
       call chain_component_center_of_mass_calc
       call build_and_calc_ensemble
       !this ends the loop over the chain.
    end do

    call calc_average_ensemble
    !$OMP parallel sections
    !$OMP section
    write (90,fmt201) l, t, g1avg
    !$OMP section
    write (91,fmt201) l, t, g2avg
    !$OMP section
    write (92,fmt201) l, t, g3avg
    !$OMP section
    write (93,fmt201) l, t, g4avg
    !$OMP section
    write (94,fmt201) l, t, g5avg
    !$OMP end parallel sections
    return


  contains
    !----------------------------------
    !      initalize
    !----------------------------------
    subroutine initalize
      g1sum=0
      g2sum=0
      g3sum=0
      g4fsum=0
      g4lsum=0
      g5fsum=0
      g5lsum=0

      g1=0
      g2=0
      g3=0
      g4f=0
      g4l=0
      g5f=0
      g5l=0
    end subroutine initalize

    !----------------------------------
    !      ZERO_VAR
    !----------------------------------
    subroutine zero_var
      cm = nw(k)/2 - 1
      i = 0

      dx = 0
      dy = 0
      dz = 0

      dxmv = 0
      dymv = 0
      dzmv = 0

      atemp = 0
    end subroutine zero_var

    !----------------------------------
    !     chain_step_down
    !----------------------------------
    subroutine chain_step_down
      !Find kth chain end and label as (x,y,z)
      !Setup calculation of mean-square displacement of chain end
      !$OMP parallel sections
      !$OMP section
      x=xx(k)
      xnwf = real(xd(k), kind=pm_dbl)
      !$OMP section
      y=yy(k)
      ynwf = real(yd(k), kind=pm_dbl)
      !$OMP section
      z=zz(k)
      znwf = real(zd(k), kind=pm_dbl)
      !$OMP end parallel sections


      !Set lab coordinates of chain end at time t=0
      IF (t==0) THEN
         x0nwf(k) = xnwf
         y0nwf(k) = ynwf
         z0nwf(k) = znwf
      END IF

      !Step down the chain in the a-direction
      atemp = a(x,y,z)
    end subroutine chain_step_down

    !----------------------------------
    !      calc_central_monomer
    !----------------------------------
    subroutine calc_central_monomer
      !Increment the delta counters dxmv, dymv, dzmv for each
      !segmental step taken. This gives the segmental position
      !in "chain" coordinates (origin at the b=13 chain end)

      !Compile distances from first monomer
      !(Sum segmental positions in chain coordinates)

      !Obtain the position of the next segment in box coordinates
      !based on present position and the value of the a-code

      dx = dx + px(atemp)
      dxmv = dxmv + dx
      xnew = xnp(atemp,x)
      x = xnew

      dy = dy + py(atemp)
      dymv = dymv + dy
      ynew = ynp(atemp,y)
      y = ynew

      dz = dz + pz(atemp)
      dzmv = dzmv + dz
      znew = znp(atemp,z)
      z = znew

      !Find central monomer position in lab coordinates
      IF (i==cm) THEN
         xcm = xd(k) + dx
         ycm = yd(k) + dy
         zcm = zd(k) + dz

         !Set lab coordinates of central monomer at t=0
         IF (t==0) THEN
            x0cm(k) = xcm
            y0cm(k) = ycm
            z0cm(k) = zcm
         END IF
      END IF

      i = i + 1
      atemp = a(x,y,z)
    end subroutine calc_central_monomer

    !----------------------------------
    !      increment_position
    !----------------------------------
    subroutine increment_position
      xnwl = xd(k) + dx
      ynwl = yd(k) + dy
      znwl = zd(k) + dz

      !Set lab coordinates of chain end at t = 0
      IF (t==0) THEN
         x0nwl(k) = xnwl
         y0nwl(k) = ynwl
         z0nwl(k) = znwl
      END IF
    end subroutine increment_position

    !----------------------------------
    !      chain_component_center_of_mass_calc
    !----------------------------------
    subroutine chain_component_center_of_mass_calc
      !Calculate components of center of mass in the chain coordinates
      !Now add first monomer position to get the position of center of mass
      !in the 'lab' reference frame

      dxcom = real(dxmv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)
      xcom = dxcom + xd(k)

      dycom = real(dymv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)
      ycom = dycom + yd(k)

      dzcom = real(dzmv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)
      zcom = dzcom + zd(k)

      !Set initial center of mass position at t=0
      IF (t==0) THEN
         x0com(k) = xcom
         y0com(k) = ycom
         z0com(k) = zcom
      END IF
    end subroutine chain_component_center_of_mass_calc

    !----------------------------------
    !      build_and_calc_ensemble
    !----------------------------------
    subroutine build_and_calc_ensemble
      !Calculate mean-squared displacement of central monomer at time t
      !$OMP parallel sections
      !$OMP section
      g1 = (xcm-x0cm(k))**2 + (ycm-y0cm(k))**2 + (zcm-z0cm(k))**2
      g1sum=g1sum+g1
      !$OMP section
      !Calculate mean-squared displacement of central monomers in
      !center of mass coordinates at time t
      g2=((xcm-xcom)**2+(ycm-ycom)**2+(zcm-zcom)**2)-2*((xcm-xcom)*(x0cm(k)-x0com(k))+ &
           (ycm-ycom)*(y0cm(k)-y0com(k))+(zcm-zcom)*(z0cm(k)-z0com(k)))+ &
           ((x0cm(k)-x0com(k))**2+(y0cm(k)-y0com(k))**2+(z0cm(k)-z0com(k))**2)
      g2sum=g2sum+g2

      !$OMP section
      !Calculate mean-squared center of mass displacement at time t
      g3=(xcom-x0com(k))**2+(ycom-y0com(k))**2+(zcom-z0com(k))**2
      g3sum=g3sum+g3

      !$OMP section
      !Calculate mean-squared displacement of chain ends at time t
      g4f=(xnwf-x0nwf(k))**2+(ynwf-y0nwf(k))**2+(znwf-z0nwf(k))**2
      g4fsum=g4fsum+g4f

      g4l=(xnwl-x0nwl(k))**2+(ynwl-y0nwl(k))**2+(znwl-z0nwl(k))**2
      g4lsum=g4lsum+g4l

      !$OMP section
      !Calculate mean-squared displacement of chain ends in center of mass coordinates at time t
      g5f=((xnwf-xcom)**2+(ynwf-ycom)**2+(znwf-zcom)**2)-2*((xnwf-xcom)*(x0nwf(k)-x0com(k))+ &
           (ynwf-ycom)*(y0nwf(k)-y0com(k))+(znwf-zcom)*(z0nwf(k)-z0com(k)))+ &
           ((x0nwf(k)-x0com(k))**2+(y0nwf(k)-y0com(k))**2+(z0nwf(k)-z0com(k))**2)

      g5fsum=g5fsum+g5f

      g5l=((xnwl-xcom)**2+(ynwl-ycom)**2+(znwl-zcom)**2)-2*((xnwl-xcom)*(x0nwl(k)-x0com(k))+ &
           (ynwl-ycom)*(y0nwl(k)-y0com(k))+(znwl-zcom)*(z0nwl(k)-z0com(k)))+((x0nwl(k)-x0com(k))**2+ &
           (y0nwl(k)-y0com(k))**2+(z0nwl(k)-z0com(k))**2)
      g5lsum=g5lsum+g5l
      !$OMP end parallel sections
    end subroutine build_and_calc_ensemble

    !----------------------------------
    !      calc_average_ensemble
    !----------------------------------
    subroutine calc_average_ensemble

      !Calculate ensemble averages
      dnk=real(nkt,kind=pm_dbl)
      g1avg=g1sum/dnk
      g2avg=g2sum/dnk
      g3avg=g3sum/dnk
      g4favg=g4fsum/dnk
      g4lavg=g4lsum/dnk
      g4avg=(g4favg+g4lavg)/two
      g5favg=g5fsum/dnk
      g5lavg=g5lsum/dnk
      g5avg=(g5favg+g5lavg)/two
    end subroutine calc_average_ensemble


  END SUBROUTINE chaindynamics


  subroutine open_chaindyn_files(dir_name)
    character(len=*) :: dir_name
    open (unit=90,file=dir_name//'g1.dat',form='formatted',status='unknown')
    open (unit=91,file=dir_name//'g2.dat',form='formatted',status='unknown')
    open (unit=92,file=dir_name//'g3.dat',form='formatted',status='unknown')
    open (unit=93,file=dir_name//'g4.dat',form='formatted',status='unknown')
    open (unit=94,file=dir_name//'g5.dat',form='formatted',status='unknown')
    return
  end subroutine open_chaindyn_files

  subroutine close_chaindyn_files
    close (unit=90)
    close (unit=91)
    close (unit=92)
    close (unit=93)
    close (unit=94)
    return
  end subroutine close_chaindyn_files

end module chaindyn
