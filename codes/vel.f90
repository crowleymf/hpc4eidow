!This subroutine is used to generate the velocity profile.
SUBROUTINE vel(xn,yn,zn,d, dir_name)

USE param
use pmtypes
IMPLICIT NONE
character(len=*) :: dir_name
INTEGER :: x,y,z,xn,yn,zn
INTEGER :: d,dmask
!
! print *,  "kind: ===> ", kind(dmask)
! print *,  "3000: " , kind(3000)
! dmask = 7
! print , "#: ", kind(dmask)
! stop

!Initialization
if (t==0) THEN
  call initalize
  lold  = 0
end if
!if maxsta loops have been reached
if ((mod(l,maxsta)==0).AND.(l/=lold)) THEN
  call update_velocity
  call calc_avg_vel
  call write_velocity_output
  lold=l
end if
if (d < 13) call make_mc_move

contains
  subroutine initalize
    dmask = 0
    open(unit=60,file=dir_name//'velocity_y.dat',status='unknown')
    close(unit=60)

    open(unit=61,file=dir_name//'velocity_x.dat',status='unknown')
    close(unit=61)

    dispy = 0
    dispx = 0
    vy    = 0.D0
    vx    = 0.D0
  end subroutine initalize

  subroutine update_velocity
    do x=1,nx
      do y=1,ny
        do z=1,nz
          vy(x,y,z)=real(dispy(x,y,z),kind=pm_dbl)/t
          vx(x,y,z)=real(dispx(x,y,z),kind=pm_dbl)/t
        end do
      end do
    end do

    do x=1,nx
      sumvy(x)=0.D0
      avgvy(x)=0.D0
      do y=1,ny
        do z=1,nz
          sumvy(x)=sumvy(x)+vy(x,y,z)
        end do
      end do
      avgvy(x)=sumvy(x)/real((ny*nz/2), kind=pm_dbl)
    end do
  end subroutine update_velocity

  subroutine calc_avg_vel
    !The average value for the velocity of a given yz-plane is the total velocity divided by
    !the number of lattice sites in a plane (ny*nz/2).
    do y=1,ny
      sumvx(y)=0.D0
      avgvx(y)=0.D0
      do x=1,nx
        do z=1,nz
          sumvx(y)=sumvx(y)+vx(x,y,z)
        end do
      end do
      avgvx(y)=sumvx(y)/real((nx*nz/2), kind=pm_dbl)
    end do
  end subroutine calc_avg_vel

  subroutine write_velocity_output
    !write out velocity values to data file
    open(unit=60,file=dir_name//'velocity_y.dat',status='old',position='append')
      write(60,*)'l=',l
      do x=1,nx
        write(60,*) x, avgvy(x)
      end do
    close(unit=60)

    open(unit=61,file=dir_name//'velocity_x.dat',status='old',position='append')
      write(61,*)'l=',l
      do y=1,ny
        write(61,*)y,avgvx(y)
      end do
    close(unit=61)
  end subroutine write_velocity_output

  subroutine make_mc_move
    dmask = ishft(1,d-1)

    dispy(xn,yn,zn) = dispy(xn,yn,zn) + ishft(iand(dmask,2384), 1-d) - ishft(iand(dmask,1696), 1-d)
    dispx(xn,yn,zn) = dispx(xn,yn,zn) + ishft(iand(dmask,1285), 1-d) - ishft(iand(dmask,2570), 1-d)
  end subroutine make_mc_move

end subroutine vel
