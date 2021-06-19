!This subroutine is used to generate the velocity profile.
SUBROUTINE vel(xn,yn,zn,d, dir_name)

USE param

IMPLICIT NONE
character(len=*) :: dir_name
INTEGER :: x,y,z,xn,yn,zn
INTEGER :: d

!Initialization
IF (t==0) THEN
 OPEN(unit=60,file=dir_name//'velocity_y.dat',status='unknown')
 CLOSE(unit=60)

 OPEN(unit=61,file=dir_name//'velocity_x.dat',status='unknown')
 CLOSE(unit=61)
 !$OMP PARALLEL DO COLLAPSE(3)
 DO x=1,nx
  DO y=1,ny
   DO z=1,nz
    dispy(x,y,z)=0
    dispx(x,y,z)=0
    vy(x,y,z)=0.D0
    vx(x,y,z)=0.D0
   END DO
  END DO
 END DO
 !$OMP END PARALLEL DO
 lold=0
END IF

!If maxsta loops have been reached, write out the updated velocities
IF ((mod(l,maxsta)==0).AND.(l/=lold)) THEN
 !$OMP PARALLEL
 !$OMP DO COLLAPSE(3) 
 DO x=1,nx
  DO y=1,ny
   DO z=1,nz
    vy(x,y,z)=dble(dispy(x,y,z))/t
    vx(x,y,z)=dble(dispx(x,y,z))/t
   END DO
  END DO
 END DO
 !$OMP END DO

 !$OMP DO
 DO x=1,nx
  sumvy(x)=0.D0
  avgvy(x)=0.D0
  DO y=1,ny
   DO z=1,nz
    sumvy(x)=sumvy(x)+vy(x,y,z)
   END DO
  END DO
  !The average value for the velocity of a given yz-plane is the total velocity divided by
  !the number of lattice sites in a plane (ny*nz/2).
  avgvy(x)=sumvy(x)/dble((ny*nz/2))
 END DO
 !$OMP END DO
 !$OMP END PARALLEL

 !$OMP DO
 DO y=1,ny
  sumvx(y)=0.D0
  avgvx(y)=0.D0
  DO x=1,nx
   DO z=1,nz
    sumvx(y)=sumvx(y)+vx(x,y,z)
   END DO
  END DO
   avgvx(y)=sumvx(y)/dble((nx*nz/2))
 END DO
 !$OMP END DO
 !Write out velocity values to data file
 !$OMP SECTIONS 
 !$OMP SECTION
 OPEN(unit=60,file=dir_name//'velocity_y.dat',status='old',position='append')
  WRITE(60,*)'l=',l
  DO x=1,nx
   WRITE(60,*)x,avgvy(x)
  END DO
 CLOSE(unit=60)
 !$OMP SECTION
 OPEN(unit=61,file=dir_name//'velocity_x.dat',status='old',position='append')
  WRITE(61,*)'l=',l
  DO y=1,ny
   WRITE(61,*)y,avgvx(y)
  END DO
 CLOSE(unit=61)
 !$OMP END SECTIONS
 lold=l
END IF

!TV moves in the -y direction => bead moves in +y direction
IF ((d==5) .OR. (d==7) .OR. (d==9) .OR. (d==12)) dispy(xn,yn,zn) = dispy(xn,yn,zn) + 1

!TV moves in the +y direction => bead moves in -y direction
IF ((d==6) .OR. (d==8) .OR. (d==10) .OR. (d==11)) dispy(xn,yn,zn) = dispy(xn,yn,zn) - 1

!TV moves in the -x direction => bead moves in +x direction
IF ((d==1) .OR. (d==3) .OR. (d==9) .OR. (d==11)) dispx(xn,yn,zn) = dispx(xn,yn,zn) + 1

!TV moves in the +x direction => bead moves in -x direction
IF ((d==2) .OR. (d==4) .OR. (d==10) .OR. (d==12)) dispx(xn,yn,zn) = dispx(xn,yn,zn) - 1

END SUBROUTINE vel
