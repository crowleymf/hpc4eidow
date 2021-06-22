!This subroutine is used to make the random direction moves in this Monte Carlo Scheme.
SUBROUTINE biasd(x,d,dir_name)

use pmtypes
USE param
use pm_subs,only:ran2

IMPLICIT NONE
character(len=*) :: dir_name
INTEGER :: x,d
!---REAL :: ran2
DOUBLE PRECISION :: dtest
Double PRECISION :: pxzx,ppyx,pmyx

pxzx=pxz(x)
ppyx=ppy(x)
pmyx=pmy(x)

iseed=10

dtest=ran2(iseed)

!Move backward
IF (dtest<=pmyx) THEN
 d=int(ran2(iseed)*4)+5
 IF (d==6) d=9
 IF (d==8) d=12
!Move laterally
ELSE IF ((dtest>pmyx).AND.(dtest<=(pmyx+pxzx))) THEN
 d=int(ran2(iseed)*4)+1
!Move forward
ELSE IF (dtest>(pmyx+pxzx)) THEN
 d=int(ran2(iseed)*4)+5
 IF (d==5) d=10
 IF (d==7) d=11
END IF

END SUBROUTINE biasd
