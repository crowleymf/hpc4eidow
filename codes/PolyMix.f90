!This is the main program for our implementation of the Cooperative Motion
!Algorithm originally developed by Tadeusz Pakula.
!This program is for bipolar shear flow of linear chains.
MODULE param
!Box sizes
INTEGER :: nx, ny, nz, nkt
!Different chain lengths
INTEGER :: nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
INTEGER :: xa,ya,za,xb,yb,zb
!Length of the respective chains
INTEGER,ALLOCATABLE :: nw(:)
!Coordinates of the first monomer in each chain
INTEGER,ALLOCATABLE :: xx(:), yy(:), zz(:)

!Loop counters
INTEGER :: i,j,k,ii,jj,kk,k0
INTEGER :: l,e,lold
DOUBLE PRECISION :: t,tu,tmk


!Total direction code number and direction code number without chain end or a single bead
INTEGER :: code, code1
PARAMETER (code=13, code1=12)
!Forward Direction Codes
INTEGER :: px(code), py(code), pz(code)
!Forward bond code
INTEGER,ALLOCATABLE :: a(:,:,:)
!Reverse bond code
INTEGER,ALLOCATABLE :: b(:,:,:)
!Chain number of the respective number
INTEGER,ALLOCATABLE :: ket(:,:,:)
!New point coordinates in derection
INTEGER,ALLOCATABLE :: xnp(:,:), ynp(:,:), znp(:,:)

CHARACTER*128 :: infile,outfile
INTEGER :: maxloops,maxsta,maxdyn,Nequil,Nauto

INTEGER :: iseed

!Variables for chaincalcs routine
INTEGER :: dx,dy,dz,dxmv,dymv,dzmv
INTEGER :: atemp, xnew,ynew, znew
INTEGER,ALLOCATABLE :: bin(:),bincount(:)

!Variables for controlling the chaindynamics call
INTEGER(4) :: tcount,value,tcountmk,valuemk

!Variables for dynamic calculation
INTEGER,ALLOCATABLE :: xd(:),yd(:),zd(:)
!Variables for vel routine
DOUBLE PRECISION,ALLOCATABLE :: ppy(:),pmy(:),pxz(:)
INTEGER,ALLOCATABLE :: dispy(:,:,:),dispx(:,:,:)
DOUBLE PRECISION,ALLOCATABLE :: vy(:,:,:),vx(:,:,:)
DOUBLE PRECISION,ALLOCATABLE :: avgvy(:),avgvx(:)
DOUBLE PRECISION,ALLOCATABLE :: sumvy(:),sumvx(:)

!Variables for autocorrelate routine
DOUBLE PRECISION,ALLOCATABLE :: ete(:)
DOUBLE PRECISION,ALLOCATABLE :: vyold(:),eteold(:)

!Variable for averaging 
DOUBLE PRECISION :: denom

END MODULE param

!Variables for chaindyanmics routine
MODULE correlation
 INTEGER,ALLOCATABLE :: x0cm(:),y0cm(:),z0cm(:)
 DOUBLE PRECISION,ALLOCATABLE :: x0com(:),y0com(:),z0com(:)
 DOUBLE PRECISION,ALLOCATABLE :: x0nwf(:),y0nwf(:),z0nwf(:) 
 DOUBLE PRECISION,ALLOCATABLE :: x0nwl(:),y0nwl(:),z0nwl(:)
END MODULE correlation
PROGRAM PolyMix

USE param

IMPLICIT NONE

!Number of different chains with different lengths
INTEGER :: nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8
!Direction code number
INTEGER :: d,dcon
!Reverse Direction Codes
INTEGER :: qx(code), qy(code), qz(code)
!Reverse direction code number
INTEGER :: pp(code)
!Variables for biasing flow
DOUBLE PRECISION :: xdiv,pzero,pmax,pnew
!Box coordinates and new point coordinates
INTEGER :: x, y, z
INTEGER :: xn, yn, zn
INTEGER :: xnw, ynw, znw
INTEGER :: pat,c,da,db,pa,pb
INTEGER :: patw,paw,pbw,kw,st,w
INTEGER :: ba,can,cbn
INTEGER :: con(code,code), cn(code1,code)
INTEGER :: cr(code1,code,code), bn(code1,code1)
INTEGER :: ds(code1),dm(code1)
REAL :: ran2



! Save Directory Name
character(len = 15):: dir_name
!Command Line Arguments Needed ==> save directory
IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
  WRITE(*,*)'ERROR, ONE COMMAND-LINE ARGUMENTS REQUIRED, STOPPING; MISSING SAVE DIRECTORY NAME'
  STOP
ENDIF
CALL GET_COMMAND_ARGUMENT(1, dir_name)   !first, read in the value
WRITE (*,*) 'Saving Code to model/runs', dir_name

!Initialization
value=10
l=0
t=0.D0
d=0
pzero=1.D0/3.D0
iseed=10
!st=0
!w=0

OPEN (unit=14,file='data.in',status='old')
 READ (14,*) infile
 READ (14,*) outfile
 READ (14,*) maxloops,maxsta,maxdyn
 READ (14,*) pmax
 READ (14,*) Nequil,Nauto
 READ (14,*) pnew
CLOSE (unit=14)

!Read FCC lattice configuration data from 'conf.in' 
OPEN(unit=15,file='conf.in',status='old')

DO i=1,13
 READ (15,*) qx(i),qy(i),qz(i),px(i),py(i),pz(i),pp(i)
END DO

DO i=1,13
 READ (15,*) (con(i,j),j=1,13)
END DO

DO i=1,12
 READ (15,*) (bn(i,j),j=1,12) 
END DO

DO i=1,12
 DO j=1,13
  READ (15,*) (cr(i,j,k),k=1,13)
 END DO
END DO

DO i=1,12
 READ (15,*) (cn(i,j),j=1,12)
END DO
CLOSE (unit=15)

!Input the created model
!OPEN(unit=9, file='model.bin', form='unformatted',access='sequential')
OPEN(unit=9, file=infile, status='old')

!Read the box dimensions and total chain number
READ (9,*) nx, ny, nz, nkt

ALLOCATE (a(nx,ny,nz),b(nx,ny,nz))
ALLOCATE (nw(nkt),ket(nx,ny,nz))
ALLOCATE (xx(nkt),yy(nkt),zz(nkt))
ALLOCATE (xd(nkt),yd(nkt),zd(nkt))
ALLOCATE (xnp(code1,nx), ynp(code1,ny), znp(code1,nz))
ALLOCATE (dispy(nx,ny,nz),dispx(nx,ny,nz))
ALLOCATE (vy(nx,ny,nz),vx(nx,ny,nz))
ALLOCATE (avgvy(nx),avgvx(ny))
ALLOCATE (sumvy(nx),sumvx(ny))

!Read the chain lengths and chain numbers
READ (9,*) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
READ (9,*) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

!Read the pointers that map new points
Do d=1,12
  READ (9,*) (xnp(d,i),i=1,nx)
  READ (9,*) (ynp(d,i),i=1,ny)
  READ (9,*) (znp(d,i),i=1,nz)
END DO

!Read bond code and chain number for every site
DO i=1,nx
 DO j=1,ny
  DO k=1,nz
   READ (9,*) ii,jj,kk,a(i,j,k),b(i,j,k),ket(i,j,k)
  END DO
 END DO
END DO

!Read the first monomer position of each chain and the corresponding chain length.
DO i=1,nkt
    READ(9,*) xx(i),yy(i),zz(i),nw(i)
END DO
CLOSE(unit=9)

PRINT *, 'data.in,conf.in and model.bin file read in'
!Define fundamental time step
tu=2.D0/(dble(nx)*dble(ny)*dble(nz))

ALLOCATE (ppy(nx),pmy(nx),pxz(nx))

!PRINT *, 'Assigning biasing value. pmax=', pmax
!!Assign the initial biasing to generate the bipolar shear flow
!DO x=1,nx
! xdiv=dble(x)/dble(nx+1)-0.5D0
! ppy(x)=pzero*(1.D0+pmax*xdiv)
! pmy(x)=pzero*(1.D0-pmax*xdiv)
! pxz(x)=pzero
!END DO

PRINT *, 'Assigning biasing value. pmax=', pmax
!Assign the initial biasing to generate the bipolar parabolic flow
DO x=1,nx
 xdiv=dble(x-1)/dble(nx-1)
 ppy(x)=pzero-pmax*(xdiv-xdiv*xdiv)
 pmy(x)=pzero+pmax*(xdiv-xdiv*xdiv)
 pxz(x)=pzero
END DO

PRINT *, 'Bipolar Velocity Profile Assigned'

!Assign the initial xd,yd,and zd array value
DO i=1,nkt
 xd(i)=xx(i)
 yd(i)=yy(i)
 zd(i)=zz(i)
END DO

!Call subroutines for initialization
CALL chaincalcs(dir_name)
CALL boxcalcs(dir_name)
CALL chaindynamics(dir_name)
CALL vel(xn,yn,zn,d, dir_name)


!Beginning of the Monte Carlo routine
!Do the Monte Carlo moves
DO WHILE (l<maxloops)
2 CONTINUE
 c=0
 k0=0

!Randomly select a point on the lattice. We keep trying for a point until we get a kink, or chain end.
 DO WHILE (c<=3 .OR. c==7)
  xn=int(ran2(iseed)*dble(nx))+1
  yn=int(ran2(iseed)*dble(ny))+1
  zn=int(ran2(iseed)*dble(nz))+1
!Access the number of the selected chain; save it as k0.
  k0=ket(xn,yn,zn)
!Get the a and b codes so that the bond angle can be determined
!from data structure con(a,b). If position is a chain "kink" or end,
!con returns a value greater than 3.

  IF (k0>0) THEN
  pa=a(xn,yn,zn)
  pb=b(xn,yn,zn)
  c=con(pa,pb)
  END IF
 END DO

!k is a label which keeps track of the chain on which the tv was first created.
!k will never be a ket number in the box as it is one more than the amount of chains that are present in the box
!e is a label which has value equal to 0, when the tv is created.
!e is reset equal to 1 when the tv returns to the initial lattice position

 k=nkt+1
 e=0
 
!c=4 represents a kink.
!Remove a chain segment and label an adjacent segment.
!Which adjacent segment to label is chosen randomly.

 IF (c==4) THEN
  da=a(xn,yn,zn)
  db=b(xn,yn,zn)
  
  xa=xnp(da,xn)
  ya=ynp(da,yn)
  za=znp(da,zn)
  
  xb=xnp(db,xn)
  yb=ynp(db,yn)
  zb=znp(db,zn)
  
  ba=bn(da,db)
  a(xb,yb,zb)=ba
  b(xa,ya,za)=pp(ba)
  
  IF (ran2(iseed)>0.5D0) THEN
   ket(xa,ya,za)=k
   dcon=da
  ELSE 
   ket(xb,yb,zb)=k
   dcon=db
  END IF

!c=5 specifies a chain end at b=13
!Shorten the chain by one and label the end of the chain
!position with a value of k=nk+1 for the chain number.

 ELSE IF (c==5) THEN
  d=a(xn,yn,zn)
  x=xnp(d,xn)
  y=ynp(d,yn)
  z=znp(d,zn)
  
  b(x,y,z)=13
  ket(x,y,z)=k
  
  dcon=d

  xd(k0)=xd(k0)+px(dcon)
  yd(k0)=yd(k0)+py(dcon)
  zd(k0)=zd(k0)+pz(dcon)
 
  xx(k0)=x
  yy(k0)=y
  zz(k0)=z
  
!c=6 specifies a chain end at a=13
!Shorten the chain by one and label the end of the chain
!position with a value of k=nk+1 for the chain number.

 ELSE IF (c==6) THEN
  d=b(xn,yn,zn)
  x=xnp(d,xn)
  y=ynp(d,yn)
  z=znp(d,zn) 
  
  a(x,y,z)=13
  ket(x,y,z)=k
  
  dcon=d

 END IF
 
 CALL biasd(xn,d, dir_name)
 
 tmk=t
 tcountmk=tcount
 valuemk=value

3 CONTINUE

 xnw=xnp(d,xn)
 
4 IF ((xnw<1).OR.(xnw>nx)) THEN
  CALL biasd(xn,d, dir_name)
  xnw=xnp(d,xn)
  GO TO 4
 END IF
 ynw=ynp(d,yn)
 znw=znp(d,zn)
  
 paw=a(xnw,ynw,znw)
 pbw=b(xnw,ynw,znw)
 patw=cr(d,paw,pbw)
  
 kw=ket(xnw,ynw,znw)
  
 IF (kw==nkt+1) THEN
  IF (c==4) THEN
   IF (ket(xa,ya,za)==nkt+1) THEN
    a(xb,yb,zb)=pp(db)
    b(xa,ya,za)=pp(da)
    ket(xa,ya,za)=k0   
   ELSE
    a(xb,yb,zb)=pp(db)
    b(xa,ya,za)=pp(da)
    ket(xb,yb,zb)=k0
   END IF
   
  ELSEIF (c==5) THEN
   b(x,y,z)=pp(dcon)
   ket(x,y,z)=k0

   xd(k0)=xd(k0)+qx(dcon)
   yd(k0)=yd(k0)+qy(dcon)
   zd(k0)=zd(k0)+qz(dcon)
 
   xx(k0)=xn
   yy(k0)=yn
   zz(k0)=zn

  ELSEIF (c==6) THEN
   a(x,y,z)=pp(dcon)
   ket(x,y,z)=k0
  END IF
  
  t=tmk
  tcount=tcountmk
  value=valuemk
  
  GO TO 2 
 
 ELSEIF (patw==1) THEN
  t=t+tu
  tcount=tcount+1
  
  IF ((mod(tcount,value)==0).AND.(tcount/=0)) THEN
   value=value*1.30-modulo(value*1.30,1E0)
  CALL chaindynamics(dir_name)
   tcount=0
  END IF
   
  CALL biasd(xn,d, dir_name)
 
  GO TO 3
 END IF
 
 CALL vel(xn,yn,zn,pp(dcon),dir_name)

!Recall e is the loop condition 

 DO WHILE (e/=1)
  x=xn
  y=yn
  z=zn
 
!Based on the "old" coordinate and the d-code, obtain the "new"
!trial position (xn,yn,zn) for the tv.  
 
  xn=xnp(d,x)
  
!If the move is through one of the walls, select a new direction.
!Repeat until the tv is not moving through the wall.

!This "IF" loop prevents steps outside the periodic
!boundary and simply makes another biased move.

!Simply put, if the TV tries to step out of the box, try again
!Since this isn't one of the "Monte Carlo Moves" time is not incremented

1 IF((xn<1).OR.(xn>nx)) THEN
   CALL biasd(x,d,dir_name)
   xn=xnp(d,x)
   GO TO 1
  END IF
  
!If walls were set up in other directions in MODLAY, include IF-THEN blocks
!equivalent to the above block for x
  
  yn=ynp(d,y)
  zn=znp(d,z)

!Find the bond directional codes (a,b) at the new point.

  pa=a(xn,yn,zn)
  pb=b(xn,yn,zn)

!Look in the data structure cr(d,a,b) to determine situational code.
!The cr data structure contains every possible combination of d, pa, and pb.

  pat=cr(d,pa,pb)


!Set k equal to the chain number at the new point.

  k=ket(xn,yn,zn)

!If the move is possible, call the velocity routine. 
  
  IF (pat/=1) THEN
   CALL vel(xn,yn,zn,d,dir_name)
  END IF
  
!START OF THE MC MOVES

!Time increments are added only upon entry of the tv into a new chain.
  
!If the tv stays within the same chain, time is not incremented.
!Sections of a chain which can move are "movable groups"; time is incremented
!when a movable group is moved, not when each segment is moved.

!Note: In each loop except in the impossible move loop (pat=1) there is
!a check to see if the TV has returned to its original position by
!checking if the chain number is greater than the total number of chains
!if so this is where the TV was created, e is now 1 and the original
!chain number is returned
 
!pat=1 is an impossible move
!This impossible move is two bonds being simultaneously stretched
!when the TV attempts to step on a kink with a bond angle of 60 degrees
!Generate a new "random" directional code.
!Increment time and set the coordinate of the tv back to the "old" value.
  
  IF (pat==1) THEN
   CALL biasd(x,d,dir_name)
   xn=x
   yn=y
   zn=z
   
   t=t+tu
   tcount=tcount+1
   
!pat=2 represents a move along the chain in a-direction.
!Set the direction of the d code coincident with a.
!Time is not incremented because this is a continuation movement;
!it is part of the movement of a movable group.  
 
  ELSE IF (pat==2) THEN
   d=pa
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    e=1
   END IF

!pat=3 represents a move along the chain in b-direction
!Set the direction of the d code coincident with b.
   
  ELSE IF (pat==3) THEN
     d=pb
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    e=1
   END IF

!pat=4 represents entry into a new chain at b=13 end.
!Set the directional code coincident with the a-direction.
!Time is incremented because a new movable group is going to move.

  ELSE IF (pat==4) THEN
   a(x,y,z)=d
   b(x,y,z)=13
   b(xn,yn,zn)=pp(d)
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    e=1
    k=k0
   END IF
   
   xd(k)=xd(k)+qx(d)
   yd(k)=yd(k)+qy(d)
   zd(k)=zd(k)+qz(d)

   xx(k)=x
   yy(k)=y
   zz(k)=z
   
   ket(x,y,z)=k
   
   d=pa
   
   t=t+tu
   tcount=tcount+1

!pat=5 represents entry into a new chain at a=13 end
!Set the directional code coincident with the b-direction.
!Time is incremented because a new movable group is going to move.
!Same logic as with b=13 but since this represents the chain end
!xx,yy,zz do not need to be updated

  ELSE IF (pat==5) THEN
   b(x,y,z)=d
   a(x,y,z)=13
   a(xn,yn,zn)=pp(d)

   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    e=1
    k=k0
   END IF
   
   ket(x,y,z)=k
   
   d=pb
   
   t=t+tu
   tcount=tcount+1

!pat=6 represents leaving a chain at the b=13 end.
!Generate a new "random" directional code.
!Time is not incremented because this is a continuation movement;
!it completes the movement of a movable group.
    
  ELSE IF (pat==6) THEN
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    e=1
    xx(k0)=xn
    yy(k0)=yn
    zz(k0)=zn
   ELSE
    
    b(x,y,z)=13
    
    xd(k)=xd(k)+qx(d)
    yd(k)=yd(k)+qy(d)
    zd(k)=zd(k)+qz(d)
    
    xx(k)=x
    yy(k)=y
    zz(k)=z
    
    CALL biasd(xn,d,dir_name)
   END IF

!pat=7 represents leaving a chain at a=13 end

  ELSE IF (pat==7) THEN
  IF (k>nkt)THEN
   ket(xn,yn,zn)=k0
   k=k0
   e=1
  ELSE
   a(x,y,z)=13
   CALL biasd(xn,d,dir_name)
  END IF

!pat=8 represents rotation of a b=13 end.
!Generate a new "random" directional code.
!The chain end flip comprises a movable group move so time is incremented
  
  ELSE IF (pat==8) THEN
   can=cn(d,pa)
   a(x,y,z)=can
   
   xa=xnp(pa,xn)
   ya=ynp(pa,yn)
   za=znp(pa,zn)
   
   b(xa,ya,za)=pp(can)
   
   b(x,y,z)=13
   
   ket(x,y,z)=k
   
   IF (k>nkt) THEN     
    xd(k0)=xd(k0)+qx(d)
    yd(k0)=yd(k0)+qy(d)
    zd(k0)=zd(k0)+qz(d) 
    
    xx(k0)=x
    yy(k0)=y
    zz(k0)=z
   ELSE
    xd(k)=xd(k)+qx(d)
    yd(k)=yd(k)+qy(d)
    zd(k)=zd(k)+qz(d)
   
    xx(k)=x
    yy(k)=y
    zz(k)=z
   END IF   
   CALL biasd(xn,d,dir_name)
   
   t=t+tu
   tcount=tcount+1

!pat=9 represents rotation of an a=13 end.
 
  ELSE IF (pat==9) THEN
   cbn=cn(d,pb)
   b(x,y,z)=cbn
   
   xb=xnp(pb,xn)
   yb=ynp(pb,yn)
   zb=znp(pb,zn)
   
   a(xb,yb,zb)=pp(cbn)
   
   a(x,y,z)=13
   
   ket(x,y,z)=k
   
   CALL biasd(xn,d,dir_name)
   
   t=t+tu
   tcount=tcount+1

!pat=10 represents entry into a chain at a "kink" (bond angle of 60 degree)
!in the b-direction. Chain bonds and the tv lie in the same plane.
!Set the direction code to be coincident with the b code at the kink.
!Time is incremented because a new movable group is going to move.
!The TV is also going to actually enter the chain and move along it or
!reverse its step and leave the chain right away.

  ELSE IF (pat==10) THEN
   can=cn(d,pa)
   a(x,y,z)=can
   b(x,y,z)=d
   
   xa=xnp(pa,xn)
   ya=ynp(pa,yn)
   za=znp(pa,zn)
   
   b(xa,ya,za)=pp(can)
   
   a(xn,yn,zn)=pp(d)
   
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    k=k0
    e=1
   END IF
   
   ket(x,y,z)=k
   
   d=pb
   
   t=t+tu
   tcount=tcount+1
   
!pat=11 represents entry into a chain at a kink in a-direction 

  ELSE IF (pat==11) THEN
   cbn=cn(d,pb)
   b(x,y,z)=cbn
   a(x,y,z)=d
   
   xb=xnp(pb,xn)
   yb=ynp(pb,yn)
   zb=znp(pb,zn)
   
   a(xb,yb,zb)=pp(cbn)
   
   b(xn,yn,zn)=pp(d)
   
   IF (k>nkt) THEN
    ket(xn,yn,zn)=k0
    k=k0
    e=1
   END IF
   
   ket(x,y,z)=k
   
   d=pa
   
   t=t+tu
   tcount=tcount+1
 
 
!pat=12 represents two-bond rotation.
!Generate a new "random" directional code.
!This "crankshaft" type move moves a movable group so time is incremented
!This move requires us to complete the triangle twice from the a and b
!side because it needs bonds from the other
!two beads. The TV does not enter the chain
  ELSE IF (pat==12) THEN
   can=cn(d,pa)
   a(x,y,z)=can
   
   xa=xnp(pa,xn)
   ya=ynp(pa,yn)
   za=znp(pa,zn)
   
   b(xa,ya,za)=pp(can)
   
   cbn=cn(d,pb)
   b(x,y,z)=cbn

   xb=xnp(pb,xn)
   yb=ynp(pb,yn)
   zb=znp(pb,zn)
   
   a(xb,yb,zb)=pp(cbn)
   ket(x,y,z)=k
   
   CALL biasd(xn,d,dir_name)
   
   t=t+tu
   tcount=tcount+1

!pat=13 represents leaving a chain at a kink in the a-direction
!Generate a new "random" directional code.
!Time is not incremented because this is a continuation movement;
!it completes the movement of a movable group.
!This move is the opposite of entering at the kink.
  
  ELSE IF (pat==13) THEN
   IF (k>nkt)THEN
    ket(xn,yn,zn)=k0
    k=k0
    e=1
   Else
    can=cn(d,pa)
    a(x,y,z)=can
   
    xa=xnp(pa,xn)
    ya=ynp(pa,yn)
    za=znp(pa,zn)
   
    b(xa,ya,za)=pp(can)
   
    CALL biasd(xn,d,dir_name)
   END IF

!pat=14 represents leaving a chain at a kink in b-direction
 
  ELSE IF (pat==14) THEN
   IF (k>nkt)THEN
    ket(xn,yn,zn)=k0
    k=k0
    e=1
   ELSE
    cbn=cn(d,pb)
    b(x,y,z)=cbn
   
    xb=xnp(pb,xn)
    yb=ynp(pb,yn)
    zb=znp(pb,zn)
   
    a(xb,yb,zb)=pp(cbn)
   
    CALL biasd(xn,d,dir_name)
   END IF


!pat=15 represents motion of a solvent bead
!Generate a new "random" directional code.

  ELSE IF (pat==15) THEN
   a(x,y,z)=13
   b(x,y,z)=13
   
   xx(k)=x
   yy(k)=y
   zz(k)=z
   
   ket(x,y,z)=k
   
   CALL biasd(xn,d,dir_name)
   
   t=t+tu
   tcount=tcount+1             
  END IF
  
  IF ((mod(tcount,value)==0).AND.(tcount/=0)) THEN
   value=value*1.30-modulo(value*1.30,1E0)
  CALL chaindynamics(dir_name)
   tcount=0
  END IF
       
 END DO
 l=l+1
 
 CALL chaincalcs(dir_name)
 
 CALL boxcalcs(dir_name)

!Bipolar shear flow update
 IF (l>=Nequil) THEN
!Update bipolar shear flow
!  DO x=1,nx
!   xdiv=dble(x)/dble(nx+1)-0.5D0
!   ppy(x)=pzero*(1.D0+pnew*xdiv)
!   pmy(x)=pzero*(1.D0-pnew*xdiv)
!   pxz(x)=pzero
!  END DO
  
DO x=1,nx
 xdiv=dble(x-1)/dble(nx-1)
 ppy(x)=pzero-pnew*(xdiv-xdiv*xdiv)
 pmy(x)=pzero+pnew*(xdiv-xdiv*xdiv)
 pxz(x)=pzero
END DO

  CALL autocorrelate(dir_name)
  
 END IF
 
!Output the undated model
 IF (mod(l,maxsta)==0) THEN
  PRINT *,'Writing model file to disk after every maxsta loops'
  PRINT *,'Total loops of:',l,'MC time:',t
  OPEN(unit=21, file=outfile, status='unknown')

!Write the box dimensions and total chain number
  WRITE (21,100) nx, ny, nz, nkt

!Write the chain lengths and chain numbers
  WRITE (21,100) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
  WRITE (21,100) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

!Write the pointers that map new points
  Do d=1,12
   WRITE (21,100) (xnp(d,i),i=1,nx)
   WRITE (21,100) (ynp(d,i),i=1,ny)
   WRITE (21,100) (znp(d,i),i=1,nz)
  END DO

!Write bond code and chain number for every site
  DO i=1,nx
   DO j=1,ny
    DO k=1,nz
     WRITE (21,100) i,j,k,a(i,j,k),b(i,j,k),ket(i,j,k)
    END DO
   END DO
  END DO

!Write the first monomer position of each chain and the corresponding chain length.
  DO i=1,nkt
   WRITE (21,100) xx(i),yy(i),zz(i),nw(i)
  END DO
  CLOSE(unit=21)
 END IF
 
END DO
! Final subroutine call
CALL vel(xn,yn,zn,d, dir_name)

CALL chaindynamics(dir_name)

PRINT *,'Preparing to write out model file to disk for the final time'

!Output the final model

OPEN(unit=21, file=outfile, status='unknown')

!Write the box dimensions and total chain number
WRITE (21,100) nx, ny, nz, nkt

!Write the chain lengths and chain numbers
WRITE (21,100) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
WRITE (21,100) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

!Write the pointers that map new points
Do d=1,12
  WRITE (21,100) (xnp(d,i),i=1,nx)
  WRITE (21,100) (ynp(d,i),i=1,ny)
  WRITE (21,100) (znp(d,i),i=1,nz)
END DO

!Write bond code and chain number for every site
DO i=1,nx
 DO j=1,ny
  DO k=1,nz
   WRITE (21,100) i,j,k,a(i,j,k),b(i,j,k),ket(i,j,k)
  END DO
 END DO
END DO

!Write the first monomer position of each chain and the corresponding chain length.
DO i=1,nkt
 WRITE (21,100) xx(i),yy(i),zz(i),nw(i)
END DO
CLOSE(unit=21)

100 FORMAT(1024(I8,x))   

END PROGRAM PolyMix

!This function is the random number generator used for all simulation.
!Function from numerical recipies      
    FUNCTION ran2(idum)
      INTEGER idum, im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, ndiv
      REAL ran2, am, eps, rnmx
      PARAMETER (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1)
      PARAMETER (ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211)
      PARAMETER (ir2=3791,ntab=32,ndiv=1+imm1/ntab)
      PARAMETER (eps=1.2E-7,rnmx=1.-eps)
      INTEGER idum2, j, k, iv(ntab), iy
      SAVE iv, iy, idum2
      DATA idum2/123456789/, iv/ntab*0/, iy/0/

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
