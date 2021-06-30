program polymix
  !This is the main program for our implementation of the Cooperative Motion
  !Algorithm originally developed by Tadeusz Pakula.
  !This program is for bipolar shear flow of linear chains.

  use pmtypes
  use param
  use pm_subs, only: get_pm_save_dir, dir_name, ran2, biasd, pm_close_all_files
  use chaindyn, only: chaindynamics
  use chaincalc, only: chaincalcs

  implicit none

  !Number of different chains with different lengths
  integer :: nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8
  !Direction code number
  integer :: d,dcon
  !Reverse Direction Codes
  integer :: qx(code), qy(code), qz(code)
  !Reverse direction code number
  integer :: pp(code)
  !Variables for biasing flow
  real(kind=pm_dbl) :: xdiv,pzero,pmax,pnew
  !Box coordinates and new point coordinates
  integer :: x, y, z
  integer :: xn, yn, zn
  integer :: xnw, ynw, znw
  integer :: pat,c,da,db,pa,pb
  integer :: patw,paw,pbw,kw,st,w
  integer :: ba,can,cbn
  integer :: con(code,code), cn(code1,code)
  integer :: cr(code1,code,code), bn(code1,code1)
  integer :: ds(code1),dm(code1)
  character(len=10) :: fmta='("PM> ",a)'

  !Initialization
  call get_pm_save_dir
  call init_vars_1  !contained subroutine
  call read_infile
  call read_conf_file
  call read_and_allocate_model
  call init_vars_2

  !Call subroutines for initialization
  call chaincalcs(dir_name)
  call boxcalcs(dir_name)
  call chaindynamics(dir_name)
  call vel(xn,yn,zn,d, dir_name)


  !================================================================A
  !Beginning of the Monte Carlo routine
  !Do the Monte Carlo moves
  !================================================================A
  mcloop: do while (l<maxloops)
     c=0
     k0=0

     !-- Randomly select a point on the lattice. 
     !-- Keep trying for a point until we get a kink, or chain end.
     l1: do while (c<=3 .or. c==7)
        xn=int(ran2(iseed)*real(nx,pm_dbl))+1
        yn=int(ran2(iseed)*real(ny,pm_dbl))+1
        zn=int(ran2(iseed)*real(nz,pm_dbl))+1
        !access the number of the selected chain; save it as k0.
        k0=ket(xn,yn,zn)
        !get the a and b codes so that the bond angle can be determined
        !from data structure con(a,b). if position is a chain "kink" or end,
        !con returns a value greater than 3.

        if (k0>0) then
           pa=a(xn,yn,zn)
           pb=b(xn,yn,zn)
           c=con(pa,pb)
        end if
     end do l1

     !k is a label which keeps track of the chain on which the tv was first created.
     !k will never be a ket number in the box as it is one more than the amount of chains that are present in the box
     !e is a label which has value equal to 0, when the tv is created.
     !e is reset equal to 1 when the tv returns to the initial lattice position

     k=nkt+1
     e=0


     cselect1: select case (c)
     !c=4 represents a kink.
     !Remove a chain segment and label an adjacent segment.
     !Which adjacent segment to label is chosen randomly.
     case (4)
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

        if (ran2(iseed)>0.5d0) then
           ket(xa,ya,za)=k
           dcon=da
        else 
           ket(xb,yb,zb)=k
           dcon=db
        end if

        !c=5 specifies a chain end at b=13
        !shorten the chain by one and label the end of the chain
        !position with a value of k=nk+1 for the chain number.

     case (5)
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
        !shorten the chain by one and label the end of the chain
        !position with a value of k=nk+1 for the chain number.

     case(6)
        d=b(xn,yn,zn)
        x=xnp(d,xn)
        y=ynp(d,yn)
        z=znp(d,zn) 

        a(x,y,z)=13
        ket(x,y,z)=k

        dcon=d

     end select cselect1

     call biasd(xn,d, dir_name)

     tmk=t
     tcountmk=tcount
     valuemk=value

!FIX if this works delete this line !!! 3    continue
     loop3: do while(1 == 1)
        xnw=xnp(d,xn)
        do while ((xnw<1).OR.(xnw>nx)) 
           call biasd(xn,d, dir_name)
           xnw=xnp(d,xn)
        end do
        ynw=ynp(d,yn)
        znw=znp(d,zn)
        paw=a(xnw,ynw,znw)
        pbw=b(xnw,ynw,znw)
        patw=cr(d,paw,pbw)
        kw=ket(xnw,ynw,znw)
        patwblock: if (patw==1 .and. kw/=nkt+1) then
           t=t+tu
           tcount=tcount+1
           if ((mod(tcount,value)==0).and.(tcount/=0)) then
              value=value*1.30-modulo(value*1.30,1e0)
              call chaindynamics(dir_name)
              tcount=0
           end if
           call biasd(xn,d, dir_name)
!FIX if this works delete this line !!go to 3
        else
           exit loop3
        end if patwblock
     end do loop3
     kwblock: if (kw==nkt+1) then
        cselect2: select case (c)
        case(4)
           if (ket(xa,ya,za)==nkt+1) then
              a(xb,yb,zb)=pp(db)
              b(xa,ya,za)=pp(da)
              ket(xa,ya,za)=k0   
           else
              a(xb,yb,zb)=pp(db)
              b(xa,ya,za)=pp(da)
              ket(xb,yb,zb)=k0
           end if

        case(5)
           b(x,y,z)=pp(dcon)
           ket(x,y,z)=k0

           xd(k0)=xd(k0)+qx(dcon)
           yd(k0)=yd(k0)+qy(dcon)
           zd(k0)=zd(k0)+qz(dcon)

           xx(k0)=xn
           yy(k0)=yn
           zz(k0)=zn

        case(6)
           a(x,y,z)=pp(dcon)
           ket(x,y,z)=k0
        end select cselect2

        t=tmk
        tcount=tcountmk
        value=valuemk
        cycle mcloop   
     end if kwblock

     call vel(xn,yn,zn,pp(dcon),dir_name)

     !Recall e is the loop condition 

     !FIX this is a huge loop. Reorganizing into parts/subroutines would be good
     eloop: do while (e/=1)
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
        
        xnloop: do while((xn<1).or.(xn>nx))
           call biasd(x,d,dir_name)
           xn=xnp(d,x)
        end do xnloop

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
        if (pat/=1) call vel(xn,yn,zn,d,dir_name)

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

        patselect: select case (pat)
        case (1)
           call biasd(x,d,dir_name)
           xn=x
           yn=y
           zn=z

           t=t+tu
           tcount=tcount+1

           !pat=2 represents a move along the chain in a-direction.
           !set the direction of the d code coincident with a.
           !time is not incremented because this is a continuation movement;
           !it is part of the movement of a movable group.  

        case(2)
           d=pa
           if (k>nkt) then
              ket(xn,yn,zn)=k0
              e=1
           end if

           !pat=3 represents a move along the chain in b-direction
           !set the direction of the d code coincident with b.

        case(3)
           d=pb
           if (k>nkt) then
              ket(xn,yn,zn)=k0
              e=1
           end if

           !pat=4 represents entry into a new chain at b=13 end.
           !set the directional code coincident with the a-direction.
           !time is incremented because a new movable group is going to move.

        case(4)
           a(x,y,z)=d
           b(x,y,z)=13
           b(xn,yn,zn)=pp(d)
           if (k>nkt) then
              ket(xn,yn,zn)=k0
              e=1
              k=k0
           end if

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
           !set the directional code coincident with the b-direction.
           !time is incremented because a new movable group is going to move.
           !same logic as with b=13 but since this represents the chain end
           !xx,yy,zz do not need to be updated

        case(5)
           b(x,y,z)=d
           a(x,y,z)=13
           a(xn,yn,zn)=pp(d)

           if (k>nkt) then
              ket(xn,yn,zn)=k0
              e=1
              k=k0
           end if

           ket(x,y,z)=k

           d=pb

           t=t+tu
           tcount=tcount+1

           !pat=6 represents leaving a chain at the b=13 end.
           !generate a new "random" directional code.
           !time is not incremented because this is a continuation movement;
           !it completes the movement of a movable group.

        case(6)
           if (k>nkt) then
              ket(xn,yn,zn)=k0
              e=1
              xx(k0)=xn
              yy(k0)=yn
              zz(k0)=zn
           else

              b(x,y,z)=13

              xd(k)=xd(k)+qx(d)
              yd(k)=yd(k)+qy(d)
              zd(k)=zd(k)+qz(d)

              xx(k)=x
              yy(k)=y
              zz(k)=z

              call biasd(xn,d,dir_name)
           end if

           !pat=7 represents leaving a chain at a=13 end

        case (7)
           if (k>nkt)then
              ket(xn,yn,zn)=k0
              k=k0
              e=1
           else
              a(x,y,z)=13
              call biasd(xn,d,dir_name)
           end if

           !pat=8 represents rotation of a b=13 end.
           !generate a new "random" directional code.
           !the chain end flip comprises a movable group move so time is incremented

        case(8)
           can=cn(d,pa)
           a(x,y,z)=can

           xa=xnp(pa,xn)
           ya=ynp(pa,yn)
           za=znp(pa,zn)

           b(xa,ya,za)=pp(can)

           b(x,y,z)=13

           ket(x,y,z)=k

           if (k>nkt) then     
              xd(k0)=xd(k0)+qx(d)
              yd(k0)=yd(k0)+qy(d)
              zd(k0)=zd(k0)+qz(d) 

              xx(k0)=x
              yy(k0)=y
              zz(k0)=z
           else
              xd(k)=xd(k)+qx(d)
              yd(k)=yd(k)+qy(d)
              zd(k)=zd(k)+qz(d)

              xx(k)=x
              yy(k)=y
              zz(k)=z
           end if
           call biasd(xn,d,dir_name)

           t=t+tu
           tcount=tcount+1

           !pat=9 represents rotation of an a=13 end.

        case(9)
           cbn=cn(d,pb)
           b(x,y,z)=cbn

           xb=xnp(pb,xn)
           yb=ynp(pb,yn)
           zb=znp(pb,zn)

           a(xb,yb,zb)=pp(cbn)

           a(x,y,z)=13

           ket(x,y,z)=k

           call biasd(xn,d,dir_name)

           t=t+tu
           tcount=tcount+1

           !pat=10 represents entry into a chain at a "kink" (bond angle of 60 degree)
           !in the b-direction. chain bonds and the tv lie in the same plane.
           !set the direction code to be coincident with the b code at the kink.
           !time is incremented because a new movable group is going to move.
           !the tv is also going to actually enter the chain and move along it or
           !reverse its step and leave the chain right away.

        case(10)
           can=cn(d,pa)
           a(x,y,z)=can
           b(x,y,z)=d

           xa=xnp(pa,xn)
           ya=ynp(pa,yn)
           za=znp(pa,zn)

           b(xa,ya,za)=pp(can)

           a(xn,yn,zn)=pp(d)

           if (k>nkt) then
              ket(xn,yn,zn)=k0
              k=k0
              e=1
           end if

           ket(x,y,z)=k

           d=pb

           t=t+tu
           tcount=tcount+1

           !pat=11 represents entry into a chain at a kink in a-direction 

        case(11)
           cbn=cn(d,pb)
           b(x,y,z)=cbn
           a(x,y,z)=d

           xb=xnp(pb,xn)
           yb=ynp(pb,yn)
           zb=znp(pb,zn)

           a(xb,yb,zb)=pp(cbn)

           b(xn,yn,zn)=pp(d)

           if (k>nkt) then
              ket(xn,yn,zn)=k0
              k=k0
              e=1
           end if

           ket(x,y,z)=k

           d=pa

           t=t+tu
           tcount=tcount+1


           !pat=12 represents two-bond rotation.
           !generate a new "random" directional code.
           !this "crankshaft" type move moves a movable group so time is incremented
           !this move requires us to complete the triangle twice from the a and b
           !side because it needs bonds from the other
           !two beads. the tv does not enter the chain
        case(12)
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

           call biasd(xn,d,dir_name)

           t=t+tu
           tcount=tcount+1

           !pat=13 represents leaving a chain at a kink in the a-direction
           !generate a new "random" directional code.
           !time is not incremented because this is a continuation movement;
           !it completes the movement of a movable group.
           !this move is the opposite of entering at the kink.

        case(13)
           if (k>nkt)then
              ket(xn,yn,zn)=k0
              k=k0
              e=1
           else
              can=cn(d,pa)
              a(x,y,z)=can

              xa=xnp(pa,xn)
              ya=ynp(pa,yn)
              za=znp(pa,zn)

              b(xa,ya,za)=pp(can)

              call biasd(xn,d,dir_name)
           end if

           !pat=14 represents leaving a chain at a kink in b-direction

        case(14)
           if (k>nkt)then
              ket(xn,yn,zn)=k0
              k=k0
              e=1
           else
              cbn=cn(d,pb)
              b(x,y,z)=cbn

              xb=xnp(pb,xn)
              yb=ynp(pb,yn)
              zb=znp(pb,zn)

              a(xb,yb,zb)=pp(cbn)

              call biasd(xn,d,dir_name)
           end if


           !pat=15 represents motion of a solvent bead
           !generate a new "random" directional code.

        case(15)
           a(x,y,z)=13
           b(x,y,z)=13

           xx(k)=x
           yy(k)=y
           zz(k)=z

           ket(x,y,z)=k

           call biasd(xn,d,dir_name)

           t=t+tu
           tcount=tcount+1             
        end select patselect

        if ((mod(tcount,value)==0).and.(tcount/=0)) then
           value=value*1.30-modulo(value*1.30,1e0)
           call chaindynamics(dir_name)
           tcount=0
        end if
     end do eloop

     l=l+1
     call chaincalcs(dir_name)
     call boxcalcs(dir_name)

     !Bipolar shear flow update
     if (l >= nequil) then
        !update bipolar shear flow

        !FIX can we delete this?
        !  do x=1,nx
        !   xdiv=real(x,kind=pm_dbl)/real(nx+1,kind=pm_dbl)-0.5d0
        !   ppy(x)=pzero*(1.d0+pnew*xdiv)
        !   pmy(x)=pzero*(1.d0-pnew*xdiv)
        !   pxz(x)=pzero
        !  end do

        do x=1,nx
           xdiv=real(x-1,kind=pm_dbl)/real(nx-1,kind=pm_dbl)
           ppy(x)=pzero-pnew*(xdiv-xdiv*xdiv)
           pmy(x)=pzero+pnew*(xdiv-xdiv*xdiv)
           pxz(x)=pzero
        end do
        call autocorrelate(dir_name)
     end if

     !output the undated model
     if (mod(l,maxsta)==0) call write_updated_model

  end do mcloop

  ! final subroutine call
  call vel(xn,yn,zn,d, dir_name)

  call chaindynamics(dir_name)

  write(outu,fmta)'Write out model file to disk for the final time'
  call write_final_model
  stop

  !==========================================================================
  ! contained subroutines
  !==========================================================================
contains

  !==========================================================================
  ! init_vars_1
  !==========================================================================
  subroutine init_vars_1
    value=10
    l=0
    t=0.D0
    d=0
    pzero=1.D0/3.D0
    iseed=10
    !st=0
    !w=0
  end subroutine init_vars_1

  !==========================================================================
  ! init_vars_2
  !==========================================================================
  subroutine init_vars_2

    write(outu,fmta)'data.in,conf.in and model.bin file read in'
    !Define fundamental time step
    tu=2.D0/(real(nx,kind=pm_dbl)*real(ny,kind=pm_dbl)*real(nz,kind=pm_dbl))

    allocate (ppy(nx),pmy(nx),pxz(nx))

    !=========================================================================
    !FIX -- can we delete this stuff?
    !FIX -- looks like an alternative way to run the program
    !FIX -- Should be made into actual code with an input parameter to turn it on
    !
    !PRINT *, 'Assigning biasing value. pmax=', pmax
    !!Assign the initial biasing to generate the bipolar shear flow
    !DO x=1,nx
    ! xdiv=real(x,kind=pm_dbl)/real(nx+1,kind=pm_dbl)-0.5D0
    ! ppy(x)=pzero*(1.D0+pmax*xdiv)
    ! pmy(x)=pzero*(1.D0-pmax*xdiv)
    ! pxz(x)=pzero
    !END DO
    !=========================================================================


    write(outu,'("PM> ",a,2x,g15.5)') 'Assigning biasing value. pmax=', pmax
    !Assign the initial biasing to generate the bipolar parabolic flow
    do x=1,nx
       xdiv=real(x-1,kind=pm_dbl)/real(nx-1,kind=pm_dbl)
       ppy(x)=pzero-pmax*(xdiv-xdiv*xdiv)
       pmy(x)=pzero+pmax*(xdiv-xdiv*xdiv)
       pxz(x)=pzero
    end do

    write(outu,fmta)'bipolar velocity profile assigned'

    !Assign the initial xd,yd,and zd array value
    xd(1:nkt)=xx(1:nkt)
    yd(1:nkt)=yy(1:nkt)
    zd(1:nkt)=zz(1:nkt)
    return
  end subroutine init_vars_2

  !==========================================================================
  ! read_infile
  !==========================================================================
  subroutine read_infile
    open (unit=14,file='data.in',status='old')
    read (14,*) infile
    read (14,*) outfile
    read (14,*) maxloops,maxsta,maxdyn
    read (14,*) pmax
    read (14,*) nequil,nauto
    read (14,*) pnew
    close (unit=14)
  end subroutine read_infile

  !==========================================================================
  ! read_conf_file
  !==========================================================================
  subroutine read_conf_file
    !Read FCC lattice configuration data from 'conf.in' 
    open(unit=15,file='conf.in',status='old')

    do i=1,13
       read (15,*) qx(i),qy(i),qz(i),px(i),py(i),pz(i),pp(i)
    end do

    do i=1,13
       read (15,*) (con(i,j),j=1,13)
    end do

    do i=1,12
       read (15,*) (bn(i,j),j=1,12) 
    end do

    do i=1,12
       do j=1,13
          read (15,*) (cr(i,j,k),k=1,13)
       end do
    end do

    do i=1,12
       read (15,*) (cn(i,j),j=1,12)
    end do
    close (unit=15)
  end subroutine read_conf_file

  !==========================================================================
  ! read_and_allocate_model 
  !==========================================================================
  subroutine read_and_allocate_model
    !Input the created model
    !OPEN(unit=9, file='model.bin', form='unformatted',access='sequential')
    OPEN(unit=9, file=infile, status='old')

    !Read the box dimensions and total chain number
    read (9,*) nx, ny, nz, nkt

    allocate (a(nx,ny,nz),b(nx,ny,nz))
    allocate (nw(nkt),ket(nx,ny,nz))
    allocate (xx(nkt),yy(nkt),zz(nkt))
    allocate (xd(nkt),yd(nkt),zd(nkt))
    allocate (xnp(code1,nx), ynp(code1,ny), znp(code1,nz))
    allocate (dispy(nx,ny,nz),dispx(nx,ny,nz))
    allocate (vy(nx,ny,nz),vx(nx,ny,nz))
    allocate (avgvy(nx),avgvx(ny))
    allocate (sumvy(nx),sumvx(ny))

    !Read the chain lengths and chain numbers
    read (9,*) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
    read (9,*) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

    !Read the pointers that map new points
    do d=1,12
       read (9,*) (xnp(d,i),i=1,nx)
       read (9,*) (ynp(d,i),i=1,ny)
       read (9,*) (znp(d,i),i=1,nz)
    end do

    !Read bond code and chain number for every site
    do i=1,nx
       do j=1,ny
          do k=1,nz
             read (9,*) ii,jj,kk,a(i,j,k),b(i,j,k),ket(i,j,k)
          end do
       end do
    end do

    !Read the first monomer position of each chain and the corresponding chain length.
    do i=1,nkt
       read(9,*) xx(i),yy(i),zz(i),nw(i)
    end do
    close(unit=9)
  end subroutine read_and_allocate_model

  !==========================================================================
  ! write_final_model
  !==========================================================================
  subroutine write_final_model
    character(len=10) :: fmtu21="(1024(I8))"

    !Output the final model
    OPEN(unit=21, file=outfile, status='unknown')

    !Write the box dimensions and total chain number
    write (21,fmtu21) nx, ny, nz, nkt

    !Write the chain lengths and chain numbers
    write (21,fmtu21) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
    write (21,fmtu21) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

    !Write the pointers that map new points
    do d=1,12
       write (21,fmtu21) (xnp(d,i),i=1,nx)
       write (21,fmtu21) (ynp(d,i),i=1,ny)
       write (21,fmtu21) (znp(d,i),i=1,nz)
    end do

    !Write bond code and chain number for every site
    WRITE (21,fmtu21) (((i,j,k,a(i,j,k),b(i,j,k),ket(i,j,k),k=1,nz),j=1,ny),i=1,nx)

    !Write the first monomer position of each chain and the corresponding chain length.
    do i=1,nkt
       write (21,fmtu21) xx(i),yy(i),zz(i),nw(i)
    end do
    close(unit=21)

    return
  end subroutine write_final_model

  !==========================================================================
  ! write_updated_model
  !==========================================================================
  subroutine write_updated_model
    character(len=40) :: fmt100="(1024(i8))"
    write(outu,fmta)'writing model file to disk after every maxsta loops'
    write(outu,'("PM> ",a,i8,a,g20.5)')'total loops of :',l,'mc time: ',t
    open(unit=21, file=outfile, status='unknown')

    !write the box dimensions and total chain number
    write (21,fmt100) nx, ny, nz, nkt

    !write the chain lengths and chain numbers
    write (21,fmt100) nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
    write (21,fmt100) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

    !write the pointers that map new points
    do d=1,12
       write (21,fmt100) (xnp(d,i),i=1,nx)
       write (21,fmt100) (ynp(d,i),i=1,ny)
       write (21,fmt100) (znp(d,i),i=1,nz)
    end do

    !write bond code and chain number for every site
    write (21,fmt100) (((i,j,k,a(i,j,k),b(i,j,k),ket(i,j,k),k=1,nz),j=1,ny),i=1,nx)

    !write the first monomer position of each chain and the corresponding chain length.
    do i=1,nkt
       write (21,fmt100) xx(i),yy(i),zz(i),nw(i)
    end do
    close(unit=21)
  end subroutine write_updated_model

end program polymix

