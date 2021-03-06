module chaincalc
  use pmtypes

  implicit none


  integer x, y, z
  integer bincheck
  integer point, loop
  real(kind=pm_dbl), allocatable:: dxcom(:), dycom(:), dzcom(:)

  real(kind=pm_dbl) parasquared, perpsquared
  real(kind=pm_dbl), allocatable:: etepara(:), eteperp(:)
  real(kind=pm_dbl) eteboxtot, eteparatot, eteperptot
  real(kind=pm_dbl) etebox, eteparabox, eteperpbox
  real(kind=pm_dbl) eteboxtemp,eteparaboxtemp, eteperpboxtemp
  real(kind=pm_dbl) eteboxsum, eteparaboxsum, eteperpboxsum
  real(kind=pm_dbl) eteboxsumdev, eteparaboxsumdev, eteperpboxsumdev
  real(kind=pm_dbl) etefinal, eteparafinal, eteperpfinal
  real(kind=pm_dbl) etedev, eteparadev, eteperpdev
  real(kind=pm_dbl) maxete

  real(kind=pm_dbl), allocatable:: etebintot(:), eteparabintot(:), eteperpbintot(:)
  real(kind=pm_dbl), allocatable:: etebin(:), eteparabin(:), eteperpbin(:)
  real(kind=pm_dbl) etebintemp, eteparabintemp, eteperpbintemp
  real(kind=pm_dbl), allocatable:: etebinsum(:), eteparabinsum(:), eteperpbinsum(:)
  real(kind=pm_dbl), allocatable:: etebinsumdev(:),eteparabinsumdev(:), eteperpbinsumdev(:)
  real(kind=pm_dbl) etebinfinal, eteparabinfinal,eteperpbinfinal
  real(kind=pm_dbl) etebindev, eteparabindev, eteperpbindev

  real(kind=pm_dbl) dxdiff, dydiff, dzdiff
  real(kind=pm_dbl) rogsqtot,rogsq
  real(kind=pm_dbl), allocatable:: rog(:)
  real(kind=pm_dbl) rogtot
  real(kind=pm_dbl) rogbox
  real(kind=pm_dbl) rogboxtemp
  real(kind=pm_dbl) rogboxsum
  real(kind=pm_dbl) rogboxsumdev
  real(kind=pm_dbl) rogfinal
  real(kind=pm_dbl) rogdev

  real(kind=pm_dbl), allocatable:: rogbintot(:)
  real(kind=pm_dbl), allocatable:: rogbin(:)
  real(kind=pm_dbl) rogbintemp
  real(kind=pm_dbl), allocatable:: rogbinsum(:)
  real(kind=pm_dbl), allocatable:: rogbinsumdev(:)
  real(kind=pm_dbl) rogbinfinal
  real(kind=pm_dbl) rogbindev

  real(kind=pm_dbl), allocatable:: thetax(:), thetay(:), thetaz(:)
  real(kind=pm_dbl) thetaxtot, thetaytot, thetaztot
  real(kind=pm_dbl) thetaxbox, thetaybox, thetazbox
  real(kind=pm_dbl) thetaxboxtemp, thetayboxtemp, thetazboxtemp
  real(kind=pm_dbl) thetaxsum, thetaysum, thetazsum
  real(kind=pm_dbl) thetaxsumdev, thetaysumdev, thetazsumdev
  real(kind=pm_dbl) thetaxfinal, thetayfinal, thetazfinal
  real(kind=pm_dbl) thetaxdev, thetaydev, thetazdev

  real(kind=pm_dbl), allocatable:: thetaxbintot(:), thetaybintot(:), thetazbintot(:)
  real(kind=pm_dbl), allocatable:: thetaxbin(:), thetaybin(:), thetazbin(:)
  real(kind=pm_dbl) thetaxbintemp, thetaybintemp, thetazbintemp
  real(kind=pm_dbl), allocatable:: thetaxbinsum(:), thetaybinsum(:), thetazbinsum(:)
  real(kind=pm_dbl), allocatable:: thetaxbinsumdev(:), thetaybinsumdev(:), thetazbinsumdev(:)
  real(kind=pm_dbl) thetaxbinfinal, thetaybinfinal, thetazbinfinal
  real(kind=pm_dbl) thetaxbindev, thetaybindev, thetazbindev

  integer, allocatable:: binnk1(:), binnk2(:), binnk3(:), binnk4(:), binnk5(:), binnk6(:), binnk7(:), binnk8(:)
  real(kind=pm_dbl), allocatable:: nk1avg(:), nk2avg(:), nk3avg(:), nk4avg(:), nk5avg(:), nk6avg(:), nk7avg(:), nk8avg(:)
  real(kind=pm_dbl) nk1temp, nk2temp, nk3temp, nk4temp, nk5temp, nk6temp, nk7temp, nk8temp
  real(kind=pm_dbl), allocatable:: nk1sum(:), nk2sum(:), nk3sum(:), nk4sum(:), nk5sum(:), nk6sum(:), nk7sum(:), nk8sum(:)
  real(kind=pm_dbl) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

  integer, allocatable:: chainend(:)
  real(kind=pm_dbl)  chainbintemp
  real(kind=pm_dbl), allocatable:: chainsum(:)
  real(kind=pm_dbl), allocatable:: chainsumdev(:)
  real(kind=pm_dbl) chainfinal
  real(kind=pm_dbl) chaindev

  character(len=100) :: fmt300="(3(A20,x))"
  character(len=100) :: fmt301="(5(A20,x))"
  character(len=100) :: fmt302="(7(A20,x))"
  character(len=100) :: fmt303="(9(A20,x))"
  character(len=100) :: fmt304="(I20,x,2(f20.14,x))"
  character(len=100) :: fmt305="(I20,x,4(f20.14,x))"
  character(len=100) :: fmt306="(I20,x,6(f20.14,x))"
  character(len=100) :: fmt307="(I20,x,8(f20.14,x))"
  character(len=100) :: fmt308="(8(A20,x))"

contains
  subroutine chaincalcs(dir_name)
    !this subroutine is used to calculate the end-to-end vector, radius of gyration,
    !order parameters and center of mass density
    use param
    implicit none
    character(len=*) :: dir_name

    dir_name = trim(adjustl(dir_name))

    !initialization
    if (t == 0) then
       allocate (ete(nkt))
       call open_chaincalcs_files(dir_name)
       call write_chaincalcs_headers
    endif

    !array allocation
    call allocate_chaincalcs_arrays
    call initalize

    do k = 1, nkt
       call iteration_init
       call calc_single_bead_data


       !calculate the maximum ete
       if(k==1) then
          maxete=ete(1)
       end if

       if (ete(k)>maxete) then
          maxete=ete(k)
       end if

    enddo
    call calc_box_avg

    !file unit 30 is the maxdyn loop and 31 is maxsta
    write(30) etebox, rogbox
    write(31) etebox, eteparabox, eteperpbox, rogbox, thetaxbox, thetaybox, thetazbox
    !calculate the binned averages
    do x = 1, nx
       denom = real(bincount(x))
       call bin_avg
       !$OMP parallel sections
       !dynamic (calculations made every maxdyn loop) temp file
       write(32) x, etebin(x), rogbin(x)
       !$OMP section
       !static (calculations made every maxsta loops) temp files
       write(33) x, etebin(x), eteparabin(x), eteperpbin(x), rogbin(x), real(chainend(x),kind=pm_dbl)

       !$OMP section
       write(34) x, nk1avg(x), nk2avg(x), nk3avg(x), nk4avg(x), nk5avg(x), nk6avg(x), nk7avg(x), nk8avg(x)
       !$OMP section
       write(35) x, thetaxbin(x), thetaybin(x), thetazbin(x)
       !$OMP end parallel sections
    enddo

    !averaging after maxdyn loops
    if (mod(l,maxdyn)==0 .and. l /= 0) then
       eteboxsum = zero
       rogboxsum = zero

       eteboxsumdev=zero
       rogboxsumdev=zero

       ALLOCATE (etebinsum(nx), rogbinsum(nx))
       ALLOCATE (etebinsumdev(nx), rogbinsumdev(nx))

       DO x = 1, nx
          etebinsum(x) = zero
          rogbinsum(x) = zero

          etebinsumdev(x)=zero
          rogbinsumdev(x)=zero
       enddo

       rewind (30)
       rewind (32)

       do loop = 1, maxdyn
          read(30) eteboxtemp, rogboxtemp
          eteboxsum = eteboxsum + eteboxtemp
          rogboxsum = rogboxsum + rogboxtemp

          do x = 1, nx
             read(32) point, etebintemp, rogbintemp
             etebinsum(point) = etebinsum(point) + etebintemp
             rogbinsum(point) = rogbinsum(point) + rogbintemp
          enddo
       enddo

       denom = real(maxdyn)

       rewind(30)
       do loop=1,maxdyn
          read(30) eteboxtemp, rogboxtemp
          eteboxsumdev = eteboxsumdev + (eteboxtemp-eteboxsum/denom)**2
          rogboxsumdev = rogboxsumdev + (rogboxtemp-rogboxsum/denom)**2
       end do
       etefinal = eteboxsum/denom
       etedev = sqrt(eteboxsumdev/(denom-one))

       rogfinal = rogboxsum/denom
       rogdev = sqrt(rogboxsumdev/(denom-one))

       write(40,fmt306) l, etefinal, etedev, rogfinal, rogdev, maxete, t


       rewind(32)
       write(42,*) 'l=', l

       do loop=1,maxdyn
          do x=1,nx
             read(32) point, etebintemp, rogbintemp
             etebinsumdev(point)=etebinsumdev(point)+(etebintemp-etebinsum(x)/denom)**2
             rogbinsumdev(point)=rogbinsumdev(point)+(rogbintemp-rogbinsum(x)/denom)**2
          end do
       end do

       do x=1,nx
          etebinfinal=etebinsum(x)/denom
          etebindev=sqrt(etebinsumdev(x)/(denom-1))
          rogbinfinal=rogbinsum(x)/denom
          rogbindev=sqrt(rogbinsumdev(x)/(denom-1))
          write(42,fmt305) x, etebinfinal, etebindev, rogbinfinal, rogbindev
       end do

       rewind(30)
       rewind(32)

       deallocate (etebinsum, rogbinsum)
       deallocate (etebinsumdev, rogbinsumdev)
    ENDIF
    !Averaging after maxsta loops
    IF (mod(l, maxsta) == 0 .AND. l /= 0) THEN
       call zero_average_init

       REWIND (31)
       REWIND (33)
       REWIND (34)
       REWIND (35)

       DO loop = 1, maxsta
          call maxsta_avg_calc
       ENDDO

       denom = real(maxsta,kind=pm_dbl)

       REWIND(31)
       DO loop=1,maxsta
          READ(31) eteboxtemp, eteparaboxtemp, eteperpboxtemp, rogboxtemp, thetaxboxtemp, thetayboxtemp, thetazboxtemp
          eteboxsumdev=eteboxsumdev+(eteboxtemp-eteboxsum/denom)**2
          eteparaboxsumdev=eteparaboxsumdev+(eteparaboxtemp-eteparaboxsum/denom)**2
          eteperpboxsumdev=eteperpboxsumdev+(eteperpboxtemp-eteperpboxsum/denom)**2
          rogboxsumdev=rogboxsumdev+(rogboxtemp-rogboxsum/denom)**2
          thetaxsumdev=thetaxsumdev+(thetaxboxtemp-thetaxsum/denom)**2
          thetaysumdev=thetaysumdev+(thetayboxtemp-thetaysum/denom)**2
          thetazsumdev=thetazsumdev+(thetazboxtemp-thetazsum/denom)**2
       END DO
       etefinal = eteboxsum/denom
       eteparafinal = eteparaboxsum/denom
       eteperpfinal = eteperpboxsum/denom

       etedev = sqrt(eteboxsumdev/(denom-one))
       eteparadev=sqrt(eteparaboxsumdev/(denom-one))
       eteperpdev=sqrt(eteperpboxsumdev/(denom-one))

       rogfinal = rogboxsum/denom
       rogdev = sqrt(rogboxsumdev/(denom-one))

       thetaxfinal = thetaxsum/denom
       thetayfinal = thetaysum/denom
       thetazfinal = thetazsum/denom

       thetaxdev = sqrt(thetaxsumdev/(denom-one))
       thetaydev = sqrt(thetaysumdev/(denom-one))
       thetazdev = sqrt(thetazsumdev/(denom-one))

       write(41,fmt307) l, etefinal, etedev, eteparafinal, eteparadev, eteperpfinal, eteperpdev, rogfinal, rogdev

       write(43,fmt306) l, thetaxfinal, thetaxdev, thetayfinal, thetaydev, thetazfinal, thetazdev

       rewind(33)
       rewind(35)

       open(unit = 44, file=dir_name//'chainends.dat', position = 'append')
       write(44,*) 'l=', l

       open(unit = 45, file=dir_name//'binchainproperties.dat', position = 'append')
       write(45,*) 'l=', l

       open(unit = 46, file=dir_name//'nkcom.dat', position = 'append')
       write(46,*) 'l=', l

       open(unit = 47, file=dir_name//'binorder.dat', position = 'append')
       write(47,*) 'l=', l

       do loop=1,maxsta
          do x=1,nx
             read(33) point, etebintemp, eteparabintemp, eteperpbintemp, rogbintemp, chainbintemp
             etebinsumdev(point) = etebinsumdev(point) + (etebintemp-etebinsum(point)/denom)**2
             eteparabinsumdev(point) = eteparabinsumdev(point) + (eteparabintemp-eteparabinsum(point)/denom)**2
             eteperpbinsumdev(point) = eteperpbinsumdev(point) + (eteperpbintemp-eteperpbinsum(point)/denom)**2
             rogbinsumdev(point) = rogbinsumdev(point) + (rogbintemp-rogbinsum(point)/denom)**2
             chainsumdev(point) = chainsumdev(point) + (chainbintemp-chainsum(point)/denom)**2

             READ(35) point, thetaxbintemp, thetaybintemp, thetazbintemp
             thetaxbinsumdev(point) = thetaxbinsumdev(point) + (thetaxbintemp-thetaxbinsum(point)/denom)**2
             thetaybinsumdev(point) = thetaybinsumdev(point) + (thetaybintemp-thetaybinsum(point)/denom)**2
             thetazbinsumdev(point) = thetazbinsumdev(point) + (thetazbintemp-thetazbinsum(point)/denom)**2
          END DO
       END DO

       DO x = 1, nx
          chainfinal = chainsum(x)/denom
          chaindev = sqrt(chainsumdev(x)/(denom-one))
          WRITE(44,fmt304) x, chainfinal, chaindev

          etebinfinal = etebinsum(x)/denom
          rogbinfinal = rogbinsum(x)/denom
          eteparabinfinal = eteparabinsum(x)/denom
          eteperpbinfinal = eteperpbinsum(x)/denom

          etebindev = sqrt(etebinsumdev(x)/(denom-one))
          eteparabindev = sqrt(eteparabinsumdev(x)/(denom-one))
          eteperpbindev = sqrt(eteperpbinsumdev(x)/(denom-one))
          rogbindev = sqrt(rogbinsumdev(x)/(denom-one))


          nk1 = nk1sum(x)/denom
          nk2 = nk2sum(x)/denom
          nk3 = nk3sum(x)/denom
          nk4 = nk4sum(x)/denom
          nk5 = nk5sum(x)/denom
          nk6 = nk6sum(x)/denom
          nk7 = nk7sum(x)/denom
          nk8 = nk8sum(x)/denom

          thetaxbinfinal = thetaxbinsum(x)/denom
          thetaybinfinal = thetaybinsum(x)/denom
          thetazbinfinal = thetazbinsum(x)/denom

          thetaxbindev = sqrt(thetaxbinsumdev(x)/(denom-one))
          thetaybindev = sqrt(thetaybinsumdev(x)/(denom-one))
          thetazbindev = sqrt(thetazbinsumdev(x)/(denom-one))
          !$OMP parallel sections
          !$OMP section
          WRITE(45,fmt307) x, etebinfinal, etebindev, eteparabinfinal, &
               eteparabindev, eteperpbinfinal, eteperpbindev, &
               rogbinfinal, rogbindev
          !$OMP section
          WRITE(46,fmt308) x, nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8
          !$OMP section
          WRITE(47,fmt306) x, thetaxbinfinal, thetaxbindev, thetaybinfinal, thetaybindev, thetazbinfinal, thetazbindev
          !$OMP end parallel sections

       ENDDO

       REWIND(31)
       REWIND(33)
       REWIND(34)
       REWIND(35)

       !$OMP parallel sections
       !$OMP section
       CLOSE (unit = 44)
       CLOSE (unit = 45)
       CLOSE (unit = 46)
       CLOSE (unit = 47)

       !$OMP section
       DEALLOCATE (etebinsum, eteparabinsum, eteperpbinsum, rogbinsum, chainsum)
       DEALLOCATE (etebinsumdev, eteparabinsumdev, eteperpbinsumdev, rogbinsumdev, chainsumdev)
       DEALLOCATE (nk1sum, nk2sum, nk3sum, nk4sum, nk5sum, nk6sum, nk7sum, nk8sum)
       DEALLOCATE (thetaxbinsum, thetaybinsum, thetazbinsum)
       DEALLOCATE (thetaxbinsumdev, thetaybinsumdev, thetazbinsumdev)
       !$OMP end parallel sections
    ENDIF

    call cleanup


  contains
    subroutine initalize
      eteboxtot = zero
      eteparatot = zero
      eteperptot = zero

      rogtot = zero

      thetaxtot = zero
      thetaytot = zero
      thetaztot = zero
      !$OMP parallel do
      DO x  = 1, nx
         chainend(x) = 0
         etebintot(x) = zero
         eteparabintot(x) = zero
         eteperpbintot(x) = zero

         rogbintot(x) = zero

         thetaxbintot(x) = zero
         thetaybintot(x) = zero
         thetazbintot(x) = zero

         bin(x) = 0
         bincount(x) = 0

         binnk1(x) = 0
         binnk2(x) = 0
         binnk3(x) = 0
         binnk4(x) = 0
         binnk5(x) = 0
         binnk6(x) = 0
         binnk7(x) = 0
         binnk8(x) = 0
      ENDDO
    end subroutine initalize

    subroutine zero_average_init
      eteboxsum = zero
      eteparaboxsum = zero
      eteperpboxsum = zero

      eteboxsumdev = zero
      eteparaboxsumdev = zero
      eteperpboxsumdev = zero

      rogboxsum = zero
      rogboxsumdev = zero

      thetaxsum = zero
      thetaysum = zero
      thetazsum = zero

      thetaxsumdev = zero
      thetaysumdev = zero
      thetazsumdev = zero

      ALLOCATE (etebinsum(nx), eteparabinsum(nx), eteperpbinsum(nx), rogbinsum(nx), chainsum(nx))
      ALLOCATE (etebinsumdev(nx), eteparabinsumdev(nx), eteperpbinsumdev(nx),rogbinsumdev(nx),chainsumdev(nx))
      ALLOCATE (nk1sum(nx), nk2sum(nx), nk3sum(nx), nk4sum(nx), nk5sum(nx), nk6sum(nx), nk7sum(nx), nk8sum(nx))
      ALLOCATE (thetaxbinsum(nx), thetaybinsum(nx), thetazbinsum(nx))
      ALLOCATE (thetaxbinsumdev(nx), thetaybinsumdev(nx), thetazbinsumdev(nx))

      DO x = 1, nx
         etebinsum(x) = zero
         eteparabinsum(x) = zero
         eteperpbinsum(x) = zero

         rogbinsum(x) = zero

         chainsum(x) = zero

         etebinsumdev(x) = zero
         eteparabinsumdev(x) = zero
         eteperpbinsumdev(x) = zero

         rogbinsumdev(x) = zero

         chainsumdev(x) = zero

         nk1sum(x) = zero
         nk2sum(x) = zero
         nk3sum(x) = zero
         nk4sum(x) = zero
         nk5sum(x) = zero
         nk6sum(x) = zero
         nk7sum(x) = zero
         nk8sum(x) = zero

         thetaxbinsum(x) = zero
         thetaybinsum(x) = zero
         thetazbinsum(x) = zero

         thetaxbinsumdev(x) = zero
         thetaybinsumdev(x) = zero
         thetazbinsumdev(x) = zero
      ENDDO
    end subroutine zero_average_init

    subroutine zero_bincount_init
      etebin(x) = zero
      eteparabin(x) = zero
      eteperpbin(x) = zero

      rogbin(x) = zero

      thetaxbin(x) = zero
      thetaybin(x) = zero
      thetazbin(x) = zero

      nk1avg(x) =  zero
      nk2avg(x) =  zero
      nk3avg(x) =  zero
      nk4avg(x) =  zero
      nk5avg(x) =  zero
      nk6avg(x) =  zero
      nk7avg(x) =  zero
      nk8avg(x) =  zero
    end subroutine zero_bincount_init

    subroutine iteration_init
      x = xx(k)
      y = yy(k)
      z = zz(k)
      !The dx, dy, dz are the components of the end-to-end vector
      dx = 0
      dy = 0
      dz = 0
      !The mv is for the center of mass
      dxmv = 0
      dymv = 0
      dzmv = 0

      atemp = a(x,y,z)

      chainend(x) = chainend(x) + 1
      DO WHILE (atemp /= 13)
         dx = dx + px(atemp)
         dy = dy + py(atemp)
         dz = dz + pz(atemp)

         dxmv = dxmv + dx
         dymv = dymv + dy
         dzmv = dzmv + dz

         xnew = xnp(atemp, x)
         ynew = ynp(atemp, y)
         znew = znp(atemp, z)

         x = xnew
         y = ynew
         z = znew

         atemp = a(x,y,z)

         IF (atemp == 13) THEN
            chainend(x) = chainend(x) + 1
         ENDIF
      ENDDO
    end subroutine iteration_init

    subroutine zero_bead
      !Single bead considerations
      dxcom(k) = zero
      dycom(k) = zero
      dzcom(k) = zero

      ete(k) = zero
      etepara(k) = zero
      eteperp(k) = zero

      thetax(k) = zero
      thetay(k) = zero
      thetaz(k) = zero
    end subroutine zero_bead

    subroutine box_average
      !The box averages are number averages
      !The totals are divided by nkt at the end
      thetaxtot = thetaxtot + (three*thetax(k)*thetax(k) - one)/two
      thetaytot = thetaytot + (three*thetay(k)*thetay(k) - one)/two
      thetaztot = thetaztot + (three*thetaz(k)*thetaz(k) - one)/two

      eteboxtot = eteboxtot + ete(k)
      eteparatot = eteparatot + etepara(k)
      eteperptot = eteperptot + eteperp(k)

      bincheck = xx(k) + int(dxcom(k))
      bin(k) = bincheck
      bincount(bin(k)) = bincount(bin(k)) + 1

      IF (nw(k) == nw1) THEN
         binnk1(bin(k)) = binnk1(bin(k)) + 1
      ELSEIF (nw(k) == nw2) THEN
         binnk2(bin(k)) = binnk2(bin(k)) + 1
      ELSEIF (nw(k) == nw3) THEN
         binnk3(bin(k)) = binnk3(bin(k)) + 1
      ELSEIF (nw(k) == nw4) THEN
         binnk4(bin(k)) = binnk4(bin(k)) + 1
      ELSEIF (nw(k) == nw5) THEN
         binnk5(bin(k)) = binnk5(bin(k)) + 1
      ELSEIF (nw(k) == nw6) THEN
         binnk6(bin(k)) = binnk6(bin(k)) + 1
      ELSEIF (nw(k) == nw7) THEN
         binnk7(bin(k)) = binnk7(bin(k)) + 1
      ELSEIF (nw(k) == nw8) THEN
         binnk8(bin(k)) = binnk8(bin(k)) + 1
      ENDIF
    end subroutine box_average

    subroutine radii_of_gyration
      dx = dx + px(atemp)
      dy = dy + py(atemp)
      dz = dz + pz(atemp)

      dxdiff = real(dx,kind=pm_dbl) - dxcom(k)
      dydiff = real(dy,kind=pm_dbl) - dycom(k)
      dzdiff = real(dz,kind=pm_dbl) - dzcom(k)

      rogsqtot = rogsqtot + dxdiff*dxdiff + dydiff*dydiff + dzdiff*dzdiff

      xnew = xnp(atemp, x)
      ynew = ynp(atemp, y)
      znew = znp(atemp, z)

      x = xnew
      y = ynew
      z = znew

      atemp = a(x,y,z)
    end subroutine radii_of_gyration

    subroutine calc_single_bead_data
      call zero_bead

      IF (nw(k) > 1) THEN
         dxcom(k) = real(dxmv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)
         dycom(k) = real(dymv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)
         dzcom(k) = real(dzmv,kind=pm_dbl)/real(nw(k),kind=pm_dbl)

         parasquared = real(dx*dx,kind=pm_dbl)
         perpsquared = real(dy*dy + dz*dz,kind=pm_dbl)

         ete(k) = sqrt(real(parasquared,kind=pm_dbl) + real(perpsquared,kind=pm_dbl))
         etepara(k) = sqrt(real(parasquared,kind=pm_dbl))
         eteperp(k) = sqrt(real(perpsquared,kind=pm_dbl))

         thetax(k) = real(dx,kind=pm_dbl)/ete(k)
         thetay(k) = real(dy,kind=pm_dbl)/ete(k)
         thetaz(k) = real(dz,kind=pm_dbl)/ete(k)
         call box_average
         !Radii of gyration
         !The center of mass used here is the distance to the first bead
         dx = 0
         dy = 0
         dz = 0

         rogsqtot = dxcom(k)*dxcom(k) + dycom(k)*dycom(k) + dzcom(k)*dzcom(k)

         x = xx(k)
         y = yy(k)
         z = zz(k)

         atemp = a(x,y,z)

         DO WHILE (atemp /= 13)
            call radii_of_gyration
         ENDDO

         rogsq = rogsqtot/real(nw(k),kind=pm_dbl)
         rog(k) = sqrt(rogsq)
         rogtot = rogtot + rog(k)

         etebintot(bin(k)) = etebintot(bin(k)) + ete(k)
         eteparabintot(bin(k)) = eteparabintot(bin(k)) + etepara(k)
         eteperpbintot(bin(k)) = eteperpbintot(bin(k)) + eteperp(k)

         rogbintot(bin(k)) = rogbintot(bin(k)) + rog(k)

         thetaxbintot(bin(k)) = thetaxbintot(bin(k)) + (three*thetax(k)*thetax(k) - one)/two
         thetaybintot(bin(k)) = thetaybintot(bin(k)) + (three*thetay(k)*thetay(k) - one)/two
         thetazbintot(bin(k)) = thetazbintot(bin(k)) + (three*thetaz(k)*thetaz(k) - one)/two
      ENDIF
    end subroutine calc_single_bead_data

    subroutine calc_box_avg
      !Calcualte the box averages
      denom = real(nkt,kind=pm_dbl)

      etebox = eteboxtot/denom
      eteparabox = eteparatot/denom
      eteperpbox = eteperptot/denom

      rogbox = rogtot/denom

      thetaxbox = thetaxtot/denom
      thetaybox = thetaytot/denom
      thetazbox = thetaztot/denom
    end subroutine calc_box_avg

    subroutine bin_avg
      IF (bincount(x) == 0) THEN
         call zero_bincount_init

      ELSE
         etebin(x) = etebintot(x)/denom
         eteparabin(x) = eteparabintot(x)/denom
         eteperpbin(x) = eteperpbintot(x)/denom

         rogbin(x) = rogbintot(x)/denom

         thetaxbin(x) = thetaxbintot(x)/denom
         thetaybin(x) = thetaybintot(x)/denom
         thetazbin(x) = thetazbintot(x)/denom

         nk1avg(x) =  real(binnk1(x),kind=pm_dbl)
         nk2avg(x) =  real(binnk2(x),kind=pm_dbl)
         nk3avg(x) =  real(binnk3(x),kind=pm_dbl)
         nk4avg(x) =  real(binnk4(x),kind=pm_dbl)
         nk5avg(x) =  real(binnk5(x),kind=pm_dbl)
         nk6avg(x) =  real(binnk6(x),kind=pm_dbl)
         nk7avg(x) =  real(binnk7(x),kind=pm_dbl)
         nk8avg(x) =  real(binnk8(x),kind=pm_dbl)
      ENDIF
    end subroutine bin_avg

    subroutine cleanup
      !Array deallocation
      DEALLOCATE (dxcom, dycom, dzcom)

      DEALLOCATE (chainend, rog)

      DEALLOCATE (bin, bincount)
      DEALLOCATE (binnk1, binnk2, binnk3, binnk4, binnk5, binnk6, binnk7, binnk8)
      DEALLOCATE (nk1avg, nk2avg, nk3avg, nk4avg, nk5avg, nk6avg, nk7avg, nk8avg)

      DEALLOCATE (etepara, eteperp)
      DEALLOCATE (etebin, eteparabin, eteperpbin, rogbin)
      DEALLOCATE (etebintot, eteparabintot, eteperpbintot)

      DEALLOCATE (thetax, thetay, thetaz)
      DEALLOCATE (rogbintot, thetaxbintot, thetaybintot, thetazbintot)
      DEALLOCATE (thetaxbin, thetaybin, thetazbin)
    end subroutine cleanup

    subroutine maxsta_avg_calc
      READ(31) eteboxtemp, eteparaboxtemp, eteperpboxtemp, rogboxtemp, thetaxboxtemp, thetayboxtemp, thetazboxtemp
      eteboxsum = eteboxsum + eteboxtemp
      eteparaboxsum = eteparaboxsum + eteparaboxtemp
      eteperpboxsum = eteperpboxsum + eteperpboxtemp

      rogboxsum = rogboxsum + rogboxtemp

      thetaxsum = thetaxsum + thetaxboxtemp
      thetaysum = thetaysum + thetayboxtemp
      thetazsum = thetazsum + thetazboxtemp

      DO x = 1, nx
         READ(33) point, etebintemp, eteparabintemp, eteperpbintemp, rogbintemp, chainbintemp
         etebinsum(point) = etebinsum(point) + etebintemp
         eteparabinsum(point) = eteparabinsum(point) + eteparabintemp
         eteperpbinsum(point) = eteperpbinsum(point) + eteperpbintemp

         rogbinsum(point) = rogbinsum(point) + rogbintemp

         chainsum(point) = chainsum(point) + chainbintemp

         READ(34) point, nk1temp, nk2temp, nk3temp, nk4temp, nk5temp, nk6temp, nk7temp, nk8temp
         nk1sum(point) = nk1sum(point) + nk1temp
         nk2sum(point) = nk2sum(point) + nk2temp
         nk3sum(point) = nk3sum(point) + nk3temp
         nk4sum(point) = nk4sum(point) + nk4temp
         nk5sum(point) = nk5sum(point) + nk5temp
         nk6sum(point) = nk6sum(point) + nk6temp
         nk7sum(point) = nk7sum(point) + nk7temp
         nk8sum(point) = nk8sum(point) + nk8temp

         READ(35) point, thetaxbintemp, thetaybintemp, thetazbintemp
         thetaxbinsum(point) = thetaxbinsum(point) + thetaxbintemp
         thetaybinsum(point) = thetaybinsum(point) + thetaybintemp
         thetazbinsum(point) = thetazbinsum(point) + thetazbintemp
      end do
    end subroutine maxsta_avg_calc


  END SUBROUTINE chaincalcs

  subroutine open_chaincalcs_files(dir_name)
    implicit none
    character(len=*) :: dir_name

    dir_name = trim(adjustl(dir_name))
    !$OMP parallel sections
    !$OMP section
    !--- Formatted files --------
    open(unit = 40, file=dir_name//'etedynamicbox.dat', form = 'formatted')
    !$OMP section
    open(unit = 41, file=dir_name//'chainproperties.dat', form = 'formatted')
    !$OMP section
    open(unit = 42, file=dir_name//'etedynamicbin.dat', form = 'formatted')
    !$OMP section
    open(unit = 43, file=dir_name//'boxorder.dat', form = 'formatted')
    !$OMP section
    open(unit = 44, file=dir_name//'chainends.dat', form = 'formatted')
    !$OMP section
    open(unit = 45, file=dir_name//'binchainproperties.dat', form = 'formatted')
    !$OMP section
    open(unit = 46, file=dir_name//'nkcom.dat', form = 'formatted')
    !$OMP section
    open(unit = 47, file=dir_name//'binorder.dat', form = 'formatted')


    !--- unformatted files --------
    !$OMP section
    open(unit = 30, file=dir_name//'tmpetedynamicbox.tmp', form = 'unformatted')
    !$OMP section
    open(unit = 31, file=dir_name//'tmpchainproperties.tmp',form='unformatted')
    !$OMP section
    open(unit = 32, file=dir_name//'tmpetedynamicbin.tmp', form = 'unformatted')
    !$OMP section
    open(unit = 33, file=dir_name//'tmpetestaticbox.tmp', form = 'unformatted')
    !$OMP section
    open(unit = 34, file=dir_name//'tmpnkcombin.tmp', form = 'unformatted')
    !$OMP section
    open(unit = 35, file=dir_name//'tmporderbin.tmp', form = 'unformatted')
    !$OMP end parallel sections
    return
  end subroutine open_chaincalcs_files

  subroutine close_chaincalcs_files
    !$OMP parallel sections
    !--- Formatted files --------
    !$OMP section
    close(unit = 40)
    !$OMP section
    close(unit = 41)
    !$OMP section
    close(unit = 42)
    !$OMP section
    close(unit = 43)
    !$OMP section
    close(unit = 44)
    !$OMP section
    close(unit = 45)
    !$OMP section
    close(unit = 46)
    !$OMP section
    close(unit = 47)

    !--- Unformatted files --------
    !$OMP section
    close(unit = 30)
    !$OMP section
    close(unit = 31)
    !$OMP section
    close(unit = 32)
    !$OMP section
    close(unit = 33)
    !$OMP section
    close(unit = 34)
    !$OMP section
    close(unit = 35)
    !$OMP end parallel sections
    return
  end subroutine close_chaincalcs_files

  subroutine write_chaincalcs_headers
    !$OMP parallel sections
    !$OMP section
    write(40,fmt302) 'l', 'ete', 'etedev', 'rog', 'rogdev', 'maxete', 't'
    !$OMP section
    write(41,fmt303) 'l', 'ete', 'etedev', 'etepara', 'eteparadev', 'eteperp', 'eteperpdev', 'rog', 'rogdev'
    !$OMP section
    write(42,fmt301) 'x', 'etebin', 'etedev', 'rogbin', 'rogdev'
    !$OMP section
    write(43,fmt302) 'l', 'thetax', 'thetaxdev', 'thetay', 'thetaydev', 'thetaz', 'thetazdev'
    !$OMP section
    write(44,fmt300) 'x', 'chainend', 'chainenddev'
    !$OMP section
    write(45,fmt303) 'x', 'etebin', 'etedbinev', 'eteparabin', 'eteparabindev', 'eteperpbin', 'eteperpbindev', 'rogbin', 'rogbindev'
    !$OMP section
    write(46,fmt308) 'x', 'nk1', 'nk2', 'nk3', 'nk4', 'nk5', 'nk6', 'nk7', 'nk8'
    !$OMP section
    write(47,fmt302) 'l', 'thetaxbin', 'thetaxbindev', 'thetaybin', 'thetaybindev', 'thetazbin', 'thetazdevbin'
    !$OMP end parallel sections
    return
  end subroutine write_chaincalcs_headers

  subroutine allocate_chaincalcs_arrays
    use param
    !$OMP parallel sections
    !$OMP section
    allocate (dxcom(nkt), dycom(nkt), dzcom(nkt))
    !$OMP section
    allocate (etepara(nkt), eteperp(nkt))
    !$OMP section
    allocate (thetax(nkt), thetay(nkt), thetaz(nkt))
    !$OMP section
    allocate (chainend(nx), rog(nkt))
    !$OMP section
    allocate (bin(nkt), bincount(nx))
    !$OMP section
    allocate (binnk1(nx), binnk2(nx), binnk3(nx), binnk4(nx), binnk5(nx), binnk6(nx), binnk7(nx), binnk8(nx))
    !$OMP section
    allocate (etebintot(nx), eteparabintot(nx), eteperpbintot(nx))
    !$OMP section
    allocate (rogbintot(nx), thetaxbintot(nx), thetaybintot(nx), thetazbintot(nx))
    !$OMP section
    allocate (etebin(nx), eteparabin(nx), eteperpbin(nx), rogbin(nx))
    !$OMP section
    allocate (thetaxbin(nx), thetaybin(nx), thetazbin(nx))
    !$OMP section
    allocate (nk1avg(nx), nk2avg(nx), nk3avg(nx), nk4avg(nx), nk5avg(nx), nk6avg(nx), nk7avg(nx), nk8avg(nx))
    !$OMP end parallel sections
    return
  end subroutine allocate_chaincalcs_arrays

end module chaincalc
