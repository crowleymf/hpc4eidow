!This subroutine is used to calculate the end-to-end vector, radius of gyration,
!order parameters and center of mass density
SUBROUTINE chaincalcs(dir_name)
USE param
use pmtypes

IMPLICIT NONE


character(len=*) :: dir_name
INTEGER x, y, z
INTEGER bincheck
INTEGER point, loop
real(kind=pm_dbl), ALLOCATABLE:: dxcom(:), dycom(:), dzcom(:)

real(kind=pm_dbl) parasquared, perpsquared
real(kind=pm_dbl), ALLOCATABLE:: etepara(:), eteperp(:)
real(kind=pm_dbl) eteboxtot, eteparatot, eteperptot
real(kind=pm_dbl) etebox, eteparabox, eteperpbox
real(kind=pm_dbl) eteboxtemp,eteparaboxtemp, eteperpboxtemp
real(kind=pm_dbl) eteboxsum, eteparaboxsum, eteperpboxsum
real(kind=pm_dbl) eteboxsumdev, eteparaboxsumdev, eteperpboxsumdev
real(kind=pm_dbl) etefinal, eteparafinal, eteperpfinal
real(kind=pm_dbl) etedev, eteparadev, eteperpdev
real(kind=pm_dbl) maxete

real(kind=pm_dbl), ALLOCATABLE:: etebintot(:), eteparabintot(:), eteperpbintot(:)
real(kind=pm_dbl), ALLOCATABLE:: etebin(:), eteparabin(:), eteperpbin(:)
real(kind=pm_dbl) etebintemp, eteparabintemp, eteperpbintemp
real(kind=pm_dbl), ALLOCATABLE:: etebinsum(:), eteparabinsum(:), eteperpbinsum(:)
real(kind=pm_dbl), ALLOCATABLE:: etebinsumdev(:),eteparabinsumdev(:), eteperpbinsumdev(:)
real(kind=pm_dbl) etebinfinal, eteparabinfinal,eteperpbinfinal
real(kind=pm_dbl) etebindev, eteparabindev, eteperpbindev

real(kind=pm_dbl) dxdiff, dydiff, dzdiff
real(kind=pm_dbl) rogsqtot,rogsq
real(kind=pm_dbl), ALLOCATABLE:: rog(:)
real(kind=pm_dbl) rogtot
real(kind=pm_dbl) rogbox
real(kind=pm_dbl) rogboxtemp
real(kind=pm_dbl) rogboxsum
real(kind=pm_dbl) rogboxsumdev
real(kind=pm_dbl) rogfinal
real(kind=pm_dbl) rogdev

real(kind=pm_dbl), ALLOCATABLE:: rogbintot(:)
real(kind=pm_dbl), ALLOCATABLE:: rogbin(:)
real(kind=pm_dbl) rogbintemp
real(kind=pm_dbl), ALLOCATABLE:: rogbinsum(:)
real(kind=pm_dbl), ALLOCATABLE:: rogbinsumdev(:)
real(kind=pm_dbl) rogbinfinal
real(kind=pm_dbl) rogbindev

real(kind=pm_dbl), ALLOCATABLE:: thetax(:), thetay(:), thetaz(:)
real(kind=pm_dbl) thetaxtot, thetaytot, thetaztot
real(kind=pm_dbl) thetaxbox, thetaybox, thetazbox
real(kind=pm_dbl) thetaxboxtemp, thetayboxtemp, thetazboxtemp
real(kind=pm_dbl) thetaxsum, thetaysum, thetazsum
real(kind=pm_dbl) thetaxsumdev, thetaysumdev, thetazsumdev
real(kind=pm_dbl) thetaxfinal, thetayfinal, thetazfinal
real(kind=pm_dbl) thetaxdev, thetaydev, thetazdev

real(kind=pm_dbl), ALLOCATABLE:: thetaxbintot(:), thetaybintot(:), thetazbintot(:)
real(kind=pm_dbl), ALLOCATABLE:: thetaxbin(:), thetaybin(:), thetazbin(:)
real(kind=pm_dbl) thetaxbintemp, thetaybintemp, thetazbintemp
real(kind=pm_dbl), ALLOCATABLE:: thetaxbinsum(:), thetaybinsum(:), thetazbinsum(:)
real(kind=pm_dbl), ALLOCATABLE:: thetaxbinsumdev(:), thetaybinsumdev(:), thetazbinsumdev(:)
real(kind=pm_dbl) thetaxbinfinal, thetaybinfinal, thetazbinfinal
real(kind=pm_dbl) thetaxbindev, thetaybindev, thetazbindev

INTEGER, ALLOCATABLE:: binnk1(:), binnk2(:), binnk3(:), binnk4(:), binnk5(:), binnk6(:), binnk7(:), binnk8(:)
real(kind=pm_dbl), ALLOCATABLE:: nk1avg(:), nk2avg(:), nk3avg(:), nk4avg(:), nk5avg(:), nk6avg(:), nk7avg(:), nk8avg(:)
real(kind=pm_dbl) nk1temp, nk2temp, nk3temp, nk4temp, nk5temp, nk6temp, nk7temp, nk8temp
real(kind=pm_dbl), ALLOCATABLE:: nk1sum(:), nk2sum(:), nk3sum(:), nk4sum(:), nk5sum(:), nk6sum(:), nk7sum(:), nk8sum(:)
real(kind=pm_dbl) nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

INTEGER, ALLOCATABLE:: chainend(:)
real(kind=pm_dbl)  chainbintemp
real(kind=pm_dbl), ALLOCATABLE:: chainsum(:)
real(kind=pm_dbl), ALLOCATABLE:: chainsumdev(:)
real(kind=pm_dbl) chainfinal
real(kind=pm_dbl) chaindev
dir_name = TRIM(ADJUSTL(dir_name))
!Initialization
IF (t == 0) THEN
 ALLOCATE (ete(nkt))

 OPEN(unit = 40, file=dir_name//'etedynamicbox.dat', form = 'formatted')
 WRITE(40,302) 'l', 'ete', 'etedev', 'rog', 'rogdev', 'maxete', 't'
 CLOSE(unit = 40)
 OPEN(unit = 41, file=dir_name//'chainproperties.dat', form = 'formatted')
 WRITE(41,303) 'l', 'ete', 'etedev', 'etepara', 'eteparadev', 'eteperp', 'eteperpdev', 'rog', 'rogdev'
 CLOSE(unit = 41)

 OPEN(unit = 42, file=dir_name//'etedynamicbin.dat', form = 'formatted')
 WRITE(42,301) 'x', 'etebin', 'etedev', 'rogbin', 'rogdev'
 CLOSE(unit = 42)

 OPEN(unit = 43, file=dir_name//'boxorder.dat', form = 'formatted')
 WRITE(43,302) 'l', 'thetax', 'thetaxdev', 'thetay', 'thetaydev', 'thetaz', 'thetazdev'
 CLOSE(unit = 43)

 OPEN(unit = 44, file=dir_name//'chainends.dat', form = 'formatted')
 WRITE(44,300) 'x', 'chainend', 'chainenddev'
 CLOSE(unit = 44)

 OPEN(unit = 45, file=dir_name//'binchainproperties.dat', form = 'formatted')
 WRITE(45,303) 'x', 'etebin', 'etedbinev', 'eteparabin', 'eteparabindev', 'eteperpbin', 'eteperpbindev', 'rogbin', 'rogbindev'
 CLOSE(unit = 45)

 OPEN(unit = 46, file=dir_name//'nkcom.dat', form = 'formatted')
 WRITE(46,308) 'x', 'nk1', 'nk2', 'nk3', 'nk4', 'nk5', 'nk6', 'nk7', 'nk8'
 CLOSE(unit = 46)

 OPEN(unit = 47, file=dir_name//'binorder.dat', form = 'formatted')
 WRITE(47,302) 'l', 'thetaxbin', 'thetaxbindev', 'thetaybin', 'thetaybindev', 'thetazbin', 'thetazdevbin'
 CLOSE(unit = 47)

 OPEN(unit = 30, file=dir_name//'tmpetedynamicbox.tmp', form = 'unformatted')
 OPEN(unit = 31, file=dir_name//'tmpchainproperties.tmp', form = 'unformatted')
 OPEN(unit = 32, file=dir_name//'tmpetedynamicbin.tmp', form = 'unformatted')
 OPEN(unit = 33, file=dir_name//'tmpetestaticbox.tmp', form = 'unformatted')
 OPEN(unit = 34, file=dir_name//'tmpnkcombin.tmp', form = 'unformatted')
 OPEN(unit = 35, file=dir_name//'tmporderbin.tmp', form = 'unformatted')

ENDIF
!Array allocation

ALLOCATE (dxcom(nkt), dycom(nkt), dzcom(nkt))
ALLOCATE (etepara(nkt), eteperp(nkt))
ALLOCATE (thetax(nkt), thetay(nkt), thetaz(nkt))
ALLOCATE (chainend(nx), rog(nkt))

ALLOCATE (bin(nkt), bincount(nx))
ALLOCATE (binnk1(nx), binnk2(nx), binnk3(nx), binnk4(nx), binnk5(nx), binnk6(nx), binnk7(nx), binnk8(nx))

ALLOCATE (etebintot(nx), eteparabintot(nx), eteperpbintot(nx))
ALLOCATE (rogbintot(nx), thetaxbintot(nx), thetaybintot(nx), thetazbintot(nx))

ALLOCATE (etebin(nx), eteparabin(nx), eteperpbin(nx), rogbin(nx))
ALLOCATE (thetaxbin(nx), thetaybin(nx), thetazbin(nx))
ALLOCATE (nk1avg(nx), nk2avg(nx), nk3avg(nx), nk4avg(nx), nk5avg(nx), nk6avg(nx), nk7avg(nx), nk8avg(nx))

call initalize

DO k = 1, nkt
 call iteration_init
 call calc_single_bead_data


!Calculate the maximum ete
IF(k==1) THEN
 maxete=ete(1)
END IF

IF (ete(k)>maxete) THEN
 maxete=ete(k)
END IF

ENDDO
call calc_box_avg
!File unit 30 is the maxdyn loop and 31 is maxsta
WRITE(30) etebox, rogbox
WRITE(31) etebox, eteparabox, eteperpbox, rogbox, thetaxbox, thetaybox, thetazbox
!Calculate the binned averages
DO x = 1, nx
 denom = real(bincount(x))
 call bin_avg

!Dynamic (calculations made every maxdyn loop) temp file
  WRITE(32) x, etebin(x), rogbin(x)
!Static (calculations made every maxsta loops) temp files
  WRITE(33) x, etebin(x), eteparabin(x), eteperpbin(x), rogbin(x), real(chainend(x),kind=pm_dbl)
  WRITE(34) x, nk1avg(x), nk2avg(x), nk3avg(x), nk4avg(x), nk5avg(x), nk6avg(x), nk7avg(x), nk8avg(x)
  WRITE(35) x, thetaxbin(x), thetaybin(x), thetazbin(x)
ENDDO

!Averaging after maxdyn loops
IF (mod(l,maxdyn)==0 .AND. l /= 0) THEN
 eteboxsum = 0.D0
 rogboxsum = 0.D0

 eteboxsumdev=0.D0
 rogboxsumdev=0.D0

 ALLOCATE (etebinsum(nx), rogbinsum(nx))
 ALLOCATE (etebinsumdev(nx), rogbinsumdev(nx))

 DO x = 1, nx
  etebinsum(x) = 0.D0
  rogbinsum(x) = 0.D0

  etebinsumdev(x)=0.D0
  rogbinsumdev(x)=0.D0
 ENDDO

 REWIND (30)
 REWIND (32)

 DO loop = 1, maxdyn
  READ(30) eteboxtemp, rogboxtemp
  eteboxsum = eteboxsum + eteboxtemp
  rogboxsum = rogboxsum + rogboxtemp

  DO x = 1, nx
   READ(32) point, etebintemp, rogbintemp
   etebinsum(point) = etebinsum(point) + etebintemp
   rogbinsum(point) = rogbinsum(point) + rogbintemp
  ENDDO
 ENDDO

 denom = real(maxdyn)

 REWIND(30)
 DO loop=1,maxdyn
  READ(30) eteboxtemp, rogboxtemp
  eteboxsumdev = eteboxsumdev + (eteboxtemp-eteboxsum/denom)**2
  rogboxsumdev = rogboxsumdev + (rogboxtemp-rogboxsum/denom)**2
 END DO
 etefinal = eteboxsum/denom
 etedev = sqrt(eteboxsumdev/(denom-1.D0))

 rogfinal = rogboxsum/denom
 rogdev = sqrt(rogboxsumdev/(denom-1.D0))

 OPEN(unit = 40, file=dir_name//'etedynamicbox.dat', position = 'append')
 WRITE(40,306) l, etefinal, etedev, rogfinal, rogdev, maxete, t
 CLOSE(unit = 40)


 REWIND(32)
 OPEN(unit = 42, file=dir_name//'etedynamicbin.dat', position = 'append')
 WRITE(42,*) 'l=', l

 DO loop=1,maxdyn
  DO x=1,nx
   READ(32) point, etebintemp, rogbintemp
   etebinsumdev(point)=etebinsumdev(point)+(etebintemp-etebinsum(x)/denom)**2
   rogbinsumdev(point)=rogbinsumdev(point)+(rogbintemp-rogbinsum(x)/denom)**2
  END DO
 END DO

 DO x=1,nx
  etebinfinal=etebinsum(x)/denom
  etebindev=sqrt(etebinsumdev(x)/(denom-1))
  rogbinfinal=rogbinsum(x)/denom
  rogbindev=sqrt(rogbinsumdev(x)/(denom-1))
  WRITE(42,305) x, etebinfinal, etebindev, rogbinfinal, rogbindev
 END DO
 CLOSE(unit = 42)

 REWIND(30)
 REWIND(32)

 DEALLOCATE (etebinsum, rogbinsum)
 DEALLOCATE (etebinsumdev, rogbinsumdev)
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

 etedev = sqrt(eteboxsumdev/(denom-1.D0))
 eteparadev=sqrt(eteparaboxsumdev/(denom-1.D0))
 eteperpdev=sqrt(eteperpboxsumdev/(denom-1.D0))

 rogfinal = rogboxsum/denom
 rogdev = sqrt(rogboxsumdev/(denom-1.D0))

 thetaxfinal = thetaxsum/denom
 thetayfinal = thetaysum/denom
 thetazfinal = thetazsum/denom

 thetaxdev = sqrt(thetaxsumdev/(denom-1.D0))
 thetaydev = sqrt(thetaysumdev/(denom-1.D0))
 thetazdev = sqrt(thetazsumdev/(denom-1.D0))

 OPEN(unit = 41, file=dir_name//'chainproperties.dat', position = 'append')
 WRITE(41,307) l, etefinal, etedev, eteparafinal, eteparadev, eteperpfinal, eteperpdev, rogfinal, rogdev
 CLOSE (unit = 41)

 OPEN(unit = 43, file=dir_name//'boxorder.dat', position = 'append')
 WRITE(43,306) l, thetaxfinal, thetaxdev, thetayfinal, thetaydev, thetazfinal, thetazdev
 CLOSE(unit = 43)

 REWIND(33)
 REWIND(35)

 OPEN(unit = 44, file=dir_name//'chainends.dat', position = 'append')
 WRITE(44,*) 'l=', l

 OPEN(unit = 45, file=dir_name//'binchainproperties.dat', position = 'append')
 WRITE(45,*) 'l=', l

 OPEN(unit = 46, file=dir_name//'nkcom.dat', position = 'append')
 WRITE(46,*) 'l=', l

 OPEN(unit = 47, file=dir_name//'binorder.dat', position = 'append')
 WRITE(47,*) 'l=', l

 DO loop=1,maxsta
  DO x=1,nx
   READ(33) point, etebintemp, eteparabintemp, eteperpbintemp, rogbintemp, chainbintemp
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
  chaindev = sqrt(chainsumdev(x)/(denom-1.D0))
  WRITE(44,304) x, chainfinal, chaindev

  etebinfinal = etebinsum(x)/denom
  rogbinfinal = rogbinsum(x)/denom
  eteparabinfinal = eteparabinsum(x)/denom
  eteperpbinfinal = eteperpbinsum(x)/denom

  etebindev = sqrt(etebinsumdev(x)/(denom-1.D0))
  eteparabindev = sqrt(eteparabinsumdev(x)/(denom-1.D0))
  eteperpbindev = sqrt(eteperpbinsumdev(x)/(denom-1.D0))
  rogbindev = sqrt(rogbinsumdev(x)/(denom-1.D0))
  WRITE(45,307) x, etebinfinal, etebindev, eteparabinfinal, eteparabindev, eteperpbinfinal, eteperpbindev, rogbinfinal, rogbindev

  nk1 = nk1sum(x)/denom
  nk2 = nk2sum(x)/denom
  nk3 = nk3sum(x)/denom
  nk4 = nk4sum(x)/denom
  nk5 = nk5sum(x)/denom
  nk6 = nk6sum(x)/denom
  nk7 = nk7sum(x)/denom
  nk8 = nk8sum(x)/denom
  WRITE(46,308) x, nk1, nk2, nk3, nk4, nk5, nk6, nk7, nk8

  thetaxbinfinal = thetaxbinsum(x)/denom
  thetaybinfinal = thetaybinsum(x)/denom
  thetazbinfinal = thetazbinsum(x)/denom

  thetaxbindev = sqrt(thetaxbinsumdev(x)/(denom-1.D0))
  thetaybindev = sqrt(thetaybinsumdev(x)/(denom-1.D0))
  thetazbindev = sqrt(thetazbinsumdev(x)/(denom-1.D0))
  WRITE(47,306) x, thetaxbinfinal, thetaxbindev, thetaybinfinal, thetaybindev, thetazbinfinal, thetazbindev

 ENDDO

 REWIND(31)
 REWIND(33)
 REWIND(34)
 REWIND(35)

 CLOSE (unit = 44)
 CLOSE (unit = 45)
 CLOSE (unit = 46)
 CLOSE (unit = 47)

 DEALLOCATE (etebinsum, eteparabinsum, eteperpbinsum, rogbinsum, chainsum)
 DEALLOCATE (etebinsumdev, eteparabinsumdev, eteperpbinsumdev, rogbinsumdev, chainsumdev)
 DEALLOCATE (nk1sum, nk2sum, nk3sum, nk4sum, nk5sum, nk6sum, nk7sum, nk8sum)
 DEALLOCATE (thetaxbinsum, thetaybinsum, thetazbinsum)
 DEALLOCATE (thetaxbinsumdev, thetaybinsumdev, thetazbinsumdev)

ENDIF

call cleanup

300 FORMAT (3(A20,x))
301 FORMAT (5(A20,x))
302 FORMAT (7(A20,x))
303 FORMAT (9(A20,x))
304 FORMAT (I20,x,2(f20.14,x))
305 FORMAT (I20,x,4(f20.14,x))
306 FORMAT (I20,x,6(f20.14,x))
307 FORMAT (I20,x,8(f20.14,x))
308 FORMAT (8(A20,x))

contains
  subroutine initalize
    eteboxtot = 0.D0
    eteparatot = 0.D0
    eteperptot = 0.D0

    rogtot = 0.D0

    thetaxtot = 0.D0
    thetaytot = 0.D0
    thetaztot = 0.D0

    DO x  = 1, nx
     chainend(x) = 0
     etebintot(x) = 0.D0
     eteparabintot(x) = 0.D0
     eteperpbintot(x) = 0.D0

     rogbintot(x) = 0.D0

     thetaxbintot(x) = 0.D0
     thetaybintot(x) = 0.D0
     thetazbintot(x) = 0.D0

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
    eteboxsum = 0.D0
    eteparaboxsum = 0.D0
    eteperpboxsum = 0.D0

    eteboxsumdev = 0.D0
    eteparaboxsumdev = 0.D0
    eteperpboxsumdev = 0.D0

    rogboxsum = 0.D0
    rogboxsumdev = 0.D0

    thetaxsum = 0.D0
    thetaysum = 0.D0
    thetazsum = 0.D0

    thetaxsumdev = 0.D0
    thetaysumdev = 0.D0
    thetazsumdev = 0.D0

    ALLOCATE (etebinsum(nx), eteparabinsum(nx), eteperpbinsum(nx), rogbinsum(nx), chainsum(nx))
    ALLOCATE (etebinsumdev(nx), eteparabinsumdev(nx), eteperpbinsumdev(nx),rogbinsumdev(nx),chainsumdev(nx))
    ALLOCATE (nk1sum(nx), nk2sum(nx), nk3sum(nx), nk4sum(nx), nk5sum(nx), nk6sum(nx), nk7sum(nx), nk8sum(nx))
    ALLOCATE (thetaxbinsum(nx), thetaybinsum(nx), thetazbinsum(nx))
    ALLOCATE (thetaxbinsumdev(nx), thetaybinsumdev(nx), thetazbinsumdev(nx))

    DO x = 1, nx
     etebinsum(x) = 0.D0
     eteparabinsum(x) = 0.D0
     eteperpbinsum(x) = 0.D0

     rogbinsum(x) = 0.D0

     chainsum(x) = 0.D0

     etebinsumdev(x) = 0.D0
     eteparabinsumdev(x) = 0.D0
     eteperpbinsumdev(x) = 0.D0

     rogbinsumdev(x) = 0.D0

     chainsumdev(x) = 0.D0

     nk1sum(x) = 0.D0
     nk2sum(x) = 0.D0
     nk3sum(x) = 0.D0
     nk4sum(x) = 0.D0
     nk5sum(x) = 0.D0
     nk6sum(x) = 0.D0
     nk7sum(x) = 0.D0
     nk8sum(x) = 0.D0

     thetaxbinsum(x) = 0.D0
     thetaybinsum(x) = 0.D0
     thetazbinsum(x) = 0.D0

     thetaxbinsumdev(x) = 0.D0
     thetaybinsumdev(x) = 0.D0
     thetazbinsumdev(x) = 0.D0
    ENDDO
  end subroutine zero_average_init

  subroutine zero_bincount_init
    etebin(x) = 0.D0
    eteparabin(x) = 0.D0
    eteperpbin(x) = 0.D0

    rogbin(x) = 0.D0

    thetaxbin(x) = 0.D0
    thetaybin(x) = 0.D0
    thetazbin(x) = 0.D0

    nk1avg(x) =  0.D0
    nk2avg(x) =  0.D0
    nk3avg(x) =  0.D0
    nk4avg(x) =  0.D0
    nk5avg(x) =  0.D0
    nk6avg(x) =  0.D0
    nk7avg(x) =  0.D0
    nk8avg(x) =  0.D0
  end subroutine

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
     dxcom(k) = 0.D0
     dycom(k) = 0.D0
     dzcom(k) = 0.D0

     ete(k) = 0.D0
     etepara(k) = 0.D0
     eteperp(k) = 0.D0

     thetax(k) = 0.D0
     thetay(k) = 0.D0
     thetaz(k) = 0.D0
  end subroutine zero_bead

  subroutine box_average
    !The box averages are number averages
    !The totals are divided by nkt at the end
      thetaxtot = thetaxtot + (3.D0*thetax(k)*thetax(k) - 1.0D0)/2.D0
      thetaytot = thetaytot + (3.D0*thetay(k)*thetay(k) - 1.0D0)/2.D0
      thetaztot = thetaztot + (3.D0*thetaz(k)*thetaz(k) - 1.0D0)/2.D0

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

      thetaxbintot(bin(k)) = thetaxbintot(bin(k)) + (3.D0*thetax(k)*thetax(k) - 1.0D0)/2.D0
      thetaybintot(bin(k)) = thetaybintot(bin(k)) + (3.D0*thetay(k)*thetay(k) - 1.0D0)/2.D0
      thetazbintot(bin(k)) = thetazbintot(bin(k)) + (3.D0*thetaz(k)*thetaz(k) - 1.0D0)/2.D0
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
