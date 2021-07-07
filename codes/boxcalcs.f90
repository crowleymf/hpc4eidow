module box_calcs
  use pmtypes
  implicit none

contains
  SUBROUTINE boxcalcs(dir_name)
    !The subroutine presented here is used to calculate the stress, lattice site density and segmental density.

    USE param

    IMPLICIT NONE
    character(len=*) :: dir_name

    INTEGER :: x,y,z
    INTEGER :: loop,length,point,boxcount,count
    real(kind=pm_dbl) :: rdsqrt,thetax,thetay,thetaz

    real(kind=pm_dbl), ALLOCATABLE :: sxy(:),sxz(:),syz(:),sxx(:),syy(:),szz(:)
    INTEGER, ALLOCATABLE :: binnw1(:),binnw2(:),binnw3(:),binnw4(:), &
                            binnw5(:),binnw6(:),binnw7(:),binnw8(:)
    real(kind=pm_dbl), ALLOCATABLE :: &
                  density1(:),density2(:),density3(:),density4(:), &
                  density5(:),density6(:),density7(:),density8(:)
    real(kind=pm_dbl), ALLOCATABLE :: &
         segdensity1(:),segdensity2(:),segdensity3(:),segdensity4(:)
    real(kind=pm_dbl), ALLOCATABLE :: &
         segdensity5(:),segdensity6(:),segdensity7(:),segdensity8(:)

    real(kind=pm_dbl) :: sxytot,sxztot,syztot
    real(kind=pm_dbl) :: sxxtot,syytot,szztot
    real(kind=pm_dbl) :: sxybox,sxzbox,syzbox
    real(kind=pm_dbl) :: sxxbox,syybox,szzbox
    real(kind=pm_dbl) :: sxybin,sxzbin,syzbin
    real(kind=pm_dbl) :: sxxbin,syybin,szzbin

    real(kind=pm_dbl) :: sxytemp,sxztemp,syztemp
    real(kind=pm_dbl) :: sxxtemp,syytemp,szztemp
    real(kind=pm_dbl) :: sxybintmp,sxzbintmp,syzbintmp
    real(kind=pm_dbl) :: sxxbintmp,syybintmp,szzbintmp
    real(kind=pm_dbl) :: density1tmp,density2tmp,density3tmp,density4tmp,density5tmp,density6tmp,density7tmp,density8tmp
    real(kind=pm_dbl) :: segdensity1tmp,segdensity2tmp,segdensity3tmp,segdensity4tmp, &
         segdensity5tmp,segdensity6tmp,segdensity7tmp,segdensity8tmp

    real(kind=pm_dbl) :: sxysumtemp,sxzsumtemp,syzsumtemp
    real(kind=pm_dbl) :: sxxsumtemp,syysumtemp,szzsumtemp
    real(kind=pm_dbl), ALLOCATABLE :: sxysumbintmp(:),sxzsumbintmp(:),syzsumbintmp(:)
    real(kind=pm_dbl), ALLOCATABLE :: sxxsumbintmp(:),syysumbintmp(:),szzsumbintmp(:)
    real(kind=pm_dbl), ALLOCATABLE :: density1sum(:),density2sum(:),density3sum(:),density4sum(:)
    real(kind=pm_dbl), ALLOCATABLE :: density5sum(:),density6sum(:),density7sum(:),density8sum(:)
    real(kind=pm_dbl), ALLOCATABLE :: segdensity1sum(:),segdensity2sum(:),segdensity3sum(:),segdensity4sum(:)
    real(kind=pm_dbl), ALLOCATABLE :: segdensity5sum(:),segdensity6sum(:),segdensity7sum(:),segdensity8sum(:)

    real(kind=pm_dbl) :: sxyfinal,sxzfinal,syzfinal
    real(kind=pm_dbl) :: sxxfinal,syyfinal,szzfinal
    real(kind=pm_dbl) :: sxybinfinal,sxzbinfinal,syzbinfinal
    real(kind=pm_dbl) :: sxxbinfinal,syybinfinal,szzbinfinal
    real(kind=pm_dbl) :: density1final,density2final,density3final,density4final
    real(kind=pm_dbl) :: density5final,density6final,density7final,density8final
    real(kind=pm_dbl) :: segdensity1final,segdensity2final,segdensity3final,segdensity4final
    real(kind=pm_dbl) :: segdensity5final,segdensity6final,segdensity7final,segdensity8final

    ALLOCATE (sxy(nx),sxz(nx),syz(nx),sxx(nx),syy(nx),szz(nx))
    ALLOCATE (binnw1(nx),binnw2(nx),binnw3(nx),binnw4(nx),binnw5(nx),binnw6(nx),binnw7(nx),binnw8(nx))
    ALLOCATE (density1(nx),density2(nx),density3(nx),density4(nx),density5(nx),density6(nx),density7(nx),density8(nx))
    ALLOCATE (segdensity1(nx),segdensity2(nx),segdensity3(nx),segdensity4(nx),segdensity5(nx), &
         segdensity6(nx),segdensity7(nx),segdensity8(nx))

    call write_col_header

    call init_setup

    DO x = 1, nx
       call init_nx_loop

       DO y = 1, ny
          IF (((mod(x,2)==1).AND. (mod(y,2)==1)).OR.((mod(x,2)==0).AND. (mod(y,2)==0))) THEN
             DO z = 1, nz, 2
                call calc_density
                call calc_stress
             ENDDO
          ELSE
             DO z = 2, nz, 2

                k = ket(x,y,z)
                atemp = a(x,y,z)
                length = nw(k)
                call density_shift
                call stress_shift
             ENDDO
          END IF
       END DO

       call normalize

       call lattice_density
       call segmental_density
       call write_density_data
    ENDDO

    denom = real(nx, kind=pm_dbl)

    sxybox = sxytot/denom
    sxzbox = sxztot/denom
    syzbox = syztot/denom
    sxxbox = sxxtot/denom
    syybox = syytot/denom
    szzbox = szztot/denom

    WRITE(51) sxybox, sxzbox, syzbox, sxxbox, syybox, szzbox

    IF (mod(l,maxsta) == 0 .AND. l /= 0) THEN
       ALLOCATE (sxysumbintmp(nx),sxzsumbintmp(nx),syzsumbintmp(nx))
       ALLOCATE (sxxsumbintmp(nx),syysumbintmp(nx),szzsumbintmp(nx))
       ALLOCATE (density1sum(nx),density2sum(nx),density3sum(nx),density4sum(nx), &
            density5sum(nx),density6sum(nx),density7sum(nx),density8sum(nx))
       ALLOCATE (segdensity1sum(nx),segdensity2sum(nx),segdensity3sum(nx),segdensity4sum(nx))
       ALLOCATE (segdensity5sum(nx),segdensity6sum(nx),segdensity7sum(nx),segdensity8sum(nx))

       REWIND(51)
       REWIND(52)
       REWIND(53)
       REWIND(54)

       call zero_out

       DO loop = 1, maxsta

          READ(51) sxytemp, sxztemp, syztemp, sxxtemp, syytemp, szztemp
          sxysumtemp = sxysumtemp + sxytemp
          sxzsumtemp = sxzsumtemp + sxztemp
          syzsumtemp = syzsumtemp + syztemp
          sxxsumtemp = sxxsumtemp + sxxtemp
          syysumtemp = syysumtemp + syytemp
          szzsumtemp = szzsumtemp + szztemp

          DO x = 1, nx

             READ(52) point, sxybintmp, sxzbintmp, syzbintmp, sxxbintmp, syybintmp, szzbintmp
             sxysumbintmp(point) = sxysumbintmp(point) + sxybintmp
             sxzsumbintmp(point) = sxzsumbintmp(point) + sxzbintmp
             syzsumbintmp(point) = syzsumbintmp(point) + syzbintmp
             sxxsumbintmp(point) = sxxsumbintmp(point) + sxxbintmp
             syysumbintmp(point) = syysumbintmp(point) + syybintmp
             szzsumbintmp(point) = szzsumbintmp(point) + szzbintmp

             READ(53) point, density1tmp, density2tmp, density3tmp, density4tmp, density5tmp, density6tmp, density7tmp, density8tmp
             density1sum(point) = density1sum(point) + density1tmp
             density2sum(point) = density2sum(point) + density2tmp
             density3sum(point) = density3sum(point) + density3tmp
             density4sum(point) = density4sum(point) + density4tmp
             density5sum(point) = density5sum(point) + density5tmp
             density6sum(point) = density6sum(point) + density6tmp
             density7sum(point) = density7sum(point) + density7tmp
             density8sum(point) = density8sum(point) + density8tmp

             READ(54) point, segdensity1tmp, segdensity2tmp, segdensity3tmp, segdensity4tmp, &
                  segdensity5tmp, segdensity6tmp, segdensity7tmp, segdensity8tmp
             segdensity1sum(point) = segdensity1sum(point) + segdensity1tmp
             segdensity2sum(point) = segdensity2sum(point) + segdensity2tmp
             segdensity3sum(point) = segdensity3sum(point) + segdensity3tmp
             segdensity4sum(point) = segdensity4sum(point) + segdensity4tmp
             segdensity5sum(point) = segdensity5sum(point) + segdensity5tmp
             segdensity6sum(point) = segdensity6sum(point) + segdensity6tmp
             segdensity7sum(point) = segdensity7sum(point) + segdensity7tmp
             segdensity8sum(point) = segdensity8sum(point) + segdensity8tmp
          ENDDO
       ENDDO

       denom = real(maxsta, kind=pm_dbl)

       sxyfinal = sxysumtemp/denom
       sxzfinal = sxzsumtemp/denom
       syzfinal = syzsumtemp/denom
       sxxfinal = sxxsumtemp/denom
       syyfinal = syysumtemp/denom
       szzfinal = szzsumtemp/denom

       REWIND(51)
       REWIND(52)
       REWIND(53)
       REWIND(54)

       OPEN(unit=55, file=dir_name//'stressbox.dat', position='append')
       WRITE(55,402) l, sxyfinal, sxzfinal, syzfinal, sxxfinal, syyfinal, szzfinal
       CLOSE(unit=55)

       OPEN(unit=56, file=dir_name//'stressbin.dat', position='append')
       WRITE(56,*) 'l=', l

       OPEN(unit=57, file=dir_name//'latticedensity.dat', position='append')
       WRITE(57,*) 'l=', l

       OPEN(unit=58, file=dir_name//'segmentdensity.dat', position='append')
       WRITE(58,*) 'l=', l, 't=', t

       DO x = 1, nx

          sxybinfinal = sxysumbintmp(x)/denom
          sxzbinfinal = sxzsumbintmp(x)/denom
          syzbinfinal = syzsumbintmp(x)/denom
          sxxbinfinal = sxxsumbintmp(x)/denom
          syybinfinal = syysumbintmp(x)/denom
          szzbinfinal = szzsumbintmp(x)/denom
          WRITE(56,402) x, sxybinfinal, sxzbinfinal, syzbinfinal, sxxbinfinal, syybinfinal, szzbinfinal

          density1final = density1sum(x)/denom
          density2final = density2sum(x)/denom
          density3final = density3sum(x)/denom
          density4final = density4sum(x)/denom
          density5final = density5sum(x)/denom
          density6final = density6sum(x)/denom
          density7final = density7sum(x)/denom
          density8final = density8sum(x)/denom
          WRITE(57,403) x, density1final, density2final, density3final, density4final, &
               density5final, density6final, density7final, density8final

          segdensity1final = segdensity1sum(x)/denom
          segdensity2final = segdensity2sum(x)/denom
          segdensity3final = segdensity3sum(x)/denom
          segdensity4final = segdensity4sum(x)/denom
          segdensity5final = segdensity5sum(x)/denom
          segdensity6final = segdensity6sum(x)/denom
          segdensity7final = segdensity7sum(x)/denom
          segdensity8final = segdensity8sum(x)/denom
          WRITE(58,403) x, segdensity1final, segdensity2final, segdensity3final, &
               segdensity4final, segdensity5final, segdensity6final, &
               segdensity7final, segdensity8final

       END DO

       CLOSE(unit = 56)
       CLOSE(unit = 57)
       CLOSE(unit = 58)

       DEALLOCATE (sxysumbintmp,sxzsumbintmp,syzsumbintmp)
       DEALLOCATE (sxxsumbintmp,syysumbintmp,szzsumbintmp)
       DEALLOCATE (density1sum,density2sum,density3sum,density4sum,density5sum,density6sum,density7sum,density8sum)
       DEALLOCATE (segdensity1sum,segdensity2sum,segdensity3sum,segdensity4sum, &
            segdensity5sum,segdensity6sum,segdensity7sum,segdensity8sum)

    END IF

    DEALLOCATE (sxy,sxz,syz,sxx,syy,szz)
    DEALLOCATE (binnw1,binnw2,binnw3,binnw4,binnw5,binnw6,binnw7,binnw8)
    DEALLOCATE (density1,density2,density3,density4,density5,density6,density7,density8)
    DEALLOCATE (segdensity1,segdensity2,segdensity3,segdensity4,segdensity5,segdensity6,segdensity7,segdensity8)

400 FORMAT (7(A20,x))
401 FORMAT (9(A20,x))
402 FORMAT (I20,x,6(f20.14,x))
403 FORMAT (I20,x,8(f20.14,x))

  contains
    subroutine write_density_data
      WRITE(52) x, sxy(x),sxz(x),syz(x),sxx(x),syy(x),szz(x)
      WRITE(53) x, density1(x),density2(x),density3(x),density4(x),density5(x),density6(x),density7(x),density8(x)
      WRITE(54) x, segdensity1(x),segdensity2(x),segdensity3(x),segdensity4(x),segdensity5(x), &
           segdensity6(x),segdensity7(x),segdensity8(x)
    end subroutine write_density_data

    subroutine write_col_header
400   FORMAT (7(A20,x))
401   FORMAT (9(A20,x))
402   FORMAT (I20,x,6(f20.14,x))
403   FORMAT (I20,x,8(f20.14,x))

      IF (t==0) THEN
         OPEN(unit=51, file=dir_name//'tmpstressbox.tmp', form='unformatted')
         OPEN(unit=52, file=dir_name//'tmpstressbin.tmp', form='unformatted')
         OPEN(unit=53, file=dir_name//'tmpdensity.tmp', form='unformatted')
         OPEN(unit=54, file=dir_name//'tmpsegdensity.tmp',form='unformatted')

         OPEN(unit=55, file=dir_name//'stressbox.dat', form='formatted')
         WRITE(55,400) 'l', 'sxy', 'sxz', 'syz', 'sxx', 'syy', 'szz'
         CLOSE(unit=55)

         OPEN(unit=56, file=dir_name//'stressbin.dat', form='formatted')
         WRITE(56,400) 'x', 'sxy', 'sxz', 'syz', 'sxx', 'syy', 'szz'
         CLOSE(unit=56)

         OPEN(unit=57, file=dir_name//'latticedensity.dat', form='formatted')
         WRITE(57,401) 'x', 'density1', 'density2', 'density3', 'density4', 'density5', 'density6', 'density7', 'density8'
         CLOSE(unit=57)

         OPEN(unit=58, file=dir_name//'segmentdensity.dat', form='formatted')
         WRITE(58,401) 'x', 'segdensity2', 'segdensity2', 'segdensity3', &
              'segdensity4', 'segdensity5', 'segdensity6', 'segdensity7', 'segdensity8'
         CLOSE(unit=58)
      ENDIF
    end subroutine write_col_header
    subroutine init_setup
      rdsqrt = sqrt(2.D0)

      boxcount = 0

      sxytot = 0.D0
      sxztot = 0.D0
      syztot = 0.D0
      sxxtot = 0.D0
      syytot = 0.D0
      szztot = 0.D0
    end subroutine init_setup

    subroutine init_nx_loop
      count = 0

      sxybin = 0.D0
      sxzbin = 0.D0
      syzbin = 0.D0
      sxxbin = 0.D0
      syybin = 0.D0
      szzbin = 0.D0

      binnw1(x) = 0
      binnw2(x) = 0
      binnw3(x) = 0
      binnw4(x) = 0
      binnw5(x) = 0
      binnw6(x) = 0
      binnw7(x) = 0
      binnw8(x) = 0
    end subroutine init_nx_loop

    subroutine calc_density
      k = ket(x,y,z)
      atemp = a(x,y,z)
      length = nw(k)
      !Density calculations
      IF (length == nw1) THEN
         binnw1(x) = binnw1(x) + 1
      ELSEIF (length == nw2) THEN
         binnw2(x) = binnw2(x) + 1
      ELSEIF (length == nw3) THEN
         binnw3(x) = binnw3(x) + 1
      ELSEIF (length == nw4) THEN
         binnw4(x) = binnw4(x) + 1
      ELSEIF (length == nw5) THEN
         binnw5(x) = binnw5(x) + 1
      ELSEIF (length == nw6) THEN
         binnw6(x) = binnw6(x) + 1
      ELSEIF (length == nw7) THEN
         binnw7(x) = binnw7(x) + 1
      ELSEIF (length == nw8) THEN
         binnw8(x) = binnw8(x) + 1
      ENDIF
    end subroutine calc_density

    subroutine calc_stress
      !Stress Calculations
      IF(atemp /= 13) THEN
         dx = px(atemp)
         dy = py(atemp)
         dz = pz(atemp)

         thetax = real(dx, kind=pm_dbl)/rdsqrt
         thetay = real(dy, kind=pm_dbl)/rdsqrt
         thetaz = real(dz, kind=pm_dbl)/rdsqrt

         sxybin = sxybin + thetax*thetay
         sxzbin = sxzbin + thetax*thetaz
         syzbin = syzbin + thetay*thetaz
         sxxbin = sxxbin + thetax*thetax
         syybin = syybin + thetay*thetay
         szzbin = szzbin + thetaz*thetaz

         count = count + 1
         boxcount = boxcount + 1
      END IF
    end subroutine calc_stress

    subroutine density_shift
      !Density calculations
      IF (length == nw1) THEN
         binnw1(x) = binnw1(x) + 1
      ELSEIF (length == nw2) THEN
         binnw2(x) = binnw2(x) + 1
      ELSEIF (length == nw3) THEN
         binnw3(x) = binnw3(x) + 1
      ELSEIF (length == nw4) THEN
         binnw4(x) = binnw4(x) + 1
      ELSEIF (length == nw5) THEN
         binnw5(x) = binnw5(x) + 1
      ELSEIF (length == nw6) THEN
         binnw6(x) = binnw6(x) + 1
      ELSEIF (length == nw7) THEN
         binnw7(x) = binnw7(x) + 1
      ELSEIF (length == nw8) THEN
         binnw8(x) = binnw8(x) + 1
      ENDIF
    end subroutine density_shift

    subroutine stress_shift
      !Stress Calculations
      IF(atemp /= 13) THEN
         dx = px(atemp)
         dy = py(atemp)
         dz = pz(atemp)

         thetax = real(dx, kind=pm_dbl)/rdsqrt
         thetay = real(dy, kind=pm_dbl)/rdsqrt
         thetaz = real(dz, kind=pm_dbl)/rdsqrt

         sxybin = sxybin + thetax*thetay
         sxzbin = sxzbin + thetax*thetaz
         syzbin = syzbin + thetay*thetaz
         sxxbin = sxxbin + thetax*thetax
         syybin = syybin + thetay*thetay
         szzbin = szzbin + thetaz*thetaz

         count = count + 1
         boxcount = boxcount + 1
      ENDIF
    end subroutine stress_shift

    subroutine normalize
      denom = real(count, kind=pm_dbl)

      sxy(x) = sxybin/denom
      sxz(x) = sxzbin/denom
      syz(x) = syzbin/denom
      sxx(x) = sxxbin/denom
      syy(x) = syybin/denom
      szz(x) = szzbin/denom

      sxytot = sxytot + sxy(x)
      sxztot = sxztot + sxz(x)
      syztot = syztot + syz(x)
      sxxtot = sxxtot + sxx(x)
      syytot = syytot + syy(x)
      szztot = szztot + szz(x)

      denom = 0.5D0*real(ny*nz, kind=pm_dbl)
    end subroutine normalize

    subroutine lattice_density

      !Lattice density
      density1(x) = real(binnw1(x), kind=pm_dbl)/denom
      density2(x) = real(binnw2(x), kind=pm_dbl)/denom
      density3(x) = real(binnw3(x), kind=pm_dbl)/denom
      density4(x) = real(binnw4(x), kind=pm_dbl)/denom
      density5(x) = real(binnw5(x), kind=pm_dbl)/denom
      density6(x) = real(binnw6(x), kind=pm_dbl)/denom
      density7(x) = real(binnw7(x), kind=pm_dbl)/denom
      density8(x) = real(binnw8(x), kind=pm_dbl)/denom

    end subroutine lattice_density

    subroutine segmental_density
      !Segmental density
      segdensity1(x) = real(binnw1(x), kind=pm_dbl)
      segdensity2(x) = real(binnw2(x), kind=pm_dbl)
      segdensity3(x) = real(binnw3(x), kind=pm_dbl)
      segdensity4(x) = real(binnw4(x), kind=pm_dbl)
      segdensity5(x) = real(binnw5(x), kind=pm_dbl)
      segdensity6(x) = real(binnw6(x), kind=pm_dbl)
      segdensity7(x) = real(binnw7(x), kind=pm_dbl)
      segdensity8(x) = real(binnw8(x), kind=pm_dbl)
    end subroutine segmental_density

    subroutine zero_out
      sxysumtemp = 0.D0
      sxzsumtemp = 0.D0
      syzsumtemp = 0.D0
      sxxsumtemp = 0.D0
      syysumtemp = 0.D0
      szzsumtemp = 0.D0

      DO x = 1, nx
         sxysumbintmp(x) = 0.D0
         sxzsumbintmp(x) = 0.D0
         syzsumbintmp(x) = 0.D0
         sxxsumbintmp(x) = 0.D0
         syysumbintmp(x) = 0.D0
         szzsumbintmp(x) = 0.D0

         density1sum(x) = 0.D0
         density2sum(x) = 0.D0
         density3sum(x) = 0.D0
         density4sum(x) = 0.D0
         density5sum(x) = 0.D0
         density6sum(x) = 0.D0
         density7sum(x) = 0.D0
         density8sum(x) = 0.D0

         segdensity1sum(x) = 0.D0
         segdensity2sum(x) = 0.D0
         segdensity3sum(x) = 0.D0
         segdensity4sum(x) = 0.D0
         segdensity5sum(x) = 0.D0
         segdensity6sum(x) = 0.D0
         segdensity7sum(x) = 0.D0
         segdensity8sum(x) = 0.D0
      ENDDO
    end subroutine zero_out

  END SUBROUTINE boxcalcs
end module box_calcs

