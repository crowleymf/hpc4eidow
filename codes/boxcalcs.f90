module box_calcs
  use pmtypes
  use param
  implicit none

  integer :: x,y,z
  integer :: loop,length,point,boxcount,count
  real(kind=pm_dbl) :: rdsqrt,thetax,thetay,thetaz

  real(kind=pm_dbl), allocatable :: sxy(:),sxz(:),syz(:),sxx(:),syy(:),szz(:)
  integer, allocatable :: binnw1(:),binnw2(:),binnw3(:),binnw4(:), &
       binnw5(:),binnw6(:),binnw7(:),binnw8(:)
  real(kind=pm_dbl), allocatable :: &
       density1(:),density2(:),density3(:),density4(:), &
       density5(:),density6(:),density7(:),density8(:)
  real(kind=pm_dbl), allocatable :: &
       segdensity1(:),segdensity2(:),segdensity3(:),segdensity4(:)
  real(kind=pm_dbl), allocatable :: &
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
  real(kind=pm_dbl), allocatable :: sxysumbintmp(:),sxzsumbintmp(:),syzsumbintmp(:)
  real(kind=pm_dbl), allocatable :: sxxsumbintmp(:),syysumbintmp(:),szzsumbintmp(:)
  real(kind=pm_dbl), allocatable :: density1sum(:),density2sum(:),density3sum(:),density4sum(:)
  real(kind=pm_dbl), allocatable :: density5sum(:),density6sum(:),density7sum(:),density8sum(:)
  real(kind=pm_dbl), allocatable :: segdensity1sum(:),segdensity2sum(:),segdensity3sum(:),segdensity4sum(:)
  real(kind=pm_dbl), allocatable :: segdensity5sum(:),segdensity6sum(:),segdensity7sum(:),segdensity8sum(:)

  real(kind=pm_dbl) :: sxyfinal,sxzfinal,syzfinal
  real(kind=pm_dbl) :: sxxfinal,syyfinal,szzfinal
  real(kind=pm_dbl) :: sxybinfinal,sxzbinfinal,syzbinfinal
  real(kind=pm_dbl) :: sxxbinfinal,syybinfinal,szzbinfinal
  real(kind=pm_dbl) :: density1final,density2final,density3final,density4final
  real(kind=pm_dbl) :: density5final,density6final,density7final,density8final
  real(kind=pm_dbl) :: segdensity1final,segdensity2final,segdensity3final,segdensity4final
  real(kind=pm_dbl) :: segdensity5final,segdensity6final,segdensity7final,segdensity8final

  character(len=30) :: fmt400="(7(a20,x))"
  character(len=30) :: fmt401="(9(a20,x))"
  character(len=30) :: fmt402="(i20,x,6(f20.14,x))"
  character(len=30) :: fmt403="(i20,x,8(f20.14,x))"


contains
  !--------------------------------------------------------
  !    BOXCALCS
  !--------------------------------------------------------
  subroutine boxcalcs
    !--- the subroutine presented here is used to calculate the stress, 
    !---        lattice site density and segmental density.
    use pm_subs,only: dir_name
    implicit none
  
    allocate (sxy(nx),sxz(nx),syz(nx),sxx(nx),syy(nx),szz(nx))
    allocate (binnw1(nx),binnw2(nx),binnw3(nx),binnw4(nx),binnw5(nx),binnw6(nx),binnw7(nx),binnw8(nx))
    allocate (density1(nx),density2(nx),density3(nx),density4(nx),density5(nx),density6(nx),density7(nx),density8(nx))
    allocate (segdensity1(nx),segdensity2(nx),segdensity3(nx),segdensity4(nx),segdensity5(nx), &
         segdensity6(nx),segdensity7(nx),segdensity8(nx))

    if ( t == 0 ) then
       call open_boxcalcs_files(dir_name)
       call write_boxcalcs_col_header
    endif
    call init_setup

    xloop1: do x = 1, nx
       call init_nx_loop
       yloop1: do y = 1, ny
          if ( ((mod(x,2)==1) .and. (mod(y,2)==1)) .or. &
               ((mod(x,2)==0) .and. (mod(y,2)==0))) then
             do z = 1, nz, 2
                call calc_density
                call calc_stress
             enddo
          else
             do z = 2, nz, 2
                k = ket(x,y,z)
                atemp = a(x,y,z)
                length = nw(k)
                call density_shift
                call stress_shift
             enddo
          end if
       end do yloop1

       call normalize
       call lattice_density
       call segmental_density
       call write_density_data
    enddo xloop1

    denom = 1.d0 / real(nx, kind=pm_dbl)

    sxybox = sxytot*denom
    sxzbox = sxztot*denom
    syzbox = syztot*denom
    sxxbox = sxxtot*denom
    syybox = syytot*denom
    szzbox = szztot*denom

    write(51) sxybox, sxzbox, syzbox, sxxbox, syybox, szzbox

    writeblock: if (mod(l,maxsta) == 0 .and. l /= 0) then
       allocate (sxysumbintmp(nx),sxzsumbintmp(nx),syzsumbintmp(nx))
       allocate (sxxsumbintmp(nx),syysumbintmp(nx),szzsumbintmp(nx))
       allocate (density1sum(nx),density2sum(nx),density3sum(nx),density4sum(nx), &
            density5sum(nx),density6sum(nx),density7sum(nx),density8sum(nx))
       allocate (segdensity1sum(nx),segdensity2sum(nx),segdensity3sum(nx),segdensity4sum(nx))
       allocate (segdensity5sum(nx),segdensity6sum(nx),segdensity7sum(nx),segdensity8sum(nx))

       rewind(51)
       rewind(52)
       rewind(53)
       rewind(54)

       call zero_out

       read1: do loop = 1, maxsta
          read(51) sxytemp, sxztemp, syztemp, sxxtemp, syytemp, szztemp
          sxysumtemp = sxysumtemp + sxytemp
          sxzsumtemp = sxzsumtemp + sxztemp
          syzsumtemp = syzsumtemp + syztemp
          sxxsumtemp = sxxsumtemp + sxxtemp
          syysumtemp = syysumtemp + syytemp
          szzsumtemp = szzsumtemp + szztemp

          readx: do x = 1, nx

             read(52) point, sxybintmp, sxzbintmp, syzbintmp, sxxbintmp, syybintmp, szzbintmp
             sxysumbintmp(point) = sxysumbintmp(point) + sxybintmp
             sxzsumbintmp(point) = sxzsumbintmp(point) + sxzbintmp
             syzsumbintmp(point) = syzsumbintmp(point) + syzbintmp
             sxxsumbintmp(point) = sxxsumbintmp(point) + sxxbintmp
             syysumbintmp(point) = syysumbintmp(point) + syybintmp
             szzsumbintmp(point) = szzsumbintmp(point) + szzbintmp

             read(53) point, density1tmp, density2tmp, density3tmp, density4tmp, density5tmp, density6tmp, density7tmp, density8tmp
             density1sum(point) = density1sum(point) + density1tmp
             density2sum(point) = density2sum(point) + density2tmp
             density3sum(point) = density3sum(point) + density3tmp
             density4sum(point) = density4sum(point) + density4tmp
             density5sum(point) = density5sum(point) + density5tmp
             density6sum(point) = density6sum(point) + density6tmp
             density7sum(point) = density7sum(point) + density7tmp
             density8sum(point) = density8sum(point) + density8tmp

             read(54) point, segdensity1tmp, segdensity2tmp, segdensity3tmp, segdensity4tmp, &
                  segdensity5tmp, segdensity6tmp, segdensity7tmp, segdensity8tmp
             segdensity1sum(point) = segdensity1sum(point) + segdensity1tmp
             segdensity2sum(point) = segdensity2sum(point) + segdensity2tmp
             segdensity3sum(point) = segdensity3sum(point) + segdensity3tmp
             segdensity4sum(point) = segdensity4sum(point) + segdensity4tmp
             segdensity5sum(point) = segdensity5sum(point) + segdensity5tmp
             segdensity6sum(point) = segdensity6sum(point) + segdensity6tmp
             segdensity7sum(point) = segdensity7sum(point) + segdensity7tmp
             segdensity8sum(point) = segdensity8sum(point) + segdensity8tmp
          enddo readx
       enddo read1

       denom = 1.d0 / real(maxsta, kind=pm_dbl)

       sxyfinal = sxysumtemp*denom
       sxzfinal = sxzsumtemp*denom
       syzfinal = syzsumtemp*denom
       sxxfinal = sxxsumtemp*denom
       syyfinal = syysumtemp*denom
       szzfinal = szzsumtemp*denom

       rewind(51)
       rewind(52)
       rewind(53)
       rewind(54)

       write(55,fmt402) l, sxyfinal, sxzfinal, syzfinal, sxxfinal, syyfinal, szzfinal
       write(56,*) 'l=', l
       write(57,*) 'l=', l
       write(58,*) 'l=', l, 't=', t

       writex: do x = 1, nx

          sxybinfinal = sxysumbintmp(x)*denom
          sxzbinfinal = sxzsumbintmp(x)*denom
          syzbinfinal = syzsumbintmp(x)*denom
          sxxbinfinal = sxxsumbintmp(x)*denom
          syybinfinal = syysumbintmp(x)*denom
          szzbinfinal = szzsumbintmp(x)*denom
          write(56,fmt402) x, sxybinfinal, sxzbinfinal, syzbinfinal, &
                              sxxbinfinal, syybinfinal, szzbinfinal

          density1final = density1sum(x)*denom
          density2final = density2sum(x)*denom
          density3final = density3sum(x)*denom
          density4final = density4sum(x)*denom
          density5final = density5sum(x)*denom
          density6final = density6sum(x)*denom
          density7final = density7sum(x)*denom
          density8final = density8sum(x)*denom
          write(57,fmt403) x, density1final, density2final, density3final, &
               density4final, &
               density5final, density6final, density7final, density8final

          segdensity1final = segdensity1sum(x)*denom
          segdensity2final = segdensity2sum(x)*denom
          segdensity3final = segdensity3sum(x)*denom
          segdensity4final = segdensity4sum(x)*denom
          segdensity5final = segdensity5sum(x)*denom
          segdensity6final = segdensity6sum(x)*denom
          segdensity7final = segdensity7sum(x)*denom
          segdensity8final = segdensity8sum(x)*denom
          write(58,fmt403) x, segdensity1final, segdensity2final, &
               segdensity3final, &
               segdensity4final, segdensity5final, segdensity6final, &
               segdensity7final, segdensity8final

       end do writex

       deallocate (sxysumbintmp,sxzsumbintmp,syzsumbintmp)
       deallocate (sxxsumbintmp,syysumbintmp,szzsumbintmp)
       deallocate (density1sum,density2sum,density3sum,density4sum,density5sum,density6sum,density7sum,density8sum)
       deallocate (segdensity1sum,segdensity2sum,segdensity3sum,segdensity4sum, &
            segdensity5sum,segdensity6sum,segdensity7sum,segdensity8sum)

    end if writeblock

    deallocate (sxy,sxz,syz,sxx,syy,szz)
    deallocate (binnw1,binnw2,binnw3,binnw4,binnw5,binnw6,binnw7,binnw8)
    deallocate (density1,density2,density3,density4,density5,density6,density7,density8)
    deallocate (segdensity1,segdensity2,segdensity3,segdensity4,segdensity5,segdensity6,segdensity7,segdensity8)
    return

  contains
    subroutine write_density_data
      write(52) x, sxy(x),sxz(x),syz(x),sxx(x),syy(x),szz(x)
      write(53) x, density1(x),density2(x),density3(x),density4(x),density5(x),density6(x),density7(x),density8(x)
      write(54) x, segdensity1(x),segdensity2(x),segdensity3(x),segdensity4(x),segdensity5(x), &
           segdensity6(x),segdensity7(x),segdensity8(x)
      return
    end subroutine write_density_data

    subroutine init_setup
      rdsqrt = 1.d0/sqrt(2.d0)

      boxcount = 0

      sxytot = 0.d0
      sxztot = 0.d0
      syztot = 0.d0
      sxxtot = 0.d0
      syytot = 0.d0
      szztot = 0.d0
      return
    end subroutine init_setup

    subroutine init_nx_loop
      count = 0

      sxybin = 0.d0
      sxzbin = 0.d0
      syzbin = 0.d0
      sxxbin = 0.d0
      syybin = 0.d0
      szzbin = 0.d0

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
      !density calculations
      if (length == nw1) then
         binnw1(x) = binnw1(x) + 1
      elseif (length == nw2) then
         binnw2(x) = binnw2(x) + 1
      elseif (length == nw3) then
         binnw3(x) = binnw3(x) + 1
      elseif (length == nw4) then
         binnw4(x) = binnw4(x) + 1
      elseif (length == nw5) then
         binnw5(x) = binnw5(x) + 1
      elseif (length == nw6) then
         binnw6(x) = binnw6(x) + 1
      elseif (length == nw7) then
         binnw7(x) = binnw7(x) + 1
      elseif (length == nw8) then
         binnw8(x) = binnw8(x) + 1
      endif
    end subroutine calc_density

    subroutine calc_stress
      !stress calculations
      if(atemp /= 13) then
         dx = px(atemp)
         dy = py(atemp)
         dz = pz(atemp)

         thetax = real(dx, kind=pm_dbl)*rdsqrt
         thetay = real(dy, kind=pm_dbl)*rdsqrt
         thetaz = real(dz, kind=pm_dbl)*rdsqrt

         sxybin = sxybin + thetax*thetay
         sxzbin = sxzbin + thetax*thetaz
         syzbin = syzbin + thetay*thetaz
         sxxbin = sxxbin + thetax*thetax
         syybin = syybin + thetay*thetay
         szzbin = szzbin + thetaz*thetaz

         count = count + 1
         boxcount = boxcount + 1
      end if
    end subroutine calc_stress

    subroutine density_shift
      !density calculations
      if (length == nw1) then
         binnw1(x) = binnw1(x) + 1
      elseif (length == nw2) then
         binnw2(x) = binnw2(x) + 1
      elseif (length == nw3) then
         binnw3(x) = binnw3(x) + 1
      elseif (length == nw4) then
         binnw4(x) = binnw4(x) + 1
      elseif (length == nw5) then
         binnw5(x) = binnw5(x) + 1
      elseif (length == nw6) then
         binnw6(x) = binnw6(x) + 1
      elseif (length == nw7) then
         binnw7(x) = binnw7(x) + 1
      elseif (length == nw8) then
         binnw8(x) = binnw8(x) + 1
      endif
    end subroutine density_shift

    subroutine stress_shift
      !stress calculations
      if(atemp /= 13) then
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
      endif
    end subroutine stress_shift

    subroutine normalize
      denom = 1.d0/real(count, kind=pm_dbl)

      sxy(x) = sxybin*denom
      sxz(x) = sxzbin*denom
      syz(x) = syzbin*denom
      sxx(x) = sxxbin*denom
      syy(x) = syybin*denom
      szz(x) = szzbin*denom

      sxytot = sxytot + sxy(x)
      sxztot = sxztot + sxz(x)
      syztot = syztot + syz(x)
      sxxtot = sxxtot + sxx(x)
      syytot = syytot + syy(x)
      szztot = szztot + szz(x)
      return
    end subroutine normalize

    subroutine lattice_density
      denom = 2.d0 / real(ny*nz, kind=pm_dbl)
      !lattice density
      density1(x) = real(binnw1(x), kind=pm_dbl)*denom
      density2(x) = real(binnw2(x), kind=pm_dbl)*denom
      density3(x) = real(binnw3(x), kind=pm_dbl)*denom
      density4(x) = real(binnw4(x), kind=pm_dbl)*denom
      density5(x) = real(binnw5(x), kind=pm_dbl)*denom
      density6(x) = real(binnw6(x), kind=pm_dbl)*denom
      density7(x) = real(binnw7(x), kind=pm_dbl)*denom
      density8(x) = real(binnw8(x), kind=pm_dbl)*denom
      return
    end subroutine lattice_density

    subroutine segmental_density
      !segmental density
      segdensity1(x) = real(binnw1(x), kind=pm_dbl)
      segdensity2(x) = real(binnw2(x), kind=pm_dbl)
      segdensity3(x) = real(binnw3(x), kind=pm_dbl)
      segdensity4(x) = real(binnw4(x), kind=pm_dbl)
      segdensity5(x) = real(binnw5(x), kind=pm_dbl)
      segdensity6(x) = real(binnw6(x), kind=pm_dbl)
      segdensity7(x) = real(binnw7(x), kind=pm_dbl)
      segdensity8(x) = real(binnw8(x), kind=pm_dbl)
      return
    end subroutine segmental_density

    subroutine zero_out
      sxysumtemp = 0.d0
      sxzsumtemp = 0.d0
      syzsumtemp = 0.d0
      sxxsumtemp = 0.d0
      syysumtemp = 0.d0
      szzsumtemp = 0.d0

      do x = 1, nx
         sxysumbintmp(x) = 0.d0
         sxzsumbintmp(x) = 0.d0
         syzsumbintmp(x) = 0.d0
         sxxsumbintmp(x) = 0.d0
         syysumbintmp(x) = 0.d0
         szzsumbintmp(x) = 0.d0

         density1sum(x) = 0.d0
         density2sum(x) = 0.d0
         density3sum(x) = 0.d0
         density4sum(x) = 0.d0
         density5sum(x) = 0.d0
         density6sum(x) = 0.d0
         density7sum(x) = 0.d0
         density8sum(x) = 0.d0

         segdensity1sum(x) = 0.d0
         segdensity2sum(x) = 0.d0
         segdensity3sum(x) = 0.d0
         segdensity4sum(x) = 0.d0
         segdensity5sum(x) = 0.d0
         segdensity6sum(x) = 0.d0
         segdensity7sum(x) = 0.d0
         segdensity8sum(x) = 0.d0
      enddo
    end subroutine zero_out

  end subroutine boxcalcs

  subroutine open_boxcalcs_files(dir_name)
    character(len=*) :: dir_name
    open(unit=51, file=dir_name//'tmpstressbox.tmp', form='unformatted')
    open(unit=52, file=dir_name//'tmpstressbin.tmp', form='unformatted')
    open(unit=53, file=dir_name//'tmpdensity.tmp', form='unformatted')
    open(unit=54, file=dir_name//'tmpsegdensity.tmp',form='unformatted')
    
    open(unit=55, file=dir_name//'stressbox.dat', form='formatted')
    open(unit=56, file=dir_name//'stressbin.dat', form='formatted')
    open(unit=57, file=dir_name//'latticedensity.dat', form='formatted')
    open(unit=58, file=dir_name//'segmentdensity.dat', form='formatted')
    return  
  end subroutine open_boxcalcs_files

  subroutine write_boxcalcs_col_header
    write(55,fmt400) 'l', 'sxy', 'sxz', 'syz', 'sxx', 'syy', 'szz'
    write(56,fmt400) 'x', 'sxy', 'sxz', 'syz', 'sxx', 'syy', 'szz'
    write(57,fmt401) 'x', 'density1', 'density2', 'density3', 'density4', 'density5', 'density6', 'density7', 'density8'
    write(58,fmt401) 'x', 'segdensity2', 'segdensity2', 'segdensity3', &
         'segdensity4', 'segdensity5', 'segdensity6', 'segdensity7', 'segdensity8'
    return
  end subroutine write_boxcalcs_col_header


end module box_calcs

