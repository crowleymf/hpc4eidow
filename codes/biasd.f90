!this subroutine is used to make the random direction moves in this monte carlo scheme.
subroutine old_biasd(x,d,dir_name)

  use pmtypes
  use param
  use pm_subs,only:ran2

  implicit none
  character(len=*) :: dir_name
  integer :: x,d
  !---real :: ran2
  double precision :: dtest
  double precision :: pxzx,ppyx,pmyx

  pxzx=pxz(x)
  ppyx=ppy(x)
  pmyx=pmy(x)

  iseed=10

  dtest=ran2(iseed)

  !move backward
  if (dtest<=pmyx) then
     d=int(ran2(iseed)*4)+5
     if (d==6) d=9
     if (d==8) d=12
     !move laterally
  else if ((dtest>pmyx).and.(dtest<=(pmyx+pxzx))) then
     d=int(ran2(iseed)*4)+1
     !move forward
  else if (dtest>(pmyx+pxzx)) then
     d=int(ran2(iseed)*4)+5
     if (d==5) d=10
     if (d==7) d=11
  end if

end subroutine old_biasd
