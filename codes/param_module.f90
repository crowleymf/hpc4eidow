module param

use pmtypes

!Box sizes
integer :: nx, ny, nz, nkt
!Different chain lengths
integer :: nw1, nw2, nw3, nw4, nw5, nw6, nw7, nw8
integer :: xa,ya,za,xb,yb,zb
!Length of the respective chains
integer,allocatable :: nw(:)
!Coordinates of the first monomer in each chain
integer,allocatable :: xx(:), yy(:), zz(:)

!Loop counters
integer :: i,j,k,ii,jj,kk,k0
integer :: l,e,lold
real(kind=pm_dbl) :: t,tu,tmk


!Total direction code number and direction code number without chain end or a single bead
integer,parameter :: code=13, code1=12

!Forward Direction Codes
integer :: px(code), py(code), pz(code)
!Forward bond code
integer,allocatable :: a(:,:,:)
!Reverse bond code
integer,allocatable :: b(:,:,:)
!Chain number of the respective number
integer,allocatable :: ket(:,:,:)
!New point coordinates in derection
integer,allocatable :: xnp(:,:), ynp(:,:), znp(:,:)

character(len=128) :: infile,outfile
integer :: maxloops,maxsta,maxdyn,Nequil,Nauto

integer :: iseed

!Variables for chaincalcs routine
integer :: dx,dy,dz,dxmv,dymv,dzmv
integer :: atemp, xnew,ynew, znew
integer,allocatable :: bin(:),bincount(:)

!Variables for controlling the chaindynamics call
integer(kind=4) :: tcount,value,tcountmk,valuemk

!Variables for dynamic calculation
integer,allocatable :: xd(:),yd(:),zd(:)

!Variables for vel routine
real(kind=pm_dbl),allocatable :: ppy(:),pmy(:),pxz(:)
integer,allocatable :: dispy(:,:,:),dispx(:,:,:)
real(kind=pm_dbl),allocatable :: vy(:,:,:),vx(:,:,:)
real(kind=pm_dbl),allocatable :: avgvy(:),avgvx(:)
real(kind=pm_dbl),allocatable :: sumvy(:),sumvx(:)

!Variables for autocorrelate routine
real(kind=pm_dbl),allocatable :: ete(:)
real(kind=pm_dbl),allocatable :: vyold(:),eteold(:)

!Variable for averaging 
real(kind=pm_dbl) :: denom

end module param
