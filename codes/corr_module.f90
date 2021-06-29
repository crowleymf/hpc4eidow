!Variables for chaindyanmics routine
module correlation
  use pmtypes
  integer,allocatable :: x0cm(:),y0cm(:),z0cm(:)
  real(kind=pm_dbl),allocatable :: x0com(:),y0com(:),z0com(:)
  real(kind=pm_dbl),allocatable :: x0nwf(:),y0nwf(:),z0nwf(:) 
  real(kind=pm_dbl),allocatable :: x0nwl(:),y0nwl(:),z0nwl(:)
end module correlation
