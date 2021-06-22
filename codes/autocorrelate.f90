SUBROUTINE autocorrelate(dir_name)
  use pmtypes
  USE param
  
  implicit none
  integer ::ix, loop
  character(len=*) :: dir_name      
  real(kind=pm_dbl) :: autovelbin, autoetechain
  real(kind=pm_dbl) :: autoetetot, autoveltot
  real(kind=pm_dbl) :: autoete, autovel
  real(kind=pm_dbl) :: autoetetmp, autoveltmp
  real(kind=pm_dbl) :: autoetesum, autovelsum
  real(kind=pm_dbl) :: autoetefinal, autovelfinal
  
  if (l==nequil) then
     allocate (eteold(nkt),vyold(nx))
     !Old mean squared end-to-end vector
     eteold(1:nkt)=ete(1:nkt)
     !Old velocity in y direction at each bin     
     vyold(1:nx)=avgvy(1:nx)
     
     open(unit=70, file=dir_name//'tmpautocorrelate.tmp', form='unformatted') 
     open(unit=71, file=dir_name//'autocorrelatebox.dat', form='formatted')
     write(71,'(i8,3e25.5)') 'l,t,autoete,autovel'
     close(unit=71)
  endif
  
  autoetetot=0.D0
  do k=1,nkt
     autoetechain=ete(k)*eteold(k)
     autoetetot=autoetetot+autoetechain
  enddo
  autoete=autoetetot/real(nkt,kind=pm_dbl)
  
  autoveltot = 0.d0
  do ix=1,nx
     autovelbin=avgvy(ix)*vyold(ix)
     autoveltot=autoveltot+autovelbin
  enddo
  autovel=autoveltot/real(nx,kind=pm_dbl)
  
  if (l==nequil) then
     open(unit=71, file=dir_name//'autocorrelatebox.dat', position='append')
     write(71,'(i8,3e25.5)') l, t, autoete, autovel 
     close(unit=71)
  endif
  
  write(70) autoete,autovel 
  
  if (mod((l-nequil),nauto)==0 .and. l/=nequil) then
     autoetesum=0.d0
     autovelsum=0.d0
     
     rewind(70)
     
     do loop=1,nauto
        read(70) autoetetmp,autoveltmp
        autoetesum=autoetesum+autoetetmp
        autovelsum=autovelsum+autoveltmp
     enddo
     
     denom=real(nauto,kind=pm_dbl)
     autoetefinal=autoetesum/denom
     autovelfinal=autovelfinal/denom
     
     open(unit=71, file=dir_name//'autocorrelatebox.dat', position='append')
     write(71,*) l,t,autoetefinal,autovelfinal 
     close (unit=71)
     
     rewind(70)  
     
  endif
  
end subroutine autocorrelate
