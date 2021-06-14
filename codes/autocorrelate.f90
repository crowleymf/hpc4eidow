SUBROUTINE autocorrelate(dir_name)

USE param
      
IMPLICIT NONE
      
INTEGER x, y, z, loop
character(len=*) :: dir_name      
DOUBLE PRECISION autovelbin, autoetechain
DOUBLE PRECISION autoetetot, autoveltot
DOUBLE PRECISION autoete, autovel
DOUBLE PRECISION autoetetmp, autoveltmp
DOUBLE PRECISION autoetesum, autovelsum
DOUBLE PRECISION autoetefinal, autovelfinal
      
IF (l==Nequil) THEN
 ALLOCATE (eteold(nkt),vyold(nx))
!Old mean squared end-to-end vector
 DO k=1,nkt
  eteold(k)=ete(k)
 ENDDO
!Old velocity in y direction at each bin     
 DO x=1,nx
  vyold(x)=avgvy(x)
 ENDDO
          
 OPEN(unit=70, file=dir_name//'tmpautocorrelate.tmp', form='unformatted') 
  
 OPEN(unit=71, file=dir_name//'autocorrelatebox.dat', form='formatted')
 WRITE(71,*) 'l,t,autoete,autovel'
 CLOSE(unit=71)
ENDIF      
      
autoetetot=0.D0
DO k=1,nkt
 autoetechain=ete(k)*eteold(k)
 autoetetot=autoetetot+autoetechain
ENDDO
autoete=autoetetot/dble(nkt)
 
autoveltot = 0.D0
DO x=1,nx
 autovelbin=avgvy(x)*vyold(x)
 autoveltot=autoveltot+autovelbin
ENDDO
autovel=autoveltot/dble(nx)
      
IF (l==Nequil) THEN
 OPEN(unit=71, file=dir_name//'autocorrelatebox.dat', position='append')
 WRITE(71,*) l, t, autoete, autovel 
 CLOSE(unit=71)
ENDIF
      
WRITE(70) autoete,autovel 
 
IF (mod((l-Nequil),Nauto)==0 .AND. l/=Nequil) THEN
 autoetesum=0.D0
 autovelsum=0.D0
        
 REWIND(70)
        
 DO loop=1,Nauto
  READ(70) autoetetmp,autoveltmp
  autoetesum=autoetesum+autoetetmp
  autovelsum=autovelsum+autoveltmp
 ENDDO
        
 denom=dble(Nauto)
 autoetefinal=autoetesum/denom
 autovelfinal=autovelfinal/denom
  
 OPEN(unit=71, file=dir_name//'autocorrelatebox.dat', position='append')
 WRITE(71,*) l,t,autoetefinal,autovelfinal 
 CLOSE (unit=71)
        
 REWIND(70)  
         
ENDIF
 
END SUBROUTINE autocorrelate
