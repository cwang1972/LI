      PROGRAM CH_2nd_1
      IMPLICIT NONE
      
      INTEGER sx,sy,syp1
      PARAMETER(sx=512,sy=512,syp1=sy+1)
      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)
      
C general variables

      INTEGER  n,m,i,j,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa



C physical variables (true solution and force terms)  

      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy),  
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy)    

      REAL*8  eng,dev,height


C time variables

      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime
      

C aux matrix variables

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)
	   
	   REAL*8  rpsi((sx+1)*(sy+1)),
     *        reng(400000+1),rdev(400000+1),rhei(400000+1)
	 
      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen,bphi

      INTEGER lenwrec,len,nmd,ndv,ipos,isign




      CALL SETUP(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     &           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     &           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,     
     &           phi,phi1,ftmp,ftmp1,
     &           aux,aux1,aux2,bux,bux1,bux2,bphi)


      ctime=0.0D+0


CCCCCCCCCCCCCCCC
C MAIN TIME LOOP
CCCCCCCCCCCCCCCC

      DO kk=1,nts

         ctime=DBLE(kk-1)*dt
         IF (MOD(kk,5) .EQ. 0) THEN
            WRITE(6,*) 'step :',kk,' time ;',ctime
         END IF


         CALL Euler(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     *           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     *           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,
     *           phi,phi1,ftmp,ftmp1,
     *           aux,aux1,aux2,bux,bux1,bux2,eng,dev,height,bphi)


         ctime=ctime + dt
		 
      WRITE(6,*) 'step :',kk-1,'energy;',eng, 
     *            'deviation:',dev

      reng(kk)=eng
      rdev(kk)=dev
	  rhei(kk)=height    


      END DO
	  
      DO j=0,n-1
        phi(m,j)=phi(0,j)
      ENDDO

      DO i=0,m
        phi(i,n)=phi(i,0)
      ENDDO

      lenwrec=100
      len=(m+1)*(n+1)
      nmd=MOD(len,lenwrec)
      ndv=(len-nmd)/lenwrec

      ipos=1
         DO i=0,m
           DO j=0,n
             rpsi(ipos)=phi(i  ,j)
             ipos=ipos+1
           ENDDO
         END DO

      WRITE(6,*) 'WRITING phi.dat '  
      WRITE(10,*) 'WRITING phi.dat '  
      OPEN(unit=50,file='phi.dat',access='direct',
     *     form='unformatted',recl=8*lenwrec,status='unknown')
      ipos=1
      DO j=1,ndv
         WRITE(50,rec=j) (rpsi(ipos+i-1),i=1,lenwrec)
         ipos=ipos+lenwrec
      END DO
      WRITE(50,rec=ndv+1) (rpsi(ipos+i-1),i=1,nmd)
      CLOSE(unit=50) 

         CALL Euler(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     *           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     *           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,
     *           phi,phi1,ftmp,ftmp1,
     *           aux,aux1,aux2,bux,bux1,bux2,eng,dev,height,bphi)

      reng(nts+1)=eng
      rdev(nts+1)=dev
	  rhei(nts+1)=height

      lenwrec=100
      len=nts+1
      nmd=MOD(len,lenwrec)
      ndv=(len-nmd)/lenwrec


      WRITE(6,*) 'WRITING eng.dat '  
      WRITE(10,*) 'WRITING eng.dat '  
      OPEN(unit=50,file='eng.dat',access='direct',
     *     form='unformatted',recl=8*lenwrec,status='unknown')
      ipos=1
      DO j=1,ndv
         WRITE(50,rec=j) (reng(ipos+i-1),i=1,lenwrec)
         ipos=ipos+lenwrec
      END DO
      WRITE(50,rec=ndv+1) (reng(ipos+i-1),i=1,nmd)
      CLOSE(unit=50) 

      WRITE(6,*) 'WRITING dev.dat '  
      WRITE(10,*) 'WRITING dev.dat '  
      OPEN(unit=50,file='dev.dat',access='direct',
     *     form='unformatted',recl=8*lenwrec,status='unknown')
      ipos=1
      DO j=1,ndv
         WRITE(50,rec=j) (rdev(ipos+i-1),i=1,lenwrec)
         ipos=ipos+lenwrec
      END DO
      WRITE(50,rec=ndv+1) (rdev(ipos+i-1),i=1,nmd)
      CLOSE(unit=50) 
	  

      STOP
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     &           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     &           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,     
     &           phi,phi1,ftmp,ftmp1,
     &           aux,aux1,aux2,bux,bux1,bux2,bphi)
      IMPLICIT NONE

      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)


      INTEGER  n,m,sx,sy,i,j,k,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa

      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy),  
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy)
      
      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime
      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen,bphi
	  
      REAL*8  aaarandom,phix,phiy,phixx,phiyy,phixy,
     *        phixxxx,phixxyy,phiyyyy,
     *        aaa1,aaa2,aaa3,aaa4 


C local variables
      REAL*8  aap1,aap2,aap3,aap4,aap5,aap6,aap7,aap8,
     *        a1p,a2p,emaxp


      m=sx
      n=sy

      Re=1000.0D+0
      nu=1.0D+0/re

      eps=0.02D+0
      aaa=2.0D+0 

C compute Chebyshev differentiation matrices

      xlen=12.8D+0
      ylen=12.8D+0
      dx=xlen/FLOAT(m)
      dy=ylen/FLOAT(n)
   
      DO i=0,m
         x(i)=FLOAT(i)*dx 
      END DO

      DO j=0,n
         y(j)=FLOAT(j)*dy
      END DO



C initialize for FFT
c     CALL VSINQI(m  ,WSAVESQX)
c     CALL VSINQI(n  ,WSAVESQY)
c     CALL VSINQI(nz ,WSAVESQZ)
c     CALL VCOSQI(m  ,WSAVECQX)
c     CALL VCOSQI(n  ,WSAVECQY)
c     CALL VCOSQI(nz ,WSAVECQZ)
c     CALL VSINTI(m-1,WSAVESX )
c     CALL VSINTI(m-1,WSAVESY )
c     CALL VSINTI(nz-1,WSAVESZ)

      CALL VRFFTI(m  ,WSAVEX)
      CALL VRFFTI(n  ,WSAVEY)


      fmx1(0)=0.0D+0
      fmx2(0)=0.0D+0
      DO j=1,m-3,2
        fmx1(j  )=-FLOAT(j+1)*pi/xlen
        fmx1(j+1)=-fmx1(j)
        fmx2(j  )=-FLOAT(j+1)*FLOAT(j+1)*pi*pi/(xlen*xlen)
        fmx2(j+1)= fmx2(j)
      ENDDO
      fmx1(m-1)=-FLOAT(m)*pi/xlen
      fmx2(m-1)=-FLOAT(m)*FLOAT(m)*pi*pi/(xlen*xlen)

      fmy1(0)=0.0D+0
      fmy2(0)=0.0D+0
      DO j=1,n-3,2
        fmy1(j  )=-FLOAT(j+1)*pi/ylen
        fmy1(j+1)=-fmy1(j)
        fmy2(j  )=-FLOAT(j+1)*FLOAT(j+1)*pi*pi/(ylen*ylen)
        fmy2(j+1)= fmy2(j)
      ENDDO
      fmy1(n-1)=-FLOAT(n)*pi/ylen
      fmy2(n-1)=-FLOAT(n)*FLOAT(n)*pi*pi/(ylen*ylen)
      
      
C time stepping 
      cfl=0.5D+0
      dt=cfl*DMIN1(dx,dy)
	  dt=0.004D+0
      final_time=10.0D+0
      nts=IDINT(final_time/dt)
      tfinal=DBLE(nts)*dt

      WRITE(6,*) ' '
      WRITE(6,*) 'tfinal, nts, dt :',tfinal,nts,dt
      WRITE(6,*) ' '


C Initial data
      
      DO j=0,n
      DO i=0,m
	  
        CALL RANDOM_NUMBER(aaarandom)
        aaarandom = 0.1D+0*(2.0D+0*aaarandom - 1.0D+0)
        phi(i,j) = aaarandom
     
        phi1(i,j)=phi (i,j)	  	  
      ENDDO
      ENDDO
	  
      bphi=0.0D+0
      DO j=0,n-1
      DO i=0,m-1
        bphi=bphi + phi(i,j)
      ENDDO
      ENDDO
      bphi=bphi/(FLOAT(m)*FLOAT(n))


      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Euler(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     *           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     *           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,
     *           phi,phi1,ftmp,ftmp1,
     *           aux,aux1,aux2,bux,bux1,bux2,eng,dev,height,bphi)
      IMPLICIT NONE

      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)

      INTEGER  n,m,sx,sy,i,j,k,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa
      
      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy), 
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy)
    
      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)


      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen,bphi
      
      REAL*8  phix,phiy,phixx,phiyy,phixy,
     *        phixxxx,phixxyy,phiyyyy,
     *        aaa1,aaa2,aaa3,aaa4
	  
	   REAL*8  eng,eng2,dev,height 
   
 
C    Taking phi_x 
      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=phi(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (n,m,aux (0,0),aux2(0,0),n+1,WSAVEX)

      DO i=0,n-1

      aux1(i,0)=0.0D+0
      DO j=1,m-3,2
        aux1(i,j  )= aux(i,j+1)*fmx1(j)

        aux1(i,j+1)= aux(i,j  )*fmx1(j+1)
      ENDDO
      aux1(i,m-1)=0.0D+0
      
      ENDDO

      CALL VRFFTB (n,m,aux1(0,0),aux2(0,0),n+1,WSAVEX)

      t1=0.0D+0
      DO j=0,n-1
      DO i=0,m-1
        t1=t1 + aux1(j,i)*aux1(j,i)
      ENDDO
      ENDDO
	  t1=t1*0.5D+0*eps*eps*dx*dy
	  
C Taking phi_y 
      DO j=0,n-1
      DO i=0,m-1
        bux(i,j)=phi(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (m,n,bux (0,0),bux2(0,0),m+1,WSAVEY)

      DO i=0,m-1

      bux1(i,0)=0.0D+0
      DO j=1,n-3,2
        bux1(i,j  )= bux(i,j+1)*fmy1(j)

        bux1(i,j+1)= bux(i,j  )*fmy1(j+1)
      ENDDO
      bux1(i,n-1)=0.0D+0
      
      ENDDO

      CALL VRFFTB (m,n,bux1(0,0),bux2(0,0),m+1,WSAVEY)
	  
	  t2=0.0D+0
      DO j=0,n-1
      DO i=0,m-1
        t2=t2 + aux1(i,j)*aux1(i,j)
      ENDDO
      ENDDO
	  t1=t1 + t2*0.5D+0*eps*eps*dx*dy
	  
C Computation of energy
      
      eng2=t1

      t2=0.0D+0
      DO j=0,n-1
      DO i=0,m-1
	    t3=phi(i,j)*phi(i,j)
	  
        t2=t2 + 0.25D+0*t3*t3 - 0.5D+0*t3
      ENDDO
      ENDDO
      eng=dx*dy*t2 + eng2
CCCCCCCCCCCCCCC

C Computation of standard deviation 
      t1=0.0D+0
      DO j=0,n-1
      DO i=0,m-1
        t2=phi(i,j) - bphi

        t1=t1 + t2*t2
      ENDDO
      ENDDO
      dev=SQRT(dx*dy*t1/(xlen*ylen))
CCCCCCCCCCCCCCCCCCCCCC   
 
      
C Compute \Delta (-(2 A + 1.5) phi + (A +0.5) phi_1)

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=-(2.0D+0*aaa + 1.5D+0)*phi (i,j) 
     *          + (       aaa + 0.5D+0)*phi1(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        bux(i,j)=aux(j,i)
      ENDDO
      ENDDO

      CALL VRFFTF (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        t1=fmx2(i) + fmy2(j)

        bux(i,j)=bux(i,j)*t1
      ENDDO
      ENDDO

      CALL VRFFTB (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=bux(i,j)
      ENDDO
      ENDDO

      CALL VRFFTB (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        ftmp(i,j) = aux(j,i)
      ENDDO
      ENDDO


C Compute \Delta^2 phi_1

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=phi1(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        bux(i,j)=aux(j,i)
      ENDDO
      ENDDO

      CALL VRFFTF (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        t1=fmx2(i) + fmy2(j)

        bux(i,j)=bux(i,j)*t1*t1
      ENDDO
      ENDDO

      CALL VRFFTB (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=bux(i,j)
      ENDDO
      ENDDO

      CALL VRFFTB (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        ftmp(i,j) = ftmp(i,j) - 0.25D+0*eps*eps*aux(j,i)
      ENDDO
      ENDDO	  

C Increment term 

      DO j=0,n-1
      DO i=0,m-1
        ftmp(i,j)=1.0D+0/dt*phi(i,j) + ftmp(i,j)
      ENDDO
      ENDDO
      
C Linear iteration     
      DO j=0,n-1
      DO i=0,m-1
        ftmp1(i,j)=ftmp(i,j)
      ENDDO
      ENDDO
      
C Initial guess
      DO j=0,n-1
      DO i=0,m-1
        bux (i,j)= phi (i,j)
        bux1(i,j)= phi1(i,j)
        
        phi (i,j)= 2.0D+0*bux(i,j) - bux1(i,j)
        
        phi1(i,j)= bux (i,j)
      ENDDO
      ENDDO
      
      
      DO kk=1,6
	 

C Update the nonlinear terms

      DO j=0,n-1
      DO i=0,m-1
	    t1= phi(i,j)*phi(i,j) + phi1(i,j)*phi1(i,j)

       bux(i,j)= 0.25D+0*t1*(phi(i,j) + phi1(i,j)) 
      ENDDO
      ENDDO
	  
	  
C Compute \Delta bux

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=bux(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        bux(i,j)=aux(j,i)
      ENDDO
      ENDDO

      CALL VRFFTF (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        t1=fmx2(i) + fmy2(j)

        bux(i,j)=bux(i,j)*t1
      ENDDO
      ENDDO

      CALL VRFFTB (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=bux(i,j)
      ENDDO
      ENDDO

      CALL VRFFTB (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        ftmp(i,j) = ftmp1(i,j) + aux(j,i)
      ENDDO
      ENDDO
	  


C Solve for phi 

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=ftmp(i,j)
      ENDDO
      ENDDO

      CALL VRFFTF (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        bux(i,j)=aux(j,i)
      ENDDO
      ENDDO

      CALL VRFFTF (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        t1=fmx2(i) + fmy2(j)
        t3=1.0D+0/dt + 0.75D+0*eps*eps*t1*t1 - aaa*t1

        bux(i,j)=bux(i,j)/t3
      ENDDO
      ENDDO

      CALL VRFFTB (m,n,bux(0,0),bux2(0,0),m+1,WSAVEY)

      DO j=0,n-1
      DO i=0,m-1
        aux(j,i)=bux(i,j)
      ENDDO
      ENDDO

      CALL VRFFTB (n,m,aux(0,0),aux2(0,0),n+1,WSAVEX)

      DO j=0,n-1
      DO i=0,m-1
        phi(i,j)=aux(j,i)
      ENDDO
      ENDDO


      ENDDO


      RETURN
      END


