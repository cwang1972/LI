      PROGRAM AC_CH_2nd_1
      IMPLICIT NONE
      
      INTEGER sx,sy,syp1
      PARAMETER(sx=640,sy=640,syp1=sy+1)
      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)
      
C general variables

      INTEGER  n,m,i,j,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa



C physical variables (true solution and force terms)  

      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy),
c    *       px ( 0:sx, 0:sy),py ( 0:sx, 0:sy),
c    *       fx ( 0:sx, 0:sy),fy ( 0:sx, 0:sy),  
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy),
     *       phie(0:sx, 0:sy)

      REAL*8 fpe1(0:sx, 0:sy),fpe2(0:sx, 0:sy),
     *       fpe3(0:sx, 0:sy)
     

      REAL*8  eng


C time variables

      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime
      

C aux matrix variables

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)

      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen



      CALL SETUP(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     &           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     &           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,     
     &           phi,phi1,phie,ftmp,ftmp1,
     &           fpe1,fpe2,fpe3,
     &           aux,aux1,aux2,bux,bux1,bux2)


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
     *           phi,phi1,phie,fpe1,fpe2,fpe3,
     *           ftmp,ftmp1,
     *           aux,aux1,aux2,bux,bux1,bux2)


         ctime=ctime + dt

         CALL CHECK(n,m,sx,sy,x,y,re,nu,eps,ctime,dt,
     *           dx,dy,WSAVEX,WSAVEY,
     *           fmx1,fmx2,fmy1,fmy2,phi,phie)


      eng=0.0D+0
      DO j=0,n-1
      DO i=0,m-1 
        eng=eng + phi(i,j)*phi(i,j) 
      ENDDO
      ENDDO
      eng=eng*dx*dy

      WRITE(6,*) 'step :',kk,'energy;', eng
    


      END DO



         CALL CHECK(n,m,sx,sy,x,y,re,nu,eps,ctime,dt,
     *           dx,dy,WSAVEX,WSAVEY,
     *           fmx1,fmx2,fmy1,fmy2,phi,phie)


      STOP
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     &           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     &           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,     
     &           phi,phi1,phie,ftmp,ftmp1,
     &           fpe1,fpe2,fpe3,
     &           aux,aux1,aux2,bux,bux1,bux2)

      IMPLICIT NONE

      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)


      INTEGER  n,m,sx,sy,i,j,k,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa

      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy), 
c    *       px ( 0:sx, 0:sy),py ( 0:sx, 0:sy),
c    *       fx ( 0:sx, 0:sy),fy ( 0:sx, 0:sy),  
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy),
     *       phie(0:sx, 0:sy)

      REAL*8 fpe1(0:sx, 0:sy),fpe2(0:sx, 0:sy),
     *       fpe3(0:sx, 0:sy) 

      
      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime
      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen

      REAL*8  phix,phiy,phixx,phiyy,phixy,
     *        phixxxx,phixxyy,phiyyyy,
     *        aaa1,aaa2,aaa3,aaa4 


C local variables
      REAL*8  aap1,aap2,aap3,aap4,aap5,aap6,aap7,aap8,
     *        a1p,a2p,emaxp


      m=sx
      n=sy

      Re=1000.0D+0
      nu=1.0D+0/re

      eps=0.1D+0
      aaa=1.0D+0 

C compute Chebyshev differentiation matrices

      xlen=1.0D+0
      ylen=1.0D+0
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
      final_time=1.0D+0
      nts=IDINT(final_time/dt)
      tfinal=DBLE(nts)*dt

      WRITE(6,*) ' '
      WRITE(6,*) 'tfinal, nts, dt :',tfinal,nts,dt
      WRITE(6,*) ' '


C Initial data
      
      DO j=0,n
      DO i=0,m
        phie(i,j) =   SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phi (i,j)=phie(i,j)

        phix       =  2.0D+0*pi 
     *              * COS(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phiy       = -2.0D+0*pi 
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * SIN(2.0D+0*pi*y(j)    )

        phixx      = -4.0D+0*pi*pi
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phiyy      = -4.0D+0*pi*pi
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phixy      = -4.0D+0*pi*pi
     *              * COS(2.0D+0*pi*x(i)    ) 
     *              * SIN(2.0D+0*pi*y(j)    )

        phixxxx    =  16.0D+0*pi*pi*pi*pi
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phiyyyy    =  16.0D+0*pi*pi*pi*pi
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )

        phixxyy    =  16.0D+0*pi*pi*pi*pi
     *              * SIN(2.0D+0*pi*x(i)    ) 
     *              * COS(2.0D+0*pi*y(j)    )


        t1=phix*phix + phiy*phiy

        fpe1 (i,j)=phixx   + phiyy

        fpe2 (i,j)= 3.0D+0*phi(i,j)*phi(i,j) 
     *                    *fpe1(i,j) 
     *            + 6.0D+0*phi(i,j)*t1

        fpe3 (i,j)=phixxxx + phiyyyy
     *          + 2.0D+0*phixxyy

        phi1(i,j)=phie(i,j)*COS(-dt)	  	  
      ENDDO
      ENDDO


      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Euler(n,m,sx,sy,x,y,re,nu,eps,aaa,xlen,ylen,
     *           final_time,tfinal,ctime,cfl,dx,dy,dt,nts,
     *           WSAVEX,WSAVEY,fmx1,fmx2,fmy1,fmy2,
     *           phi,phi1,phie,fpe1,fpe2,fpe3,
     *           ftmp,ftmp1,
     *           aux,aux1,aux2,bux,bux1,bux2)

      IMPLICIT NONE

      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)

      INTEGER  n,m,sx,sy,i,j,k,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps,aaa
      
      REAL*8 phi( 0:sx, 0:sy),phi1(0:sx, 0:sy),
c    *       px ( 0:sx, 0:sy),py ( 0:sx, 0:sy),
c    *       fx ( 0:sx, 0:sy),fy ( 0:sx, 0:sy),  
     *       ftmp(0:sx, 0:sy),ftmp1(0:sx,0:sy),
     *       phie(0:sx, 0:sy)

      REAL*8 fpe1(0:sx, 0:sy),fpe2(0:sx, 0:sy),
     *       fpe3(0:sx, 0:sy)
     
      INTEGER nts
      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime

      REAL*8 aux(0:sy,0:sx),aux1(0:sy,0:sx),aux2(0:sy,0:sx), 
     *       bux(0:sx,0:sy),bux1(0:sx,0:sy),bux2(0:sx,0:sy)

      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)


      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen
      
      REAL*8  phix,phiy,phixx,phiyy,phixy,
     *        phixxxx,phixxyy,phiyyyy,
     *        aaa1,aaa2,aaa3,aaa4 
   
      
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



C Force terms 
      DO j=0,n-1
      DO i=0,m-1 
        t4= COS(ctime+0.5D+0*dt)
        t5=-SIN(ctime+0.5D+0*dt)

        ftmp(i,j) = ftmp(i,j) 
     *           +  fpe1(i,j)*t4
     *           -  fpe2(i,j)*t4*t4*t4
     *           +  eps*eps*fpe3(i,j)*t4
     *           + t5*phie(i,j)
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
      
      
      DO kk=1,8
	 

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


      SUBROUTINE CHECK(n,m,sx,sy,x,y,re,nu,eps,ctime,dt,
     *           dx,dy,WSAVEX,WSAVEY,
     *           fmx1,fmx2,fmy1,fmy2,phi,phie)
      IMPLICIT NONE

      REAL*8 pi
      PARAMETER(pi=3.14159265358979323846264338327950288D+0)

      INTEGER  n,m,sx,sy,i,j,k,kk
      REAL*8   x(0:sx),y(0:sy),nu,re,eps

      REAL*8 phi (0:sx,0:sy), 
     *       phie(0:sx,0:sy)


      REAL*8 final_time,tfinal,cfl,dt
      REAL*8 ctime
      REAL*8 WSAVEX(2*sx+15),WSAVEY(2*sy+15),
     *       fmx1(0:sx),fmx2(0:sx),fmy1(0:sy),fmy2(0:sy)

      REAL*8  dx,dy,t1,t2,t3,t4,t5,t6,xlen,ylen


C local variables
      REAL*8  aap1,aap2,aap3,aap4,aap5,aap6,aap7,aap8,
     *        a1p,a2p,emaxp

      write(6,*) 'The time is', ctime

      
      a1p=0.0D+0
      a2p=0.0D+0
      emaxp=0.0D+0

      do j=0,n-1
      do i=0,m-1
        aap1=abs(phi(i,j) - phie(i,j)*cos(ctime))

	a1p=a1p + dx*dy*aap1
        a2p=a2p + dx*dy*aap1*aap1
        
        if(emaxp.le.aap1) emaxp=aap1
      enddo
      enddo

      
      write(6,*) 'The max norm for error for the phase variable'
      write(6,*) emaxp, emaxp/dt**2
       
      write(6,*) 'The L1 norm for error for the phase variable'
      write(6,*) a1p, a1p/dt**2
  
      a2p=sqrt(a2p)
 
      write(6,*) 'The L2 norm for error for the phase variable'
      write(6,*) a2p, a2p/dt**2



      RETURN
      END


   
