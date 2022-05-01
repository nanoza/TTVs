c-----------------------------------------------------------------------------
c   Integrates orbits and calculates transits from tstart to tend
c-----------------------------------------------------------------------------

	subroutine ntrans_rmvs3(mstar,rmin,rmax,npl,mpl,rpl,aflag,
     &	        apl,epl,ipl,opl,vpl,lpl,t0,tstart,tend,dt,
     &          xobs,yobs,zobs,xaux,yaux,zaux,ntrans,jtrans,
     &          nt_,tc_,bpar_,vsky_)

	include 'swift.inc'

        real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
        real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
	
	real*8 mass(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

	integer istat(NTPMAX,NSTAT),i1st,aflag
	integer nbod,ntp,nleft
        real*8 rstat(NTPMAX,NSTATR),r2,rmin,rmax
	real*8 t0,dt,t

c...    Added for transits
	integer npl,ialpha
	real*8 mstar,mpl(NPLMAX),apl(NPLMAX),epl(NPLMAX),ipl(NPLMAX)
        real*8 opl(NPLMAX),vpl(NPLMAX),lpl(NPLMAX),capm,omega
	real*8 rcrit(NPLMAX),rpl(NPLMAX)
	integer ntrans,i,j,jtpl(NPLMAX),jrtpl(NPLMAX)
	real*8 tstart,tend
	real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
	real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
     	real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
	real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)  
	real*8 xobs,yobs,zobs,xaux,yaux,zaux,ex,ey,ez
	integer it,nt_(NPLMAX)
        integer MT
	parameter (MT=10000)
	real*8 tc_(NPLMAX,MT),bpar_(NPLMAX,MT)
        real*8 vsky_(NPLMAX,MT)	
	real*8 sprod,xproj,yproj,xsky,ysky,vxsky,vysky
        real*8 tc,tcnew,tcold,bpar,vsky,zsky
	real*8 gauss,gmsun,dpi,gm
	integer k,jr,jrad(NPLMAX)
	real*8 xr,yr,zr,vxr,vyr,vzr,r2crit,r2critp
	real*8 q,period,dtmin,dtmax

	integer jtrans(NPLMAX)    
	real*8 sema(NPLMAX),ecc(NPLMAX),inc(NPLMAX)
	real*8 capom(NPLMAX),varpi(NPLMAX),lambda(NPLMAX)
        real*8 radius(NPLMAX),massr(NPLMAX)

c     Constants and initializations
         gauss = 0.01720209895d0
         gmsun = gauss*gauss	! units are Julian days and AU
         dpi = acos(-1.d0)

c       Initializations
	 ntp = 0
	 nleft = ntp
	 j2rp2 = 0.d0
	 j4rp4 = 0.d0
	 do i=1,ntrans
	    nt_(i) = 0                       
	 end do
	
c...  Construct X vector on the sky plane
	ex = yaux*zobs-zaux*yobs
	ey = zaux*xobs-xaux*zobs
	ez = xaux*yobs-yaux*xobs

c...    Shift indices
	nbod = npl+1
	do j=2,nbod
	   massr(j)=mpl(j-1)
	   sema(j)=apl(j-1)
	   ecc(j)=epl(j-1)
	   inc(j)=ipl(j-1)
	   capom(j)=opl(j-1)
	   varpi(j)=vpl(j-1)
	   lambda(j)=lpl(j-1)
	   radius(j)=rpl(j-1)
	end do
	do i=1,ntrans
	   jtpl(i)=jtrans(i)+1
	end do

c       Reorder planets in increasing radial distance (needed for Jacobi)
	do j=2,nbod 
	   jrad(j)=2
	   do k=2,nbod
	      if(sema(k).lt.sema(j)) jrad(j)=jrad(j)+1   ! radial index
	   end do
	end do

c       Transformation to Cartesian coordinates
	mass(1) = mstar*gmsun
	xh(1)=0.d0
	yh(1)=0.d0
	zh(1)=0.d0
	vxh(1)=0.d0
	vyh(1)=0.d0
	vzh(1)=0.d0
	do j=2,nbod
	   jr = jrad(j)
	   mass(jr) = massr(j)*mstar*gmsun    ! <-this was corrected...
	   gm =  mass(1)+mass(jr)
	   capm = lambda(j)-varpi(j)
	   capm = mod(capm,2.d0*dpi)
	   omega = varpi(j)-capom(j)
	   omega = mod(omega,2.d0*dpi)
	   if(ecc(j).lt.1) ialpha = -1
	   if(ecc(j).eq.1) ialpha = 0
	   if(ecc(j).gt.1) ialpha = 1	   
	   call orbel_el2xv(gm,ialpha,sema(j),ecc(j),inc(j),capom(j),
     &          omega,capm,xh(jr),yh(jr),zh(jr),vxh(jr),vyh(jr),vzh(jr))
           rcrit(jr) = radius(j)
	end do

c       Re-link transiting planets after re-ordering 
	do i=1,ntrans
	   jrtpl(i) = jrad(jtpl(i))
	end do

c	Step control
	dtmin = 0.1d0
	dtmax = dt
	do j=2,nbod
	   gm = mass(1)+mass(j)
	   period = 2.d0*dpi*sqrt(sema(j)**3/gm)
	   if(dtmax.gt.0.1d0*period) then 
	      dtmax = 0.1d0*period
	   end if
	end do	
	if(dt.gt.dtmax) dt=0.9d0*dtmax ! to avoid changing step too often
	if(dt.lt.dtmin) dt=dtmin       ! do not accept dt too small

c-----------------------Integration loops------------------------
	i1st = 0
	t = t0

! step backward even if tstart=t0 (to be sure that we do not miss 1st transit)
	if(tstart.le.t0) then  

c... Integrate backward
	do j=2,nbod
	   vxh(j)=-vxh(j)
	   vyh(j)=-vyh(j)
	   vzh(j)=-vzh(j)
	end do
	
	do while (t.le.-tstart)

c... Check on planet removal during the next timestep
	do j=2,nbod-1
	   do k=j+1,nbod
	      xr = xh(k)-xh(j)
	      yr = yh(k)-yh(j)
	      zr = zh(k)-zh(j)
	      vxr = vxh(k)-vxh(j)
	      vyr = vyh(k)-vyh(j)
	      vzr = vzh(k)-vzh(j)
	      r2crit = (rcrit(k)+rcrit(j))**2
	      r2critp = r2crit
	      aflag = 0
	      call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &                         r2crit,r2critp,aflag)
	      if(aflag.ne.0) return
	   end do
	end do

c       Check if planet too close or too far
	do j=2,nbod
	   r2 = xh(j)*xh(j)+yh(j)*yh(j)+zh(j)*zh(j)
	   if(r2.lt.rmin*rmin) then
	      aflag = 2
	      return
	   end if
	   if(r2.gt.rmax*rmax) then
	      aflag = 3
	      return
	   end if
	end do

	call rmvs3_step(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &	        xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &          vzht,istat,rstat,dt)
	   t = t + dt
	end do
 
	t = -(t-t0)
	do j=2,nbod
	   vxh(j)=-vxh(j)
	   vyh(j)=-vyh(j)
	   vzh(j)=-vzh(j)
	end do
	end if     ! for backward stepping

c       Integrate forward 
	do while (t.lt.tend)

c... Check on planet removal during the next timestep
	   do j=2,nbod-1
	      do k=j+1,nbod
		 xr = xh(k)-xh(j)
		 yr = yh(k)-yh(j)
		 zr = zh(k)-zh(j)
		 vxr = vxh(k)-vxh(j)
		 vyr = vyh(k)-vyh(j)
		 vzr = vzh(k)-vzh(j)
		 r2crit = (rcrit(k)+rcrit(j))**2
		 r2critp = r2crit
		 aflag = 0
		 call rmvs_chk_ind(xr,yr,zr,vxr,vyr,vzr,dt,
     &                         r2crit,r2critp,aflag)
		 if(aflag.ne.0) return
	      end do
	   end do

c       Check if planet too close or too far
	do j=2,nbod
	   r2 = xh(j)*xh(j)+yh(j)*yh(j)+zh(j)*zh(j)
	   if(r2.lt.rmin*rmin) then
	      aflag = 2
	      return
	   end if
	   if(r2.gt.rmax*rmax) then
	      aflag = 3
	      return
	   end if
	end do

c       Step forward
	   do i=1,ntrans
	      j = jrtpl(i)
	      xbeg(i) = xh(j)
	      ybeg(i) = yh(j)
	      zbeg(i) = zh(j)
	      vxbeg(i) = vxh(j)
	      vybeg(i) = vyh(j)
	      vzbeg(i) = vzh(j)
	   end do
	   
	   call rmvs3_step(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &	        xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &          vzht,istat,rstat,dt)
	   t = t + dt

	   do i=1,ntrans
	      j = jrtpl(i)
	      xend(i) = xh(j)
	      yend(i) = yh(j)
	      zend(i) = zh(j)
	      vxend(i) = vxh(j)
	      vyend(i) = vyh(j)
	      vzend(i) = vzh(j)
	   end do

c..     Check if transits of any planet occurred in the past time step
	   do i=1,ntrans
   
c...    Project beg position vector on the sky plane
	      xsky = xbeg(i)*ex+ybeg(i)*ey+zbeg(i)*ez
	      ysky = xbeg(i)*xaux+ybeg(i)*yaux+zbeg(i)*zaux
	      vxsky = vxbeg(i)*ex+vybeg(i)*ey+vzbeg(i)*ez
	      vysky = vxbeg(i)*xaux+vybeg(i)*yaux+vzbeg(i)*zaux
	      tcold = -(xsky*vxsky+ysky*vysky)/(vxsky*vxsky+vysky*vysky)

c...    Project end position vector on the sky plane
	      xsky = xend(i)*ex+yend(i)*ey+zend(i)*ez
	      ysky = xend(i)*xaux+yend(i)*yaux+zend(i)*zaux
	      vxsky = vxend(i)*ex+vyend(i)*ey+vzend(i)*ez
	      vysky = vxend(i)*xaux+vyend(i)*yaux+vzend(i)*zaux
c..     Estimate the time of conjunction 
	      tcnew = -(xsky*vxsky+ysky*vysky)/(vxsky*vxsky+vysky*vysky)

c..     Find if transit or occultation occurred in the last step
	      if(tcold.ge.0.d0.and.tcnew.lt.0.d0) then 

c..     If so, use the tc estimate that is smaller (avoids most problems for e~1)
	      if(abs(tcnew).lt.abs(tcold)) then
		 tc = tcnew
	      else
		 tc = tcold - dt
	      end if
	      if(tc.ge.-dt.and.tc.lt.0.d0) then

c..     Interpolate and iterate to the time of conjunction
		 call rmvs3_interp_transit(ntrans,xbeg,ybeg,zbeg,vxbeg,vybeg,
     &		      vzbeg,xend,yend,zend,vxend,vyend,vzend,mass,dt,
     &                xobs,yobs,zobs,xaux,yaux,zaux,ex,ey,ez,
     &                jrtpl,i,tc,bpar,vsky,zsky)

c..     Record transit into returned variables
		 if(zsky.gt.0.d0) then   ! yes transits, no occultations
		    nt_(i) = nt_(i) + 1
		    tc_(i,nt_(i)) = t + tc
		    bpar_(i,nt_(i)) = bpar
		    vsky_(i,nt_(i)) = vsky
		 end if
	      end if
	   end if 
	   end do
	enddo

	return
        end			! ntrans_rmvs3.f
c---------------------------------------------------------------------

c*************************************************************************
c                            RMVS3_INTERP_transit.F
c*************************************************************************
c This subroutine interpolates between two kepler orbits.
c For outer region only
c
c             Input:
c                 nbod                ==>  number of massive bodies 
c                                          (int scalar)
c                 xbeg,ybeg,zbeg      ==>  initial planet position in helio 
c                                            (real arrays)
c                 vxbeg,vybeg,vzbeg   ==>  initial planet vel in helio 
c                                            (real arrays)
c                 xend,yend,zend      ==>  final planet position in helio 
c                                            (real arrays)
c                 vxend,vyend,vzend   ==>  final planet position in helio 
c                                            (real arrays)
c                 dt                   ==>  time step (real sclar)
c                 msun                 ==>  mass of sun (real sclar)
c             Output:
c                 xtmp,ytmp,ztmp      ==>  position of planet wrt time 
c                                          for outer region
c                                            (2d real arrays)
c                 vxtmp,vytmp,vztmp   ==>  velocoty of planet wrt time 
c                                          for outer region
c                                            (2d real arrays)
c
c Remarks: Based on rmvs2_interp_o 
c Authors:  Hal Levison 
c Date:    7/10/96
c Last revision: 

      subroutine rmvs3_interp_transit(nbod,xbeg,ybeg,zbeg,vxbeg,vybeg,
     &     vzbeg,xend,yend,zend,vxend,vyend,vzend,mass,dt,
     &     xobs,yobs,zobs,xaux,yaux,zaux,ex,ey,ez,
     &     jtpl,id,tc,bpar,vsky,zsky)

      include 'swift.inc'
      include 'rmvs.inc'

c...  Inputs Only: 
      integer nbod,jtpl(NPLMAX),id
      real*8 dt
      real*8 mass(NPLMAX)
      real*8 xbeg(NPLMAX),ybeg(NPLMAX),zbeg(NPLMAX)
      real*8 vxbeg(NPLMAX),vybeg(NPLMAX),vzbeg(NPLMAX)
      real*8 xend(NPLMAX),yend(NPLMAX),zend(NPLMAX)
      real*8 vxend(NPLMAX),vyend(NPLMAX),vzend(NPLMAX)
      real*8 xaux,yaux,zaux,ex,ey,ez,xobs,yobs,zobs

c...  Input and output:
      real*8 tc     ! relative to the current step

c...  Output:
      real*8 bpar,vsky

c...  Internals
      integer i,iflg,iter,miter
      real*8 xc2,yc2,zc2
      real*8 vxc2,vyc2,vzc2
      real*8 xc1,yc1,zc1
      real*8 vxc1,vyc1,vzc1
      real*8 dt1,dt2,frac,onemf,gm
      real*8 xtr,ytr,ztr,vxtr,vytr,vztr
      real*8 xprod,xsky,ysky,zsky,xproj,yproj,vxsky,vysky
      real*8 tcorr,sprod 	

c     Interpolate to get transit parameters 
	i=id

c...  Loop till tcorr small enough or iter exceeds miter
	iter = 0
	miter = 5

 100	continue

	dt1 = dt + tc		! forward timestep 
	dt2 = tc		! backward timestep

	xc1 = xbeg(i)
	yc1 = ybeg(i)
	zc1 = zbeg(i)
	vxc1 = vxbeg(i)
	vyc1 = vybeg(i)
	vzc1 = vzbeg(i)
	   
	xc2 = xend(i)
	yc2 = yend(i)
	zc2 = zend(i)
	vxc2 = vxend(i)
	vyc2 = vyend(i)
	vzc2 = vzend(i)

c...Propagate forward
	gm=mass(1)+mass(jtpl(i))
	call drift_one(gm,xc1,yc1,zc1,vxc1,vyc1,vzc1,dt1,iflg)
	if(iflg.ne.0) then
c	   write(*,*)'iflg'
c       This is indication that a planet was lost
	endif

c... Propagate backward 
	call drift_one(gm,xc2,yc2,zc2,vxc2,vyc2,vzc2,dt2,iflg)
	if(iflg.ne.0) then
c	   write(*,*)'iflg'
c	This is indication that a planet was lost
	endif
 
c...    Interpolate
	frac = (dt+tc)/dt
        
	onemf = 1.0d0 - frac
	xtr = onemf*xc1 + frac*xc2
	ytr = onemf*yc1 + frac*yc2
	ztr = onemf*zc1 + frac*zc2
	vxtr = onemf*vxc1 + frac*vxc2
	vytr = onemf*vyc1 + frac*vyc2
	vztr = onemf*vzc1 + frac*vzc2

c...    Project position vector on the sky plane
	xsky = xtr*ex+ytr*ey+ztr*ez
	ysky = xtr*xaux+ytr*yaux+ztr*zaux
	vxsky = vxtr*ex+vytr*ey+vztr*ez
	vysky = vxtr*xaux+vytr*yaux+vztr*zaux
	      
c.. Estimate the time of conjunction relative to previous guess 
	tcorr = -(xsky*vxsky+ysky*vysky)/(vxsky*vxsky+vysky*vysky)
	tc = tc + tcorr

c... Go back if tcorr too big	   
	iter = iter+1
	if(abs(tcorr).gt.1.d0/86400.d0.and.iter.lt.miter) goto 100

c... Calculate 
        bpar = sqrt(xsky*xsky+ysky*ysky)
	vsky = sqrt(vxsky*vxsky+vysky*vysky)
	zsky = xobs*xtr+yobs*ytr+zobs*ztr

	return
	end			! rmvs3_interp.f

c-----------------------------------------------------------------------

c*************************************************************************
c                        HELIO_DRIFT.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in helio coord 
c                                    (real arrays)
c                 vxb,vyb,vzb   ==>  initial position in bary coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xh,yh,zh      ==>  final position in helio coord 
c                                       (real arrays)
c                 vxb,vyb,vzb   ==>  final position in bary coord 
c                                       (real arrays)
c
c Remarks:  Based on drift.f
c Authors:  Hal Levison 
c Date:    11/14/96
c Last revision: 1/8/97  for symba

      subroutine helio_drift(nbod,mass,xh,yh,zh,vxb,vyb,vzb,dt)	

      include 'swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt

c...  Inputs and Outputs:
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxb(nbod),vyb(nbod),vzb(nbod)

c...  Internals:
      integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

      do j = 2,nbod
         if(mass(j).ne.0.0d0) then
            call drift_one(mass(1),xh(j),yh(j),zh(j),
     &           vxb(j),vyb(j),vzb(j),dt,iflg)
            if(iflg.ne.0) then
c	       write(*,*)'iflg'
c           This is indication that a planet was lost
            endif
         endif
      enddo

      return
      end
c--------------------------------------------------------------------------

c*************************************************************************
c                        DRIFT.F
c*************************************************************************
c This subroutine loops thorugh the particles and calls the danby routine
c
c             Input:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xj,yj,zj      ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 vxj,vyj,vzj   ==>  initial position in jacobi coord 
c                                    (real arrays)
c                 dt            ==>  time step
c             Output:
c                 xj,yj,zj      ==>  final position in jacobi coord 
c                                       (real arrays)
c                 vxj,vyj,vzj   ==>  final position in jacobi coord 
c                                       (real arrays)
c
c Authors:  Hal Levison 
c Date:    2/12/93
c Last revision: 9/5/94

      subroutine drift(nbod,mass,xj,yj,zj,vxj,vyj,vzj,dt)	

      include 'swift.inc'

c...  Inputs Only: 
      integer nbod
      real*8 mass(nbod),dt

c...  Inputs and Outputs:
      real*8 xj(nbod),yj(nbod),zj(nbod)
      real*8 vxj(nbod),vyj(nbod),vzj(nbod)

c...  Internals:
	real*8 etajm1,etaj,mu
	integer j,iflg

c----
c...  Executable code 

c Take a drift forward dth

	etajm1 = mass(1)
	do j = 2,nbod
	   etaj = etajm1 + mass(j)
	   mu = mass(1)*etaj/etajm1
	   call drift_one(mu,xj(j),yj(j),zj(j),
     &             vxj(j),vyj(j),vzj(j),dt,iflg)
           if(iflg.ne.0) then
c	      write(*,*)'iflg'
c       This is an indication that a planet was lost
           endif
	   etajm1 = etaj
	enddo

	return
	end
c--------------------------------------------------------------------------

	subroutine fit1(nt,xx,yy,dyy,ttv)
	implicit NONE

	integer MT
	parameter (MT=10000)
	integer it,nt
	real*8 ttv(MT),xx(MT),yy(MT),dyy(MT)
	real*8 s,sx,sy,sxy,sxx,x,y
	real*8 a1,b1,a2,b2,siga,sigb,delta,sigy2,c2

c..     Fit ttv data by a very straight line
	s = 0.d0
	sx = 0.d0
	sy = 0.d0
	sxx = 0.d0
	sxy = 0.d0     

	do it=1,nt
	   x = xx(it)
	   y = yy(it)         
	   sigy2 = dyy(it)*dyy(it)
	   s = s + 1.d0/sigy2
	   sx = sx + x/sigy2
	   sy = sy + y/sigy2
	   sxx = sxx + x*x/sigy2
	   sxy = sxy + x*y/sigy2
	end do
	
	delta = S*Sxx - Sx*Sx 
	a1 = (Sxx*Sy - Sx*Sxy)/delta 
	b1 = (S*Sxy - Sx*Sy)/delta 
c	siga = sqrt(Sxx/delta) 
c	sigb = sqrt(S/delta)

	do it=1,nt
	   ttv(it) = yy(it)-(a1+b1*xx(it))
	enddo

!	write(*,*)a1,b1
!	stop


	return
	end

	subroutine fit2(nt,xx,yy,dyy,tdv)
	implicit NONE

	integer MT
	parameter (MT=10000)
	integer it,nt
	real*8 tdv(MT),xx(MT),yy(MT),dyy(MT)
	real*8 s,sy,y
	real*8 a1,b1,a2,b2,siga,sigb,delta,sigy2,c2

c..     Fit and subtract average from the tdv data 
	s = 0.d0
	sy = 0.d0
	do it=1,nt
	   y = yy(it)         
	   sigy2 = dyy(it)*dyy(it)
	   s = s + 1.d0/sigy2
	   sy = sy + y/sigy2
	end do 
	a1 = sy/s

	do it=1,nt
	   tdv(it) = yy(it)-a1
	enddo

	return
	end


***********************************************************************
c                    ORBEL_FGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_fget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
*           Cel. Mech. ".  Quartic convergence from Danby's book.
***********************************************************************

	real*8 function orbel_fget(e,capn)

      include 'swift.inc'

c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

c----
c...  Executable code 

c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. 

c  begin with a guess proposed by Danby	
	if( capn .lt. 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
	  call orbel_schget(x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
c	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) RETURN
	  x = orbel_fget
	enddo	

c	write(*,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end   ! orbel_fget
c------------------------------------------------------------------
