PROGRAM main

	use params
	use nestwrapper
      
	implicit none
	
	INTEGER :: i,j,it,it2
        REAL(8) :: tref,buf
        REAL(8),DIMENSION(MT) :: xx,yy,dyy,yyfit

      	!no parameters to wrap around
      	nest_pWrap = 0
        nest_pWrap(3) = 1   
        nest_pWrap(7) = 1 
        nest_pWrap(8) = 1
        nest_pWrap(10) = 1

        sflag = 1

        ! OBTAIN PRIORS
        !--------------------------------------------------------------------------
        ! Open up jammi_in.txt and read file
!        open(1,FILE='example_LRVplan/jammi_in.jam',FORM='FORMATTED',STATUS='UNKNOWN')
        open(1,FILE='priors.in',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,sdim
          read(1,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
          write(*,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
        END DO
        close(1)

        ! Mode 0 override
        DO j=1,sdim
          IF( Rflag(j) .EQ. 0 ) THEN
            Rmin(j) = Rsol(j) !- 1.0D-9
            Rmax(j) = Rsol(j) !+ 1.0D-9
          END IF
        END DO
        
        ! Write out
        DO j=1,sdim
        write(*,*) Rmin(j),Rmax(j)
        END DO
        !--------------------------------------------------------------------------


! Read TTV/TDV data
!        IF( priflag .EQ. 1 ) THEN
           ntrans = 1
           jtrans(1)=1
           nt(1) = 100
           tref = 54955.d0

           i = 1
           open(1,FILE='times.txt',STATUS='old') 
           do it=1,nt(i)
              read(1,*)cycl(i,it),tobs(i,it),dtobs(i,it)
              cycl(i,it)=cycl(i,it) + 1
             tobs(i,it)=tobs(i,it)-tref
!              tobs(i,it)=tobs(i,it)/60.d0/24.d0
!              dtobs(i,it)=dtobs(i,it)/60.d0/24.d0
           end do
           close(1)

!           open(1,FILE='TDVs.txt',STATUS='old') 
!           do it=1,nt(i)
!              read(1,*)cycl(i,it),tobs2(i,it),dtobs2(i,it)
!              cycl(i,it)=cycl(i,it) + 1
!           end do
!           close(1)

! Define the original variables
        if(sflag.eq.1.or.sflag.eq.3) then
           do it=1,nt(1)
              fp(it) = tobs(1,it)
              sigfp(it) = dtobs(1,it)
              fpwei(it) = fp(it)/(sigfp(it)**2)
              sigfpwei(it) = 1.0D0/(sigfp(it)**2)
              logsigfp(it) = DLOG(sigfp(it)**2)
           end do
        end if
        if(sflag.eq.3) then
           do it=1,nt(1)
              it2 = nt(1)+it
              fp(it2) = tobs2(1,it)
              sigfp(it2) = dtobs2(1,it)
              fpwei(it2) = fp(it2)/(sigfp(it2)**2)
              sigfpwei(it2) = 1.0D0/(sigfp(it2)**2)
              logsigfp(it2) = DLOG(sigfp(it2)**2)
           end do
        end if
        if(sflag.eq.2) then 
           do it=1,nt(1)
              fp(it) = tobs2(1,it)
              sigfp(it) = dtobs2(1,it)
              fpwei(it) = fp(it)/(sigfp(it)**2)
              sigfpwei(it) = 1.0D0/(sigfp(it)**2)
              logsigfp(it) = DLOG(sigfp(it)**2)
           end do
        end if

! Read planetary system template      
	 open(1,file='planet.in',status='old')
   	 read(1,*)mstar,rstar,rminr,rmaxr,rhostar
	 read(1,*)npl
	 do j=1,npl
	    read(1,*)mpl(j),rpl(j)
            read(1,*)apl(j),epl(j),ipl(j),opl(j),vpl(j),lpl(j)
	    ipl(j) = ipl(j)*pi/1.8d2
	    opl(j) = opl(j)*pi/1.8d2
	    vpl(j) = vpl(j)*pi/1.8d2
	    lpl(j) = lpl(j)*pi/1.8d2
	 end do
	 close(1)

!  rstar = mstar/(rhostar/1.408d0)
!  rstar = rstar**(1.d0/3.d0)

! Read transit parameters
	 open(1,file='time.in',status='old')
	 read(1,*)t0,tstart,tend,dt0
!	 read(8848,*)rmin,rmax
!	 read(8848,*)xobs,yobs,zobs 
!        read(8848,*)xaux,yaux,zaux
	 close(1)

        ! Read-in the seriesP.dat
!        IF( priflag .EQ. 1 ) THEN
!          write(*,*) 'Reading in primary data...'
!          open(unit=30301,file='example_LRVplan/seriesP.jam')
!          DO i = 1,nplen
!            read(30301,*) tp(i),fp(i),sigfp(i),epochp_seq(i)
!            fpwei(i) = fp(i)/(sigfp(i)**2)
!            sigfpwei(i) = 1.0D0/(sigfp(i)**2)
!            logsigfp(i) = DLOG(sigfp(i)**2)
!          END DO
!          close(30301)
!        END IF

        ! Read-in the seriesR.dat
        IF( rvflag .EQ. 1 ) THEN
          write(*,*) 'Reading in rv data...'
          open(unit=30302,file='seriesR.dat')
          DO i = 1,nrlen
            read(30302,*) tr(i),rv(i),sigrv(i),rvmode(i)
          END DO
          close(30302)
        END IF

      	call nest_Sample
END
