MODULE like

use params
use transitmod

implicit none
      
contains
      
!=======================================================================

SUBROUTINE slikelihood(R,slhood)
         
	implicit none
      
	REAL(8), DIMENSION(nest_nPar) :: R
  REAL(8) :: slhood, loglike
	INTEGER :: i

        ! Scaling of the parameters from hypercube
	DO i = 1, sdim
          IF( Rflag(i) .EQ. 0 .OR. Rflag(i) .EQ. 2 ) THEN ! Modes 0 and 2
            ! Uniform: Rmax = max, Rmin = min
	    R(i) = Rmin(i) + (Rmax(i)-Rmin(i))*R(i)
          ELSE IF( Rflag(i) .EQ. 3 ) THEN ! Mode 3
            ! Gaussian: Rmax = mean, Rmin = stdev
            R(i) = Rmax(i) + roottwo*Rmin(i)*inverf(-1.0D0+2.0D0*R(i))
          ELSE IF( Rflag(i) .EQ. 4 ) THEN ! Mode 4
            ! Jeffrey's: Rmax = max, Rmin = min
            R(i) = ( Rmax(i)**R(i) )*( Rmin(i)**(1.0D0-R(i)) )
          ELSE IF( Rflag(i) .EQ. 5 ) THEN ! Mode 5
            ! Modified Jefrey's: Rmax = max, Rmin = inflection point
            R(i) = -( Rmin(i)**(1.0D0-R(i)) )*( Rmin(i)**R(i) - ( Rmin(i)+Rmax(i) )**R(i) )
          END IF
	END DO

        ! Call transit to get chi^2
        call models(R,loglike,0)

	slhood = loglike
 
!        write(*,*)'slhood=',slhood  
!        stop

        !!write(*,*) 'R(1) = ',R(1),' yielding chi2 = ',chi2

END SUBROUTINE slikelihood
      
!=======================================================================

!! ======================================================================
SUBROUTINE models(Rvec,loglike,showpri)

!implicit none

 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 INTEGER :: showpri, showrv
! REAL(8), DIMENSION(nplen) :: resP
 REAL(8), DIMENSION(MT) :: resP
 REAL(8) :: loglike, loglikeP, loglikeR
 INTEGER :: i, j
 REAL(8) :: omegab, eccb, logPsib, mu, sigmaneg, sigmapos

 ! === Call transit to primary transit ===

!call transit(Rvec,nplen,resP,&
!              showpri,tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,&
!              NresamP,integ_bigP,&
!              .FALSE.)

 call transit(Rvec,resP)

 loglikeP = 0.0D0
 if(sflag.eq.1.or.sflag.eq.2) then
    DO i=1,nt(1)
       loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
    END DO
    loglikeP = -0.5D0*nt(1)*LogTwoPi - 0.5D0*loglikeP
 end if
 if(sflag.eq.3) then
    DO i=1,2*nt(1)
       loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
    END DO
    loglikeP = -0.5D0*2.0d0*nt(1)*LogTwoPi - 0.5D0*loglikeP
 end if
 
 ! Add on photoeccentric penalty
 ! I *think* that omega_b = varpi_b - Omega_b and that Omega_b = 1.5pi by definition
 omegab = (Rvec(7) - 270.0D0)*radian
 eccb = Rvec(5)
 logPsib = DLOG( ( (1.0D0+eccb*DSIN(omegab))**3 ) / ( (DSQRT(1.0D0-eccb**2))**3 ) )
 mu = -0.3555572027034898D0
 sigmaneg = 0.12157839230062002D0
 sigmapos = 0.32574170730200136D0
 IF( logPsib .LE. mu ) THEN
	 loglikeR = (logPsib - mu)**2 / (2.0D0*sigmaneg**2)
	 loglikeR = DSQRT(2.0/pi)*DEXP(-loglikeR)
	 loglikeR = DLOG( loglikeR/(sigmaneg+sigmapos) )
 ELSE
	 loglikeR = (logPsib - mu)**2 / (2.0D0*sigmapos**2)
	 loglikeR = DSQRT(2.0/pi)*DEXP(-loglikeR)
	 loglikeR = DLOG( loglikeR/(sigmaneg+sigmapos) )
 END IF

 ! === Sum up all loglikes (product of likes) ===
 loglike = loglikeP + loglikeR

END SUBROUTINE models
! ======================================================================

! ======================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ======================================================================

END MODULE like

