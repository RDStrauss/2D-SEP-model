   program Frankel
      
   
  IMPLICIT NONE
  
  REAL, PARAMETER :: PI = 3.14159, time_max = 20.
  INTEGER, PARAMETER :: N = 121, M = 108, NK = 49
  REAL :: delta_r, delta_phi, R(N), PHI(M), delta_t, time, left_lim, right_lim
  REAL :: f_old(N,M,NK), f_new(N,M,NK), CFL_coeff, limiter, f_holder(N,M,NK), K_rphi(N,M,NK), MU(NK)
  REAL :: A(N,M,NK), B(N,M,NK), C(N,M,NK), D(N,M,NK), max_a, max_b, max_c, max_d, V_r(N,M,NK), V_phi(N,M,NK), left_deriv, right_deriv
  INTEGER :: i,j,k, indeks_r, indeks_phi
  REAL :: delta_mu, f_omni, F(N,M,NK), max_f, E(N,M,NK), G(N,M,NK), max_g, time_printer, anis
  REAL :: time_begin, time_end 
 
   REAL :: injection_broadness = 15./180.*PI, acceleration_time = 0.1, escape_time = 1. 
   REAL :: maximum_intensity(M), maximum_time(M), maximum_anisotropy(M), onset_time(M), background_intensity = 0.01
   INTEGER :: onset_time_flag(M)
   
  CALL omp_set_num_threads(6)
 !---------------------------------------------------
 ! Files to write to
 
 OPEN(100,file='r_grid.txt',status='unknown')
 OPEN(200,file='phi_grid.txt',status='unknown')
 OPEN(300,file='pitch_angle_dependence.txt',status='unknown')
 OPEN(400,file='output.txt',status='unknown')
 OPEN(500,file='mu_grid.txt',status='unknown')
 OPEN(600,file='coeffs_mu.txt',status='unknown')
 OPEN(700,file='coeffs_r.txt',status='unknown') 
 OPEN(800,file='omni_time_best_connection.txt',status='unknown') 
 OPEN(801,file='omni_time_earth.txt',status='unknown') 
 OPEN(802,file='omni_time_stereo_A.txt',status='unknown') 
 OPEN(803,file='omni_time_stereo_B.txt',status='unknown') 
 OPEN(900,file='omni_phi.txt',status='unknown')
 OPEN(999,file='maximums.txt',status='unknown')
 !---------------------------------------------------
  CALL CPU_TIME ( time_begin )
 
 ! TODO: Fix radia grid to align with Samira
 
 !Set up the grids
 delta_r = (1.2 - 0.05)/(N - 1)
 delta_phi = 2.*PI/(M)
 delta_mu = 2./NK

   MU(1) = -1. + Delta_mu/2.
  
     WRITE(500,"(1(ES18.8))") MU(1)
     
 DO k = 2, NK
 
  MU(k) = MU(k - 1) + Delta_mu
  
    WRITE(500,"(1(ES18.8))") MU(k)
    
 END DO
 
 R(1) = 0.05
 
  WRITE(100,"(1(ES18.8))") R(1)
 
 DO i = 2, N
 
  R(i) = R(i - 1) + delta_r
 
  WRITE(100,"(1(ES18.8))") R(i)
 
 END DO
 
 PHI(1) = -PI
 
 WRITE(200,"(1(ES18.8))") PHI(1)
 
 DO j = 2, M
 
  PHI(j) = PHI(j - 1) + delta_phi
 
 WRITE(200,"(1(ES18.8))") PHI(j)
 
 END DO
 !--------------------------------------------------- 
 ! Set up initial conditions and Coeffs

 DO j = 1, M
 
    onset_time_flag(j) = 0
    maximum_intensity(j) = 0.
    maximum_time(j) = 0.
    maximum_anisotropy(j) = 0.
    onset_time(j) = time_max   ! Peuter-baar, maar nie nul nie. 
 
 END DO
 
 
 
 DO i = 1, N
 
  DO j = 1, M
  
    DO k = 1, NK
    
      f_new(i,j,k) = 0.0001
      f_old(i,j,k) = 0.0001
   
   END DO
   
  END DO
  
END DO


CALL DEFINE_COEFFCIENTS(N,M,NK,R,A,B,C,D,delta_r,V_r,V_phi,K_rphi,MU,E,F,G,delta_mu,delta_phi)

 !---------------------------------------------------
 ! Set up CFL conditions
 
 CFL_coeff = 0.9
 
 max_a = MAXVAL(A)
 max_b = MAXVAL(B)
 max_c = MAXVAL(C)
 max_d = MAXVAL(D)
 max_f = MAXVAL(F)
 max_g = MAXVAL(ABS(G))
 
 WRITE(*,*) 'max_g', max_g
 
 delta_t = CFL_coeff*MIN(ABS(0.5*delta_r*delta_r/max_b),ABS(0.5*delta_phi*delta_phi/max_d),ABS(0.5*delta_r/MAXVAL(V_r)),ABS(0.5*delta_phi/MAXVAL(V_phi)),ABS(0.5*delta_mu*delta_mu/max_f),ABS(0.5*delta_mu/max_g))
 
 WRITE(*,*) 'Delta t: ', delta_t
 !---------------------------------------------------
 ! Iterate in time
 
 time = 0.
 time_printer = 0.
 
 DO WHILE (time.LT.time_max)
   
  time = time + delta_t
  
 !Time dependent boundary conditions
  DO k = 1, NK
 
 DO j = 1, M
  
 !This is the wird constant injection form WOlfgang, injecting only for first 10 minutes
  f_old(1,j,k) = 0.
  
 ! IF (time.LT.0.1667) THEN
  
 !f_old(1,j,k) = 1000.*exp(-(PHI(j) - pi/2.)*(PHI(j) - pi/2.)/2./injection_broadness/injection_broadness)
 
 ! ENDIF
 
 ! This is the normnal Reid injection in time
  f_old(1,j,k) = 1000.*exp(-(PHI(j) - pi/2.)*(PHI(j) - pi/2.)/2./injection_broadness/injection_broadness)/time*exp(-acceleration_time/time - time/escape_time)
  
 f_new(1,j,k) = f_old(1,j,k)
 
 END DO
 
 END DO
 
!---------------------------------------------------
! Do the r integration - convection

    !$omp parallel do private(i,j,left_lim,right_lim,limiter)
DO k = 1, NK

  DO j = 1, M
  
  !Do the convertion to conserved flux; F = r*r*f
  
  DO i = 1, N
  
    f_old(i,j,k) = f_old(i,j,k)*R(i)*R(i)
  
  END DO
  
      DO i = 2, N - 1
      
	IF (V_r(i,j,k).GT.0.) THEN
	
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_r*(V_r(i,j,k)*f_old(i,j,k) - V_r(i-1,j,k)*f_old(i-1,j,k))
 
	ELSE
	
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_r*(V_r(i+1,j,k)*f_old(i+1,j,k) - V_r(i,j,k)*f_old(i,j,k))
	
	END IF
	
      END DO
 

      DO i = 2, N - 1
      
	IF (V_r(i,j,k).GT.0.) THEN
      
      	  left_lim = (V_r(i,j,k)*f_new(i,j,k) - V_r(i-1,j,k)*f_new(i-1,j,k))/2.
	  right_lim = (V_r(i+1,j,k)*f_new(i+1,j,k) - V_r(i,j,k)*f_new(i,j,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

          
            !limiter = 0.
  
IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 
	    f_holder(i,j,k) = V_r(i,j,k)*f_new(i,j,k) + limiter
	    
	  ELSE
	  
      	  left_lim = (V_r(i-1,j,k)*f_new(i-1,j,k) - V_r(i,j,k)*f_new(i,j,k))/2.
	  right_lim = (V_r(i,j,k)*f_new(i,j,k) - V_r(i+1,j,k)*f_new(i+1,j,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

          
            !limiter = 0.
  
IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
	  
	  f_holder(i,j,k) = V_r(i,j,k)*f_new(i,j,k) + limiter
	  
	  	  
	  ENDIF
	  
	  f_holder(1,j,k) = V_r(1,j,k)*f_new(1,j,k)
 
	END DO
	
	DO i = 2, N - 1
 
	  IF (V_r(i,j,k).GT.0.) THEN
 
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_r*(f_holder(i,j,k) - f_holder(i-1,j,k))
	
	  ELSE
	  
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_r*(f_holder(i+1,j,k) - f_holder(i,j,k))
	  
	  ENDIF
	
	END DO
	

!Boundary equations

    
    
!    IF (V_r(N,k).GT.0.) THEN
    
!      f_new(N,j,k) = f_old(N,j,k) + Delta_t/Delta_r*f_old(N-1,j,k)*V_r(N-1,k) - Delta_t/Delta_r*f_old(N,j,k)*V_r(N,k)    
    
!    ELSE
    
!      f_new(N,j,k) = f_old(N,j,k) + Delta_t/Delta_r*f_old(N-1,j,k)*V_r(N-1,k)    
    
!    ENDIF
       
  END DO ! End do for loop over j
  
!Update f_new and convert back to f =F/r/r

  DO i = 1, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)/R(i)/R(i)
    
    END DO
    
  END DO
  
END DO
 !---------------------------------------------------
! Do the phi integration - convection

!$omp parallel do private(i,j,left_lim,right_lim,limiter)
DO k = 1, NK

  DO i = 2, N - 1
  
  IF (V_phi(i,j,k).LT.0.) THEN

    DO j = 2, M - 1

    	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_phi*(V_phi(i,j+1,k)*f_old(i,j+1,k) - V_phi(i,j,k)*f_old(i,j,k))
    
    END DO
    
    DO j = 2, M - 1
    
       	  left_lim = (f_new(i,j-1,k) - f_new(i,j,k))/2.
	  right_lim = (f_new(i,j,k) - f_new(i,j+1,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

          
            !limiter = 0.
  
IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 
	    f_holder(i,j,k) = f_new(i,j,k) + limiter
 
    END DO
    
    
    !Boundary equations
    f_new(i,1,k) = f_old(i,1,k) - delta_t/2./delta_phi*(V_phi(i,2,k)*f_old(i,2,k) - V_phi(i,1,k)*f_old(i,1,k))/2.
    
          left_lim = (f_new(i,M,k) - f_new(i,1,k))/2.
	  right_lim = (f_new(i,1,k) - f_new(i,2,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
    IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
    IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
    
    	    f_holder(i,1,k) = f_new(i,1,k) + limiter
    	    
    f_new(i,M,k) = f_old(i,M,k) - delta_t/2./delta_phi*(V_phi(i,1,k)*f_old(i,1,k) - V_phi(i,M,k)*f_old(i,M,k))/2.    
    
          left_lim = (f_new(i,M-1,k) - f_new(i,M,k))/2.
	  right_lim = (f_new(i,M,k) - f_new(i,1,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
    IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
    IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
    
    	    f_holder(i,M,k) = f_new(i,M,k) + limiter
    
           
    DO j = 2, M - 1
    
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_phi*(V_phi(i,j+1,k)*f_holder(i,j+1,k) - V_phi(i,j,k)*f_holder(i,j,k))   
        
    END DO
    
    
    f_new(i,1,k) = f_old(i,1,k) - delta_t/delta_phi*(V_phi(i,2,k)*f_holder(i,2,k) - V_phi(i,1,k)*f_holder(i,1,k))
    f_new(i,M,k) = f_old(i,M,k) - delta_t/delta_phi*(V_phi(i,1,k)*f_holder(i,1,k) - V_phi(i,M,k)*f_holder(i,M,k))

      
  ELSE
  
     DO j = 2, M - 1

    	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_phi*(V_phi(i,j,k)*f_old(i,j,k) - V_phi(i,j-1,k)*f_old(i,j-1,k))
    
    END DO

    	  f_new(i,1,k) = f_old(i,1,k) - delta_t/2./delta_phi*(V_phi(i,1,k)*f_old(i,1,k) - V_phi(i,M,k)*f_old(i,M,k))
    	  f_new(i,M,k) = f_old(i,M,k) - delta_t/2./delta_phi*(V_phi(i,M,k)*f_old(i,M,k) - V_phi(i,M-1,k)*f_old(i,M-1,k))
    
    DO j = 2, M - 1
    
       	  left_lim = (f_new(i,j,k) - f_new(i,j-1,k))/2.
	  right_lim = (f_new(i,j+1,k) - f_new(i,j,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 
	    f_holder(i,j,k) = f_new(i,j,k) + limiter
 
    END DO

           left_lim = (f_new(i,1,k) - f_new(i,M,k))/2.
	  right_lim = (f_new(i,2,k) - f_new(i,1,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 
	    f_holder(i,1,k) = f_new(i,1,k) + limiter

	  left_lim = (f_new(i,M,k) - f_new(i,M-1,k))/2.
	  right_lim = (f_new(i,1,k) - f_new(i,M,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)

IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 
	    f_holder(i,M,k) = f_new(i,M,k) + limiter
    
    DO j = 2, M - 1
 
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_phi*(V_phi(i,j,k)*f_holder(i,j,k) - V_phi(i,j-1,k)*f_holder(i,j-1,k))   
    
    
    END DO
    
    
!Boundary equations
      
    f_new(i,1,k) = f_old(i,1,k) - delta_t/delta_phi*(V_phi(i,1,k)*f_holder(i,1,k) - V_phi(i,M,k)*f_holder(i,M,k))
    f_new(i,M,k) = f_old(i,M,k) - delta_t/delta_phi*(V_phi(i,M,k)*f_holder(i,M,k) - V_phi(i,M-1,k)*f_holder(i,M-1,k)) 
  
  
  ENDIF
    
  END DO

  END DO
  
!Update f_new

DO k = 1, Nk

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
    
    END DO
    
  END DO

END DO
!---------------------------------------------------
! Do the r integration - diffusion

!$omp parallel do private(i,j,left_deriv,right_deriv)
DO k = 1, NK

  DO j = 1, M
  
      DO i = 2, N - 1
      
	IF ((j.NE.1).AND.(j.NE.M)) THEN
      
      left_deriv = (f_old(i-1,j+1,k) - f_old(i-1,j-1,k))/2./delta_phi
      right_deriv = (f_old(i+1,j+1,k) - f_old(i+1,j-1,k))/2./delta_phi 
      
	ENDIF

	IF (j.EQ.1) THEN
      
      left_deriv = (f_old(i-1,2,k) - f_old(i-1,M,k))/2./delta_phi
      right_deriv = (f_old(i+1,2,k) - f_old(i+1,M,k))/2./delta_phi 
      
	ENDIF
	
	IF (j.EQ.M) THEN
      
      left_deriv = (f_old(i-1,1,k) - f_old(i-1,M-1,k))/2./delta_phi
      right_deriv = (f_old(i+1,1,k) - f_old(i+1,M-1,k))/2./delta_phi 
      
	ENDIF
	
	f_new(i,j,k) = f_old(i,j,k) + A(i,j,k)*delta_t/delta_r/2.*(f_old(i + 1,j,k) - f_old(i - 1,j,k)) + B(i,j,k)*delta_t/delta_r/delta_r*(f_old(i + 1,j,k) - 2.*f_old(i,j,k) + f_old(i - 1,j,k)) + K_rphi(i,j,k)/R(i)*delta_t/2./delta_r*(right_deriv - left_deriv)
     	      
      END DO
      
!Boundary equations
      
!      f_new(1,j) = f_old(1,j) - Delta_t/Delta_r*(-(B(1) + B(2))/2.*(f_old(2,j) - f_old(1,j))/delta_r)
!      f_new(N,j,k) = f_old(N,j,k) + Delta_t/Delta_r*(-(B(N,k) + B(N-1,k))/2.*(f_old(N,j,k) - f_old(N-1,j,k))/delta_r) - Delta_t/Delta_r*B(N,k)*f_old(N,j,k)/delta_r 

  END DO
  
!Update f_new

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
    
    END DO
    
  END DO
 
 END DO

!---------------------------------------------------
! Do the phi integration - diffusion

!$omp parallel do private(i,j,left_deriv,right_deriv)
DO k = 1, NK

  DO i = 2, N - 1
  
      DO j = 2, M - 1
      	
      left_deriv = (f_old(i+1,j-1,k) - f_old(i-1,j-1,k))/2./delta_r
      right_deriv = (f_old(i+1,j+1,k) - f_old(i-1,j+1,k))/2./delta_r 
      	
	  f_new(i,j,k) = f_old(i,j,k) + C(i,j,k)*delta_t/delta_phi/2.*(f_old(i,j + 1,k) - f_old(i,j - 1,k)) + D(i,j,k)*delta_t/delta_phi/delta_phi*(f_old(i,j + 1,k) - 2.*f_old(i,j,k) + f_old(i,j - 1,k)) + K_rphi(i,j,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
	      
      END DO
      
!Boundary equations
      
      left_deriv = (f_old(i+1,M,k) - f_old(i-1,M,k))/2./delta_r
      right_deriv = (f_old(i+1,2,k) - f_old(i-1,2,k))/2./delta_r 
      
    f_new(i,1,k) = f_old(i,1,k) + C(i,1,k)*delta_t/delta_phi/2.*(f_old(i,2,k) - f_old(i,M,k)) + D(i,1,k)*delta_t/delta_phi/delta_phi*(f_old(i,2,k) - 2.*f_old(i,1,k) + f_old(i,M,k)) + K_rphi(i,1,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
    
    
    left_deriv = (f_old(i+1,M-1,k) - f_old(i-1,M-1,k))/2./delta_r
    right_deriv = (f_old(i+1,1,k) - f_old(i-1,1,k))/2./delta_r 
      
    
    f_new(i,M,k) = f_old(i,M,k) + C(i,M,k)*delta_t/delta_phi/2.*(f_old(i,1,k) - f_old(i,M-1,k)) + D(i,M,k)*delta_t/delta_phi/delta_phi*(f_old(i,1,k) - 2.*f_old(i,M,k) + f_old(i,M-1,k)) + K_rphi(i,M,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
      
  END DO
  
!Update f_new

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
    
    END DO
    
  END DO

END DO
!---------------------------------------------------
! Mu diffusion

!$omp parallel do private(j,k)
DO i = 1, N

  DO j = 1, M
  
    DO k = 2, NK - 1
    
      f_new(i,j,k) = f_old(i,j,k) + E(i,j,k)*Delta_t/2./Delta_mu*(f_old(i,j,k+1) - f_old(i,j,k-1)) + F(i,j,k)*Delta_t/Delta_mu/Delta_mu*(f_old(i,j,k+1)- 2.*f_old(i,j,k) + f_old(i,j,k-1))
    
    END DO

!At the boundaries

  f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*(-(F(i,j,1)+F(i,j,2))/2.*(f_old(i,j,2) - f_old(i,j,1))/Delta_mu)
  f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu*(-(F(i,j,NK)+F(i,j,NK-1))/2.*(f_old(i,j,NK) - f_old(i,j,NK-1))/Delta_mu)
    
   END DO
   
  END DO

!Update f_new
DO k = 1, NK

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
    
    END DO
    
  END DO

END DO
!---------------------------------------------------
! Mu convection

!$omp parallel do private(j,k,left_lim,right_lim,limiter)
DO i = 2, N

  DO j = 1, M
   
    DO k = 2, NK - 1

      IF (G(i,j,k).GE.0.) THEN
    
      f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(G(i,j,k)*f_old(i,j,k) - G(i,j,k-1)*f_old(i,j,k-1))/2.
    
     ELSE
     
     f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(G(i,j,k+1)*f_old(i,j,k+1) - G(i,j,k)*f_old(i,j,k))/2.
     
     END IF
    
    END DO
    
    
    IF (G(i,j,1).GE.0.) THEN
    
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu/2.*f_old(i,j,1)*G(i,j,1)
    
    ELSE
    
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*(G(i,j,2)*f_old(i,j,2) - G(i,j,1)*f_old(i,j,1))/2.
    
    END IF
    
    IF (G(i,j,Nk).GE.0.) THEN    
    
     f_new(i,j,Nk) = f_old(i,j,Nk) - Delta_t/Delta_mu*(G(i,j,Nk)*f_old(i,j,Nk) - G(i,j,Nk-1)*f_old(i,j,Nk-1))/2.   
    
    ELSE
        
    f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu/2.*f_old(i,j,NK)*G(i,j,NK)  
    
    END IF
    
        
    DO k = 2, NK - 1
    
    IF (G(i,j,k).GE.0.) THEN
    
      left_lim = (G(i,j,k)*f_new(i,j,k) - G(i,j,k-1)*f_new(i,j,k-1))/2.
	  right_lim = (G(i,j,k+1)*f_new(i,j,k+1) - G(i,j,k)*f_new(i,j,k))/2.

    ELSE
    
      left_lim = (G(i,j,k-1)*f_new(i,j,k-1) - G(i,j,k)*f_new(i,j,k))/2.
	  right_lim = (G(i,j,k)*f_new(i,j,k) - G(i,j,k+1)*f_new(i,j,k+1))/2.
    
    END IF
	  
	  limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
	  
	  IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN
	  
	  IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs are different
    
  
    f_holder(i,j,k) = G(i,j,k)*f_new(i,j,k) + limiter
    
    END DO
    
    f_holder(i,j,1) = G(i,j,1)*f_new(i,j,1)
    f_holder(i,j,NK) = G(i,j,NK)*f_new(i,j,NK)
    
    DO k = 2, NK - 1
    
        IF (G(i,j,k).GE.0.) THEN
    
      	f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(f_holder(i,j,k) - f_holder(i,j,k-1))
      	
      	ELSE
      	
       	f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(f_holder(i,j,k+1) - f_holder(i,j,k))
       	
      	END IF
    
    END DO
    
     
    IF (G(i,j,1).GE.0.) THEN 
    
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*f_holder(i,j,1)
 
  ELSE
  
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*(f_holder(i,j,2) - f_holder(i,j,1))
  
  END IF

    IF (G(i,j,Nk).GE.0.) THEN 
    
    f_new(i,j,Nk) = f_old(i,j,Nk) - Delta_t/Delta_mu*(f_holder(i,j,Nk) - f_holder(i,j,Nk-1))
 
  ELSE
  
    f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu*f_holder(i,j,NK)
  
  END IF
 
  
 
   END DO
   
 END DO

 
!Update f_new
DO k = 1, NK

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
          
    END DO
    
  END DO

END DO


! Write omni-directional intensity at r=?, phi=? vs time
! OPEN(801,file='omni_time_earth.txt',status='unknown') 
! OPEN(802,file='omni_time_stereo_A.txt',status='unknown') 
! OPEN(803,file='omni_time_stereo_B.txt',status='unknown') 


time_printer = time_printer + delta_t
IF (time_printer.GE.0.01) THEN

  time_printer = 0.
 
 
  ! This should correspond to best magnetic connection at Earth
  !-------------------------------------------------------------------------------------
  indeks_phi = (33./180.*PI - phi(1) )/delta_phi + 1.
  indeks_r = (1. - R(1))/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)/PI*180., 'best connection'
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(800,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni
       
   ! This should correspond to azimuth of Earth
  !-------------------------------------------------------------------------------------
  indeks_phi = (33./180.*PI - phi(1) + 82./180.*PI)/delta_phi + 1.
  indeks_r = (1. - R(1))/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)/PI*180., 'Earth'
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(801,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

   ! This should correspond to azimuth of STEREO A
  !-------------------------------------------------------------------------------------
  indeks_phi = (33./180.*PI - phi(1) + 120./180.*PI)/delta_phi + 1.
  indeks_r = (1. - R(1))/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)/PI*180., 'STEREO A'
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(802,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

   ! This should correspond to azimuth of STEREO B
  !-------------------------------------------------------------------------------------
  indeks_phi = (33./180.*PI - phi(1) + -21./180.*PI)/delta_phi + 1.
  indeks_r = (1. - R(1))/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)/PI*180., 'STEREO B'
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(803,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni
       
  !-------------------------------------------------------------------------------------       
       
  ENDIF

!========================================================
! Write peak intensity vs phi_grid

  indeks_r = (1 - R(1) + delta_r)/delta_r + 1.
  
  i = indeks_r
  
  DO j = 1, M
  
  anis = 0.
  f_omni = 0.
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)*delta_mu/2.
    anis = anis + mu(k)*f_old(i,j,k)*delta_mu/2.
 
  END DO ! For k
  
  IF (f_omni.GT.maximum_intensity(j)) THEN
  
    maximum_intensity(j) = f_omni
    maximum_time(j) = time
    
  ENDIF
  
  anis = 3.*anis/f_omni
  
  IF ((anis.GT.maximum_anisotropy(j)).AND.(onset_time_flag(j).EQ.1)) THEN
  
    maximum_anisotropy(j) = anis
    
  ENDIF  
  
  
  IF ((onset_time_flag(j).EQ.0).AND.(f_omni.GT.background_intensity)) THEN 
  
    onset_time(j) = time
    onset_time_flag(j) = 1
  
  END IF
  
  
  END DO ! for j



!---------------------------------------------------
 
 END DO ! While time loop
 
  DO j = 1, M
 
        WRITE(999,"(5(ES18.8))") phi(j), maximum_time(j), maximum_intensity(j), maximum_anisotropy(j), onset_time(j)
 
 END DO
 
 !---------------------------------------------------
 !Write omni directional intensity vs r,phi at last time step.
   DO i = 1, N
  
    DO j = 1, M
    
    f_omni = 0.
    
      DO k = 1, NK
      
	f_omni = f_omni + f_old(i,j,k)
      
      END DO
    
       WRITE(400,"(1(ES18.8))") f_omni*delta_mu/2.
    
    END DO
    
  END DO


!Write pitch-angle dependence as r=?, phi=? at last time step 
  i = indeks_r

  j = indeks_phi
  
    DO k = 1, NK
  
    WRITE(300,"(4(ES18.8))") r(i), phi(j), mu(k), f_old(i,j,k)
  
    END DO

!Write the phi dependence of omni-directional intensity at r=?

  i = indeks_r
  
     DO j = 1, M
    
    f_omni = 0.
    
      DO k = 1, NK
      
	f_omni = f_omni + f_old(i,j,k)
      
      END DO
    
       WRITE(900,"(3(ES18.8))") r(i), phi(j), f_omni*delta_mu/2.
    
    END DO 
    
      CALL CPU_TIME ( time_end )
    
   PRINT *, 'Time of operation was ', &
time_end - time_begin, ' seconds'     

 
 CLOSE(100)
 CLOSE(200)
 CLOSE(300)
 CLOSE(400)
 CLOSE(500)
 CLOSE(600)
 CLOSE(700)
 CLOSE(800)
 CLOSE(801)
 CLOSE(802)
 CLOSE(803)
 CLOSE(900)
CLOSE(999)
 
 END
 
 !----------------------------------------------------
 
 SUBROUTINE DEFINE_COEFFCIENTS(N,M,NK,R,A,B,C,D,delta_r,V_r,V_phi,K_rphi,MU,E,F,G,delta_mu,delta_phi)

 IMPLICIT NONE
 
 REAL :: R(N), A(N,M,NK), B(N,M,NK), C(N,M,NK), D(N,M,NK), delta_r, delta_mu, delta_phi
 INTEGER :: i, j, k, N, M, NK
 REAL, PARAMETER :: PI = 3.14159
 REAL :: cospsi(N,M), sinpsi(N,M), Kappa_perp(N,M,NK), V(NK)
 REAL :: K_rr(N,M,NK), K_phiphi(N,M,NK), K_rphi(N,M,NK), K_rr_dr(N,M,NK), K_phiphi_dr(N,M,NK), K_rphi_dr(N,M,NK)
 REAL :: K_rphi_dphi(N,M,NK), K_phiphi_dphi(N,M,NK)
 REAL :: V_r(N,M,NK), V_phi(N,M,NK), MU(NK), L(N,M), speed, F(N,M,NK), E(N,M,NK), G(N,M,NK), integral
 REAL :: energy, speed_of_light, rigidity, beta, electron_rest_energy
 REAL :: lambda, alpha
 
 
 lambda = 0.08 !0.08 was the Droege default value
 alpha = 0.13 !0.13 is the Droege default value

 energy = 0.100 ! in MeV
 speed_of_light = 7.2 !in units of AU/hr
 electron_rest_energy = 0.51 ! in MeV
 rigidity = SQRT(energy*(energy + 2.*electron_rest_energy)) ! in MV
 beta = rigidity/(energy + electron_rest_energy)
 speed = beta*speed_of_light

 CALL READ_FILES(sinpsi, cospsi, L, N, M)


 DO i = 1, N
 
 DO j = 1, M
 
  DO k = 1, NK
  
    ! F = D_mumu; here just pitch-angle dependence
    F(i,j,k) = (1. - mu(k)*mu(k))*((ABS(MU(k)))**(1.67 - 1.) + 0.5)
  
  END DO
 
 END DO
 
 END DO
 
 integral = 0.
 
  DO k = 1, NK
    
    i = 1
    j = 1
    
    ! To nomalize mean free path to lambda
    integral = integral + (1. - mu(k)*mu(k))*(1. - mu(k)*mu(k))/F(i,j,k)*delta_mu
  
  END DO 
 
 
 DO i = 1, N
 
 DO j = 1, M
 
  DO k = 1, NK
  
    ! Here the radial dependce; see Droege et al. papers
    F(i,j,k) = F(i,j,k)/lambda*integral*3.*speed/8.*cospsi(i,j)*cospsi(i,j)
    !F(i,j,k) = F(i,j,k)/lambda*integral*3.*speed/8.  
  END DO
  
  END DO
  
 END DO
 
 
 DO i = 1, N
 
 DO j = 1, M
 
  DO k = 2, NK - 1
  
    ! E = dD_mumu/dmu
    E(i,j,k) = (F(i,j,k+1) - F(i,j,k-1))/2./delta_mu
  
  END DO

    E(i,j,1) = (F(i,j,2) - F(i,j,1))/delta_mu
    E(i,j,NK) = (F(i,j,NK) - F(i,j,NK - 1))/delta_mu
  
 END DO
 
 END DO

 
 
 DO i = 1, N

  DO j = 1, M
 
  DO k = 1, NK
    
    ! G = focusing term
    G(i,j,k) = speed/2./L(i,j)*(1. - mu(k)*mu(k))
  
!  WRITE(*,*) G(i,j,k), L(i,j)
  
  END DO
  
 END DO
 
 END DO
 
 
 DO k = 1, NK
 
 ! V = parallel speed in HMF-aligned coordinates
 V(k) = speed*MU(k)
 
 END DO
 
 
 ! Perpendicular coefficient
 DO i = 1, N
 
 DO j = 1, M
 
   DO k = 1, NK

   ! D_perp; see Droege et al. paper
    Kappa_perp(i,j,k) = speed/3.*alpha*lambda*R(i)*R(i)*cospsi(i,j)/cospsi(i,j)/cospsi(i,j)*SQRT(1. - mu(k)*mu(k))
    !Kappa_perp(i,j,k) = speed/3.*alpha*lambda*(4/PI)*SQRT(1. - mu(k)*mu(k))  
    END DO
 
 END DO
 
 END DO
 
  
! Diffusion tensor in spiral coordinates
  DO i = 1, N
  
  DO j = 1, M
  
   DO k = 1, NK
  
    K_rr(i,j,k) = Kappa_perp(i,j,k)*sinpsi(i,j)*sinpsi(i,j)
    K_phiphi(i,j,k) = Kappa_perp(i,j,k)*cospsi(i,j)*cospsi(i,j)
    K_rphi(i,j,k) = Kappa_perp(i,j,k)*cospsi(i,j)*sinpsi(i,j)
    
   END DO
    
  END DO
  
  END DO
  
 !The derivatives of the tensor
 DO j = 1, M
 
 DO k = 1, NK
 
  DO i = 2, N - 1
  
    K_rr_dr(i,j,k) = (K_rr(i + 1,j,k) - K_rr(i - 1,j,k))/2./delta_r
    K_phiphi_dr(i,j,k) = (K_phiphi(i + 1,j,k) - K_phiphi(i - 1,j,k))/2./delta_r
    K_rphi_dr(i,j,k) = (K_rphi(i + 1,j,k) - K_rphi(i - 1,j,k))/2./delta_r
    
  END DO
 
 
  K_rr_dr(1,j,k) = (K_rr(2,j,k) - K_rr(1,j,k))/delta_r
  K_rr_dr(N,j,k) = (K_rr(N,j,k) - K_rr(N - 1,j,k))/delta_r
  
  K_phiphi_dr(1,j,k) = (K_phiphi(2,j,k) - K_phiphi(1,j,k))/delta_r
  K_phiphi_dr(N,j,k) = (K_phiphi(N,j,k) - K_phiphi(N - 1,j,k))/delta_r
  
  K_rphi_dr(1,j,k) = (K_rphi(2,j,k) - K_rphi(1,j,k))/delta_r
  K_rphi_dr(N,j,k) = (K_rphi(N,j,k) - K_rphi(N - 1,j,k))/delta_r
 
 END DO
 
 END DO
 
 
 ! Phi derivatives 
 
DO i = 1, N

  DO k = 1, NK
 
    DO j = 2, M - 1
 
     K_rphi_dphi(i,j,k) = (K_rphi(i,j + 1,k) - K_rphi(i, j - 1,k))/2./delta_phi
     K_phiphi_dphi(i,j,k) = (K_phiphi(i,j + 1,k) - K_phiphi(i, j - 1,k))/2./delta_phi
 
    END DO
 
     K_rphi_dphi(i,1,k) = (K_rphi(i,2,k) - K_rphi(i,1,k))/delta_phi
     K_rphi_dphi(i,M,k) = (K_rphi(i,M,k) - K_rphi(i,M - 1,k))/delta_phi
     
     K_phiphi_dphi(i,1,k) = (K_phiphi(i,2,k) - K_phiphi(i,1,k))/delta_phi
     K_phiphi_dphi(i,M,k) = (K_phiphi(i,M,k) - K_phiphi(i,M - 1,k))/delta_phi
  
 END DO
 
END DO
 
 
! Define numerical coeffs
DO k = 1, NK

DO j = 1, M

 DO i = 1, N
 
    ! Add phi derivative terms
    A(i,j,k) = 2.*K_rr(i,j,k)/R(i) + K_rr_dr(i,j,k) + K_rphi_dphi(i,j,k)/R(i)
    B(i,j,k) = K_rr(i,j,k)
    C(i,j,k) = K_rphi(i,j,k)/R(i)/R(i) + K_rphi_dr(i,j,k)/R(i) + K_phiphi_dphi(i,j,k)/R(i)/R(i)
    D(i,j,k) = K_phiphi(i,j,k)/R(i)/R(i)

    ! V_r, V_Phi are the V_|| components in spherical coordinates
    V_r(i,j,k) = V(k)*cospsi(i,j)
    V_phi(i,j,k) = -V(k)*sinpsi(i,j)/R(i)
  
!  WRITE(*,*) V_r(i), V_phi(i), sinpsi(i), sinpsi(i)*sinpsi(i)+cospsi(i)*cospsi(i)
  
END DO

END DO

END DO


!Write some stuff to a file

DO i = 1, N

    j = 50

  WRITE(700,"(4(ES18.8))") R(i), sinpsi(i,j), cospsi(i,j), L(i,j)

END DO

i = N

DO k = 1, NK

    j = 50

  WRITE(600,"(4(ES18.8))") MU(k), F(i,j,k), E(i,j,k), Kappa_perp(i,j,k)/speed*3.

END DO

RETURN

END
!----------------------------------------------------
!===================================================================
	SUBROUTINE READ_FILES(sinpsi, cospsi, L, N, M)
! Read in heliospheric quantities

	IMPLICIT NONE

	INTEGER :: N, M, i, j
	REAL :: sinpsi(N,M), cospsi(N,M), L(N,M), placeholder

	OPEN(unit = 2, file = "Parker_CosPsi_108_Phi.dat")


	DO i = 1, N, 1

        DO j = 1, M, 1
	
            READ(2,*) placeholder

            cospsi(i,j) = placeholder
	
        END DO

	END DO

    CLOSE(2)

	OPEN(unit = 2, file = "Parker_SinPsi_108_Phi.dat")

	
	DO i = 1, N, 1
	
        DO j = 1, M, 1
	
	        READ(2,*) placeholder

            sinpsi(i,j) = placeholder
	
        END DO

	END DO
	
	WRITE(*,*) 'File read in . . .'

	CLOSE(2)

	OPEN(unit = 2, file = "Parker_Focus_108.dat")
	
	DO i = 1, N, 1
	
        DO j = 1, M, 1
	
	        READ(2,*) placeholder

	        !IF (placeholder.GE.1.) placeholder = placeholder + 2.
	        !IF (placeholder.LT.1.) placeholder = placeholder - 2.
	        
            L(i,j) = 1./placeholder
	
        END DO

	END DO
	
	WRITE(*,*) 'File read in . . .'

	CLOSE(2)	
	
	RETURN
	END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
