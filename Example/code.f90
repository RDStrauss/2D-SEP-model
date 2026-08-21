! Compile flags: $ gfortran -fopenmp -O2 Version_1.f90
! On most system need to increase stack size first, ulimit -s unlimited
  program AaP
  
  
  use omp_lib
   
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N = 100, M = 100, NK = 49
  INTEGER :: i, j, k, indeks_r, indeks_phi, onset_time_flag(M), species
  
  REAL(KIND=8), PARAMETER :: PI = 3.14159, time_max = 5.
  REAL(KIND=8), PARAMETER :: minute = 0.016666667
  REAL(KIND=8) :: delta_r, delta_phi, R(N), PHI(M), delta_t, time, left_lim, right_lim
  REAL(KIND=8) :: f_old(N,M,NK), f_new(N,M,NK), CFL_coeff, limiter, f_holder(N,M,NK), K_rphi(N,NK), MU(NK)
  REAL(KIND=8) :: A(N,NK), B(N,NK), C(N,NK), D(N,NK), max_a, max_b, max_c, max_d, V_r(N,NK), V_phi(N,NK), left_deriv, right_deriv
  REAL(KIND=8) :: delta_mu, f_omni, F(N,NK), max_f, E(N,NK), G(N,NK), max_g, time_printer, anis
  REAL(KIND=8) :: time_begin, time_end, energy
  REAL(KIND=8) :: maximum_intensity(M), maximum_time(M), maximum_anisotropy(M), onset_time(M), background_intensity = 0.01
  REAL(KIND=8) :: acceleration_time = 1.0, escape_time = 2.0 
  REAL(KIND=8) :: injection_broadness = 25./180.*PI, phi_source = 90./180*PI
  REAL(KIND=8) :: Earth_phi, Earth_r, Mars_phi, Mars_r, STA_phi, STA_r, STB_phi, STB_r
  REAL(KIND=8) :: Bepi_phi, Bepi_r, PSP_phi, PSP_r, SolO_phi, SolO_r
  
  ! species = a way to switch between electrons (species = 1) and protons (species = 2)  
  species = 1
  
  ! energy = energy of the mono-energetic SEP distribution in MeV
  energy = 0.085 ! in MeV
 !---------------------------------------------------
 ! Azimuthal coordinate system: We use an azimuthal coordinate system where the domain \phi \in [-Pi, Pi],
 ! similar to the Heliographic Inertial Coordinate System (HGI), Earth is fixed at Earth_phi = 0. Towards the
 ! 'west' (i.e. direction of solar rotation) is positive, i.e. a source at W60deg is phi_source = 60./180*PI
 ! Check for indexing issues when spacecraft are located close to the edges of the azimuthal domain, i.e. 
 ! near phi = +- Pi
 
 Earth_phi = 33./180.*PI + 82./180.*PI !0/180.*PI
 Earth_r = 1.0
 
 Mars_phi = 169.0638746/180.*PI
 Mars_r = 1.610973075
 
 STA_phi = 33./180.*PI + 120./180.*PI ! 65./180.*PI
 STA_r = 1.0

 STB_phi = 33./180.*PI - 21./180.*PI !-71./180.*PI
 STB_r = 1.0
 
 Bepi_phi = 89.94839563/180.*PI
 Bepi_r = 0.411289483
 
 PSP_phi = -53.64955116/180.*PI
 PSP_r = 0.623189194
 
 SolO_phi = -2.819722862/180.*PI
 SolO_r = 0.803862906
 
 !---------------------------------------------------
 ! Files to write to
 
 OPEN(100,file='r_grid.txt',status='unknown')
 OPEN(200,file='phi_grid.txt',status='unknown')
 OPEN(300,file='pitch_angle_dependence.txt',status='unknown')
 OPEN(400,file='output.txt',status='unknown')
 OPEN(500,file='mu_grid.txt',status='unknown')
 OPEN(600,file='coeffs_mu.txt',status='unknown')
 OPEN(700,file='coeffs_r.txt',status='unknown') 
 OPEN(900,file='omni_phi.txt',status='unknown')
 OPEN(999,file='maximums.txt',status='unknown')

 OPEN(800,file='Earth_omni_time.txt',status='unknown') 
 OPEN(801,file='Mars_omni_time.txt',status='unknown')
 OPEN(802,file='Boundary_omni_time.txt',status='unknown') 
 OPEN(803,file='STA_omni_time.txt',status='unknown') 
 OPEN(804,file='Bepi_omni_time.txt',status='unknown') 
 OPEN(805,file='PSP_omni_time.txt',status='unknown') 
 OPEN(806,file='SolO_omni_time.txt',status='unknown') 
 OPEN(807,file='STB_omni_time.txt',status='unknown') 
 !---------------------------------------------------
  CALL CPU_TIME(time_begin)
 
    CALL omp_set_num_threads(4)
 
 !Set up the grids
 R(1) = 0.05
 delta_r = (3. - R(1))/(REAL(N) - 1)
 delta_phi = 2.*PI/REAL(M)
 delta_mu = 2./REAL(NK)

   MU(1) = -1. + Delta_mu/2.
  
     WRITE(500,"(1(ES18.8))") MU(1)
     
 DO k = 2, NK
 
  MU(k) = MU(k - 1) + Delta_mu
  
    WRITE(500,"(1(ES18.8))") MU(k)
    
 END DO
 
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
  onset_time(j) = 750.
 
 END DO
 
 DO i = 1, N
 
  DO j = 1, M
  
    DO k = 1, NK
    
      f_new(i,j,k) = 0.0001
      f_old(i,j,k) = 0.0001
   
   END DO
   
  END DO
  
END DO


CALL DEFINE_COEFFCIENTS(N,NK,R,A,B,C,D,delta_r,V_r,V_phi,K_rphi,MU,E,F,G,delta_mu,species,energy)

 !---------------------------------------------------
 ! Set up CFL conditions
 
 CFL_coeff = 0.9
 
 max_a = MAXVAL(A)
 max_b = MAXVAL(B)
 max_c = MAXVAL(C)
 max_d = MAXVAL(D)
 max_f = MAXVAL(F)
 max_g = MAXVAL(G)
 
 delta_t = CFL_coeff*MIN(ABS(0.5*delta_r*delta_r/max_b),ABS(0.5*delta_phi*delta_phi/max_d),ABS(0.5*delta_r/MAXVAL(V_r)), &
                        ABS(0.5*delta_phi/MAXVAL(V_phi)),ABS(0.5*delta_mu*delta_mu/max_f),ABS(0.5*delta_mu/max_g))
 !---------------------------------------------------
 ! Iterate in time

WRITE(*,*) '---------------------------------------------------'
WRITE(*,*) 'Start iteration....'
WRITE(*,*) '---------------------------------------------------'
 
 time = 0.
 time_printer = 0.
 

 DO WHILE (time.LT.time_max)

  time = time + delta_t
  
 !Time dependent boundary condition
  DO k = 1, NK
   DO j = 1, M
 
  f_old(1,j,k) = exp(-(PHI(j) - phi_source)*(PHI(j) - phi_source)/2./injection_broadness/injection_broadness) &
  /time*exp(-acceleration_time/time - time/escape_time)*Delta_t*1000000. + 0.0001
  f_new(1,j,k) = f_old(1,j,k)  
  
 END DO
 
END DO

!---------------------------------------------------
! Do the r integration - convection

    !$omp parallel do private(i,j,left_lim,right_lim,limiter)
DO k = 1, NK

  DO j = 1, M
  
  !Do the conversion to conserved flux; F = r*r*f
  
  DO i = 1, N
  
    f_old(i,j,k) = f_old(i,j,k)*R(i)*R(i)
  
  END DO
  
      DO i = 2, N - 1
      
	IF (V_r(i,k).GT.0.) THEN
	
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_r*(V_r(i,k)*f_old(i,j,k) - V_r(i-1,k)*f_old(i-1,j,k))
 
	ELSE
	
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/2./delta_r*(V_r(i+1,k)*f_old(i+1,j,k) - V_r(i,k)*f_old(i,j,k))
	
	END IF
	
      END DO
 
      DO i = 2, N - 1
      
	IF (V_r(i,k).GT.0.) THEN
      
      	  left_lim = (V_r(i,k)*f_new(i,j,k) - V_r(i-1,k)*f_new(i-1,j,k))/2.
	  right_lim = (V_r(i+1,k)*f_new(i+1,j,k) - V_r(i,k)*f_new(i,j,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
            !limiter = 0.
  
IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different

 	    f_holder(i,j,k) = V_r(i,k)*f_new(i,j,k) + limiter
	    
	  ELSE
	  
      	  left_lim = (V_r(i-1,k)*f_new(i-1,j,k) - V_r(i,k)*f_new(i,j,k))/2.
	  right_lim = (V_r(i,k)*f_new(i,j,k) - V_r(i+1,k)*f_new(i+1,j,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
            !limiter = 0.
  
IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
	  
	  f_holder(i,j,k) = V_r(i,k)*f_new(i,j,k) + limiter
	  	  
	  ENDIF
	  
	  f_holder(1,j,k) = V_r(1,k)*f_new(1,j,k)
	  f_holder(N,j,k) = V_r(N,k)*f_new(N,j,k)
 
	END DO
	
	DO i = 2, N - 1
 
	  IF (V_r(i,k).GT.0.) THEN
 
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_r*(f_holder(i,j,k) - f_holder(i-1,j,k))
	
	  ELSE
	  
	  f_new(i,j,k) = f_old(i,j,k) - delta_t/delta_r*(f_holder(i+1,j,k) - f_holder(i,j,k))
	  
	  ENDIF
	
	END DO
       
  END DO
  
!Update f_new and convert back to f =F/r/r

  DO i = 1, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)/R(i)/R(i)
    
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
	
	f_new(i,j,k) = f_old(i,j,k) + A(i,k)*delta_t/delta_r/2.*(f_old(i + 1,j,k) - f_old(i - 1,j,k)) + B(i,k)*delta_t/delta_r/delta_r& 
	*(f_old(i + 1,j,k) - 2.*f_old(i,j,k) + f_old(i - 1,j,k)) + K_rphi(i,k)/R(i)*delta_t/2./delta_r*(right_deriv - left_deriv)
     	      
      END DO

  END DO
  
!Update f_new

  DO i = 2, N
  
    DO j = 1, M
    
      f_old(i,j,k) = f_new(i,j,k)
    
    END DO
    
  END DO
 
 END DO
 !---------------------------------------------------
! Do the phi integration - convection

!$omp parallel do private(i,j,left_lim,right_lim,limiter)
DO k = 1, NK

  DO i = 2, N - 1
  
  IF (V_phi(i,k).LT.0.) THEN

    DO j = 2, M - 1

    	  f_new(i,j,k) = f_old(i,j,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,j+1,k) - f_old(i,j,k))
    
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
    f_new(i,1,k) = f_old(i,1,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,2,k) - f_old(i,1,k))
    
          left_lim = (f_new(i,M,k) - f_new(i,1,k))/2.
	  right_lim = (f_new(i,1,k) - f_new(i,2,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
    IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
    IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
    
    	    f_holder(i,1,k) = f_new(i,1,k) + limiter
    	    
    
    f_new(i,M,k) = f_old(i,M,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,1,k) - f_old(i,M,k)) 
    
          left_lim = (f_new(i,M-1,k) - f_new(i,M,k))/2.
	  right_lim = (f_new(i,M,k) - f_new(i,1,k))/2.
          limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
          
    IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
    IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
    
    	    f_holder(i,M,k) = f_new(i,M,k) + limiter
    
    DO j = 2, M - 1
 
	  f_new(i,j,k) = f_old(i,j,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,j+1,k) - f_holder(i,j,k))   
    
    END DO
    
    f_new(i,1,k) = f_old(i,1,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,2,k) - f_holder(i,1,k))
    f_new(i,M,k) = f_old(i,M,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,1,k) - f_holder(i,M,k))

  ELSE
  
     DO j = 2, M - 1

    	  f_new(i,j,k) = f_old(i,j,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,j,k) - f_old(i,j-1,k))
    
    END DO

    	  f_new(i,1,k) = f_old(i,1,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,1,k) - f_old(i,M,k))
    	  f_new(i,M,k) = f_old(i,M,k) - V_phi(i,k)*delta_t/2./delta_phi*(f_old(i,M,k) - f_old(i,M-1,k))
    
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
 
	  f_new(i,j,k) = f_old(i,j,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,j,k) - f_holder(i,j-1,k))   
    
    END DO
    
!Boundary equations
      
    f_new(i,1,k) = f_old(i,1,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,1,k) - f_holder(i,M,k))
    f_new(i,M,k) = f_old(i,M,k) - V_phi(i,k)*delta_t/delta_phi*(f_holder(i,M,k) - f_holder(i,M-1,k)) 
  
  ENDIF
    
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
      	
	  f_new(i,j,k) = f_old(i,j,k) + C(i,k)*delta_t/delta_phi/2.*(f_old(i,j + 1,k) - f_old(i,j - 1,k)) + & 
	  D(i,k)*delta_t/delta_phi/delta_phi*(f_old(i,j + 1,k) - 2.*f_old(i,j,k) + f_old(i,j - 1,k)) & 
	  + K_rphi(i,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
	      
      END DO
      
!Boundary equations
      
      left_deriv = (f_old(i+1,M,k) - f_old(i-1,M,k))/2./delta_r
      right_deriv = (f_old(i+1,2,k) - f_old(i-1,2,k))/2./delta_r 
      
      
    f_new(i,1,k) = f_old(i,1,k) + C(i,k)*delta_t/delta_phi/2.*(f_old(i,2,k) - f_old(i,M,k)) + D(i,k)*delta_t/delta_phi/delta_phi* & 
    (f_old(i,2,k) - 2.*f_old(i,1,k) + f_old(i,M,k)) + K_rphi(i,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
    
    
    left_deriv = (f_old(i+1,M-1,k) - f_old(i-1,M-1,k))/2./delta_r
    right_deriv = (f_old(i+1,1,k) - f_old(i-1,1,k))/2./delta_r 
      
    
     f_new(i,M,k) = f_old(i,M,k) + C(i,k)*delta_t/delta_phi/2.*(f_old(i,1,k) - f_old(i,M-1,k)) + & 
     D(i,k)*delta_t/delta_phi/delta_phi*(f_old(i,1,k) - 2.*f_old(i,M,k) + f_old(i,M-1,k)) + & 
     K_rphi(i,k)/R(i)*delta_t/2./delta_phi*(right_deriv - left_deriv)
      
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
DO i = 2, N

  DO j = 1, M
  
    DO k = 2, NK - 1
    
      f_new(i,j,k) = f_old(i,j,k) + E(i,k)*Delta_t/2./Delta_mu*(f_old(i,j,k+1) - f_old(i,j,k-1)) & 
      + F(i,k)*Delta_t/Delta_mu/Delta_mu*(f_old(i,j,k+1)- 2.*f_old(i,j,k) + f_old(i,j,k-1))
    
    END DO

!At the boundaries

  f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*(-(F(i,1)+F(i,2))/2.*(f_old(i,j,2) - f_old(i,j,1))/Delta_mu)
  f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu*(-(F(i,NK)+F(i,NK-1))/2.*(f_old(i,j,NK) - f_old(i,j,NK-1))/Delta_mu)
    
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

      f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(G(i,k)*f_old(i,j,k) - G(i,k-1)*f_old(i,j,k-1))/2.
    
    END DO
    
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu/2.*f_old(i,j,1)*G(i,1)
    f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu/2.*f_old(i,j,NK-1)*G(i,NK-1)    
    
    DO k = 2, NK - 1
    
    	  left_lim = (G(i,k)*f_new(i,j,k) - G(i,k-1)*f_new(i,j,k-1))/2.
	  right_lim = (G(i,k+1)*f_new(i,j,k+1) - G(i,k)*f_new(i,j,k))/2.
	  
	  limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
	  
	  IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN
	  
	  IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs are different
    
    f_holder(i,j,k) = G(i,k)*f_new(i,j,k) + limiter
    
    END DO
    
    f_holder(i,j,1) = G(i,1)*f_new(i,j,1)
    f_holder(i,j,NK) = G(i,NK)*f_new(i,j,NK)
    
    DO k = 2, NK - 1
    
      	f_new(i,j,k) = f_old(i,j,k) - Delta_t/Delta_mu*(f_holder(i,j,k) - f_holder(i,j,k-1))
    
    END DO
    
    f_new(i,j,1) = f_old(i,j,1) - Delta_t/Delta_mu*f_holder(i,j,1)
    f_new(i,j,NK) = f_old(i,j,NK) + Delta_t/Delta_mu*f_holder(i,j,NK-1)
    
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

time_printer = time_printer + delta_t

!---------------------------------------------------------------------------------------
IF (time_printer.GE.minute) THEN

  time_printer = time_printer - minute
  
 !-------------------------------------------------------------------------------------
 !Earth
 WRITE(*,*) 'Earth'
  indeks_phi = (Earth_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (Earth_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(800,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
 !Mars
  WRITE(*,*) 'Mars'
  indeks_phi = (Mars_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (Mars_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(801,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni
       
!-------------------------------------------------------------------------------------
!STA
 WRITE(*,*) 'STA '
  indeks_phi = (STA_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (STA_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(803,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
!STA
 WRITE(*,*) 'STB '
  indeks_phi = (STB_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (STB_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(807,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
!Bepi
 WRITE(*,*) 'Bepi '
  indeks_phi = (Bepi_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (Bepi_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(804,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
!PSP
 WRITE(*,*) 'PSP '
  indeks_phi = (PSP_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (PSP_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(805,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
!SolO
 WRITE(*,*) 'SolO'
  indeks_phi = (SolO_phi - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = (SolO_r - R(1) + delta_r)/delta_r + 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(806,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni

!-------------------------------------------------------------------------------------
!Boundary
  WRITE(*,*) 'Boundary'
  indeks_phi = (phi_source - phi(1) + delta_phi)/delta_phi + 1.
  indeks_r = 1.
  
  anis = 0.
  f_omni = 0.
  
  i = indeks_r
  j = indeks_phi
  
  WRITE(*,*) 'Time:', time
  WRITE(*,*) indeks_r, R(i), indeks_phi, phi(j)
  
  DO k = 1, NK
  
    f_omni = f_omni + f_old(i,j,k)
    anis = anis + mu(k)*f_old(i,j,k)
  
  END DO
  
       WRITE(802,"(5(ES18.8))") R(i), phi(j), time, f_omni*delta_mu/2., 3.*anis/f_omni
!-------------------------------------------------------------------------------------
END IF

!========================================================
! Write peak intensity vs phi_grid

  indeks_r = (1. - R(1) + delta_r)/delta_r + 1.
  
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
 CLOSE(801)
 CLOSE(802)
 CLOSE(803)
 CLOSE(804)
 CLOSE(805)
 CLOSE(806)
 CLOSE(807)
 CLOSE(800)
 CLOSE(900)
 CLOSE(999)
 
 END
 
 !----------------------------------------------------
 SUBROUTINE DEFINE_COEFFCIENTS(N,NK,R,A,B,C,D,delta_r,V_r,V_phi,K_rphi,MU,E,F,G,delta_mu,species,energy)

 IMPLICIT NONE
 
 INTEGER :: i, j, k, N, NK, species
 REAL(KIND=8) :: R(N), A(N,NK), B(N,NK), C(N,NK), D(N,NK), delta_r, delta_mu
 REAL(KIND=8) :: Omega, V_sw, tanpsi, cospsi(N), sinpsi(N), Kappa_perp(N,NK), V(NK)
 REAL(KIND=8) :: K_rr(N,NK), K_phiphi(N,NK), K_rphi(N,NK), K_rr_dr(N,NK), K_phiphi_dr(N,NK), K_rphi_dr(N,NK)
 REAL(KIND=8) :: V_r(N,NK), V_phi(N,NK), MU(NK), L(N), energy, F(N,NK), E(N,NK), G(N,NK), integral
 REAL(KIND=8) :: lambda_para(N), speed, lambda_perp(N)
 REAL(KIND=8) :: B_mag(N), B_0, rest_energy, rigidity, lambda(N), beta
 REAL(KIND=8), PARAMETER :: speed_of_light = 7.2093d0 !in units of AU/hr
 REAL(KIND=8), PARAMETER :: PI = 3.14159
 !------------------------------------------------------- 
  IF (species.EQ.1) THEN
 
  rest_energy = 0.51 !Electron rest mass in MeV
 
 ELSE
 
  rest_energy = 938.27 !Proton rest mass in MeV
 
 ENDIF
 
 rigidity = SQRT(energy*(energy + 2.*rest_energy)) ! in MV
 beta = rigidity/(energy + rest_energy)
 speed = beta*speed_of_light
 
 Omega = 2.*PI/(25.38*24.) !/hr
 V_sw = 400./1.5d8*(60*60) !AU/hr
 
   !------------------------------------------
   ! Transformation angles and magn field
  DO i = 1, N
  
    tanpsi = (R(i) - 0.05)*Omega/V_sw ! no units, of course
    cospsi(i) = 1./SQRT(1. + tanpsi*tanpsi)
    sinpsi(i) = SQRT(1. - cospsi(i)*cospsi(i))
    
	B_0 = 5./SQRT(1. + Omega/V_sw*Omega/V_sw) ! reference value at B_0 in nT
    
	B_mag(i) =  B_0/R(i)/R(i)*SQRT(1. + tanpsi*tanpsi) ! magnetic field magnitude in nT
    	
  END DO
!------------------------------------------
! Focussing length

DO i = 1, N
 
  L(i) = R(i)/(2. + omega*omega/V_sw/V_sw*R(i)*R(i))*(1. + omega*omega/V_sw/V_sw*R(i)*R(i))**(3./2.) ! no units
  
END DO

!------------------------------------------
! Parallel mean-free-path

DO i = 1, N
 
 lambda(i) = 0.1
 lambda_perp(i) = 0.025*lambda(i)
  
 END DO


DO i = 1, N
 
  DO k = 1, NK
  
    F(i,k) = (1. - mu(k)*mu(k))*((ABS(MU(k)))**(1.67 - 1.) + 0.5)

  END DO
  
 END DO
 
 integral = 0.
 
  DO k = 1, NK
    
    i = 1
    
    integral = integral + (1. - mu(k)*mu(k))*(1. - mu(k)*mu(k))/F(i,k)*delta_mu
  
  END DO 
 
 DO i = 1, N
 
  DO k = 1, NK
  
    F(i,k) = F(i,k)*integral*3.*speed/8./lambda(i)
  
  END DO
  
 END DO

!------------------------------------------
 !Calculate lambda_para for visualization
 
 DO i = 1, N
 
    integral = 0.
 
  DO k = 1, NK
    
    integral = integral + (1. - mu(k)*mu(k))*(1. - mu(k)*mu(k))/F(i,k)
   
  END DO

 integral = integral*ABS(mu(1) - mu(2))
 lambda_para(i) = 3.*speed/8.*integral
 
END DO
  !------------------------------------------ 
 ! Perpendicular coefficient
 
 DO i = 1, N
 
   DO k = 1, NK

    Kappa_perp(i,k) = speed/3.*lambda_perp(i)*ABS(MU(k))*2.
  
  END DO
  
 END DO
  
  !Calculate lambda_perp for visualization
 
 DO i = 1, N
 
 integral = 0.
 
  DO k = 1, NK
  
    integral = integral + Kappa_perp(i,k)
  
  END DO
  
  lambda_perp(i) = 0.5*integral*ABS(mu(1) - mu(2))*3./speed
 
END DO
   !------------------------------------------
 ! Derivative of D_mu-mu
 
 DO i = 1, N
 
  DO k = 2, NK - 1
  
    E(i,k) = (F(i,k+1) - F(i,k-1))/2./delta_mu
  
  END DO

    E(i,1) = (F(i,2) - F(i,1))/delta_mu
    E(i,NK) = (F(i,NK) - F(i,NK - 1))/delta_mu
  
 END DO
 
 ! Focusing term
 DO i = 1, N
 
  DO k = 1, NK
  
    G(i,k) = speed/2./L(i)*(1. - mu(k)*mu(k))
  
  END DO
  
 END DO
 
 ! Advection term
 DO k = 1, NK
 
 V(k) = speed*MU(k)
 
 END DO
   !------------------------------------------
! Diffusion tensor in spiral coordinates

  DO i = 1, N
  
   DO k = 1, NK
  
    K_rr(i,k) = Kappa_perp(i,k)*sinpsi(i)*sinpsi(i)
    K_phiphi(i,k) = Kappa_perp(i,k)*cospsi(i)*cospsi(i)
    K_rphi(i,k) = Kappa_perp(i,k)*cospsi(i)*sinpsi(i)
    
   END DO
    
  END DO
  
 !The derivatives of the tensor
 
 DO k = 1, NK
 
  DO i = 2, N - 1
  
    K_rr_dr(i,k) = (K_rr(i + 1,k) - K_rr(i - 1,k))/2./delta_r
    K_phiphi_dr(i,k) = (K_phiphi(i + 1,k) - K_phiphi(i - 1,k))/2./delta_r
    K_rphi_dr(i,k) = (K_rphi(i + 1,k) - K_rphi(i - 1,k))/2./delta_r
    
  END DO
 
 
  K_rr_dr(1,k) = (K_rr(2,k) - K_rr(1,k))/delta_r
  K_rr_dr(N,k) = (K_rr(N,k) - K_rr(N - 1,k))/delta_r
  
  K_phiphi_dr(1,k) = (K_phiphi(2,k) - K_phiphi(1,k))/delta_r
  K_phiphi_dr(N,k) = (K_phiphi(N,k) - K_phiphi(N - 1,k))/delta_r
  
  K_rphi_dr(1,k) = (K_rphi(2,k) - K_rphi(1,k))/delta_r
  K_rphi_dr(N,k) = (K_rphi(N,k) - K_rphi(N - 1,k))/delta_r
 
 END DO
 
! Define numerical coeffs
DO k = 1, NK

 DO i = 1, N
 
    A(i,k) = 2.*K_rr(i,k)/R(i) + K_rr_dr(i,k)
    B(i,k) = K_rr(i,k)
    C(i,k) = K_rphi(i,k)/R(i)/R(i) + K_rphi_dr(i,k)/R(i)
    D(i,k) = K_phiphi(i,k)/R(i)/R(i)
  
    V_r(i,k) = V(k)*cospsi(i) + V_sw
    V_phi(i,k) = -V(k)*sinpsi(i)/R(i)

END DO

END DO

!Write some stuff to a file

DO i = 1, N

  WRITE(700,"(7(ES18.8))") R(i), sinpsi(i), cospsi(i), L(i), lambda_para(i), lambda_perp(i), B_mag(i)

END DO

  i = (1. - R(1) + delta_r)/delta_r + 1.

DO k = 1, NK

  WRITE(600,"(4(ES18.8))") MU(k), F(i,k), E(i,k), Kappa_perp(i,k)

END DO

WRITE(*,*) '---------------------------------------------------'
WRITE(*,*) 'Print stuff to file...'
WRITE(*,*) '---------------------------------------------------'

RETURN

END
