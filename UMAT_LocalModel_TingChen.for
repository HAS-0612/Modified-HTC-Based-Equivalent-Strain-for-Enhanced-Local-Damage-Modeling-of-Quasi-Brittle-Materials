!=======================================================================
!	Hsieh-Ting-Chen equivalent strain
!	For 2D and 3D analyses
!	 Developer: Hemam Amarjit Singh
!	Supervisor: Dr. Rimen Jamatia
! 	© Indian Institute of Technology Jammu
!=======================================================================
! 	 Properties
!       props(1) = Youngs Modulus
!       props(2) = Poisson's Ratio
!       props(3) = Compressive strength to Tensile strength Ratio
!       props(4) = Damage Initiation threshold
!       props(5) = Tensile Strength
!       props(6) = Fracture Energy
!  
! 	 The state variables are stored as:
!		stateNew(*,1) = kappa
!		stateNew(*,2) = Damage			
!		stateNew(*,3) =	coeff_damage	                                        
!=======================================================================
C=======================================================================
		SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,ELEM_COORDS,AREA)

C-----------------------------------------------------------------------
C   INCLUDE AND VARIABLE DECLARATIONS
C-----------------------------------------------------------------------
		INCLUDE 'ABA_PARAM.INC'

		CHARACTER*80 CMNAME
		DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(2,4),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
	 
		DOUBLE PRECISION CC(NTENS,NTENS),delta(NTENS),cos3theta, theta,MaxPS,
     $ II1,JJ2,JJ3,II1_d(NTENS),JJ2_d(NTENS),NODEE(1),f_t(1),G_f(1),JJ3_d(NTENS),
     $ AA,e_tilde1,kappa,kappa_n,damage,ee,gg,coeff_damage,coeff_ee(NTENS),
     $ a1,a2,a4(NTENS),Betaa(1),ELEM_COORDS(2,4),MaxPS_d(NTENS),ConstTERM,
     $ STRESS_EQ(NTENS),CCC(4,4),X(4),Y(4),ee_dII1,ee_dJJ2,ee_dMaxPS,costheta_d(NTENS),
     $ dk_dee(1);

C   MATERIAL PROPERTIES
		DOUBLE PRECISION E,NU,k,kappa_o,ft,Gf
C   PARAMETERS AND CONSTANTS
		INTEGER i, j, a, b       
		
		PARAMETER(Zero=0d0, One=1.0d0, Two=2.0d0, Three=3.0d0, 
     $ Six=6.0d0, TOL=1d-25, Alpha=1.0d0, Beta=305.0d0, Ne=3576.0d0,
     $ Eta=5.0d0) 

C-----------------------------------------------------------------------
C   EXTRACT PROPERTIES FROM INPUT
C-----------------------------------------------------------------------
        E=PROPS(1)          ! Young's modulus
        NU=PROPS(2)         ! Poisson's ratio
		k=PROPS(3)          ! Ratio f_t/f_c 
		kappa_o=props(4)    ! Damage initiation threshold
		ft = props(5)
		Gf = props(6)
		! print*, 'TERM = ', One
		ConstA=7.6924
		ConstB=6.3689
		ConstC=0.8252
		ConstD=5.2412
C-----------------------------------------------------------------------
C   INITIALIZATION
C-----------------------------------------------------------------------
		CC=0d0
		delta=0d0
		DDSDDE=0d0

C	TOTAL STRAIN 
		DO a=1,NTENS		
			STRAN(a)=STRAN(a)+DSTRAN(a)	
		END DO

C-----------------------------------------------------------------------
C   ELASTIC STIFFNESS MATRIX AND INVARIANTS (PLANE STRAIN)
C-----------------------------------------------------------------------
		IF (NDI==3 .and. NSHR==1) THEN
		! Elastic stiffness CC       
			DO i=1,NDI
            	DO j=1,NDI
               		CC(i,j)=NU
            	END DO
        	END DO

        	DO i=1,NDI
				CC(i,i)=1.d0-NU
        	END DO

			CC(4,4)=(One-Two*NU)/Two                   		
			CC=CC*E/((One+NU)*(One-Two*NU))

		! Kronecker delta
			delta(1)=One
        	delta(2)=One
        	delta(3)=One
        	delta(4)=Zero

		! Invariants II1&JJ2
			II1=0d0
        	DO a=1,NTENS
				II1=II1+STRAN(a)*delta(a)
			END DO
		JJ2 = STRAN(1)**2/three + STRAN(2)**2/three
     &           - STRAN(1)*STRAN(2)/three + (STRAN(4)/two)**2
		JJ3 = (2.d0/27.d0)*(STRAN(1)**3 + STRAN(2)**3 + STRAN(3)**3)
     &     	   - (1.d0/9.d0)*(STRAN(1)**2*(STRAN(2)+STRAN(3)) 
     &              + STRAN(2)**2*(STRAN(1)+STRAN(3)) 
     &              + STRAN(3)**2*(STRAN(1)+STRAN(2)))
     &         + (4.d0/9.d0)*STRAN(1)*STRAN(2)*STRAN(3)
     &         - (1.d0/3.d0)*((2.d0*STRAN(3) - STRAN(1) - STRAN(2))*(STRAN(4)/two)**2)
		IF (JJ2 .LT. 1D-40) then
		   cos3theta = 0.d0
		ELSE
		   cos3theta = 1.5d0 * sqrt(three) * (JJ3 / (JJ2**1.5d0))
		ENDIF
		IF (cos3theta .ge. 1.0d0) cos3theta = 1.0d0
		IF (cos3theta .le. -1.0d0) cos3theta = -1.0d0
		theta = (one/three) * acos(cos3theta)
		MaxPS = (II1/three) + (sqrt(4.d0/3.d0)*sqrt(JJ2)*cos(theta))
		IF (MaxPS .lt. 0.d0) MaxPS = 0.d0

			II1_d(1)=One
			II1_d(2)=One
			II1_d(3)=One
			II1_d(4)=Zero
			
			JJ2_d(1)=STRAN(1)-II1*II1_d(1)/3.0d0
			JJ2_d(2)=STRAN(2)-II1*II1_d(2)/3.0d0
			JJ2_d(3)=STRAN(3)-II1*II1_d(3)/3.0d0
			JJ2_d(4)=(STRAN(4)/two)-II1*II1_d(4)/3.0d0
			
			JJ3_d(1) = (2d0/9d0)*STRAN(1)*STRAN(1)-
     &         - (2d0/9d0)*STRAN(1)*(STRAN(2)+STRAN(3)) 
     &         - (1d0/9d0)*(STRAN(2)*STRAN(2)+STRAN(3)*STRAN(3)) 
     &         + (4d0/9d0)*STRAN(2)*STRAN(3)
     &         + (1d0/12d0)*STRAN(4)
			JJ3_d(2) = (2d0/9d0)*STRAN(2)*STRAN(2) 
     &         - (2d0/9d0)*STRAN(2)*(STRAN(1)+STRAN(3)) 
     &         - (1d0/9d0)*(STRAN(1)*STRAN(1)+STRAN(3)*STRAN(3)) 
     &         + (4d0/9d0)*STRAN(1)*STRAN(3)
     &         + (1d0/12d0)*STRAN(4)
			JJ3_d(3) = (2d0/9d0)*STRAN(3)*STRAN(3) 
     &         - (2d0/9d0)*STRAN(3)*(STRAN(1)+STRAN(2)) 
     &         - (1d0/9d0)*(STRAN(1)*STRAN(1)+STRAN(2)*STRAN(2)) 
     &         + (4d0/9d0)*STRAN(1)*STRAN(2)
     &         - (1d0/6d0)*STRAN(4)
			JJ3_d(4) = -(STRAN(4)/6d0)*(2d0*STRAN(3) - STRAN(1) - STRAN(2))
			

			costheta_d(1)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(1))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(1)))
			costheta_d(2)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(2))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(2)))
			costheta_d(3)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(3))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(3)))
			costheta_d(4)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(4))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(4)))
			
			MaxPS_d(1)=(II1_d(1)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(1))+(sqrt(Two*JJ2/Three)*costheta_d(1))
			MaxPS_d(2)=(II1_d(2)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(2))+(sqrt(Two*JJ2/Three)*costheta_d(2))
			MaxPS_d(3)=(II1_d(3)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(3))+(sqrt(Two*JJ2/Three)*costheta_d(3))
			MaxPS_d(4)=(II1_d(4)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(4))+(sqrt(Two*JJ2/Three)*costheta_d(4))
C
C-----------------------------------------------------------------------
C   EQUIVALENT STRAIN (MODIFIED Hsieh-Ting-Chen)
C-----------------------------------------------------------------------
		ConstTERM = (ConstB*sqrt(JJ2)/(One+NU))+(ConstD*II1/(One-Two*NU))+(ConstC*MaxPS/E)
		AA=((ConstTERM)**Two)+(4.0*ConstA*JJ2/((One+NU)**Two)) 
		ee=((ConstTERM)+(sqrt(AA)))/(Two*k)			
			ee_dII1=(ConstD/(Two*k*(One-Two*NU)))*(One+(ConstTERM/sqrt(AA)))
			ee_dJJ2=(ConstB/(Two*k*(One+NU)*sqrt(JJ2)))*(One+(ConstTERM/sqrt(AA)))
			ee_dMaxPS=(ConstC/(Two*k*E))*(One+(ConstTERM/sqrt(AA)))
			
		! derivative of eq. strain w.r.t strain tensor
			IF (ee .lt. 1d-14) THEN
			a4(1)=0d0
			a4(2)=0d0
			a4(3)=0d0
			a4(4)=0d0
			ELSE
			a4(1)=(ee_dII1*II1_d(1))+(ee_dJJ2*JJ2_d(1))+(ee_dMaxPS*MaxPS_d(1))
			a4(2)=(ee_dII1*II1_d(2))+(ee_dJJ2*JJ2_d(2))+(ee_dMaxPS*MaxPS_d(2))
			a4(3)=(ee_dII1*II1_d(3))+(ee_dJJ2*JJ2_d(3))+(ee_dMaxPS*MaxPS_d(3))
			a4(4)=(ee_dII1*II1_d(4))+(ee_dJJ2*JJ2_d(4))+(ee_dMaxPS*MaxPS_d(4))						
			ENDIF
		ENDIF
C-----------------------------------------------------------------------
C   ELASTIC STIFFNESS MATRIX AND INVARIANTS (3D)
C-----------------------------------------------------------------------		
		IF (NDI==3 .and. NSHR==3) THEN
		! Elastic stiffness CC      
			CC(1,1)=1.d0-NU
			CC(1,2)=NU
			CC(1,3)=NU
			CC(2,1)=NU
			CC(2,2)=1.d0-NU
			CC(2,3)=NU
			CC(3,1)=NU
			CC(3,2)=NU
			CC(3,3)=1.d0-NU
			CC(4,4)=(1.d0-2.d0*NU)/Two
			CC(5,5)=(1.d0-2.d0*NU)/Two
			CC(6,6)=(1.d0-2.d0*NU)/Two
			CC=CC*E/((1.d0+NU)*(1.d0-2.d0*NU))	

		! Kronecker delta
			delta(1)=One
			delta(2)=One
			delta(3)=One
			delta(4)=Zero		
			delta(5)=Zero
			delta(6)=Zero		
			
		! Invariants II1&JJ2		
			II1=0d0
			DO a=1,NTENS
				II1 = II1+STRAN(a)*delta(a)
			END DO

            JJ2 = (STRAN(1)**2 + STRAN(2)**2 + STRAN(3)**2
     &         - STRAN(1)*STRAN(2) - STRAN(2)*STRAN(3) 
     &         - STRAN(3)*STRAN(1))/three
     &         + (one/4.d0)*(STRAN(4)**2 + STRAN(5)**2 + STRAN(6)**2)
			JJ3 = (2.d0/27.d0)*(STRAN(1)**3 + STRAN(2)**3 + STRAN(3)**3)
     &     	   - (1.d0/9.d0)*(STRAN(1)**2*(STRAN(2)+STRAN(3)) 
     &              + STRAN(2)**2*(STRAN(1)+STRAN(3)) 
     &              + STRAN(3)**2*(STRAN(1)+STRAN(2)))
     &         + (4.d0/9.d0)*STRAN(1)*STRAN(2)*STRAN(3)
     &         - (1.d0/3.d0)*((2.d0*STRAN(1) - STRAN(2) - STRAN(3))*(STRAN(6)/two)**2
     &              + (2.d0*STRAN(2) - STRAN(1) - STRAN(3))*(STRAN(5)/two)**2
     &              + (2.d0*STRAN(3) - STRAN(1) - STRAN(2))*(STRAN(4)/two)**2)
     &         + 2.d0*(STRAN(4)/two)*(STRAN(5)/two)*(STRAN(6)/two)

            if (JJ2 .LT. 1D-40) then
               cos3theta = 0.d0
            else
               cos3theta = 1.5d0 * sqrt(three) * (JJ3 / (JJ2**1.5d0))
            end if
            if (cos3theta .ge. 1.0d0) cos3theta = 1.0d0
            if (cos3theta .le. -1.0d0) cos3theta = -1.0d0
            theta = (one/three) * acos(cos3theta)
            MaxPS = (II1/three) + (sqrt(4.d0/3.d0)*sqrt(JJ2)*cos(theta))
            if (MaxPS .lt. 0.d0) MaxPS = 0.d0

			II1_d(1)=One
			II1_d(2)=One
			II1_d(3)=One
			II1_d(4)=Zero
			II1_d(5)=Zero
			II1_d(6)=Zero
			
			JJ2_d(1)=STRAN(1)-II1*II1_d(1)/3.0d0
			JJ2_d(2)=STRAN(2)-II1*II1_d(2)/3.0d0
			JJ2_d(3)=STRAN(3)-II1*II1_d(3)/3.0d0
			JJ2_d(4)=(STRAN(4)/two)-II1*II1_d(4)/3.0d0
			JJ2_d(5)=(STRAN(5)/two)-II1*II1_d(5)/3.0d0
			JJ2_d(6)=(STRAN(6)/two)-II1*II1_d(6)/3.0d0	

			JJ3_d(1) = (2d0/9d0)*STRAN(1)*STRAN(1) 
     &         - (2d0/9d0)*STRAN(1)*(STRAN(2)+STRAN(3)) 
     &         - (1d0/9d0)*(STRAN(2)*STRAN(2)+STRAN(3)*STRAN(3)) 
     &         + (4d0/9d0)*STRAN(2)*STRAN(3)
     &         - (1d0/6d0)*STRAN(6)*STRAN(6) + (1d0/12d0)*(STRAN(5)*STRAN(5) + STRAN(4)*STRAN(4))
			JJ3_d(2) = (2d0/9d0)*STRAN(2)*STRAN(2) 
     &         - (2d0/9d0)*STRAN(2)*(STRAN(1)+STRAN(3)) 
     &         - (1d0/9d0)*(STRAN(1)*STRAN(1)+STRAN(3)*STRAN(3)) 
     &         + (4d0/9d0)*STRAN(1)*STRAN(3)
     &         - (1d0/6d0)*STRAN(5)*STRAN(5) + (1d0/12d0)*(STRAN(6)*STRAN(6) + STRAN(4)*STRAN(4))
			JJ3_d(3) = (2d0/9d0)*STRAN(3)*STRAN(3) 
     &         - (2d0/9d0)*STRAN(3)*(STRAN(1)+STRAN(2)) 
     &         - (1d0/9d0)*(STRAN(1)*STRAN(1)+STRAN(2)*STRAN(2)) 
     &         + (4d0/9d0)*STRAN(1)*STRAN(2)
     &         - (1d0/6d0)*STRAN(4)*STRAN(4) + (1d0/12d0)*(STRAN(5)*STRAN(5) + STRAN(6)*STRAN(6))
			JJ3_d(4) = -(STRAN(4)/6d0)*(2d0*STRAN(3) - STRAN(1) - STRAN(2)) + (STRAN(5)*STRAN(6))/4d0
			JJ3_d(5) = -(STRAN(5)/6d0)*(2d0*STRAN(2) - STRAN(1) - STRAN(3)) + (STRAN(4)*STRAN(6))/4d0
			JJ3_d(6) = -(STRAN(6)/6d0)*(2d0*STRAN(1) - STRAN(2) - STRAN(3)) + (STRAN(4)*STRAN(5))/4d0
	
			costheta_d(1)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(1))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(1)))
			costheta_d(2)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(2))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(2)))
			costheta_d(3)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(3))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(3)))
			costheta_d(4)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(4))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(4)))
			costheta_d(5)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(5))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(5)))
			costheta_d(6)=(1.5*sqrt(Three)/JJ2**3)*((sqrt(JJ2**3)*JJ3_d(6))-(1.5*JJ3*sqrt(JJ2)*JJ2_d(6)))
			
			MaxPS_d(1)=(II1_d(1)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(1))+(sqrt(Two*JJ2/Three)*costheta_d(1))
			MaxPS_d(2)=(II1_d(2)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(2))+(sqrt(Two*JJ2/Three)*costheta_d(2))
			MaxPS_d(3)=(II1_d(3)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(3))+(sqrt(Two*JJ2/Three)*costheta_d(3))
			MaxPS_d(4)=(II1_d(4)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(4))+(sqrt(Two*JJ2/Three)*costheta_d(4))
			MaxPS_d(5)=(II1_d(5)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(5))+(sqrt(Two*JJ2/Three)*costheta_d(5))
			MaxPS_d(6)=(II1_d(6)/Three)+((cos(theta)/(sqrt(6.*JJ2)))*JJ2_d(6))+(sqrt(Two*JJ2/Three)*costheta_d(6))
C
C-----------------------------------------------------------------------
C   EQUIVALENT STRAIN (MODIFIED Hsieh-Ting-Chen)
C-----------------------------------------------------------------------
		ConstTERM = (ConstB*sqrt(JJ2)/(One+NU))+(ConstD*II1/(One-Two*NU))+(ConstC*MaxPS/E)
		AA=((ConstTERM)**Two)+(4.0*ConstA*JJ2/((One+NU)**Two)) 
		ee=((ConstTERM)+(sqrt(AA)))/(Two*k)			
			ee_dII1=(ConstD/(Two*k*(One-Two*NU)))*(One+(ConstTERM/sqrt(AA)))
			ee_dJJ2=(ConstB/(Two*k*(One+NU)*sqrt(JJ2)))*(One+(ConstTERM/sqrt(AA)))
			ee_dMaxPS=(ConstC/(Two*k*E))*(One+(ConstTERM/sqrt(AA)))
			
		! derivative of eq. strain w.r.t strain tensor
			IF (ee .lt. 1d-14) THEN
			a4(1)=0d0
			a4(2)=0d0
			a4(3)=0d0
			a4(4)=0d0
			a4(5)=0d0
			a4(6)=0d0
			ELSE
			a4(1)=(ee_dII1*II1_d(1))+(ee_dJJ2*JJ2_d(1))+(ee_dMaxPS*MaxPS_d(1))
			a4(2)=(ee_dII1*II1_d(2))+(ee_dJJ2*JJ2_d(2))+(ee_dMaxPS*MaxPS_d(2))
			a4(3)=(ee_dII1*II1_d(3))+(ee_dJJ2*JJ2_d(3))+(ee_dMaxPS*MaxPS_d(3))
			a4(4)=(ee_dII1*II1_d(4))+(ee_dJJ2*JJ2_d(4))+(ee_dMaxPS*MaxPS_d(4))
			a4(5)=(ee_dII1*II1_d(5))+(ee_dJJ2*JJ2_d(5))+(ee_dMaxPS*MaxPS_d(5))
			a4(6)=(ee_dII1*II1_d(6))+(ee_dJJ2*JJ2_d(6))+(ee_dMaxPS*MaxPS_d(6))		
			ENDIF
		ENDIF

C   HISTORY VARIABLE kappa_n
		IF (Time(2) .eq. 0d0) THEN
		kappa_n=kappa_o
		ELSE
		kappa_n=STATEV(1)
		ENDIF
		STATEV(1)=kappa_n

C   UPDATE HISTORY PARAMETER kappa
        IF (kappa_n .lt. ee) THEN
			kappa=ee        
			STATEV(1)=ee
			dk_dee(1)=One
        ELSE
			kappa=kappa_n
			dk_dee(1)=Zero
	    ENDIF


C-----------------------------------------------------------------------
C   DAMAGE VARIABLE AND DERIVATIVES
C-----------------------------------------------------------------------
			NODEE=CELENT
		Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_0 * ft * NODEE)
        damage=One-(kappa_o/kappa)*((One-Alpha)+Alpha*(exp(-Betaa(1)*(kappa-kappa_o))))		

		coeff_damage=(kappa_o/(kappa**Two))*(One-Alpha
     $ 	+ Alpha*(exp(-Betaa(1)*(kappa-kappa_o))))+(kappa_o/kappa)*(Alpha*
     $ (exp(-Betaa(1)*(kappa-kappa_o)))*Betaa(1))

C-----------------------------------------------------------------------
C   STRESS UPDATE (WITH DAMAGE)
C-----------------------------------------------------------------------
		STRESS_EQ=0d0		! eq. stress (without damage)	
		STRESS=0d0 			! stress with damage
		DO i=1, NTENS
			DO j=1, NTENS
				STRESS_EQ(i)=STRESS_EQ(i)+CC(i,j)*STRAN(j)
				STRESS(i)=(One-damage)*STRESS_EQ(i)
			END DO
		END DO


C-----------------------------------------------------------------------
C   COMPUTE DDSDDE (TANGENT MODULUS)
C-----------------------------------------------------------------------

	    CCC = 0d0
		DO i = 1, NTENS
			DO j = 1, NTENS
				CCC(i,j)=STRESS_EQ(i)*a4(j)
			END DO
		END DO

		DO i=1, NTENS
			DO j=1, NTENS
				DDSDDE(i,j)=(1-damage)*CC(i,j)-coeff_damage*dk_dee(1)*CCC(i,j)
			END DO
		END DO	

C	STORE
		!STATEV(1)=kappa_n
		STATEV(2)=damage
		STATEV(3)=coeff_damage
		RETURN
    	END
