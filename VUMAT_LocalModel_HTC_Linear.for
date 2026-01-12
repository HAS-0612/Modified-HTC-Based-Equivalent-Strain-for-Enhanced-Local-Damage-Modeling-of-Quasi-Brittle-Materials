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
!		  #For 3D				 #For Plane Strain		 #For Plane Stress					
!		stateNew(*,3) = E11		stateNew(*,3) = E11 	stateNew(*,3) = E11                 
!		stateNew(*,4) = E22		stateNew(*,4) = E22 	stateNew(*,4) = E22                 
!		stateNew(*,5) = E33		stateNew(*,5) = E33 	stateNew(*,5) = E12                 
!		stateNew(*,6) = E12		stateNew(*,6) = E12 			            		        
!		stateNew(*,7) = E23				                                        
!		stateNew(*,8) = E13				                                        
!=======================================================================
      subroutine vumat(
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      character*80 cmname
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*), charLength(nblock),
     2  strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr),
     3  tempOld(nblock), stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+2*nshr), fieldOld(nblock,nfieldv),
     5  stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev),
     6  enerInternOld(nblock), enerInelasOld(nblock),
     7  tempNew(nblock), stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+2*nshr), fieldNew(nblock,nfieldv),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)

      double precision II1, JJ2, JJ3, cos3theta, theta, MaxPS
      double precision ConstA, ConstD, ConstC, ConstB, ConstTERM, AA
      double precision kappa, kappa_n, damage, ee, Betaa, kappa_o
      double precision totalStrain(6), s(6), E33, E33total
      double precision E, NU, k, ft, Gf, NODEE, temp_work, W1
      double precision trace, alamda, twomu
      integer i, j

      parameter( zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0,
     1  third = one/three, half = .5d0, twoThirds = two/three,
     2  threeHalfs = 1.5d0 )

      E        = props(1)
      NU       = props(2)	  
      k        = props(3)
      kappa_o  = props(4)
      ft       = props(5)
      Gf       = props(6)
	  W1	   = Gf/ft
		ConstA=7.6924
		ConstB=6.3689
		ConstC=0.8252
		ConstD=5.2412
!		ConstA=6.30933
!		ConstB=5.72027
!		ConstC=2.08138
!		ConstD=4.40571
	  Alpha = one
      twomu  = E / ( one + NU )
      alamda = twomu * NU / (one - two*NU)
C
      do 100 i = 1, nblock
C ---- 2D Plane Stress (Voigt: 11,22,12) ----
      if ((ndir.eq.2) .and. (nshr.eq.1)) then
         if (stepTime .eq. 0.d0) then
            trace  = strainInc(i,1) + strainInc(i,2)
            E33 = -alamda * trace / (alamda + twomu)

            stressNew(i,1) = stressOld(i,1) + alamda*(trace + E33) + twomu*strainInc(i,1)
            stressNew(i,2) = stressOld(i,2) + alamda*(trace + E33) + twomu*strainInc(i,2)
            stressNew(i,3) = stressOld(i,3) + twomu*strainInc(i,3)

            kappa  = kappa_o
            damage = 0.d0
            totalStrain(1) = strainInc(i,1)
            totalStrain(2) = strainInc(i,2)
            totalStrain(3) = strainInc(i,3)
            E33total = E33
         else
            kappa_n = stateOld(i,1)
            totalStrain(1) = stateOld(i,3) + strainInc(i,1)
            totalStrain(2) = stateOld(i,4) + strainInc(i,2)
            totalStrain(3) = stateOld(i,5) + strainInc(i,3)

            trace  = totalStrain(1) + totalStrain(2)
            E33total = -alamda * trace / (alamda + twomu)

            II1 = totalStrain(1) + totalStrain(2) + E33total

            JJ2 = totalStrain(1)**2/three + totalStrain(2)**2/three
     &            - totalStrain(1)*totalStrain(2)/three + (totalStrain(3))**2

            JJ3 = (2.d0/27.d0)*(totalStrain(1)**3 + totalStrain(2)**3 + E33total**3)
     &         - (1.d0/9.d0)*(totalStrain(1)**2*(totalStrain(2)+E33total) 
     &         + totalStrain(2)**2*(totalStrain(1)+E33total) 
     &         + E33total**2*(totalStrain(1)+totalStrain(2)))
     &         + (4.d0/9.d0)*totalStrain(1)*totalStrain(2)*E33total
     &         - (1.d0/3.d0)*((2.d0*E33total - totalStrain(1) - totalStrain(2))*totalStrain(3)**2)

            if (JJ2 .LT. 1D-40) then
               cos3theta = 0.d0
            else
               cos3theta = 1.5d0 * sqrt(three) * (JJ3 / (JJ2**1.5d0))
            end if

            if (cos3theta .gt. 1.d0) cos3theta = 1.d0
            if (cos3theta .lt. -1.d0) cos3theta = -1.d0

            theta = one/three * acos(cos3theta)
            MaxPS = (II1/three) + sqrt(4.d0/3.d0)*sqrt(JJ2)*cos(theta)
            if (MaxPS .lt. 0.d0) MaxPS = 0.d0

            ConstTERM = (ConstB*sqrt(JJ2)/(one+NU))
     &                  +(ConstD*II1/(one-two*NU))+(ConstC*MaxPS/E)
            AA = ConstTERM**2 + 4.d0*ConstA*JJ2/((one+NU)**2)
            ee = (ConstTERM + sqrt(AA)) / (two*k)

            if (kappa_n .lt. ee) then
               kappa = ee
            else
               kappa = kappa_n
            end if

            if (kappa .gt. kappa_o) then
               NODEE = charLength(i)
               Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_o * ft * NODEE)
				damage = One - (kappa_o / kappa)*(One-0.5*Betaa*(kappa-kappa_o))
			   if (damage .ge. one) then
					damage = 0.99
			   end if			   
					
            else
               damage = 0.d0
            end if

            s(1) = (one-damage) * (alamda*(trace + E33total) + twomu*totalStrain(1))
            s(2) = (one-damage) * (alamda*(trace + E33total) + twomu*totalStrain(2))
            s(3) = (one-damage) * twomu*totalStrain(3)

            do j = 1,3
               stressNew(i,j) = s(j)
            end do
         end if

         temp_work = 0.5d0*(
     & (stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     & + (stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     & + 2.0d0*(stressOld(i,3)+stressNew(i,3))*strainInc(i,3))

         stateNew(i,1) = kappa
         stateNew(i,2) = damage
         stateNew(i,3) = totalStrain(1)
         stateNew(i,4) = totalStrain(2)
         stateNew(i,5) = totalStrain(3)

         enerInternNew(i) = enerInternOld(i) + temp_work/density(i)
         enerInelasNew(i) = enerInelasOld(i) + (enerInternOld(i) - enerInternNew(i))
      end if

C ---- 2D Plane Strain (Voigt: 11,22,33,12) ----
      if ((ndir.eq.3) .and. (nshr.eq.1)) then
      if (stepTime .eq. 0.d0) then
         trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
         stressNew(i,1) = stressOld(i,1) + alamda*trace + twomu*strainInc(i,1)
         stressNew(i,2) = stressOld(i,2) + alamda*trace + twomu*strainInc(i,2)
         stressNew(i,3) = stressOld(i,3) + alamda*trace + twomu*strainInc(i,3)
         stressNew(i,4) = stressOld(i,4) + twomu*strainInc(i,4)
         kappa  = kappa_o
         damage = 0.d0
         do j=1,4
            totalStrain(j) = strainInc(i,j)
         enddo
      else
         kappa_n = stateOld(i,1)
		do j=1,4
		   totalStrain(j) = stateOld(i,2+j) + strainInc(i,j)
		enddo
		trace  = totalStrain(1) + totalStrain(2) + totalStrain(3)
		II1 = trace
		JJ2 = totalStrain(1)**2/three + totalStrain(2)**2/three
     &           - totalStrain(1)*totalStrain(2)/three + (totalStrain(4))**2
		JJ3 = (2.d0/27.d0)*(totalStrain(1)**3 + totalStrain(2)**3 + totalStrain(3)**3)
     &     	   - (1.d0/9.d0)*(totalStrain(1)**2*(totalStrain(2)+totalStrain(3)) 
     &              + totalStrain(2)**2*(totalStrain(1)+totalStrain(3)) 
     &              + totalStrain(3)**2*(totalStrain(1)+totalStrain(2)))
     &         + (4.d0/9.d0)*totalStrain(1)*totalStrain(2)*totalStrain(3)
     &         - (1.d0/3.d0)*((2.d0*totalStrain(3) - totalStrain(1) - totalStrain(2))*totalStrain(4)**2)
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
		ConstTERM = (ConstB*sqrt(JJ2)/(one+NU))
     &                +(ConstD*II1/(one-two*NU))+(ConstC*MaxPS/E)
		AA = ConstTERM**2 + 4.d0*ConstA*JJ2/((one+NU)**2)
		ee = (ConstTERM + sqrt(AA)) / (two*k)
		if (kappa_n .lt. ee) then
		   kappa = ee
		else
		   kappa = kappa_n
		end if
		if (kappa .gt. kappa_o) then
		   NODEE = charLength(i)
!		   Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_o * ft * NODEE)
		   Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_o * ft * NODEE)
			damage = One - (kappa_o / kappa)*(One-0.5*Betaa*(kappa-kappa_o))
		   if (damage .ge. one) then
				damage = 0.99
		   end if
		else
		   NODEE = charLength(i)
		   Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_o * ft * NODEE)
		   damage = 0.d0
		endif
            s(1) = (one-damage) * (alamda*trace + twomu*totalStrain(1))
            s(2) = (one-damage) * (alamda*trace + twomu*totalStrain(2))
            s(3) = (one-damage) * (alamda*trace + twomu*totalStrain(3))
            s(4) = (one-damage) * twomu*totalStrain(4)
            do j=1,4
               stressNew(i,j) = s(j)
            enddo
         endif
		! Internal and dissipated energy
         temp_work = 0.5d0*(
     1 (stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     2 + (stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     3 + (stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     4 + 2.0d0*(stressOld(i,4)+stressNew(i,4))*strainInc(i,4) )

         stateNew(i,1) = kappa
         stateNew(i,2) = damage
         do j=1,4
            stateNew(i,2+j) = totalStrain(j)
         enddo
		 stateNew(i,7) = Betaa
		 stateNew(i,8) = NODEE

         enerInternNew(i) = enerInternOld(i) + temp_work/density(i)
         enerInelasNew(i) = enerInelasOld(i) + (enerInternOld(i) - enerInternNew(i))
		endif
C ---- 3D (Voigt: 11,22,33,12,23,31)  ----
      if ((ndir.eq.3) .and. (nshr.eq.3)) then
      if (stepTime .eq. 0.d0) then
         trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
         stressNew(i,1) = stressOld(i,1) + alamda*trace + twomu*strainInc(i,1)
         stressNew(i,2) = stressOld(i,2) + alamda*trace + twomu*strainInc(i,2)
         stressNew(i,3) = stressOld(i,3) + alamda*trace + twomu*strainInc(i,3)
         stressNew(i,4) = stressOld(i,4) + twomu*strainInc(i,4)
		 stressNew(i,5) = stressOld(i,5) + twomu*strainInc(i,5)
		 stressNew(i,6) = stressOld(i,6) + twomu*strainInc(i,6)
         kappa  = kappa_o
         damage = 0.d0
      else
         kappa_n = stateOld(i,1)
         do j=1,6
            totalStrain(j) = stateOld(i,2+j) + strainInc(i,j)
         enddo
            trace  = totalStrain(1) + totalStrain(2) + totalStrain(3)
            II1 = trace
            JJ2 = (totalStrain(1)**2 + totalStrain(2)**2 + totalStrain(3)**2
     &         - totalStrain(1)*totalStrain(2) - totalStrain(2)*totalStrain(3) 
     &         - totalStrain(3)*totalStrain(1))/three
     &         + (totalStrain(4)**2 + totalStrain(5)**2 + totalStrain(6)**2)
			JJ3 = (2.d0/27.d0)*(totalStrain(1)**3 + totalStrain(2)**3 + totalStrain(3)**3)
     &     	   - (1.d0/9.d0)*(totalStrain(1)**2*(totalStrain(2)+totalStrain(3)) 
     &              + totalStrain(2)**2*(totalStrain(1)+totalStrain(3)) 
     &              + totalStrain(3)**2*(totalStrain(1)+totalStrain(2)))
     &         + (4.d0/9.d0)*totalStrain(1)*totalStrain(2)*totalStrain(3)
     &         - (1.d0/3.d0)*((2.d0*totalStrain(1) - totalStrain(2) - totalStrain(3))*totalStrain(6)**2
     &              + (2.d0*totalStrain(2) - totalStrain(1) - totalStrain(3))*totalStrain(5)**2
     &              + (2.d0*totalStrain(3) - totalStrain(1) - totalStrain(2))*totalStrain(4)**2)
     &         + 2.d0*totalStrain(4)*totalStrain(5)*totalStrain(6)

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
            ConstTERM = (ConstB*sqrt(JJ2)/(one+NU))
     &                +(ConstD*II1/(one-two*NU))+(ConstC*MaxPS/E)
            AA = ConstTERM**2 + 4.d0*ConstA*JJ2/((one+NU)**2)
            ee = (ConstTERM + sqrt(AA)) / (two*k)
            if (kappa_n .lt. ee) then
               kappa = ee
            else
               kappa = kappa_n
            end if
            if (kappa .gt. kappa_o) then
               NODEE = charLength(i)
               Betaa = E * kappa_o * NODEE / (Gf - 0.5d0 * kappa_o * ft * NODEE)
				damage = One - (kappa_o / kappa)*(One-0.5*Betaa*(kappa-kappa_o))
			   if (damage .ge. one) then
					damage = 0.99
			   end if	
            else
               damage = 0.d0
            endif
            s(1) = (one-damage) * (alamda*trace + twomu*totalStrain(1))
            s(2) = (one-damage) * (alamda*trace + twomu*totalStrain(2))
            s(3) = (one-damage) * (alamda*trace + twomu*totalStrain(3))
            s(4) = (one-damage) * twomu*totalStrain(4)
            s(5) = (one-damage) * twomu*totalStrain(5)
            s(6) = (one-damage) * twomu*totalStrain(6)
            do j=1,6
               stressNew(i,j) = s(j)
            enddo
         endif

         temp_work = 0.5d0*(
     1 (stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     2 + (stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     3 + (stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     4 + 2.d0*(stressOld(i,4)+stressNew(i,4))*strainInc(i,4)
     5 + 2.d0*(stressOld(i,5)+stressNew(i,5))*strainInc(i,5)
     6 + 2.d0*(stressOld(i,6)+stressNew(i,6))*strainInc(i,6) )

         stateNew(i,1) = kappa
         stateNew(i,2) = damage
         do j=1,6
            stateNew(i,2+j) = totalStrain(j)
         enddo

         enerInternNew(i) = enerInternOld(i) + temp_work/density(i)
         enerInelasNew(i) = enerInelasOld(i) + (enerInternOld(i) - enerInternNew(i))
      endif
 100  continue
      return
      end
