! =====================================
! Unsplit 2D nodal algorithm for tracer transport
! Uses tensor product Lagrange polynomial basis to simulate 2D tracer transport equations with variable windspeeds
!
! Dependencies:
! 	netCDF
!	LAPACK
!	nDGmod.f90 ; coeff_update.f90
! By: Devin Light ; Apr. 2014
! =====================================
PROGRAM EXECUTE
    USE nDGmod
    USE netCDF

    IMPLICIT NONE
    INTEGER :: startRes,polyOrder,whichTest,testEnd,nTest,ierr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: testsVec
    LOGICAL :: transient,DEBUG,doModalComparison,doTimeTest,doZSMaxCFL
    REAL(KIND=8) :: muMAX

    DEBUG = .FALSE.
    doModalComparison = .FALSE.
    doTimeTest = .FALSE.
    doZSMaxCFL = .TRUE.

    polyOrder = 4
	  startRes = 12

    IF(doZSMaxCFL) THEN
      SELECT CASE(polyOrder)
        CASE(2)
          muMAX = 0.450D0
!          muMAX  = 0.167D0
        CASE(3)
!          muMAX = 0.255D0
          muMAX  = 0.167D0
!          muMAX  = 0.083D0

!          muMAX = 0.13D0 ! modal CFL max
        CASE(4)
         muMAX = 0.168D0
!         muMAX = 0.083D0
!         muMAX  = 0.050D0

!        muMAX = 0.089D0 ! modal CFL max
        CASE(5)
!         muMAX = 0.120D0
          muMAX = 0.083D0
!         muMAX  = 0.033D0

!          muMAX = 0.066D0 ! modal CFL max
        CASE(6)
!        muMAX = 0.0910D0
         muMAX  = 0.024D0
        CASE(7)
!        muMAX = 0.0725D0
         muMAX  = 0.018D0
        CASE(8)
!        muMAX = 0.0589D0
         muMAX  = 0.014D0
        CASE(9)
          muMAX = 0.0490D0
!         muMAX  = 0.011D0
      END SELECT
    ELSE
      muMAX = 0.707D0
!      muMAX = 1.0D0
    ENDIF ! doZSMaxCFL
!    muMAX = muMAX*0.9*.707
     muMAX = muMAX*0.9
    testEnd = 1
    ALLOCATE(testsVec(1:testEnd),STAT=ierr)
    testsVec = (/ 5 /)

    write(*,*) '======================================================'
    write(*,*) '             BEGINNING RUN OF NODAL TESTS             '
    write(*,'(A27,F7.4)') 'muMAX=',muMAX
    WRITE(*,'(A,I1,A,I3,A)') ' N= ',polyOrder,' Uses a total of',(polyOrder+1)**2,' Lagrange basis polynomials per element'
    write(*,*) '======================================================'

    DO nTest=1,testEnd
      whichTest = testsVec(nTest)
      write(*,*) '======'
      SELECT CASE(whichTest)
        CASE(100)
          write(*,*) 'TEST 0: Consistency test'
          transient = .TRUE.
        CASE(1)
          write(*,*) 'TEST 1: Uniform advection (u=v=1)'
        	transient = .FALSE.
        CASE(2)
        	write(*,*) 'TEST 2: Solid body rotation of cylinder'
          transient = .FALSE.
        CASE(3)
        	write(*,*) 'TEST 3: Solid body rotation of cylinder (modified for frank)'
          transient = .FALSE.
        CASE(5)
          write(*,*) 'TEST 5: LeVeque Cosbell Deformation Test'
          transient = .TRUE.
        CASE(6)
        	write(*,*) 'TEST 6: LeVeque Smoother Cosbell Deformation Test'
        	transient = .TRUE.
        CASE(7)
        	write(*,*) 'TEST 7: Slotted cylinder deformation'
          transient = .TRUE.
        CASE(8)
          write(*,*) 'TEST 8: Constant Adv. of Gaussian Bump'
          transient = .FALSE.
        CASE(9)
          write(*,*) 'TEST 9: Constant Adv. of cosinebell'
          transient = .FALSE.
      END SELECT
    	write(*,*) '======'
    	CALL test2d_nodal(whichTest,startRes,startRes,2,3,2,muMAX)
    ENDDO
    DEALLOCATE(testsVec,STAT=ierr)

CONTAINS
	SUBROUTINE test2d_nodal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
    USE testParameters
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
		REAL(KIND=8), INTENT(IN) :: maxcfl

		! Local variables
		INTEGER, DIMENSION(10) :: tmp_method
    REAL(KIND=8), DIMENSION(nlevel) :: e1, e2, ei
		REAL(KIND=8) :: cnvg1, cnvg2, cnvgi,cons
		INTEGER :: nmethod, nmethod_final,imethod,ierr,nstep,nout
		INTEGER :: nQuadNodes,norder,p,gqOrder,nZSNodes

		LOGICAL :: dozshulimit,doStrangSplit,oddstep

		CHARACTER(len=80) :: cdf_out,outdir

		INTEGER :: nex,ney,nxiplot,netaplot
		REAL(KIND=8) :: dxel,dyel,tfinal, tmp_umax, tmp_vmax, dxm, dym,dt, time,calculatedMu
    DOUBLE PRECISION :: mux, muy
    REAL(KIND=8), DIMENSION(1:2) :: xEdge,yEdge
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: lagrangeDeriv,lagGaussVal,C0,C,tmpArray,tmpErr,lagValsZS
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xQuad, yQuad
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: gllNodes,gllWeights,gaussNodes,gaussWeights,x_elcent,y_elcent,&
                                               xplot,yplot,xiplot,etaplot,nodeSpacing,lambda,quadZSNodes,quadZSWeights
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: A,A0,& ! Coefficent array
													     u0,v0 ! Velocites within each element


		INTEGER :: i,j,l,m,k,s,t,n
		REAL(KIND=4), DIMENSION(2) :: tstart,tend
		REAL(KIND=4) :: t0,tf,t1,t2

		REAL(KIND=8) :: tmp_qmax,tmp_qmin

		if(nlevel.lt.1) STOP 'nlev should be at least 1 in test2d_modal'

		nmethod_final = 3
		tmp_method = 0
    tmp_method(1) = 1
    !tmp_method(1) = 2
		tmp_method(2) = 9
    tmp_method(3) = 10
    tmp_method(4) = 3
    tmp_method(5) = 5
    tmp_method(6) = 8

		DO nmethod=1,nmethod_final
			imethod = tmp_method(nmethod)

      doZShuLimit = .FALSE.
      doStrangSplit = .FALSE.
      doFluxMod = .FALSE.
      doPosLim = .FALSE.
      limitingMeth = -1

      write(*,*) ''
      write(*,*) '*******'

			SELECT CASE(imethod)
				CASE(1)
				  WRITE(*,*) '2D Nodal, Unsplit, No limiting'
				  norder = polyOrder
				  nQuadNodes = norder !2*(norder+1)
          gqOrder = 1
          nZSnodes = 1
          !outdir = '_ndgunlim/'
          write(outdir,'(A,I1,A)') '_ndgunlim/n',norder,'/'
				CASE(2)
				  WRITE(*,*) '2D Nodal, Unsplit, Zhang and Shu Limiting'
				  dozshulimit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 1
				  norder = polyOrder
          nQuadNodes = norder
				  gqOrder = CEILING((polyOrder+1)/2D0 )-1
          nZSNodes = CEILING((polyOrder+3)/2D0 )-1
          !outdir = '_ndgzhshu/n5/pRefine/rescale/'
          write(outdir,'(A,I1,A)') '_ndgzhshu/n',norder,'/'
          WRITE(*,'(A,I1,A)') ' Using ',gqOrder+1,' points for gauss quadrature nodes'

        CASE(3)
          write(*,*) '2D Nodal, Strang Split, No limiting'
          doStrangSplit = .TRUE.
          norder = polyOrder
          nQuadNodes = norder
          gqOrder = 1
          nZSnodes = 1
          !outdir = '_ndgsplun/'
          WRITE(outdir,'(A,I1,A)') '_ndgsplun/n',norder,'/'
        CASE(4)
          write(*,*) '2D Nodal, Strang Split, Zhang and Shu Limiting'
          dozshulimit = .TRUE.
          doStrangSplit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 1
          nOrder = polyOrder
          nQuadNodes = nOrder
          gqOrder = 1
          nZSNodes = CEILING((polyOrder+3)/2D0 )-1
!          nZSNodes = nQuadNodes
!outdir = '_ndgsplzs/n5/pRefine/rescale/'
          WRITE(outdir,'(A,I1,A)') '_ndgsplzs/n',norder,'/'

        CASE(5)
          WRITE(*,*) '2D Nodal, Strang Split, FCT + TMAR'
          doFluxMod = .TRUE.
          doStrangSplit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 2
          nOrder = polyOrder
          nQuadNodes = nOrder
          gqOrder = 1
          nZSNodes = nQuadNodes
          !outdir = '_ndgsplfc/'
          write(outdir,'(A,I1,A)') '_ndgsplfc/n',norder,'/'
        CASE(6)
          WRITE(*,*) '2D Nodal, Unsplit, MA Truncation (using ZS thm)'
          dozshulimit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 2
          norder = polyOrder
          nQuadNodes = norder
          gqOrder = CEILING((polyOrder+1)/2D0 )-1
          nZSNodes = CEILING((polyOrder+3)/2D0 )-1
          !write(outdir,'(A,I1,A)') '_matrunc/n',norder,'/meanLimitingComparison/tmarMin/'
          !write(outdir,'(A,I1,A)') '_matrunc/n',norder,'/meanLimitingComparison/tmarFull/'
          write(outdir,'(A,I1,A)') '_matrunc/n',nOrder,'/'
          !WRITE(*,'(A,I1,A)') '     Using ',gqOrder+1,' points for gauss quadrature nodes'
          !WRITE(*,'(A,I1,A)') '     Using ',nZSNodes+1,' GLL nodes for positivity rescaling.'
        CASE(7)
          WRITE(*,*) '2D Nodal, Unsplit, MA Truncation + mean limiting at ZSmin'
          dozshulimit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 3
          norder = polyOrder
          nQuadNodes = norder
          gqOrder = CEILING((polyOrder+1)/2D0 )-1
          nZSNodes = CEILING((polyOrder+3)/2D0 )-1
          !outdir =
          write(outdir,'(A,I1,A)') '_matrunc/n',norder,'/meanLimitingComparison/zsMin/'
          WRITE(*,'(A,I1,A)') ' Using ',gqOrder+1,' points for gauss quadrature nodes'

        CASE(8)
          WRITE(*,*) '2D Nodal, Split, Lambda Limiting + TMAR each stage'
          doFluxMod = .TRUE.
          doStrangSplit = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 3
          nOrder = polyOrder
          nQuadNodes = nOrder
          gqOrder = 1
          nZSNodes = nQuadNodes
          !outdir = '_ndgspllm/'
          write(outdir,'(A,I1,A)') '_ndgspllm/n',norder,'/'

        CASE(9)
          WRITE(*,*) '2D Nodal, Unsplit, Lambda Limiting + TMAR each stage'
          doFluxMod = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 5
          nOrder = polyOrder
          nQuadNodes = nOrder
          gqOrder = 1
          nZSNodes = nQuadNodes
          write(outdir,'(A,I1,A)') '_ndguslam/n',norder,'/'

        CASE(10)
          WRITE(*,*) '2D Nodal, Unsplit, sliced Lambda Limiting + TMAR'
          doFluxMod = .TRUE.
          doPosLim = .TRUE.
          limitingMeth = 6
          nOrder = polyOrder
          nQuadNodes = nOrder
          gqOrder = 1
          nZSNodes = nQuadNodes
          write(outdir,'(A,I1,A)') '_ndguslam/n',norder,'/'

			END SELECT
      IF(limitingMeth == -1) THEN
        WRITE(*,*) 'No Limiting Active'
      ELSE
        WRITE(*,'(A,I1)') 'limitingMeth = ',limitingMeth
        WRITE(*,'(A,I1,A)') 'Using ',nZSNodes+1,' GLL nodes for positivity limiting.'
      ENDIF
      WRITE(*,*) 'Output directory is: ',outdir
      write(*,*) '*******'

			! Initialize quadrature weights and nodes (only done once per call)
			ALLOCATE(gllnodes(0:nQuadNodes),gllweights(0:nQuadNodes),gaussNodes(0:gqOrder),gaussWeights(0:gqOrder),&
                lagrangeDeriv(0:norder,0:nQuadNodes),tmpArray(0:norder,0:norder),nodeSpacing(0:nQuadNodes-1),&
                lagGaussVal(0:norder,0:gqOrder),lambda(0:nQuadNodes),STAT=ierr)

      ALLOCATE(quadZSNodes(0:nZSNodes),quadZSWeights(0:nZSNodes),lagValsZS(0:norder,0:nZSNodes),STAT=ierr)

  		! --  Compute Gauss-Lobatto quadrature nodes and weights
      CALL gllquad_nodes(nQuadNodes,gllNodes)
      CALL gllquad_weights(nQuadNodes,gllNodes,gllWeights)
      CALL gllquad_nodes(nZSNodes,quadZSNodes)
      CALL gllquad_weights(nZSNodes,quadZSNodes,quadZSWeights)
      CALL gaussquad_nodes(gqOrder+1,gaussNodes)
      CALL gaussquad_weights(gqOrder+1,gaussNodes,gaussWeights)

      zsMinWeight = MINVAL(quadZSWeights)

      nodeSpacing = gllNodes(1:nQuadNodes)-gllNodes(0:nQuadNodes-1)

			nxiplot = norder+1
			netaplot = norder+1

			ALLOCATE(xiplot(1:nxiplot), etaplot(1:netaplot), STAT=ierr)

			xiplot = gllNodes(:)
			etaplot = gllNodes(:)

      ! Fill matrix of basis derivatives at GLL nodes (note that this assumes lagrangeDeriv is square!)
      CALL Dmat(nQuadNodes,gllNodes,lagrangeDeriv)

      ! Fill maxtris of basis polynomials evaluated at Gauss quadrature nodes (used in Zhang and Shu Limiter)
      CALL baryWeights(lambda,gllnodes,nQuadNodes)
      DO i=0,norder
        DO j=0,gqOrder
          lagGaussVal(i,j) = lagrange(gaussNodes(j),i,nQuadNodes,gllNodes,lambda)
        ENDDO !j
      ENDDO !i

      lagValsZS = 0D0
      ! -- Evaluate Basis polynomials at quad nodes for Zhang & Shu positivity limiting (edge nodes fixed at -1 and +1)
      DO i=0,norder
        DO j=0,nZSNodes
          lagValsZS(i,j) = lagrange(quadZSNodes(j),i,nQuadNodes,gllNodes,lambda)
        ENDDO!l
      ENDDO!k

      xEdge(1) = 0D0
      xEdge(2) = 1D0
      if(ntest .eq. 11 .or. ntest .eq. 10) xEdge(1) = -1D0
      yEdge = xEdge

      write(*,*) 'Domain: [',xEdge(1),',',xEdge(2),'].'
      IF(doTimeTest .or. DEBUG) THEN
        write(*,*) 'Warning: overwritting number of outputs!'
      ENDIF

			DO p=1,nlevel
        magCount = 0

        CALL CPU_TIME(t0)

				nex = nex0*nscale**(p-1)
				ney = ney0*nscale**(p-1)
				dxel = (xEdge(2)-xEdge(1))/DBLE(nex)
				dyel = (yEdge(2)-yEdge(1))/DBLE(ney)

				ALLOCATE(x_elcent(1:nex),y_elcent(1:ney),A(0:norder,0:norder,1:nex,1:ney),A0(0:norder,0:norder,1:nex,1:ney),&
                 u0(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), v0(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), &
			           C0(1:nex,1:ney), C(1:nex,1:ney), xplot(1:nxiplot*nex),yplot(1:netaplot*ney),&
                 xQuad(0:nQuadNodes,1:nex),yQuad(0:nQuadNodes,1:ney), STAT=ierr)

				! Elements are ordered row-wise within the domain

				! A(l,m,i,j) : Coefficent array ; (i,j) = element, (l,m) = which Lagrange polynomial (q_ij = SUM(SUM(A(i,j,l,m)*L_l(xi)*L_m(eta))))
        ! This is also the solution at the interpolating GLL points
				! u(l,m,i,j), v(l,m,i,j) : Horizontal/vertical vel. array ; (i,j) = element, (l,m) = horiz,vertical location

				! Initialize x- and y- grids and xi- and eta- plotting grids
				x_elcent(1) = xEdge(1)+dxel/2D0
				DO i=2,nex
					x_elcent(i) = x_elcent(i-1)+dxel
				ENDDO !i

				y_elcent(1) = yEdge(1)+dyel/2D0
				DO i=2,ney
					y_elcent(i) = y_elcent(i-1)+dyel
				ENDDO !i

				DO i=1,nex
					xplot(1+(i-1)*nxiplot:i*nxiplot) = x_elcent(i)+(dxel/2D0)*xiplot(:)
				ENDDO !i

				DO i=1,ney
					yplot(1+(i-1)*netaplot:i*netaplot) = y_elcent(i)+(dyel/2D0)*etaplot(:)
				ENDDO !i


				! Initialize A, u, and v
        CALL init2d(ntest,nex,ney,nQuadNodes,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,xEdge,yEdge,cdf_out,&
                    tfinal,xQuad,yQuad)

        IF(DEBUG) THEN
            write(*,*) '-------- Warning DEBUG = .TRUE. --------'
        ENDIF ! DEBUG
        A0 = A
				cdf_out = TRIM(outdir) // TRIM(cdf_out)

				! Store element averages for conservation estimation
        DO i=1,nex
          DO j=1,ney
            DO l=0,norder
              tmpArray(l,:) = gllWeights(l)*gllWeights(:)*A0(l,:,i,j)
            ENDDO!l
            C0(i,j) = 0.25D0*dxel*dyel*SUM(tmpArray)
          ENDDO !j
        ENDDO !i
        !write(*,*) 'Minimum element mean = ',MINVAL(C0)

				! Set up timestep
        IF(doModalComparison .or. doZSMaxCFL) THEN
          dxm = dxel
          dym = dyel
        ELSE
    			dxm = dxel*MINVAL(nodeSpacing)/2D0
          dym = dyel*MINVAL(nodeSpacing)/2D0
        ENDIF

				tmp_umax = MAXVAL(abs(u0))
				tmp_vmax = MAXVAL(abs(v0))

        IF(noutput .eq. -1) THEN
		      nstep = CEILING((tfinal/maxcfl)*(maxval(sqrt(u0**2 + v0**2))/min(dxm,dym)))
			    nout = nstep
        ELSE IF(doModalComparison .or. doStrangSplit) THEN
          nstep = noutput*CEILING( (tfinal/maxcfl)*MAX(tmp_umax/dxm,tmp_vmax/dym)/DBLE(noutput) )
          nout = noutput
        ELSE IF(doZSMaxCFL) THEN
          nstep = noutput*CEILING( (tfinal/maxcfl)*(tmp_umax/dxm + tmp_vmax/dym)/DBLE(noutput) )
          nout = noutput
        ELSE
      		nstep = noutput*CEILING( (tfinal/maxcfl)*(maxval(sqrt(u0**2 + v0**2))/min(dxm,dym))/DBLE(noutput) )
      		nout = noutput
        ENDIF !noutput

        IF(doTimeTest .or. DEBUG) THEN
          nout = noutput
          nstep = 800*nscale**(p-1)
        ENDIF

				dt = tfinal/DBLE(nstep)

        ! mu1 and mu2 are used in Zhang and Shu (2010) mean positivity limiting
        mux = tmp_umax*dt/dxel
        muy = tmp_vmax*dt/dyel

        mu1 = mux/(mux+muy)
        mu2 = muy/(mux+muy)

!        calculatedMu = maxval(sqrt(u0**2 + v0**2))*dt/min(dxm,dym)
!                calculatedMu = maxval(sqrt(u0**2 + v0**2))*dt/min(dxel,dyel)
!        write(*,'(A,E10.4,A,E10.4,A,E10.4)') '  mu used = ',calculatedMu, &
!                                     '  mu2 used = ',maxval(sqrt(u0**2 + v0**2))*dt/min(dxel,dyel), &
!                                     '  mux + muy = ', mux+muy
        IF(doZSMaxCFL) THEN
          calculatedMu = mux+muy
        ELSE
          calculatedMu = maxval(sqrt(u0**2 + v0**2))*dt/min(dxel,dyel)
        ENDIF
				IF(p .eq. 1) THEN ! Set up netCDF file
          write(*,*) 'Maximum velocity: |u| = ',maxval(abs(u0)),maxval(abs(v0))!maxval(abs(sqrt(u0**2+v0**2)))
					CALL output2d(A0,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,nxiplot,netaplot,&
                        tfinal,calculatedMu,cdf_out,nout,-1)
				ENDIF

        ! Set up variables for this value of p ; Write x, y, and initial conditions
				CALL output2d(A0,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,nxiplot,netaplot,0D0,calculatedMu,cdf_out,p,0)

				! Time integration
				tmp_qmax = MAXVAL(A0)
				tmp_qmin = MINVAL(A0)

				time = 0D0

        IF(doStrangSplit) THEN
          ! For Strang Split nodal methods, use strangSplitUpdate and nDGsweep to update nodal values
          oddstep = .TRUE.
          DO n=1,nstep
            CALL CPU_TIME(t1)
            CALL strangSplitUpdate(A,u0,v0,gllNodes,gllWeights,time,lagrangeDeriv,&
                         dt,dxel,dyel,nOrder,nQuadNodes,nex,ney,oddstep,transient,&
                         dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
            CALL CPU_TIME(t2)
            time = time + dt

      	    IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN
  			     ! Write output
             CALL output2d(A,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,&
                            nxiplot,netaplot,time,calculatedMu,cdf_out,p,2)
            ENDIF

            ! Track solution max and min throughout integration
            tmp_qmax = MAX(tmp_qmax,MAXVAL(A))
            tmp_qmin = MIN(tmp_qmin,MINVAL(A))

            oddstep = .NOT. oddstep

          ENDDO !n

        ELSE
        ! Use coeff_update to update 2d elements
		      DO n=1,nstep
            CALL CPU_TIME(t1)
            CALL coeff_update(A,u0,v0,gllNodes,gllWeights,gaussWeights,lagrangeDeriv,time,dt,dxel,dyel,nex,ney,&
                              norder,nQuadNodes,gqOrder,lagGaussVal,nZSnodes,lagValsZS,&
                              dozshulimit,transient,doZSMaxCFL)
            CALL CPU_TIME(t2)
      			time = time + dt

      			IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN
      				! Write output
      				CALL output2d(A,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,&
                            nxiplot,netaplot,time,calculatedMu,cdf_out,p,2)
      			ENDIF

            ! Track solution max and min throughout integration
      			tmp_qmax = MAX(tmp_qmax,MAXVAL(A))
      			tmp_qmin = MIN(tmp_qmin,MINVAL(A))

      		ENDDO !n
        ENDIF ! doStrangSplit
        CALL CPU_TIME(tf)
        tf = tf - t0

        ! Store element averages for conservation estimation
        DO i=1,nex
          DO j=1,ney
            DO l=0,norder
              tmpArray(l,:) = gllWeights(l)*gllWeights(:)*A(l,:,i,j)
            ENDDO!l
            C(i,j) = 0.25D0*dxel*dyel*SUM(tmpArray)
          ENDDO !j
        ENDDO !i
        cons = SUM(C-C0)/DBLE(nex*ney)

        ! Compute errors
        ALLOCATE(tmpErr(1:nex,1:ney),STAT=ierr)
        DO i=1,nex
          DO j=1,ney
            DO l=0,norder
              tmpArray(l,:) = gllWeights(l)*gllWeights(:)*ABS(A(l,:,i,j)-A0(l,:,i,j))
            ENDDO!l
            tmpErr(i,j) = SUM(tmpArray)
          ENDDO !j
        ENDDO !i
        e1(p) = 0.25D0*dxel*dyel*SUM(tmpErr)

        DO i=1,nex
          DO j=1,ney
            DO l=0,norder
              tmpArray(l,:) = gllWeights(l)*gllWeights(:)*(ABS(A(l,:,i,j)-A0(l,:,i,j)))**2
            ENDDO!l
            tmpErr(i,j) = SUM(tmpArray)
          ENDDO !j
        ENDDO !i
        e2(p) = sqrt(0.25D0*dxel*dyel*SUM(tmpErr))

        ei(p) = MAXVAL(ABS(A(:,:,:,:) - A0(:,:,:,:)))


  			if (p.eq.1) then
      	   write(UNIT=6,FMT='(A131)') &
'   nex    ney     E1       E2         Einf      convergence rate  overshoot  undershoot   cons        cputime  nstep   tf   cfl   '
        	cnvg1 = 0.d0
        	cnvg2 = 0.d0
        	cnvgi = 0.d0
        else
      		cnvg1 = -log(e1(p)/e1(p-1))/log(dble(nscale))
      		cnvg2 = -log(e2(p)/e2(p-1))/log(dble(nscale))
      		cnvgi = -log(ei(p)/ei(p-1))/log(dble(nscale))
        end if
        write(*,990) nex,ney, e1(p), e2(p), ei(p), &
                cnvg1, cnvg2, cnvgi, &
                tmp_qmax-MAXVAL(A0), &
                MINVAL(A0)-tmp_qmin, &
                cons, tf, nstep,tfinal,calculatedMu,magCount

				IF(p .eq. nlevel) THEN
          ! Close netCDF files and write cputime vector
          CALL output2d(A,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,nxiplot,netaplot,time,&
                        calculatedMu,cdf_out,p,1)
				ENDIF
				DEALLOCATE(A,A0,x_elcent,y_elcent,xplot,yplot,u0,v0,tmpErr,xQuad,yQuad,C,C0,STAT=ierr)
			ENDDO
        		DEALLOCATE(gllnodes,gllweights,gaussnodes,gaussweights,lagrangeDeriv,nodeSpacing,lagGaussVal,lambda, &
                       xiplot,etaplot,tmpArray, tmpErr, STAT=ierr)
            DEALLOCATE(quadZSNodes,quadZSWeights,lagValsZS,STAT=ierr)

		ENDDO


990    format(2i6,3e12.4,3f5.2,3e12.4,f8.2,i8,f8.2,e12.4,i6)

	END SUBROUTINE test2d_nodal

  SUBROUTINE strangSplitUpdate(A,u0,v0,gllNodes,gllWeights,time,lagrangeDeriv,&
                               dt,dxel,dyel,nOrder,nQuadNodes,nex,ney,oddstep,transient,&
                               dozshulimit,nZSnodes,quadZSWeights,lagValsZS)

  ! =====================================================================================================
  ! strangSplitUpdate is responsible for selecting which slice of subcell volumes is sent to mDGsweep for update to time
  ! level tn+1 following a Strang splitting.
  ! For Strang splitting:
  !   - Each slice is updated
  !   - Odd steps: x-slices are updated first (horizontal advection) then y-slices are updated (vertical advection)
  !   - Even steps: y-slices are updated first then x-slices are updated (vertical advection)
  ! =====================================================================================================

    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nOrder,nQuadNodes,nex,ney,nZSnodes
    REAL(KIND=8), INTENT(IN) :: dt,dxel,dyel,time
    REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), INTENT(IN) :: u0,v0
    REAL(KIND=8), DIMENSION(0:norder,0:nQuadNodes), INTENT(IN) :: lagrangeDeriv
    REAL(KIND=8), DIMENSION(0:nQuadNodes), INTENT(IN) :: gllNodes,gllWeights
    REAL(KIND=8), DIMENSION(0:nZSnodes), INTENT(IN) :: quadZSWeights
    REAL(KIND=8), DIMENSION(0:nOrder,0:nZSnodes), INTENT(IN) :: lagValsZS
    LOGICAL, INTENT(IN) :: oddstep,dozshulimit,transient
    ! Outputs
    REAL(KIND=8), DIMENSION(0:nOrder,0:nOrder,1:nex,1:ney), INTENT(INOUT) :: A
    ! Local variables
    INTEGER :: i,j,k,whichEl,whichLvl,totalLvls
    REAL(KIND=8), DIMENSION(0:nOrder,1:nex) :: A1dx
    REAL(KIND=8), DIMENSION(0:nOrder,1:ney) :: A1dy
    REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex) :: u1dx
    REAL(KIND=8), DIMENSION(0:nQuadNodes,1:ney) :: v1dy

    IF(oddstep) THEN
      ! ===================================
      ! Perform sweeps in x-direction first
      ! ===================================
      totalLvls = ney*(nOrder+1)
      DO i=0,totalLvls-1
        whichEl = i/(nOrder+1) + 1
        whichLvl = MOD(i,nOrder+1)
        A1dx(:,1:nex) = A(:,whichLvl,1:nex,whichEl)
        u1dx(:,1:nex) = u0(:,whichLvl,1:nex,whichEl)

        CALL nDGsweep(A1dx,nex,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u1dx,lagrangeDeriv,time,dt,transient,&
                      dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
        ! Update solution
        A(:,whichLvl,1:nex,whichEl) = A1dx(:,1:nex)
      ENDDO !i

      totalLvls = nex*(nOrder+1)
      DO i=0,totalLvls-1
        whichEl = i/(nOrder+1) + 1
        whichLvl = MOD(i,nOrder+1)
        A1dy(:,1:ney) = A(whichLvl,:,whichEl,1:ney)
        v1dy(:,1:ney) = v0(whichLvl,:,whichEl,1:ney)

        CALL nDGsweep(A1dy,ney,dyel,nOrder,nQuadNodes,gllNodes,gllWeights,v1dy,lagrangeDeriv,time,dt,transient,&
                      dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
        ! Update solution
        A(whichLvl,:,whichEl,1:ney) = A1dy(:,1:ney)
      ENDDO !i
    ELSE
      ! ===================================
      ! Perform sweeps in y-direction first
      ! ===================================
      totalLvls = nex*(nOrder+1)
      DO i=0,totalLvls-1
        whichEl = i/(nOrder+1) + 1
        whichLvl = MOD(i,nOrder+1)
        A1dy(:,1:ney) = A(whichLvl,:,whichEl,1:ney)
        v1dy(:,1:ney) = v0(whichLvl,:,whichEl,1:ney)

        CALL nDGsweep(A1dy,ney,dyel,nOrder,nQuadNodes,gllNodes,gllWeights,v1dy,lagrangeDeriv,time,dt,transient,&
                      dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
        ! Update solution
        A(whichLvl,:,whichEl,1:ney) = A1dy(:,1:ney)
      ENDDO !i

      totalLvls = ney*(nOrder+1)
      DO i=0,totalLvls-1
        whichEl = i/(nOrder+1) + 1
        whichLvl = MOD(i,nOrder+1)
        A1dx(:,1:nex) = A(:,whichLvl,1:nex,whichEl)
        u1dx(:,1:nex) = u0(:,whichLvl,1:nex,whichEl)

        CALL nDGsweep(A1dx,nex,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u1dx,lagrangeDeriv,time,dt,transient,&
                      dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
        ! Update solution
        A(:,whichLvl,1:nex,whichEl) = A1dx(:,1:nex)
      ENDDO !i
    ENDIF !oddstep
  END SUBROUTINE strangSplitUpdate

  SUBROUTINE init2d(ntest,nex,ney,nQuadNodes,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,xEdge,yEdge,cdf_out,tfinal,&
                      xQuad,yQuad)
  ! ----
  !  Computes initial conditons for coefficient matrix (by simple evaluation of IC function) and initial velocities (via a streamfunction)
  !  Also sets the final time and output file name
  ! ----
      IMPLICIT NONE
      !Inputs
      INTEGER, INTENT(IN) :: ntest,nex,ney,nQuadNodes,norder
      REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: x_elcent
      REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: y_elcent
      REAL(KIND=8), DIMENSION(0:nQuadNodes) :: gllNodes
      REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: xEdge,yEdge
      !Outputs
      REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(OUT) :: A,u0,v0
      CHARACTER(LEN=*), INTENT(OUT) :: cdf_out
      REAL(KIND=8), INTENT(OUT) :: tfinal
      !Local Variables
      REAL(KIND=8), DIMENSION(0:norder,1:nex), INTENT(OUT) :: xQuad
      REAL(KIND=8), DIMENSION(0:norder,1:ney), INTENT(OUT) :: yQuad
      REAL(KIND=8), DIMENSION(0:norder,0:norder,0:1,1:nex,1:ney) :: psiu,psiv
      REAL(KIND=8), DIMENSION(0:norder,0:norder) :: r
      REAL(KIND=8), DIMENSION(1:2) :: domainCenter
      REAL(KIND=8) :: PI,dxmin,dymin,dxel,dyel,xWidth,yWidth,xc,yc,spd
      INTEGER :: i,j,k,l

	    PI = DACOS(-1D0)
      spd = 2D0*PI
      IF(ntest .eq. 11) spd = 1D0

      xWidth = xEdge(2)-xEdge(1)
      yWidth = yEdge(2)-yEdge(1)
      domainCenter(1) = xEdge(1)+xWidth/2D0
      domainCenter(2) = yEdge(1)+yWidth/2D0

      dxel = x_elcent(2) - x_elcent(1)
      dyel = y_elcent(2) - y_elcent(1)

	    ! Minimum internode spacing, mapped to physical domain
      dxmin = (dxel/2D0)*MINVAL(gllNodes(1:nQuadNodes)-gllNodes(0:nQuadNodes-1))
      dymin = (dyel/2D0)*MINVAL(gllNodes(1:nQuadNodes)-gllNodes(0:nQuadNodes-1))

      ! Quadrature locations, mapped to physical domain
      DO i=1,nex
        xQuad(:,i) = x_elcent(i) + (dxel/2d0)*gllNodes(:)
      ENDDO!i
      DO j=1,ney
        yQuad(:,j) = y_elcent(j) + (dyel/2d0)*gllNodes(:)
      ENDDO!j

      ! Fill streamfunction array
      SELECT CASE(ntest)
      CASE(1,8:10) ! uniform velocity u = v = 1
        DO j=1,ney
          DO i=1,nex
            DO k=0,norder
              psiu(k,:,1,i,j) = (yQuad(:,j) + dymin/2d0) - xQuad(k,i)
              psiu(k,:,0,i,j) = (yQuad(:,j) - dymin/2d0) - xQuad(k,i)
!              psiu(k,:,1,i,j) = 0D0
!              psiu(k,:,0,i,j) = 0D0

              psiv(k,:,1,i,j) = yQuad(:,j) - (xQuad(k,i) + dxmin/2d0)
              psiv(k,:,0,i,j) = yQuad(:,j) - (xQuad(k,i) - dxmin/2d0)
!              psiv(k,:,1,i,j) = 0D0
!              psiv(k,:,0,i,j) = 0D0

            ENDDO !k
          ENDDO !j
        ENDDO !i
      CASE(2:4,11) ! Solid body rotation
        ! pi*( (x-0.5)**2 + (y-0.5)**2 )
        DO j=1,ney
          DO i=1,nex
            DO k=0,norder
           psiu(k,:,1,i,j) = spd*0.5D0*( (xQuad(k,i)-domainCenter(1))**2 + (yQuad(:,j) + dymin/2d0 -domainCenter(2))**2 )
           psiu(k,:,0,i,j) = spd*0.5D0*( (xQuad(k,i)-domainCenter(1))**2 + (yQuad(:,j) - dymin/2d0 -domainCenter(2))**2 )

           psiv(k,:,1,i,j) = spd*0.5D0*( (xQuad(k,i) + dxmin/2d0-domainCenter(1))**2 + (yQuad(:,j) -domainCenter(2))**2 )
           psiv(k,:,0,i,j) = spd*0.5D0*( (xQuad(k,i) - dymin/2d0-domainCenter(1))**2 + (yQuad(:,j) -domainCenter(2))**2 )
            ENDDO !k
          ENDDO !i
        ENDDO !j
      CASE(5:7,100) ! Leveque deformation flow
			!(1/pi)*sin(pi*x(i))**2 * sin(pi*y(j))**2
        DO j=1,ney
          DO i=1,nex
            DO k=0,norder
            psiu(k,:,1,i,j) = (1D0/PI)*(DSIN(PI*xQuad(k,i))**2 * DSIN(PI*(yQuad(:,j)+dymin/2d0))**2)
            psiu(k,:,0,i,j) = (1D0/PI)*(DSIN(PI*xQuad(k,i))**2 * DSIN(PI*(yQuad(:,j)-dymin/2d0))**2)

            psiv(k,:,1,i,j) = (1D0/PI)*(DSIN(PI*(xQuad(k,i)+dxmin/2d0))**2 * DSIN(PI*yQuad(:,j))**2)
            psiv(k,:,0,i,j) = (1D0/PI)*(DSIN(PI*(xQuad(k,i)-dxmin/2d0))**2 * DSIN(PI*yQuad(:,j))**2)
            ENDDO !k
          ENDDO !i
        ENDDO !j
      END SELECT !ntest

      ! Compute velocity from stream functions
      DO j=1,ney
        DO i=1,nex
          u0(:,:,i,j) = (psiu(:,:,1,i,j)-psiu(:,:,0,i,j) )/dymin
          v0(:,:,i,j) = -(psiv(:,:,1,i,j)-psiv(:,:,0,i,j) )/dxmin
        ENDDO !j
      ENDDO !i

      ! Initialize coefficent array
      SELECT CASE(ntest)
      CASE(1) ! uniform advection of a sine wave
        cdf_out = 'dg2d_sine_adv.nc'
        tfinal = 10D0
        DO j=1,ney
          DO i=1,nex
            DO k=0,norder
! A(i,j,k,:) = DSIN(2D0*PI*xQuad(i,k))*DSIN(2D0*PI*yQuad(j,:))
              A(k,:,i,j) = 1D0+DSIN(2D0*PI*xQuad(k,i))*DSIN(2D0*PI*yQuad(:,j))
            ENDDO !k
          ENDDO!i
        ENDDO!j

        CASE(2) ! solid body rotation of a cylinder
          cdf_out = 'dg2d_rot_cylinder.nc'
!tfinal = 2D0*PI
          tfinal = 1D0
          A = 0D0
          DO j=1,ney
            DO i=1,nex
              DO k=0,norder
                r(k,:) = SQRT((xQuad(i,k)-0.3d0)**2 + (yQuad(:,j)-0.3d0)**2)
              ENDDO !k
              WHERE(r .lt. 0.125D0)
                A(:,:,i,j) = 1D0
              END WHERE
            ENDDO!i
          ENDDO!j

          CASE(3,10:11) ! solid body rotation of a cylinder (comparison to frank's code)
            cdf_out = 'dg2d_rot_cylinder_modified.nc'
            tfinal = 1D0
            A = 0D0
            xc = xEdge(1)+xWidth/4D0
            yc = yEdge(1)+yWidth/2D0

            DO j=1,ney
              DO i=1,nex
                DO k=0,norder
                  r(k,:) = SQRT((xQuad(k,i)-xc)**2 + (yQuad(j,:)-yc)**2)
                ENDDO !k
                WHERE(r .lt. 0.25D0)
                  A(:,:,i,j) = 1D0
                END WHERE
              ENDDO!i
            ENDDO!j

          CASE(5,9) ! standard cosbell deformation
    				cdf_out = 'dg2d_def_cosbell.nc'
    				tfinal = 5D0
            A = 0D0
            DO j=1,ney
              DO i=1,nex
                DO k=0,norder ! Fill distance array for (i,j) element
                  r(k,:) = 4D0*SQRT( (xQuad(k,i)-0.25D0)**2 + (yQuad(:,j)-0.25D0)**2 )
                ENDDO !k
                WHERE(r .lt. 1D0)
!                            A(i,j,:,:) = (0.5d0*(1.d0 + DCOS(PI*r(:,:))))
                  A(:,:,i,j) = (0.25d0*(1.d0 + DCOS(PI*r(:,:)))**2)
                END WHERE
              ENDDO !i
            ENDDO !j

          CASE(6) ! smoother cosbell deformation
            cdf_out = 'dg2d_smth_cosbell.nc'
            tfinal = 5D0
            A = 0D0
            DO j=1,ney
              DO i=1,nex
                DO k=0,norder ! Fill distance array for (i,j) element
                  r(k,:) = 3D0*SQRT( (xQuad(k,i)-0.4D0)**2 + (yQuad(:,j)-0.4D0)**2 )
                ENDDO !k
                WHERE(r .lt. 1D0)
                  A(:,:,i,j) = (0.5d0*(1.d0 + DCOS(PI*r(:,:))))**3
                END WHERE
              ENDDO !i
            ENDDO !j
          CASE(7) ! slotted cylinder ICs
            cdf_out = 'dg2d_def_cyl.nc'
            tfinal = 5D0
            A = 0D0
            xc = 0.25D0
            yc = 0.5D0
            DO j=1,ney
              DO i=1,nex
                DO k=0,nOrder
                  r(k,:) = SQRT((xQuad(k,i)-xc)**2 + (yQuad(j,:)-yc)**2)
                ENDDO !k
                WHERE(r .lt. 0.15D0)
                  A(:,:,i,j) = 1D0
                END WHERE

                DO k=0,nOrder
                  DO l=0,nOrder
                    IF(ABS(xQuad(k,i)-xc) .lt. 0.025D0 .AND. yQuad(l,j) .gt.(yc-0.0625D0)) THEN
                      A(k,l,i,j) = 0D0
                    ENDIF
                  ENDDO !l
                ENDDO !k
              ENDDO !i
            ENDDO !j
          CASE(8) ! gaussian bump
            cdf_out = 'dg2d_gauss_bump.nc'
            tfinal = 1D0
            DO j=1,ney
              DO i=1,nex
                DO k=0,norder
                  A(k,:,i,j) = EXP(- 2*(xQuad(k,i)**2 + yQuad(:,j)**2 )/(8D0**2) )
                ENDDO !k
              ENDDO !j
            ENDDO !i
      END SELECT !ntest

      if(ntest .eq. 10) then
          tfinal = xWidth*10D0*tfinal*5
      elseif(ntest .eq. 11) then
          tfinal = 3D0*2D0*PI
      elseif(ntest .eq. 6) then
          tfinal = 1D0*tfinal
      endif

  END SUBROUTINE init2d

	SUBROUTINE output2d(A,x,y,gllWeights,gllNodes,nex,ney,norder,nQuadNodes,nxiplot,netaplot,tval_in,mu,cdf_out,ilvl,stat)
		IMPLICIT NONE

		! Inputs
		INTEGER, INTENT(IN) :: norder,nQuadNodes,nex,ney,nxiplot,netaplot,stat,ilvl
		CHARACTER(len=*), INTENT(IN) :: cdf_out
		REAL(KIND=8), INTENT(IN) :: tval_in,mu
		REAL(KIND=8), DIMENSION(1:nex*nxiplot), INTENT(IN) :: x
		REAL(KIND=8), DIMENSION(1:ney*netaplot), INTENT(IN) :: y
		REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(IN) :: A
    REAL(KIND=8), DIMENSION(0:nQuadNodes), INTENT(IN) :: gllWeights,gllNodes

		! Outputs

		! Local variables
		INTEGER :: cdfid ! ID for netCDF file
		INTEGER, PARAMETER :: NDIMS = 3
		INTEGER :: ierr
    INTEGER :: idq,idt,idx,idy,dimids(NDIMS),idweight,idnode,idmu
    INTEGER :: x_dimid, y_dimid, t_dimid,node_dimid
		INTEGER, DIMENSION(1:NDIMS) :: start, count
		CHARACTER(len=8) :: nxname,xname,nyname,yname,qname,muname

		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tmp
		REAL(KIND=8), DIMENSION(1:nex*nxiplot) :: temp
		INTEGER :: i,j,l,m,nxout,nyout,ylvl,internalLvl

    SAVE cdfid, idq, t_dimid, start, count

    nxout = nex*nxiplot
    nyout = ney*netaplot

		IF(stat .eq. -1) THEN
			! Create netCDF file and time variables
			ierr = NF90_CREATE(TRIM(cdf_out),NF90_CLOBBER,cdfid)

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, "nt", ilvl+1, t_dimid)
      ierr = NF90_DEF_DIM(cdfid, "nnodes", nQuadNodes+1, node_dimid)

			ierr = NF90_DEF_VAR(cdfid, "qweights",NF90_FLOAT, node_dimid, idweight)
			ierr = NF90_DEF_VAR(cdfid, "qnodes",NF90_FLOAT, node_dimid, idnode)
			ierr = NF90_DEF_VAR(cdfid, "time", NF90_FLOAT, t_dimid,idt)

			ierr = NF90_ENDDEF(cdfid)

			! Calculate time at output levels (note ilvl=noutput)
			ALLOCATE(tmp(1:ilvl+1), STAT=ierr)
			DO i=0,ilvl
				tmp(i+1) = DBLE(i)*tval_in/DBLE(ilvl)
			ENDDO

			! Write t values
			ierr = NF90_PUT_VAR(cdfid,idt,tmp)
			ierr = NF90_PUT_VAR(cdfid,idweight,gllWeights)
      ierr = NF90_PUT_VAR(cdfid,idnode,gllNodes)

			DEALLOCATE(tmp, STAT=ierr)

			RETURN

		ELSEIF(stat .eq. 0) THEN
			! Create dimensions and variables for this level of runs (ilvl = p)
			start = 1
			count = 1

			! Define names of variables
			WRITE(nxname,'(a2,i1)') 'nx',ilvl
			WRITE(nyname,'(a2,i1)') 'ny',ilvl
			WRITE(xname, '(a1,i1)') 'x',ilvl
			WRITE(yname, '(a1,i1)') 'y',ilvl
			WRITE(qname, '(a1,i1)') 'Q',ilvl
      WRITE(muname, '(a2,i1)') 'mu',ilvl

			ierr = NF90_REDEF(cdfid)

			ierr = NF90_DEF_DIM(cdfid, TRIM(nxname), nxout, x_dimid)
			ierr = NF90_DEF_DIM(cdfid, TRIM(nyname), nyout, y_dimid)

			dimids(1) = x_dimid
			dimids(2) = y_dimid
			dimids(3) = t_dimid

			ierr = NF90_DEF_VAR(cdfid, TRIM(qname),NF90_FLOAT,dimids,idq)
			ierr = NF90_DEF_VAR(cdfid, TRIM(xname),NF90_FLOAT,x_dimid,idx)
			ierr = NF90_DEF_VAR(cdfid, TRIM(yname),NF90_FLOAT,y_dimid,idy)
      ierr = NF90_DEF_VAR(cdfid, TRIM(muname),NF90_FLOAT,idmu)

			ierr = NF90_enddef(cdfid)

			! Write x and y values
			ierr = NF90_PUT_VAR(cdfid, idx, x)
			ierr = NF90_PUT_VAR(cdfid, idy, y)
      ierr = NF90_PUT_VAR(cdfid,idmu,mu)

			start(3) = 1

		ELSEIF(stat .eq. 1) THEN
			ierr = NF90_CLOSE(cdfid)
			RETURN
		ENDIF

		! Write out concentration field
		count(1) = nxout
    DO ylvl=1,nyout
      start(2) = ylvl
      j = 1 + (ylvl-1)/netaplot
      internalLvl = mod(ylvl-1,netaplot)
      DO i=1,nex
        temp(1+(i-1)*nxiplot:i*nxiplot) = A(:,internalLvl,i,j)
      ENDDO!i
      ierr = NF90_PUT_VAR(cdfid,idq,temp,start,count)
    ENDDO!j

		! Increment t level
		start(3) = start(3) + 1

	END SUBROUTINE output2d

END PROGRAM EXECUTE
