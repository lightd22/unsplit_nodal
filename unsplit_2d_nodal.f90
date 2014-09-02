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
	INTEGER :: start_res,polyOrder
	LOGICAL :: transient,DEBUG,doModalComparison
	REAL(KIND=8) :: muMAX

    DEBUG = .TRUE.
    doModalComparison = .FALSE.
    polyOrder = 4
	start_res = 8

!	write(*,*) '======'
!	write(*,*) 'TEST 0: Uniform field, deformation flow'
!	write(*,*) '======'
!	transient = .TRUE.
!	CALL test2d_nodal(100,start_res,start_res,2,3,20,0.01D0)

	write(*,*) '======'
	write(*,*) 'TEST 1: Uniform advection (u=v=1)'
	write(*,*) '======'
	transient = .FALSE.
    muMAX = 0.690D0
	CALL test2d_nodal(1,start_res,start_res,2,4,1,muMAX) !1D0/(2D0*4D0-1D0) !0.3D0/sqrt(2d0)
 !  CALL test2d_nodal(10,start_res,start_res,2,1,-1,0.970D0) !0.970D0
!   CALL test2d_nodal(10,start_res,start_res,2,2,40,0.500D0) !0.970D0

	write(*,*) '======'
	write(*,*) 'TEST 2: Smooth cosbell deformation'
	write(*,*) '======'
	transient = .TRUE.
    muMAX = 0.869D0
    write(*,*) 'muMAX=',muMAX
!	CALL test2d_nodal(6,start_res,start_res,2,4,1,muMAX) !1D0/(2D0*4D0-1D0)

!	write(*,*) '======'
!	write(*,*) 'TEST 3: Standard cosbell deformation'
!	write(*,*) '======'
!	transient = .TRUE.
!	CALL test2d_nodal(5,start_res,start_res,2,3,20,0.01D0) !1D0/(2D0*4D0-1D0)

!	write(*,*) '======'
!	write(*,*) 'TEST 4: Solid body rotation of cylinder'
!	write(*,*) '======'
!	transient = .FALSE.
!    start_res = 30
!	CALL test2d_nodal(2,start_res,start_res,2,3,20,0.75D0) !0.08D0

!	write(*,*) '======'
!	write(*,*) 'TEST 5: Square wave deformation'
!	write(*,*) '======'
!	transient = .TRUE.
!	CALL test2d_nodal(7,start_res,start_res,2,3,20,1D0/(2D0*4D0-1D0)) !1D0/(2D0*4D0-1D0)

!	write(*,*) '======'
!	write(*,*) 'TEST 6: Solid body rotation of cylinder (modified for frank)'
!	write(*,*) '======'
!	transient = .FALSE.
!	CALL test2d_nodal(3,start_res,start_res,2,3,nout,0.758D0) !1D0/(2D0*4D0-1D0)
!    muMAX = 0.5364D0!0.759D0/sqrt(2D0)
!    muMAX = 1.05D0!0.759D0/sqrt(2D0)
!	CALL test2d_nodal(11,start_res,start_res,2,3,50,muMAX)


CONTAINS
	SUBROUTINE test2d_nodal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
		REAL(KIND=8), INTENT(IN) :: maxcfl

		! Local variables
		INTEGER, DIMENSION(10) :: tmp_method
	    REAL(KIND=8), DIMENSION(nlevel) :: e1, e2, ei
		REAL(KIND=8) :: cnvg1, cnvg2, cnvgi,cons
		INTEGER :: nmethod, nmethod_final,imethod,ierr,nstep,nout
		INTEGER :: dgorder, norder,p

		LOGICAL :: dozshulimit

		CHARACTER(len=40) :: cdf_out

		INTEGER :: nex,ney,nxiplot,netaplot
		REAL(KIND=8) :: dxel,dyel,tfinal, tmp_umax, tmp_vmax, dxm, dym,dt, time,calculatedMu
        REAL(KIND=8), DIMENSION(1:2) :: xEdge,yEdge
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Leg,lagrangeDeriv,lagGaussVal,C0,C,tmpArray,tmpErr
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xQuad, yQuad
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: gllNodes,gllWeights,gaussNodes,gaussWeights,x_elcent,y_elcent,&
                                                   xplot,yplot,xiplot,etaplot,nodeSpacing,lambda
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: A,A0,& ! Coefficent array
													     u0,v0 ! Velocites within each element  
                                                
														
		INTEGER :: i,j,l,m,k,s,t,n

		CHARACTER(len=9) :: outdir
		REAL(KIND=4), DIMENSION(2) :: tstart,tend
		REAL(KIND=4) :: t0,tf

		REAL(KIND=8) :: tmp_qmax,tmp_qmin

		if(nlevel.lt.1) STOP 'nlev should be at least 1 in test2d_modal'

		nmethod_final = 1
		tmp_method = 0
		tmp_method(1) = 1
        tmp_method(2) = 2

		DO nmethod=1,nmethod_final
			imethod = tmp_method(nmethod)

			SELECT CASE(imethod)
				CASE(1)
				  WRITE(*,*) '2D Nodal, Unsplit, No limiting'
				  dozshulimit = .FALSE.
				  outdir = 'ndgunlim/'
				  norder = polyOrder
				  dgorder = norder !2*(norder+1)
				  WRITE(*,*) 'N=',norder,'Uses a total of',(norder+1)**2,'Lagrange basis polynomials'
				CASE(2)
				  WRITE(*,*) '2D Nodal, Unsplit, Zhang and Shu Limiting'
				  dozshulimit = .TRUE.
				  outdir = 'ndgzhshu/'
				  norder = polyOrder
				  dgorder = norder
				  WRITE(*,*) 'N=',norder,'Uses a total of',(norder+1)**2,'Lagrange basis polynomials'
			END SELECT

			! Initialize quadrature weights and nodes (only done once per call)
			ALLOCATE(gllnodes(0:dgorder),gllweights(0:dgorder),gaussNodes(0:dgorder),gaussWeights(0:dgorder),&
                    lagrangeDeriv(0:norder,0:dgorder),tmpArray(0:norder,0:norder),nodeSpacing(0:dgorder-1),&
                    lagGaussVal(0:norder,0:dgorder),lambda(0:dgorder),STAT=ierr)

            CALL gllquad_nodes(dgorder,gllNodes)
            CALL gllquad_weights(dgorder,gllNodes,gllWeights)
            CALL gaussquad_nodes(dgorder+1,gaussNodes)
            CALL gaussquad_weights(dgorder+1,gaussNodes,gaussWeights)

            nodeSpacing = gllNodes(1:dgorder)-gllNodes(0:dgorder-1)

			nxiplot = norder+1
			netaplot = norder+1

			ALLOCATE(xiplot(1:nxiplot), etaplot(1:netaplot), STAT=ierr)

			xiplot = gllNodes(:)
			etaplot = gllNodes(:)

            ! Fill matrix of basis derivatives at GLL nodes (note that this assumes lagrangeDeriv is square!)    
            CALL Dmat(dgorder,gllNodes,lagrangeDeriv)

            ! Fill maxtris of basis polynomials evaluated at Gauss quadrature nodes (used in Zhang and Shu Limiter)
            CALL baryWeights(lambda,gllnodes,dgorder)
            DO i=0,norder
                DO j=0,dgorder
                    lagGaussVal(i,j) = lagrange(gaussNodes(j),i,dgorder,gllNodes,lambda)
                ENDDO !j
            ENDDO !i

            xEdge(1) = 0D0
            xEdge(2) = 1D0
            if(ntest .eq. 11 .or. ntest .eq. 10) xEdge(1) = -1D0
            yEdge = xEdge

            write(*,*) 'Domain is: [',xEdge(1),',',xEdge(2),'].'
            write(*,*) 'Warning: Not all tests have been implemented for non-[0,1] square domain!'

			DO p=1,nlevel
				
				t0 = etime(tstart)

				nex = nex0*nscale**(p-1)
				ney = ney0*nscale**(p-1)
				dxel = (xEdge(2)-xEdge(1))/DBLE(nex)
				dyel = (yEdge(2)-yEdge(1))/DBLE(ney)


				ALLOCATE(x_elcent(1:nex),y_elcent(1:ney),A(1:nex,1:ney,0:norder,0:norder),A0(1:nex,1:ney,0:norder,0:norder),&
                         u0(1:nex,1:ney,0:dgorder,0:dgorder), v0(1:nex,1:ney,0:dgorder,0:dgorder), &
						 C0(1:nex,1:ney), C(1:nex,1:ney), xplot(1:nxiplot*nex),yplot(1:netaplot*ney),&
                         xQuad(1:nex,0:dgorder),yQuad(1:ney,0:dgorder), STAT=ierr)

				! Elements are ordered row-wise within the domain

				! A(i,j,l,m) : Coefficent array ; (i,j) = element, (l,m) = which Lagrange polynomial (q_ij = SUM(SUM(A(i,j,l,m)*L_l(xi)*L_m(eta))))
                ! This is also the solution at the interpolating GLL points
				! u(i,j,l,m), v(i,j,l,m) : Horizontal/vertical vel. array ; (i,j) = element, (l,m) = horiz,vertical location
				
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
                CALL init2d(ntest,nex,ney,dgorder,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,gllWeights,xEdge,yEdge,cdf_out,&
                            tfinal,xQuad,yQuad)

                A0 = A
				cdf_out = outdir // cdf_out

				! Store element averages for conservation estimation
                DO i=1,nex
                    DO j=1,ney
                        DO l=0,norder
                            tmpArray(l,:) = gllWeights(l)*gllWeights(:)*A0(i,j,l,:)
                        ENDDO!l
                        C0(i,j) = 0.25D0*dxel*dyel*SUM(tmpArray)
                    ENDDO !j
                ENDDO !i

				! Set up timestep
                IF(doModalComparison) THEN
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
                ELSE IF(doModalComparison) THEN
                    nstep = noutput*CEILING( (tfinal/maxcfl)*MAX(tmp_umax/dxm,tmp_vmax/dym)/DBLE(noutput) )
                    nout = noutput
                ELSE IF(DEBUG .and. p.gt.1) THEN
                    write(*,*) 'Debugging..'
                    nstep = nstep*2
                ELSE
                		nstep = noutput*CEILING( (tfinal/maxcfl)*(maxval(sqrt(u0**2 + v0**2))/min(dxm,dym))/DBLE(noutput) )
                		nout = noutput
                ENDIF !noutput

				dt = tfinal/DBLE(nstep)
                calculatedMu = maxval(sqrt(u0**2 + v0**2))*dt/min(dxm,dym)
!                calculatedMu = maxval(sqrt(u0**2 + v0**2))*dt/min(dxel,dyel)
                write(*,*) 'mu used =',calculatedMu
                !calculatedMu = (tmp_umax/dxm + tmp_vmax/dym)*dt


				IF(p .eq. 1) THEN ! Set up netCDF file
                    write(*,*) 'Maximum velocity: |u| = ',maxval(abs(sqrt(u0**2+v0**2)))
					CALL output2d(A0,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,dgorder,nxiplot,netaplot,tfinal,calculatedMu,cdf_out,nout,-1)
				ENDIF

                ! Set up variables for this value of p ; Write x, y, and initial conditions
				CALL output2d(A0,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,dgorder,nxiplot,netaplot,0D0,calculatedMu,cdf_out,p,0)

				! Time integration
				tmp_qmax = MAXVAL(A0)
				tmp_qmin = MINVAL(A0)

				time = 0D0
				DO n=1,nstep

!                    A0 = A
                    CALL coeff_update(A,u0,v0,gllNodes,gllWeights,gaussNodes,lagrangeDeriv,time,dt,dxel,dyel,nex,ney,&
                                      norder,dgorder,lagGaussVal,dozshulimit,transient)
 !                   write(*,*) 'CHANGE IN A',maxval(abs(A-A0))
					time = time + dt
	
					IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN
						! Write output
						CALL output2d(A,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,dgorder,nxiplot,netaplot,time,calculatedMu,cdf_out,p,2)
					ENDIF
					
					tmp_qmax = MAX(tmp_qmax,MAXVAL(A))
					tmp_qmin = MIN(tmp_qmin,MINVAL(A))

				ENDDO !n
                tf = etime(tend) - t0

                ! Store element averages for conservation estimation 
                DO i=1,nex
                    DO j=1,ney
                        DO l=0,norder
                            tmpArray(l,:) = gllWeights(l)*gllWeights(:)*A(i,j,l,:)
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
                            tmpArray(l,:) = gllWeights(l)*gllWeights(:)*ABS(A(i,j,l,:)-A0(i,j,l,:))
                        ENDDO!l
                        tmpErr(i,j) = SUM(tmpArray)
                    ENDDO !j
                ENDDO !i                
                e1(p) = 0.25D0*dxel*dyel*SUM(tmpErr)

                DO i=1,nex
                    DO j=1,ney
                        DO l=0,norder
                            tmpArray(l,:) = gllWeights(l)*gllWeights(:)*(ABS(A(i,j,l,:)-A0(i,j,l,:)))**2
                        ENDDO!l
                        tmpErr(i,j) = SUM(tmpArray)
                    ENDDO !j
                ENDDO !i                
                e2(p) = sqrt(0.25D0*dxel*dyel*SUM(tmpErr))

                ei(p) = MAXVAL(ABS(A(:,:,:,:) - A0(:,:,:,:)))


        			if (p.eq.1) then
		        	write(UNIT=6,FMT='(A125)') &
'   nex    ney     E1       E2         Einf      convergence rate  overshoot  undershoot   cons        cputime  nstep   tf   '
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
                    cons, tf, nstep,tfinal

				IF(p .eq. nlevel) THEN
                    CALL output2d(A,xplot,yplot,gllWeights,gllNodes,nex,ney,norder,dgorder,nxiplot,netaplot,time,&
                                  calculatedMu,cdf_out,p,1) ! Close netCDF files
				ENDIF
				DEALLOCATE(A,A0,x_elcent,y_elcent,xplot,yplot,u0,v0,tmpErr,xQuad,yQuad,C,C0,STAT=ierr)
			ENDDO
		ENDDO
		DEALLOCATE(gllnodes,gllweights,gaussnodes,gaussweights,xiplot,etaplot,tmpArray, tmpErr, STAT=ierr)

990    format(2i6,3e12.4,3f5.2,3e12.4,f8.2,i8,f8.2)

	END SUBROUTINE test2d_nodal


    SUBROUTINE init2d(ntest,nex,ney,dgorder,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,gllWeights,xEdge,yEdge,cdf_out,tfinal,&
                        xQuad,yQuad)
    ! ----
    !  Computes initial conditons for coefficient matrix (by simple evaluation of IC function) and initial velocities (via a streamfunction)
    !  Also sets the final time and output file name
    ! ----
        IMPLICIT NONE
        !Inputs
        INTEGER, INTENT(IN) :: ntest,nex,ney,dgorder,norder
        REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: x_elcent
        REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: y_elcent
        REAL(KIND=8), DIMENSION(0:dgorder) :: gllNodes,gllWeights
        REAL(KIND=8), DIMENSION(1:2), INTENT(IN) :: xEdge,yEdge
        !Outputs
        REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(OUT) :: A,u0,v0
        CHARACTER(LEN=40), INTENT(OUT) :: cdf_out
        REAL(KIND=8), INTENT(OUT) :: tfinal
        !Local Variables
        REAL(KIND=8), DIMENSION(1:nex,0:norder), INTENT(OUT) :: xQuad
        REAL(KIND=8), DIMENSION(1:ney,0:norder), INTENT(OUT) :: yQuad
        REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder,0:1) :: psiu,psiv
        REAL(KIND=8), DIMENSION(0:norder,0:norder) :: r
        REAL(KIND=8), DIMENSION(1:2) :: domainCenter
        REAL(KIND=8) :: PI,dxmin,dymin,dxel,dyel,xWidth,yWidth,xc,yc,spd
        INTEGER :: i,j,k

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
		dxmin = (dxel/2D0)*MINVAL(gllNodes(1:dgorder)-gllNodes(0:dgorder-1))
		dymin = (dyel/2D0)*MINVAL(gllNodes(1:dgorder)-gllNodes(0:dgorder-1))

        ! Quadrature locations, mapped to physical domain
        DO i=1,nex
            xQuad(i,:) = x_elcent(i) + (dxel/2d0)*gllNodes(:)
        ENDDO!i
        DO j=1,ney
            yQuad(j,:) = y_elcent(j) + (dyel/2d0)*gllNodes(:)
        ENDDO!j

        ! Fill streamfunction array
        SELECT CASE(ntest)
            CASE(1,10) ! uniform velocity u = v = 1
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                        psiu(i,j,k,:,1) = (yQuad(j,:) + dymin/2d0) - xQuad(i,k)
                        psiu(i,j,k,:,0) = (yQuad(j,:) - dymin/2d0) - xQuad(i,k)

                        psiv(i,j,k,:,1) = yQuad(j,:) - (xQuad(i,k) + dxmin/2d0)
                        psiv(i,j,k,:,0) = yQuad(j,:) - (xQuad(i,k) - dxmin/2d0)
                        ENDDO !k
                    ENDDO !j
                ENDDO !i
            CASE(2:4,11) ! Solid body rotation
                ! pi*( (x-0.5)**2 + (y-0.5)**2 )
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
!                        psiu(i,j,k,:,1) = PI*( (xQuad(i,k)-domainCenter(1))**2 + (yQuad(j,:) + dymin/2d0 -domainCenter(2))**2 )
!                        psiu(i,j,k,:,0) = PI*( (xQuad(i,k)-domainCenter(1))**2 + (yQuad(j,:) - dymin/2d0 -domainCenter(2))**2 ) 
                        
!                        psiv(i,j,k,:,1) = PI*( (xQuad(i,k) + dxmin/2d0-domainCenter(1))**2 + (yQuad(j,:) -domainCenter(2))**2 )
!                        psiv(i,j,k,:,0) = PI*( (xQuad(i,k) - dymin/2d0-domainCenter(1))**2 + (yQuad(j,:) -domainCenter(2))**2 ) 

                     psiu(i,j,k,:,1) = spd*0.5D0*( (xQuad(i,k)-domainCenter(1))**2 + (yQuad(j,:) + dymin/2d0 -domainCenter(2))**2 )
                     psiu(i,j,k,:,0) = spd*0.5D0*( (xQuad(i,k)-domainCenter(1))**2 + (yQuad(j,:) - dymin/2d0 -domainCenter(2))**2 ) 
                        
                     psiv(i,j,k,:,1) = spd*0.5D0*( (xQuad(i,k) + dxmin/2d0-domainCenter(1))**2 + (yQuad(j,:) -domainCenter(2))**2 )
                     psiv(i,j,k,:,0) = spd*0.5D0*( (xQuad(i,k) - dymin/2d0-domainCenter(1))**2 + (yQuad(j,:) -domainCenter(2))**2 ) 

                        ENDDO !k
                    ENDDO !j
                ENDDO !i
            CASE(5:7,100) ! Leveque deformation flow
				!(1/pi)*sin(pi*x(i))**2 * sin(pi*y(j))**2
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                        psiu(i,j,k,:,1) = (1D0/PI)*(DSIN(PI*xQuad(i,k))**2 * DSIN(PI*(yQuad(j,:)+dymin/2d0))**2)
                        psiu(i,j,k,:,0) = (1D0/PI)*(DSIN(PI*xQuad(i,k))**2 * DSIN(PI*(yQuad(j,:)-dymin/2d0))**2)

                        psiv(i,j,k,:,1) = (1D0/PI)*(DSIN(PI*(xQuad(i,k)+dxmin/2d0))**2 * DSIN(PI*yQuad(j,:))**2)
                        psiv(i,j,k,:,0) = (1D0/PI)*(DSIN(PI*(xQuad(i,k)-dxmin/2d0))**2 * DSIN(PI*yQuad(j,:))**2)
                        ENDDO !k
                    ENDDO !j
                ENDDO !i
        END SELECT !ntest

		! Compute velocity from stream functions
        DO i=1,nex
            DO j=1,ney
                u0(i,j,:,:) = (psiu(i,j,:,:,1)-psiu(i,j,:,:,0) )/dymin
                v0(i,j,:,:) = -(psiv(i,j,:,:,1)-psiv(i,j,:,:,0) )/dxmin
            ENDDO !j
        ENDDO !i

        ! Initialize coefficent array
        SELECT CASE(ntest)
            CASE(1) ! uniform advection of a sine wave
                cdf_out = 'dg2d_sine_adv.nc'
                tfinal = 1D0
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                            A(i,j,k,:) = DSIN(2D0*PI*xQuad(i,k))*DSIN(2D0*PI*yQuad(j,:))
                        ENDDO !k
                    ENDDO!j
                ENDDO!i

            CASE(2) ! solid body rotation of a cylinder
                cdf_out = 'dg2d_rot_cylinder.nc'
                tfinal = 2D0*PI
                A = 0D0
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                            r(k,:) = SQRT((xQuad(i,k)-0.25d0)**2 + (yQuad(j,:)-0.5d0)**2)
                        ENDDO !k
                        WHERE(r .lt. 0.125D0)
                            A(i,j,:,:) = 1D0
                        END WHERE
                    ENDDO!j
                ENDDO!i

            CASE(3,10:11) ! solid body rotation of a cylinder (comparison to frank's code)
                cdf_out = 'dg2d_rot_cylinder_modified.nc'
                tfinal = 1D0
                A = 0D0
                xc = xEdge(1)+xWidth/4D0
                yc = yEdge(1)+yWidth/2D0

                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                            r(k,:) = SQRT((xQuad(i,k)-xc)**2 + (yQuad(j,:)-yc)**2)
                        ENDDO !k
                        WHERE(r .lt. 0.25D0)
                            A(i,j,:,:) = 1D0
                        END WHERE
                    ENDDO!j
                ENDDO!i

            CASE(5) ! standard cosbell deformation
				cdf_out = 'dg2d_def_cosbell.nc'
				tfinal = 5D0
                A = 0D0
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder ! Fill distance array for (i,j) element
                            r(k,:) = 4D0*SQRT( (xQuad(i,k)-0.25D0)**2 + (yQuad(j,:)-0.25D0)**2 )
                        ENDDO !k
                        WHERE(r .lt. 1D0)                            
                            A(i,j,:,:) = (0.5d0*(1.d0 + DCOS(PI*r(:,:))))
                        END WHERE
                    ENDDO !j
                ENDDO !i

            CASE(6) ! smoother cosbell deformation
				cdf_out = 'dg2d_smth_cosbell.nc'
				tfinal = 5D0
                A = 0D0
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder ! Fill distance array for (i,j) element
                            r(k,:) = 3D0*SQRT( (xQuad(i,k)-0.4D0)**2 + (yQuad(j,:)-0.4D0)**2 )
                        ENDDO !k
                        WHERE(r .lt. 1D0)                            
                            A(i,j,:,:) = (0.5d0*(1.d0 + DCOS(PI*r(:,:))))**3
                        END WHERE
                    ENDDO !j
                ENDDO !i
        END SELECT !ntest

        if(ntest .eq. 10) then
            tfinal = xWidth*10D0*tfinal*5
        elseif(ntest .eq. 11) then
            tfinal = 3D0*2D0*PI
        elseif(ntest .eq. 6) then
            tfinal = 5D0*tfinal
        endif

    END SUBROUTINE init2d

	SUBROUTINE output2d(A,x,y,gllWeights,gllNodes,nex,ney,norder,dgorder,nxiplot,netaplot,tval_in,mu,cdf_out,ilvl,stat)
		IMPLICIT NONE

		! Inputs
		INTEGER, INTENT(IN) :: norder,dgorder,nex,ney,nxiplot,netaplot,stat,ilvl
		CHARACTER(len=40), INTENT(IN) :: cdf_out
		REAL(KIND=8), INTENT(IN) :: tval_in,mu
		REAL(KIND=8), DIMENSION(1:nex*nxiplot), INTENT(IN) :: x
		REAL(KIND=8), DIMENSION(1:ney*netaplot), INTENT(IN) :: y
		REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(IN) :: A
        REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: gllWeights,gllNodes
		
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
            ierr = NF90_DEF_DIM(cdfid, "nnodes", dgorder+1, node_dimid)

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
                temp(1+(i-1)*nxiplot:i*nxiplot) = A(i,j,:,internalLvl)
            ENDDO!i
            ierr = NF90_PUT_VAR(cdfid,idq,temp,start,count)
        ENDDO!j
		
		! Increment t level 
		start(3) = start(3) + 1 

	END SUBROUTINE output2d

END PROGRAM EXECUTE
