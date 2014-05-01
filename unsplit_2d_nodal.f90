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
	INTEGER :: ntest,start_res
	LOGICAL :: transient

	write(*,*) '======'
	write(*,*) 'TEST 0: Uniform field, deformation flow'
	write(*,*) '======'
	transient = .TRUE.
	start_res = 8
!	CALL test2d_nodal(100,start_res,start_res,2,3,20,0.01D0)

	write(*,*) '======'
	write(*,*) 'TEST 1: Uniform advection (u=v=1)'
	write(*,*) '======'
	transient = .FALSE.
	start_res = 8 ! Number of elements in each direction
	CALL test2d_nodal(1,start_res,start_res,2,2,20,0.05D0) !1D0/(2D0*4D0-1D0) !0.3D0/sqrt(2d0)

	write(*,*) '======'
	write(*,*) 'TEST 2: Smooth cosbell deformation'
	write(*,*) '======'
	transient = .TRUE.
!	CALL test2d_nodal(6,start_res,start_res,2,2,20,0.01D0) !1D0/(2D0*4D0-1D0)

	write(*,*) '======'
	write(*,*) 'TEST 3: Standard cosbell deformation'
	write(*,*) '======'
	transient = .TRUE.
!	CALL test2d_nodal(5,start_res,start_res,2,3,20,1D0/(2D0*4D0-1D0))

	write(*,*) '======'
	write(*,*) 'TEST 4: Square wave deformation'
	write(*,*) '======'
	transient = .TRUE.
!	CALL test2d_nodal(7,start_res,start_res,2,3,20,1D0/(2D0*4D0-1D0)) !1D0/(2D0*4D0-1D0)

CONTAINS
	SUBROUTINE test2d_nodal(ntest,nex0,ney0,nscale,nlevel,noutput,maxcfl)
		IMPLICIT NONE
		! Inputs
		INTEGER, INTENT(IN) :: ntest,nex0,ney0,nscale,nlevel,noutput
		REAL(KIND=8), INTENT(IN) :: maxcfl

		! Local variables
		INTEGER, DIMENSION(10) :: tmp_method
	    REAL(KIND=8), DIMENSION(nlevel) :: e1, e2, ei
		REAL(KIND=8) :: cnvg1, cnvg2, cnvgi
		INTEGER :: nmethod, nmethod_final,imethod,ierr,nstep,nout
		INTEGER :: dgorder, norder,p

		LOGICAL :: dozshulimit

		CHARACTER(len=40) :: cdf_out

		INTEGER :: nex,ney,nxiplot,netaplot
		REAL(KIND=8) :: dxel,dyel,tfinal, tmp_umax, tmp_vmax, dxm, dym,dt, time
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Leg, lagrangeDeriv,C0,C
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: gllNodes,gllWeights,gaussNodes,gaussWeights,x_elcent,y_elcent,&
                                                   xplot,yplot,xiplot,etaplot
		REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: A,A0,& ! Coefficent array
													     u0,v0 ! Velocites within each element  
														
		INTEGER :: i,j,l,m,k,s,t,n

		CHARACTER(len=9) :: outdir
		REAL(KIND=4), DIMENSION(2) :: tstart,tend
		REAL(KIND=4) :: t0,tf

		REAL(KIND=8) :: tmp_qmax,tmp_qmin

		if(nlevel.lt.1) STOP 'nlev should be at least 1 in test2d_modal'

		nmethod_final = 2
		tmp_method = 0
		tmp_method(1) = 1
        tmp_method(2) = 2

		DO nmethod=1,nmethod_final
			imethod = tmp_method(nmethod)

			SELECT CASE(imethod)
				CASE(1)
				  WRITE(*,*) '2D Nodal, Unsplit, No limiting'
				  dozshulimit = .FALSE.
				  outdir = 'mdgunlim/'
				  norder = 4
				  dgorder = norder !2*(norder+1)
				  WRITE(*,*) 'N=',norder,'Uses a total of',(norder+1)**2,'Legendre basis polynomials'
				CASE(2)
				  WRITE(*,*) '2D Nodal, Unsplit, Zhang and Shu Limiting'
				  dozshulimit = .TRUE.
				  outdir = 'mdgzhshu/'
				  norder = 4
				  dgorder = norder
				  WRITE(*,*) 'N=',norder,'Uses a total of',(norder+1)**2,'Legendre basis polynomials'
			END SELECT

			! Initialize quadrature weights and nodes (only done once per call)
			ALLOCATE(gllnodes(0:dgorder),gllweights(0:dgorder),gaussNodes(0:dgorder),gaussWeights(0:dgorder),&
                    lagrangeDeriv(0:norder,0:dgorder),STAT=ierr)

            CALL gllquad_nodes(dgorder+1,gllNodes)
            CALL gllquad_weights(dgorder+1,gllNodes,gllWeights)
            CALL gaussquad_nodes(dgorder+1,gaussNodes)
            CALL gaussquad_weights(dgorder+1,gaussNodes,gaussWeights)

			nxiplot = norder+1
			netaplot = norder+1

			ALLOCATE(xiplot(1:nxiplot), etaplot(1:netaplot), STAT=ierr)

			xiplot = gllNodes(:)
			etaplot = gllNodes(:)

            ! Fill matrix of basis derivatives at GLL nodes (note that this assumes lagrangeDeriv is square!)    
            CALL Dmat(dgorder,gllNodes,lagrangeDeriv)

			DO p=1,nlevel
				
				t0 = etime(tstart)

				nex = nex0*nscale**(p-1)
				ney = ney0*nscale**(p-1)
				dxel = 1D0/DBLE(nex)
				dyel = 1D0/DBLE(ney)

				ALLOCATE(x_elcent(1:nex),y_elcent(1:ney),A(1:nex,1:ney,0:norder,0:norder),A0(1:nex,1:ney,0:norder,0:norder),&
                         u0(1:nex,1:ney,0:dgorder,0:dgorder), v0(1:nex,1:ney,0:dgorder,0:dgorder),&
						 C0(1:nex,1:ney), C(1:nex,1:ney), xplot(1:nxiplot*nex),yplot(1:netaplot*ney), STAT=ierr)

				! Elements are ordered row-wise within the domain

				! A(i,j,l,m) : Coefficent array ; (i,j) = element, (l,m) = which Lagrange polynomial (q_ij = SUM(SUM(A(i,j,l,m)*L_l(xi)*L_m(eta))))
                ! This is also the solution at the interpolating GLL points
				! u(i,j,l,m), v(i,j,l,m) : Horizontal/vertical vel. array ; (i,j) = element, (l,m) = horiz,vertical location
				
				! Initialize x- and y- grids and xi- and eta- plotting grids
				x_elcent(1) = dxel/2D0
				DO i=2,nex
					x_elcent(i) = x_elcent(i-1)+dxel
				ENDDO !i

				y_elcent(1) = dyel/2D0
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

                CALL init2d(ntest,nex,ney,dgorder,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,gllWeights,cdf_out,tfinal)
                write(*,*) MAXVAL(v0),MINVAL(v0)
                A0 = A
				cdf_out = outdir // cdf_out

				! Store element averages for conservation estimation

				! Set up timestep
				dxm = dxel
				dym = dyel

				tmp_umax = MAXVAL(u0)
				tmp_vmax = MAXVAL(v0)

				IF(noutput .eq. -1) THEN
					nstep = CEILING( (tfinal/maxcfl)*(tmp_umax + tmp_vmax)/(dxm) )
					nout = nstep
				ELSE
					nstep = noutput*CEILING( (tfinal/maxcfl)*(sqrt(tmp_umax**2 + tmp_vmax**2))/(dxm)/DBLE(noutput) )
					nout = noutput
				ENDIF

				dt = tfinal/DBLE(nstep)

				IF(p .eq. 1) THEN ! Set up netCDF file
!					CALL output2d(q0,xplot,yplot,nxplot,nyplot,tfinal,cdf_out,nout,-1)
				ENDIF

!				CALL output2d(q0,xplot,yplot,nxplot,nyplot,0D0,cdf_out,p,0) ! Set up variables for this value of p ; Write x, y, and initial conditions

				! Time integration
				tmp_qmax = MAXVAL(A0)
				tmp_qmin = MINVAL(A0)

				time = 0D0
!				DO n=1,nstep

!					CALL coeff_update(q,A,u_tmp,v_tmp,uedge_tmp,vedge_tmp,qnodes,qweights,Leg,dLL,LL,L_xi_plot,L_eta_plot,dxel,dyel,& 
!									  dt,dgorder,norder,nxplot,nyplot,nex,ney,nxiplot,netaplot,transient,time,dozshulimit)

					! Store element averages for conservation estimation (for modal DG these are just the 0th order coeffs)

					time = time + dt
	
					IF((MOD(n,nstep/nout).eq.0).OR.(n.eq.nstep)) THEN
						! Write output
!						CALL output2d(q,xplot,yplot,nxplot,nyplot,time,cdf_out,p,2)
					ENDIF
					
					tmp_qmax = MAX(tmp_qmax,MAXVAL(A))
					tmp_qmin = MIN(tmp_qmin,MINVAL(A))

!				ENDDO !nstep

				IF(p .eq. nlevel) THEN
!					CALL output2d(q,xplot,yplot,nxplot,nyplot,tfinal,cdf_out,p,1) ! Close netCDF files
				ENDIF
				DEALLOCATE(A,A0,x_elcent,y_elcent,xplot,yplot,u0,v0, STAT=ierr)
			ENDDO
		ENDDO
		DEALLOCATE(gllnodes,gllweights,gaussnodes,gaussweights,xiplot,etaplot, STAT=ierr)

990    format(2i6,3e12.4,3f5.2,3e12.4,f8.2,i8)

	END SUBROUTINE test2d_nodal
    SUBROUTINE init2d(ntest,nex,ney,dgorder,norder,A,u0,v0,x_elcent,y_elcent,gllNodes,gllWeights,cdf_out,tfinal)
        IMPLICIT NONE
        !Inputs
        INTEGER, INTENT(IN) :: ntest,nex,ney,dgorder,norder
        REAL(KIND=8), DIMENSION(1:nex), INTENT(IN) :: x_elcent
        REAL(KIND=8), DIMENSION(1:ney), INTENT(IN) :: y_elcent
        REAL(KIND=8), DIMENSION(0:dgorder) :: gllNodes,gllWeights
        !Outputs
        REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(OUT) :: A,u0,v0
        CHARACTER(LEN=40), INTENT(OUT) :: cdf_out
        REAL(KIND=8), INTENT(OUT) :: tfinal
        !Local Variables
        REAL(KIND=8), DIMENSION(1:nex,0:norder) :: xQuad
        REAL(KIND=8), DIMENSION(1:ney,0:norder) :: yQuad
        REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder,0:1) :: psiu,psiv
        REAL(KIND=8) :: PI,dxmin,dymin,dxel,dyel
        INTEGER :: i,j,k

		PI = DACOS(-1D0)

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
            CASE(1) ! uniform velocity u = v = 1
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
            CASE(1)
                cdf_out = 'dg2d_sine_adv.nc'
                tfinal = 1D0
                DO i=1,nex
                    DO j=1,ney
                        DO k=0,norder
                            A(i,j,k,:) = DSIN(2D0*PI*xQuad(i,k))*DSIN(2D0*PI*yQuad(j,:))
                        ENDDO !k
                    ENDDO!j
                ENDDO!i
        END SELECT !ntest

    END SUBROUTINE init2d

END PROGRAM EXECUTE
