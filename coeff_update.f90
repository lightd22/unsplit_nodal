! ==========================================
!
! Coefficent update for 2D Unsplit Nodal DG
! Transport equation with cartesian gridding
!
! Updates coefficents A(i,j,l,m) for each element
! See paper: Nair et al (2005) "A Discontinuous Galerkin Transport Scheme.." for a similar idea (modal impelementation)
!
! By: Devin Light ; April 2014
! ==========================================

SUBROUTINE coeff_update(A,u0,v0,gllNodes,gllWeights,gqWeights,lagrangeDeriv,time,dt,dxel,dyel,nex,ney,norder,nQuadNodes,&
                        gqOrder,lagGaussVal,nZSnodes,lagValsZS,dozshulimit,transient,doZSMaxCFL)
  IMPLICIT NONE

	! External functions
	REAL(KIND=8), EXTERNAL :: dadt ! RHS function for evolution ODE for expansion coefficent
	REAL(KIND=8), EXTERNAL :: vel_update ! Velocity update function

  ! Inputs
  INTEGER, INTENT(IN) :: nex,ney,norder,nQuadNodes,gqOrder,nZSnodes
  REAL(KIND=8), INTENT(IN) :: dxel,dyel,dt,time
  REAL(KIND=8), DIMENSION(0:nQuadNodes), INTENT(IN) :: gllNodes,gllWeights
  DOUBLE PRECISION, DIMENSION(0:gqOrder), INTENT(IN) :: gqWeights
  REAL(KIND=8), DIMENSION(0:norder,0:nQuadNodes), INTENT(IN) :: lagrangeDeriv
  REAL(KIND=8), DIMENSION(0:norder,0:gqOrder), INTENT(IN) :: lagGaussVal
  REAL(KIND=8), DIMENSION(0:norder,0:nZSNodes), INTENT(IN) :: lagValsZS
	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), INTENT(IN) :: u0,v0
  LOGICAL, INTENT(IN) :: dozshulimit,transient,doZSMaxCFL

  ! Outputs
  REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(INOUT) :: A

  ! Local Variables
  INTEGER :: i,j,l,m,stage,k,whichElement
  REAL(KIND=8) :: coef1,sentFhat,sentGhat
  REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney) :: A1,A2
  REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney) :: quadVals,u,v
  REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes) :: currElemQuad,currElemU,currElemV
  REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,0:nex+1,1:ney) :: edgeValsEW
  REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,1:nex,0:ney+1) :: edgeValsNS
  REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nex,1:ney) :: Fhat
  REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex,0:ney) :: Ghat
  REAL(KIND=8) :: error
  DOUBLE PRECISION, DIMENSION(1:nex,1:ney) :: elemAvg

	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: tmpArray

	REAL(KIND=4) :: t0,tf,t1

  INTERFACE
    SUBROUTINE limitMeanPositivity(coeffs,elemAvgs,lagGaussVal,gqWeights,gllWeights,&
                        nex,ney,nOrder,gqOrder)
        ! Modifies approximating polynomial so that element mean value remains non-negative
        ! after forward update according to Zhang and Shu (2010) Thm 5.2
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nex,ney,nOrder,gqOrder
        DOUBLE PRECISION, DIMENSION(0:gqOrder), INTENT(IN) :: gqWeights
        DOUBLE PRECISION, DIMENSION(0:nOrder), INTENT(IN) :: gllWeights
        DOUBLE PRECISION, DIMENSION(0:norder,0:gqOrder), INTENT(IN) :: lagGaussVal
        DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvgs
        ! Outputs
        DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder,1:nex,1:ney), INTENT(INOUT) :: coeffs

    END SUBROUTINE limitMeanPositivity

    SUBROUTINE limitNodePositivity(coeffs,elemAvgs,gllWeights,nex,ney,nOrder)
      ! Modifies approximating polynomial so that nodal values are non-negative for output
      USE testParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nex,ney,nOrder
      DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvgs
      DOUBLE PRECISION, DIMENSION(0:nOrder), INTENT(IN) :: gllWeights
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder,1:nex,1:ney), INTENT(INOUT) :: coeffs
    END SUBROUTINE limitNodePositivity

    SUBROUTINE computeAverages(elemAvgs,gllWeights,coeffs,nex,ney,nOrder,nQuad)
      USE testParameters
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nex,ney,nOrder,nQuad
      DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder,1:nex,1:ney), INTENT(IN) :: coeffs
      DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: gllWeights
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(OUT) :: elemAvgs
    END SUBROUTINE computeAverages
  END INTERFACE


	! ######################
	! Time step using Strong-Stability Preserving RK3
	! ######################
	coef1 = dt/dxel/dyel
  u = u0
  v = v0

  A1 = A

  DO stage=1,3

   ! Update velocities (if necessary)
   IF(transient) THEN
    SELECT CASE(stage)
      CASE(1)
      		u = u0*vel_update(time)
      		v = v0*vel_update(time)
      CASE(2)
        	u = u0*vel_update(time+dt)
        	v = v0*vel_update(time+dt)
      CASE(3)
      		u = u0*vel_update(time+dt/2d0)
      		v = v0*vel_update(time+dt/2d0)
     END SELECT
    ENDIF

    IF(dozshulimit) THEN
      ! Check element averages
      CALL computeAverages(elemAvg,gllWeights,A1,nex,ney,nOrder,nQuadNodes)
      IF(MINVAL(elemAvg) .lt. 0d0) THEN
        write(*,*) 'WARNING-- ELEMENT MEAN IS NEGATIVE BEFORE FWD STEP',stage
        write(*,*) minval(elemAvg)
      ENDIF
      IF(stage>1) THEN
!        CALL CPU_TIME(t0)
!        CALL limitMeanPositivity(A1,elemAvg,lagGaussVal,gqWeights,gllWeights,&
!                                nex,ney,nOrder,gqOrder)
        CALL limitNodePositivity(A1,elemAvg,gllWeights,nex,ney,nOrder)
!        CALL CPU_TIME(tf)
!        tf = tf - t0
!        write(*,*) 'tend=',tf
      ENDIF
    ENDIF !dozshulimit

    ! Update fluxes
    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,A1,nQuadNodes,norder,nex,ney)
    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuadNodes,norder,nex,ney)

    ! Forward step of SSPRK3
  	DO j=1,ney
    	DO i=1,nex
        currElemQuad = quadVals(:,:,i,j)
        currElemU = u(:,:,i,j)
        currElemV = v(:,:,i,j)
    		DO l=0,norder
      		DO m=0,norder
        		A2(l,m,i,j) = A1(l,m,i,j)+ &
    				coef1*dadt(i,j,l,m,currElemQuad,Fhat,Ghat,currElemU,currElemV,gllWeights,&
                       lagrangeDeriv,dxel,dyel,nQuadNodes,norder,nex,ney)
        	ENDDO !m
        ENDDO!l
      ENDDO !i
    ENDDO!j

		SELECT CASE(stage)
		CASE(1)
			A1 = A2
		CASE(2)
			A1 = 0.75D0*A + 0.25D0*A2
		CASE(3)
			A1 = A/3d0 + 2D0*A2/3D0
		END SELECT
    ENDDO !stage
    A = A1

    IF(dozshulimit) THEN
      ! Compute element average via quadrature at next time level
      CALL computeAverages(elemAvg,gllWeights,A,nex,ney,nOrder,nQuadNodes)
!      IF(minval(elemAvg) .lt. 0d0) THEN
!        write(*,*) 'WARNING-- ELEMENT MEAN IS NEGATIVE BEFORE LIMITING NODES'
!        write(*,*) minval(elemAvg)
        !STOP
!      ENDIF
      CALL limitNodePositivity(A,elemAvg,gllWeights,nex,ney,nOrder)
    ENDIF !dozshulimit

END SUBROUTINE coeff_update

REAL(KIND=8) FUNCTION dadt(i,j,l,m,quadVals,Fhat,Ghat,uIn,vIn,gllWeight,lagrangeDeriv,dxel,dyel,nQuadNodes,norder,nex,ney)
! Computes dadt for element (i,j) for use in SSPRK3 update step
! Does this in 3 steps: dadt = 2*(A+B+C) where
!	1) A+B = Int(Int( dyel*dPs/dxi*Pt*F+dxel*dPt/deta*Ps*G )) -- Interior contribution to change
!	2) C = -dxel*Int( Ps*(Ghat(eta=1) - (-1)**t*Ghat(eta=-1)) ) -- Contribution through North/South faces
!	3) D = -dyel*Int( Pt*(Fhat(xi=1) - (-1)**s*Fhat(xi=-1)) ) -- Contribution through East/West faces
!
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: i,j,l,m,nQuadNodes,norder,nex,ney
	REAL(KIND=8), INTENT(IN) :: dxel,dyel
 	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes), INTENT(IN) :: uIn,vIn,quadVals
	REAL(KIND=8), DIMENSION(0:norder,0:nQuadNodes), INTENT(IN) :: lagrangeDeriv
	REAL(KIND=8), DIMENSION(0:nQuadNodes), INTENT(IN) :: gllWeight
	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nex,1:ney), INTENT(IN) :: Fhat
	REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex,0:ney), INTENT(IN) :: Ghat

	! Local variables
	REAL(KIND=8) :: A,B,C,D

	! Step 1: Compute A & B
    A = (dyel/gllWeight(l))*SUM( gllWeight(:)*uIn(:,m)*quadVals(:,m)*lagrangeDeriv(l,:) )
    B = (dxel/gllWeight(m))*SUM( gllWeight(:)*vIn(l,:)*quadVals(l,:)*lagrangeDeriv(m,:) )

	! Step 2: Compute NS Flux Contribution
    C = 0D0
    IF(m .eq. nQuadNodes) THEN
        C = C + Ghat(l,i,j)
    ELSEIF( m .eq. 0) THEN
        C = C - Ghat(l,i,j-1)
    ENDIF

	C = (-1D0*dxel/gllWeight(m))*C

	! Step 3: Compute EW Flux Contribution
    D = 0D0
    IF(l .eq. nQuadNodes) THEN
        D = D + Fhat(m,i,j)
    ELSEIF( l .eq. 0) THEN
        D = D - Fhat(m,i-1,j)
    ENDIF

    D = (-1D0*dyel/gllWeight(l))*D

	dadt = (2D0)*(A+B+C+D)

END FUNCTION dadt

SUBROUTINE numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuadNodes,norder,nex,ney)
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: nQuadNodes,norder,nex,ney
  REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney) :: u, v
	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,0:nex+1,1:ney), INTENT(IN) :: edgeValsEW
	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,1:nex,0:ney+1), INTENT(IN) :: edgeValsNS

	! Outputs
	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nex,1:ney), INTENT(OUT) :: Fhat
	REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex,0:ney), INTENT(OUT) :: Ghat

	! Local Variables
	INTEGER :: i,j,k,l,m,s,t
	REAL(KIND=8), DIMENSION(0:nQuadNodes) :: u_tmp, v_tmp
	INTEGER :: which_el,which_sign

	DO j=1,ney
		DO i=1,nex

			u_tmp = u(nQuadNodes,:,i,j)
			v_tmp = v(:,nQuadNodes,i,j)

			DO k=0,nQuadNodes
				! Compute numerical flux through North edges (Ghat) (There are nQuadNodes+1 nodes per element)
				which_sign = 1-(1-INT(SIGN(1D0,v_tmp(k))) )/2
				which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2
				Ghat(k,i,j) = v_tmp(k)*edgeValsNS(which_sign,k,i,which_el)

				! Compute numerical flux through East edges (Fhat) (There are nQuadNodes+1 nodes per element)
				which_sign = 1-(1-INT(SIGN(1D0,u_tmp(k))) )/2
				which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
				Fhat(k,i,j) = u_tmp(k)*edgeValsEW(which_sign,k,which_el,j)
			ENDDO !k
		ENDDO !i
	ENDDO !j

    j = 0
    DO i=1,nex
			v_tmp = v(:,0,i,j+1)
			DO k=0,nQuadNodes
  			! Compute numerical flux through North edges (Ghat) (There are nQuadNodes+1 nodes per element)
  			which_sign = 1-(1-INT(SIGN(1D0,v_tmp(k))) )/2
  			which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2
    		Ghat(k,i,j) = v_tmp(k)*edgeValsNS(which_sign,k,i,which_el)
      ENDDO!k
    ENDDO!i

    i=0
    DO j=1,ney
      u_tmp = u(0,:,i+1,j)
      DO k=0,nQuadNodes
  			! Compute numerical flux through East edges (Fhat) (There are nQuadNodes+1 nodes per element)
  			which_sign = 1-(1-INT(SIGN(1D0,u_tmp(k))) )/2
  			which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
  			Fhat(k,i,j) = u_tmp(k)*edgeValsEW(which_sign,k,which_el,j)
	    ENDDO !k
    ENDDO !j

END SUBROUTINE numFlux

SUBROUTINE evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuadNodes,norder,nex,ney)
! Evaluate current expansion at interior quadrature locations and quadrature locations along EW and NS edges of each element
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: nQuadNodes,norder,nex,ney
	REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(IN) :: coeffs

	! Outputs
	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), INTENT(OUT) :: quadVals
	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,0:nex+1,1:ney), INTENT(OUT) :: edgeValsEW
	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,1:nex,0:ney+1), INTENT(OUT) :: edgeValsNS

	! Local Variables
	INTEGER :: i,j

	DO j=1,nex
		DO i=1,nex
      ! Evaluate expansion at interior GLL quad points
      quadVals(:,:,i,j) = coeffs(:,:,i,j)

	    ! Evaluate expansion at NS edges
			edgeValsNS(1,:,i,j) = coeffs(:,norder,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:),2))
			edgeValsNS(0,:,i,j) = coeffs(:,0,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:)*tmp2,2))

			! Evaluate expansion at EW edges
			edgeValsEW(1,:,i,j) = coeffs(norder,:,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:),1))
			edgeValsEW(0,:,i,j) = coeffs(0,:,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:)*tmp1,1))
		ENDDO !j
	ENDDO !i

	! Extend edge values periodically
  edgeValsEW(:,:,0,:) = edgeValsEW(:,:,nex,:)
  edgeValsEW(:,:,nex+1,:) = edgeValsEW(:,:,1,:)

  edgeValsNS(:,:,:,0) = edgeValsNS(:,:,:,ney)
  edgeValsNS(:,:,:,ney+1) = edgeValsNS(:,:,:,1)

END SUBROUTINE evalExpansion

REAL(KIND=8) FUNCTION vel_update(t)
		IMPLICIT NONE
		REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8) :: pi
    REAL(KIND=8), parameter :: t_period = 5.d0

    pi = DACOS(-1D0)
    vel_update = DCOS(pi*t/t_period)
END FUNCTION vel_update
