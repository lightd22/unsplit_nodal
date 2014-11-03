!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light Oct 2014
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(A,nelem,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u0,lagDeriv,time,dt,transient)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! External Functions
	REAL(KIND=DOUBLE), EXTERNAL :: split_dadt ! RHS function for evolution ODE for kth expansion coefficent
	REAL(KIND=8), EXTERNAL :: split_vel_update ! Velocity update function

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nQuadNodes,nOrder
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt,time
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes), INTENT(IN):: gllNodes,gllWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,0:nQuadNodes), INTENT(IN) :: lagDeriv
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nQuadNodes), INTENT(IN) :: u0
    LOGICAL, INTENT(IN) :: transient

    ! Outputs
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nOrder), INTENT(INOUT) :: A 

    ! Local Variables
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nOrder) :: AFwd,AStar
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nQuadNodes):: u,quadVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem,0:nQuadNodes) :: utmp
    REAL(KIND=DOUBLE), DIMENSION(0:nelem+1,0:1) :: edgeVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flx
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes) :: currLagDeriv,currElemQuad,currElemVel

    INTEGER :: k,j,stage,stat,l

    AStar = A
    u = u0
    ! Do SSPRK3 Update
!    DO stage=1,3
    DO stage=1,3
        ! Update velocities (if necessary) 
        IF(transient) THEN
            SELECT CASE(stage)
                CASE(1)
                		u = u0*split_vel_update(time)
                CASE(2)
	                	u = u0*split_vel_update(time+dt)
                CASE(3)
                		u = u0*split_vel_update(time+dt/2d0)
            END SELECT
        ENDIF ! transient
        utmp(1:nelem,:) = u
        utmp(0,:) = u(nelem,:)

        CALL split_evalExpansion(quadVals,edgeVals,AStar,nelem,nQuadNodes,nOrder)
        CALL split_numFlux(flx,edgeVals,utmp,nelem,nQuadNodes)

        ! Forward Step
        DO j=1,nelem
            currElemQuad = quadVals(j,:)
            currElemVel = utmp(j,:)
            DO k=0,nOrder
                currLagDeriv = lagDeriv(k,:)
                AFwd(j,k) = AStar(j,k) + (dt/dxel)*split_dadt(currElemQuad,flx,currElemVel,gllWeights,&
                                                              currLagDeriv,k,j,nelem,nQuadNodes)
            ENDDO !k
        ENDDO ! j

		SELECT CASE(stage)
		CASE(1)
			AStar = AFwd
		CASE(2)
			AStar = 0.75D0*A + 0.25D0*AFwd
		CASE(3)
			AStar = A/3d0 + 2D0*AFwd/3D0
		END SELECT
       
    ENDDO !stage
    A = AStar

END SUBROUTINE nDGsweep

SUBROUTINE split_evalExpansion(quadVals,edgeVals,AIn,nelem,nNodes,nOrder)
! Evaluate ansatz solution at quadrature nodes and edges
! Note:  Assumes basis is Lagrange interpolating polynomials at quad nodes -> phi_j (xi_k) = qIn(k,j)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,nOrder
    REAL(KIND=8), DIMENSION(1:nelem,0:nOrder), INTENT(INOUT) :: AIn

    ! Outputs
    REAL(KIND=8), DIMENSION(1:nelem,0:nNodes), INTENT(OUT) :: quadVals
    REAL(KIND=8), DIMENSION(0:nelem+1,0:1), INTENT(OUT) :: edgeVals

    ! Local Variables

    ! Ansatz value at quad locations for the jth element is just coefficient value
    quadVals = AIn
    
    edgeVals(1:nelem,0) = AIn(:,0) ! left edge value is just left-most coefficient
    edgeVals(1:nelem,1) = AIn(:,nNodes) ! right edge value is just right-most coefficient

	! Extend edgeVals periodically
	edgeVals(0,:) = edgeVals(nelem,:)
	edgeVals(nelem+1,:) = edgeVals(1,:)

END SUBROUTINE split_evalExpansion

SUBROUTINE split_numFlux(flx,edgeVals,uin,nelem,nNodes)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,nNodes
	REAL(KIND=DOUBLE), DIMENSION(0:nelem+1,0:1), INTENT(IN) :: edgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nelem,0:nNodes), INTENT(IN) :: uin

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flx

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flx(j) = 0.5D0*edgeVals(j,1)*(uin(j,nNodes)+DABS(uin(j,nNodes)))+0.5D0*edgeVals(j+1,0)*(uin(j,nNodes)-DABS(uin(j,nNodes)))
	ENDDO
END SUBROUTINE split_numFlux

REAL(KIND=8) FUNCTION split_dadt(quadVals,flx,u,qWeights,lagDeriv,k,j,nelem,nNodes)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,k,j
    REAL(KIND=8), DIMENSION(0:nNodes) :: quadVals,u,qWeights,lagDeriv
    REAL(KIND=8), DIMENSION(0:nelem) :: flx

    split_dadt = SUM( u(:)*quadVals(:)*lagDeriv(:)*qWeights(:) )

    IF( k .eq. 0) THEN
        split_dadt = split_dadt + flx(j-1)
    ENDIF
    IF( k .eq. nNodes) THEN
        split_dadt = split_dadt - flx(j)
    ENDIF
    split_dadt = 2D0*split_dadt/qWeights(k)

END FUNCTION split_dadt

SUBROUTINE split_polyMod(qIn,qNodes,qWeights,nelem,nNodes,nQuad,lagVals,stat)
! Rescales DG polynomial around element averages for positivity, based on Zhang and Shu (2010)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
    REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: qWeights,qNodes
    REAL(KIND=8), DIMENSION(0:nNodes,0:nQuad), INTENT(IN) :: lagVals
    ! Local Variables
    INTEGER :: j,l,stat
    REAL(KIND=8), DIMENSION(0:nQuad) :: qVals
    REAL(KIND=8) :: avgVal,theta,valMin,valMax,eps

    eps = epsilon(1d0)
!    eps = 0D0
    ! First check incoming element averages
    IF(stat == 1) THEN
        DO j=1,nelem
            ! -- First evaluate DG polynomial at quad locations
            ! Note: This is not necessary if the quad nodes used are same nodes that the basis is interpolating --
            !       In this case the polynomial coefficients are the nodal values. Two exceptions are the edge
            !       values (-1 and 1) which are always part of a GLL integration. These may be read from coefficients directly
            qVals(0) = qIn(0,j)
            qVals(nQuad) = qIn(nNodes,j)
            DO l=1,nQuad-1
                qVals(l) = SUM(qIn(:,j)*lagVals(:,l))
            ENDDO !l
    
            avgVal = 0.5D0*SUM( qWeights(:)*qVals(:) )
!            valMin = MINVAL(qVals(:))-eps
             valMin = MIN( MINVAL(qIn(1:nNodes-1,j)) , MINVAL(qVals(:)) ) - eps

            ! -- Compute rescale factor
            theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )

            qIn(:,j) = theta*(qIn(:,j)-avgVal) + avgVal

            IF(avgVal .lt. 0D0) THEN
                qVals(0) = qIn(0,j)
                qVals(nQuad) = qIn(nNodes,j)
                DO l=1,nQuad-1
                    qVals(l) = SUM(qIn(:,j)*lagVals(:,l))
                ENDDO !l

                write(*,'(A,E10.4)') '   Element AVG in Z&S Limiter is negative!! Negative average value = ',avgVal
                write(*,'(A,E10.4)') '   Minimum value = ',valMin
                write(*,'(A,F5.3)') '   Rescaling factor theta = ',theta
                write(*,'(A,E10.4)') '   New Average = ', 0.5D0*SUM( qWeights(:)*qVals(:) )
                write(*,'(A,E10.4)') '   New minimum value = ', MINVAL(qVals(:))
!                STOP
            ENDIF

        ENDDO !j
    ELSE
        DO j=1,nelem
            ! (Optional) Use nodes themselves to do limiting, rather than evaluating polynomial multiple times
            avgVal = 0.5D0*SUM( qWeights(:)*qIn(:,j) )
            valMin = MINVAL(qIn(:,j))-eps

            IF(avgVal .lt. 0D0) THEN
                write(*,*) 'Element AVG in Z&S Limiter is negative!! Negative average value = ',avgVal
            ENDIF

            ! -- Compute rescale factor
            theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )
    !       theta = MIN( abs(avgVal/(valMin-avgVal)),abs((1D0-avgVal)/(valMax-avgVal)),1D0 )

            ! -- Rescale polynomial
            qIn(:,j) = theta*(qIn(:,j)-avgVal) + avgVal
        ENDDO !j
    ENDIF !stat
    
END SUBROUTINE split_polyMod

REAL(KIND=8) FUNCTION split_vel_update(t)
		IMPLICIT NONE
		REAL(KIND=8), INTENT(IN) :: t
	    REAL(KIND=8) :: pi
	    REAL(KIND=8), parameter :: t_period = 5.d0

	    pi = DACOS(-1D0)
	    split_vel_update = DCOS(pi*t/t_period)
END FUNCTION split_vel_update
