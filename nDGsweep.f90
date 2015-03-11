!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light Oct 2014
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(A,nelem,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u0,lagDeriv,time,dt,transient,&
                    dozshulimit,dofctlimit,nZSnodes,quadZSWeights,lagValsZS)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! External Functions
	REAL(KIND=DOUBLE), EXTERNAL :: split_dadt ! RHS function for evolution ODE for kth expansion coefficent
	REAL(KIND=8), EXTERNAL :: split_vel_update ! Velocity update function

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nQuadNodes,nOrder,nZSNodes
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt,time
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes), INTENT(IN):: gllNodes,gllWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,0:nQuadNodes), INTENT(IN) :: lagDeriv
    REAL(KIND=DOUBLE), DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,0:nZSNodes), INTENT(IN) :: lagValsZS
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nQuadNodes), INTENT(IN) :: u0
    LOGICAL, INTENT(IN) :: transient,dozshulimit,dofctlimit

    ! Outputs
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nOrder), INTENT(INOUT) :: A

    ! Local Variables
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nOrder) :: AFwd,AStar
    REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nQuadNodes):: u,quadVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem,0:nQuadNodes) :: utmp
    REAL(KIND=DOUBLE), DIMENSION(0:nelem+1,0:1) :: edgeVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flx
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes) :: currLagDeriv,currElemQuad,currElemVel

    REAL(KIND=DOUBLE) :: elemAvg

    INTEGER :: k,j,stage,stat,l

    u = u0

    IF(dozshulimit) THEN
        IF(nZSnodes .eq. nQuadNodes) THEN
            CALL split_polyMod(A,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,0)
        ELSE
            CALL split_polyMod(A,quadZSWeights,nelem,nOrder,nZSNodes,lagValsZS,1)
        ENDIF
    ENDIF
    AStar = A

    ! Do SSPRK3 Update
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

        IF(dofctlimit) THEN ! Use FCT adjustment for element avg nonnegativity
          CALL FLUXCOR(Astar,Astar,flx,gllWeights,dxel,dt,nelem,nOrder,1)
          !CALL FLUXCOR(Astar,Astar,flx,gllWeights,dxel,dt,nelem,nOrder,stage)
        ENDIF !dofctlimit

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

      IF(dozshulimit) THEN
        IF(nZSnodes .eq. nQuadNodes) THEN
!          CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,0)
         CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,2)
        ELSE
          CALL split_polyMod(Astar,quadZSWeights,nelem,nOrder,nZSNodes,lagValsZS,1)
        ENDIF
      ENDIF

!        IF(dozshulimit) THEN
        ! Compute element average via quadrature at next time level
!            DO j=1,nelem
!                elemAvg = 0.5D0*SUM(gllWeights(:)*Astar(j,:))
!                IF(elemAvg .lt. 0d0) THEN
!                write(*,'(A,I1)') 'WARNING-- ELEMENT MASS IS NEGATIVE AFTER STAGE',stage
!                write(*,'(A,E10.4)') '  Minimum value = ', minval(Astar(j,:))
!                    STOP
!                ENDIF
!            ENDDO !j
!        ENDIF !dozshulimit

    ENDDO !stage
    IF(dofctlimit) THEN
      ! For fct positivity limiting, only adjust polynomial at the very END of a timestep
      !CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,0)
      CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,2)
    ENDIF !dofctlimit
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
    ELSEIF( k .eq. nNodes) THEN
        split_dadt = split_dadt - flx(j)
    ENDIF
    split_dadt = 2D0*split_dadt/qWeights(k)

END FUNCTION split_dadt

SUBROUTINE split_polyMod(AIn,qWeights,nelem,nOrder,nQuad,lagVals,stat)
! Rescales DG polynomial around element averages for positivity, based on Zhang and Shu (2010)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nOrder,nQuad
    REAL(KIND=8), DIMENSION(1:nelem,0:nOrder), INTENT(INOUT) :: AIn
    REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    REAL(KIND=8), DIMENSION(0:nOrder,0:nQuad), INTENT(IN) :: lagVals
    ! Local Variables
    INTEGER :: j,l,stat
    REAL(KIND=8), DIMENSION(0:nQuad) :: qVals
    REAL(KIND=8) :: avgVal,theta,valMin,valMax,eps,Mp,Mt

    eps = epsilon(1d0)
    ! First check incoming element averages
    IF(stat == 0) THEN
      DO j=1,nelem
          ! Use nodes themselves to do limiting, rather than evaluating polynomial multiple times
          avgVal = 0.5D0*SUM( qWeights(:)*AIn(j,:) )
          valMin = MINVAL(AIn(j,:))-eps

!          IF(avgVal .lt. 0D0) THEN
!              write(*,*) 'Element AVG in Z&S Limiter is negative!! Negative average value = ',avgVal
!          ENDIF

          ! -- Compute rescale factor
          theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )
  !       theta = MIN( abs(avgVal/(valMin-avgVal)),abs((1D0-avgVal)/(valMax-avgVal)),1D0 )

          ! -- Rescale polynomial
          AIn(j,:) = theta*(AIn(j,:)-avgVal) + avgVal
      ENDDO !j
    ELSE IF(stat == 1) THEN
      DO j=1,nelem
          ! --  Evaluate DG polynomial at quad locations
          ! Note: This isn't necessary if the quad nodes used here are same nodes that the basis is interpolating at --
          !       In this case the polynomial coefficients are the nodal values. Two exceptions are the edge
          !       values (-1 and 1) which are always part of a GLL integration. These may be read from Ain directly
          qVals(0) = AIn(j,0)
          qVals(nQuad) = AIn(j,nOrder)
          DO l=1,nQuad-1
              qVals(l) = SUM(AIn(j,:)*lagVals(:,l))
          ENDDO !l

          avgVal = 0.5D0*SUM( qWeights(:)*qVals(:) )
          valMin = MIN( MINVAL(AIn(j,:)) , MINVAL(qVals(:)) ) - eps

          ! -- Compute rescale factor
          theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )

          AIn(j,:) = theta*(AIn(j,:)-avgVal) + avgVal

      ENDDO !j
    ELSE IF(stat==2) THEN
        ! ===============================================================================================
        ! ALTERNATIVE: Replace linear rescaling with node truncation + mass redistribution
        ! ===============================================================================================
        DO j=1,nelem
            Mp = 0D0
            Mt = 0D0

            DO l=0,nOrder
                Mt = Mt + qWeights(l)*AIn(j,l)
                AIn(j,l) = MAX(0D0,AIn(j,l)) ! Zero out negative nodes
                Mp = Mp + qWeights(l)*AIn(j,l)

            ENDDO !l
            theta = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
            AIn(j,:) = theta*AIn(j,:) ! Reduce remaining positive nodes by reduction factor
        ENDDO !j
    ENDIF !stat

END SUBROUTINE split_polyMod

SUBROUTINE FLUXCOR(Acur,Apre,flx,qweights,dxel,dt,nelem,nOrder,substep)
	! Computes flux reductions factors to prevent total mass within each element from going negative
	! Returns modified flx array.
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nOrder, nelem,substep
	REAL(KIND=DOUBLE), DIMENSION(1:nelem,0:nOrder), INTENT(IN) :: Acur,Apre
  REAL(KIND=DOUBLE), DIMENSION(0:nOrder), INTENT(IN) :: qweights
	REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
	! -- Outputs
  REAL(KIND=DOUBLE), DIMENSION(0:nelem), INTENT(INOUT) :: flx
	! -- Local variables
	REAL(KIND=DOUBLE) :: Pj,Qj,eps,avgValj
  REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flxcor
	REAL(KIND=DOUBLE), DIMENSION(0:nelem+1) :: R ! Reduction ratio for outward fluxes so that element j has non-negative values (1D0 indicates no limiting needed)
	INTEGER :: j,k

	eps = 1D-6 ! Small parameter used to prevent division by 0

	DO j=1,nelem
		! Compute maximum allowable flux out of element j based on which substep of ssprk3 we are on
		SELECT CASE(substep)
			CASE(1)
        avgValj = 0.5D0*SUM(qweights(:)*Acur(j,:))
!				Qj = (dxel/dt)*avgValj

			CASE(2)
        avgValj = 0.5D0*(SUM(qweights(:)*Acur(j,:))+3D0*SUM(qweights(:)*Apre(j,:)))
!				Qj = (dxel/dt)*avgValj

			CASE(3)
        avgValj = 0.5D0*(2D0*SUM(qweights(:)*Acur(j,:))+SUM(qweights(:)*Apre(j,:)))
!				Qj = (dxel/(2D0*dt))*avgValj

		END SELECT
    Qj = (dxel/dt)*avgValj
    Qj = MAX(Qj,0D0)

    !		IF(Qj .lt. 0D0) THEN
    !			write(*,*) 'Qj < 0! j=',j,'Qj=',Qj
    !			write(*,*) Apre(0,j),Acur(0,j)
    !		END IF
		! Compute actual flux out of element j
		Pj = MAX(0D0,flx(j)) - MIN(0D0,flx(j-1)) + eps

		! Compute reduction ratio
		R(j) = MIN(1D0,Qj/Pj)
	END DO
	! Periodicity
	R(0) = R(nelem)
	R(nelem+1) = R(1)


	! Compute flux corection factors
	DO j=0,nelem
		! If flux at right edge is negative, use limiting ratio in element to the right of current one
		! (since that is where we are pulling mass from)
		flxcor(j) = R(j) - 0.5D0*(1D0-INT(SIGN(1D0,flx(j))))*(R(j)-R(j+1))
    flx(j)=flxcor(j)*flx(j)
	END DO

END SUBROUTINE FLUXCOR

REAL(KIND=8) FUNCTION split_vel_update(t)
		IMPLICIT NONE
		REAL(KIND=8), INTENT(IN) :: t
	    REAL(KIND=8) :: pi
	    REAL(KIND=8), parameter :: t_period = 5.d0

	    pi = DACOS(-1D0)
	    split_vel_update = DCOS(pi*t/t_period)
END FUNCTION split_vel_update
