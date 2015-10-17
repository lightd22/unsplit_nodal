!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light Oct 2014
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(A,nelem,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u0,lagDeriv,time,dt,transient,&
                    dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
    USE testParameters
    USE splitFunctions
    USE splitEvalFlux
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nQuadNodes,nOrder,nZSNodes
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt,time
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes), INTENT(IN):: gllNodes,gllWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,0:nQuadNodes), INTENT(IN) :: lagDeriv
    REAL(KIND=DOUBLE), DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,0:nZSNodes), INTENT(IN) :: lagValsZS
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes,1:nelem), INTENT(IN) :: u0
    LOGICAL, INTENT(IN) :: transient,dozshulimit

    ! Outputs
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,1:nelem), INTENT(INOUT) :: A

    ! Local Variables
    REAL(KIND=DOUBLE), DIMENSION(0:nOrder,1:nelem) :: AFwd,AStar
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes,1:nelem):: u,quadVals
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes,0:nelem) :: utmp
    REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: edgeVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flx
    REAL(KIND=DOUBLE), DIMENSION(0:nQuadNodes) :: currLagDeriv,currElemQuad,currElemVel
    LOGICAL :: stopFlag = .FALSE.

    DOUBLE PRECISION, DIMENSION(1:nelem) :: elemAvg

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
    CALL computeAverages(elemAvg,AStar,gllWeights,nelem,nOrder,nQuadNodes)

    ! Do SSPRK3 Update
    DO stage=1,3
      ! Update velocities (if necessary)
      IF(transient) THEN
        SELECT CASE(stage)
          CASE(1)
        		u = u0*vel_update(time)
          CASE(2)
          	u = u0*vel_update(time+dt)
          CASE(3)
        		u = u0*vel_update(time+dt/2d0)
        END SELECT
      ENDIF ! transient
      utmp(:,1:nelem) = u
      utmp(:,0) = u(:,nelem)

      CALL evalFluxes(flx,quadVals,elemAvg,AStar,utmp,gllWeights,nelem,nOrder,nQuadNodes,dt,dxel)

!      CALL split_evalExpansion(quadVals,edgeVals,AStar,nelem,nQuadNodes,nOrder)
!      CALL split_numFlux(flx,edgeVals,utmp,nelem,nQuadNodes)
!      IF(dofctlimit) THEN ! Use FCT adjustment for element avg nonnegativity
!        CALL FLUXCOR(Astar,Astar,flx,gllWeights,dxel,dt,nelem,nOrder,1)
!        CALL FLUXCOR(Astar,Astar,flx,gllWeights,dxel,dt,nelem,nOrder,stage)
!      ENDIF !dofctlimit

      ! Forward Step
      DO j=1,nelem
        currElemQuad = quadVals(:,j)
        currElemVel = utmp(:,j)
        DO k=0,nOrder
          currLagDeriv = lagDeriv(k,:)
          AFwd(k,j) = AStar(k,j) + (dt/dxel)*dadt(currElemQuad,flx,currElemVel,gllWeights,&
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

      IF(doPosLim) THEN
        CALL computeAverages(elemAvg,AStar,gllWeights,nelem,nOrder,nQuadNodes)
        DO j=1,nelem
          IF(elemAvg(j).lt.0D0) THEN
            write(*,'(A,I2,A,I2)') 'Average of element',j,' is negative after stage',stage
            write(*,'(A,D10.4)') 'Problem average=',elemAvg(j)
            stopFlag = .TRUE.
          ENDIF
        ENDDO!j
        IF(limitingMeth == 3) THEN
          CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,2)
        ENDIF
      ENDIF
      IF(dozshulimit) THEN
        IF(nZSnodes .eq. nQuadNodes) THEN
          CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,0)
!         CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,2)
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

      IF(stopFlag) THEN
        write(*,*) '********* CRITICAL ERROR *********'
        write(*,*) 'S T O P P I N G...'
        write(*,*) '********* CRITICAL ERROR *********'
        STOP
      ENDIF

    ENDDO !stage
    IF(doFluxMod) THEN
      ! For fct positivity limiting, only adjust polynomial at the very END of a timestep
      !CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,0)
      CALL split_polyMod(Astar,gllWeights,nelem,nOrder,nQuadNodes,lagValsZS,2)
    ENDIF !dofctlimit
    A = AStar

END SUBROUTINE nDGsweep
