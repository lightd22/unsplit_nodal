!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light Oct 2014
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(A,nelem,dxel,nOrder,nQuadNodes,gllNodes,gllWeights,u0,lagDeriv,time,dt,transient,&
                    dozshulimit,nZSnodes,quadZSWeights,lagValsZS)
    USE testParameters
    USE splitFunctions
    USE splitEvalFlux
    USE splitPositivityLimit
    IMPLICIT NONE

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nQuadNodes,nOrder,nZSNodes
    DOUBLE PRECISION, INTENT(IN) :: dxel,dt,time
    DOUBLE PRECISION, DIMENSION(0:nQuadNodes), INTENT(IN):: gllNodes,gllWeights
    DOUBLE PRECISION, DIMENSION(0:nOrder,0:nQuadNodes), INTENT(IN) :: lagDeriv
    DOUBLE PRECISION, DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSWeights
    DOUBLE PRECISION, DIMENSION(0:nOrder,0:nZSNodes), INTENT(IN) :: lagValsZS
    DOUBLE PRECISION, DIMENSION(0:nQuadNodes,1:nelem), INTENT(IN) :: u0
    LOGICAL, INTENT(IN) :: transient,dozshulimit

    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:nOrder,1:nelem), INTENT(INOUT) :: A

    ! Local Variables
    DOUBLE PRECISION, DIMENSION(0:nOrder,1:nelem) :: AFwd,AStar
    DOUBLE PRECISION, DIMENSION(0:nQuadNodes,1:nelem):: u,quadVals
    DOUBLE PRECISION, DIMENSION(0:nQuadNodes,0:nelem) :: utmp
    DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1) :: edgeVals
    DOUBLE PRECISION, DIMENSION(0:nelem) :: flx
    DOUBLE PRECISION, DIMENSION(0:nQuadNodes) :: currLagDeriv,currElemQuad,currElemVel
    LOGICAL :: stopFlag = .FALSE.

    DOUBLE PRECISION, DIMENSION(1:nelem) :: elemAvg

    INTEGER :: k,j,stage,stat,l

    u = u0

    AStar = A
    IF(doPosLim) THEN
      CALL computeAverages(elemAvg,AStar,gllWeights,nelem,nOrder,nQuadNodes)
      IF(MINVAL(elemAvg) .lt. 0D0) THEN
        write(*,*) '*** WARNING::: Incoming element means are negative before first stage!'
        stopFlag = .TRUE.
      ENDIF
    ENDIF

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

      IF(doPosLim .and. limitingMeth == 1) THEN
        CALL limitMeanPositivity(AStar,elemAvg,quadZSWeights,nelem,nOrder,nZSNodes)
      ENDIF

      CALL evalFluxes(flx,quadVals,elemAvg,AStar,utmp,gllWeights,nelem,nOrder,nQuadNodes,dt,dxel)

      ! Forward Step
      DO j=1,nelem
        currElemQuad = quadVals(:,j)
        currElemVel = utmp(:,j)
        DO k=0,nOrder
          currLagDeriv = lagDeriv(k,:)
          AFwd(k,j) = AStar(k,j) + (dt/dxel)*dadt(currElemQuad,flx,currElemVel,gllWeights,&
                                                        currLagDeriv,k,j,nelem,nQuadNodes)
        ENDDO !k
      ENDDO !j

  		SELECT CASE(stage)
  		CASE(1)
  			AStar = AFwd
  		CASE(2)
  			AStar = 0.75D0*A + 0.25D0*AFwd
  		CASE(3)
  			AStar = A/3d0 + 2D0*AFwd/3D0
  		END SELECT

      IF(doPosLim .and. stage .lt. 3) THEN
        CALL computeAverages(elemAvg,AStar,gllWeights,nelem,nOrder,nQuadNodes)
        DO j=1,nelem
          IF(elemAvg(j).lt.0D0) THEN
            write(*,'(A,I2,A,I2)') 'Average of element',j,' is negative after stage',stage
            write(*,'(A,D10.4)') 'Problem average=',elemAvg(j)
            stopFlag = .TRUE.
          ENDIF
        ENDDO!j
!        IF(limitingMeth == 1) THEN
!          CALL limitMeanPositivity(AStar,elemAvg,quadZSWeights,nelem,nOrder,nZSNodes)
!        ELSE IF(limitingMeth == 3) THEN
!          CALL limitNodePositivity(AStar,elemAvg,gllWeights,nelem,nOrder,nQuadNodes)
!        END IF
        IF(limitingMeth == 3) THEN
          CALL limitNodePositivity(AStar,elemAvg,gllWeights,nelem,nOrder,nQuadNodes)
        ENDIF

      ENDIF

      IF(stopFlag) THEN
        write(*,*) '********* CRITICAL ERROR *********'
        write(*,*) 'S T O P P I N G...'
        write(*,*) '********* CRITICAL ERROR *********'
        STOP
      ENDIF

    ENDDO !stage
    IF(doPosLim) THEN
      ! Adjust polynomial at every node for positive definiteness at end of time step
      CALL computeAverages(elemAvg,AStar,gllWeights,nelem,nOrder,nQuadNodes)
      CALL limitNodePositivity(AStar,elemAvg,gllWeights,nelem,nOrder,nQuadNodes)
    ENDIF !doPosLim
    A = AStar

END SUBROUTINE nDGsweep
