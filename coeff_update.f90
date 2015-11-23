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
  USE unspFunctions
  USE unspEvalFlux
  USE unspPositivityLimit
  USE testParameters
  IMPLICIT NONE

  ! Inputs
  INTEGER, INTENT(IN) :: nex,ney,norder,nQuadNodes,gqOrder,nZSnodes
  DOUBLE PRECISION, INTENT(IN) :: dxel,dyel,dt,time
  DOUBLE PRECISION, DIMENSION(0:nQuadNodes), INTENT(IN) :: gllNodes,gllWeights
  DOUBLE PRECISION, DIMENSION(0:gqOrder), INTENT(IN) :: gqWeights
  DOUBLE PRECISION, DIMENSION(0:norder,0:nQuadNodes), INTENT(IN) :: lagrangeDeriv
  DOUBLE PRECISION, DIMENSION(0:norder,0:gqOrder), INTENT(IN) :: lagGaussVal
  DOUBLE PRECISION, DIMENSION(0:norder,0:nZSNodes), INTENT(IN) :: lagValsZS
	DOUBLE PRECISION, DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), INTENT(IN) :: u0,v0
  LOGICAL, INTENT(IN) :: dozshulimit,transient,doZSMaxCFL

  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(INOUT) :: A

  ! Local Variables
  INTEGER :: i,j,l,m,stage
  DOUBLE PRECISION :: coef1
  DOUBLE PRECISION, DIMENSION(0:norder,0:norder,1:nex,1:ney) :: A1,A2
  DOUBLE PRECISION, DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney) :: quadVals,u,v
  DOUBLE PRECISION, DIMENSION(0:nQuadNodes,0:nQuadNodes) :: currElemQuad,currElemU,currElemV
  DOUBLE PRECISION, DIMENSION(0:nQuadNodes,0:nex,1:ney) :: Fhat
  DOUBLE PRECISION, DIMENSION(0:nQuadNodes,1:nex,0:ney) :: Ghat
  DOUBLE PRECISION, DIMENSION(1:nex,1:ney) :: elemAvg
  LOGICAL :: stopFlag = .FALSE.

	! ######################
	! Time step using Strong-Stability Preserving RK3
	! ######################
	coef1 = dt/dxel/dyel
  u = u0
  v = v0

  A1 = A

  IF(doPosLim) THEN
    CALL computeAverages(elemAvg,gllWeights,A1,nex,ney,nOrder,nQuadNodes)
    IF(MINVAL(elemAvg) .lt. 0D0) THEN
      write(*,*) '*** WARNING::: Incoming element means are negative before first stage!'
      stopFlag = .TRUE.
    ENDIF
    if(minval(A1) .lt. 0D0) then
      write(*,*) '*** WARNING::: Incoming nodal values are negative before first stage!'
      stopFlag = .TRUE.
    endif
  ENDIF

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

    ! Update fluxes
     CALL evalFluxes(Fhat,Ghat,A1,quadVals,elemAvg,u,v,gllWeights,nex,ney,nOrder,nQuadNodes,dt,dxel,dyel)

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

    IF(doPosLim .and. stage < 3) THEN
      ! Check element averages
      CALL computeAverages(elemAvg,gllWeights,A1,nex,ney,nOrder,nQuadNodes)
      IF(MINVAL(elemAvg) .lt. 0d0) THEN
        write(*,*) 'WARNING-- ELEMENT MEAN IS NEGATIVE AFTER FWD STEP',stage
        write(*,*) minval(elemAvg)
        stopFlag = .TRUE.
      ENDIF

      IF(doZShuLimit) THEN
        CALL limitMeanPositivity(A1,elemAvg,lagGaussVal,gqWeights,gllWeights,&
                                nex,ney,nOrder,gqOrder)
      ENDIF
      IF(eachStageNodeLim) THEN
        CALL limitNodePositivity(A1,elemAvg,gllWeights,nex,ney,nOrder)
      ENDIF

    ENDIF

    IF(stopFlag) THEN
      write(*,*) '********* CRITICAL ERROR *********'
      write(*,*) 'S T O P P I N G...'
      write(*,*) '********* CRITICAL ERROR *********'
      STOP
    ENDIF

    ENDDO !stage
    A = A1

    IF(doPosLim) THEN
      ! Compute element average via quadrature at next time level
      CALL computeAverages(elemAvg,gllWeights,A,nex,ney,nOrder,nQuadNodes)
      IF(minval(elemAvg) .lt. 0d0) THEN
        write(*,*) 'WARNING-- ELEMENT MEAN IS NEGATIVE BEFORE LIMITING NODES'
        write(*,*) minval(elemAvg)
      ENDIF
      CALL limitNodePositivity(A,elemAvg,gllWeights,nex,ney,nOrder)
    ENDIF !doposlim
END SUBROUTINE coeff_update
