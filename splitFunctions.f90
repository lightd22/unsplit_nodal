MODULE splitFunctions
  IMPLICIT NONE
  PUBLIC :: numFlux,evalExpansion,dadt,computeAverages,vel_update
  CONTAINS
    SUBROUTINE numFlux(flx,edgeVals,uin,nelem,nQuad)
    	IMPLICIT NONE
    	! -- Inputs
    	INTEGER, INTENT(IN) :: nelem,nQuad
    	DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
    	DOUBLE PRECISION, DIMENSION(0:nQuad,0:nelem), INTENT(IN) :: uin
    	! -- Outputs
    	DOUBLE PRECISION,DIMENSION(0:nelem), INTENT(OUT) :: flx
    	! -- Local variables
    	INTEGER :: j

    	DO j=0,nelem
    		flx(j) = 0.5D0*edgeVals(1,j)*(uin(nQuad,j)+DABS(uin(nQuad,j)))+0.5D0*edgeVals(0,j+1)*(uin(nQuad,j)-DABS(uin(nQuad,j)))
    	ENDDO
    END SUBROUTINE numFlux

    SUBROUTINE evalExpansion(quadVals,edgeVals,AIn,nelem,nNodes,nQuad)
    ! Evaluate ansatz solution at quadrature nodes and edges
    ! Note:  Assumes basis is Lagrange interpolating polynomials at quad nodes -> phi_j (xi_k) = qIn(k,j)
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
      DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: AIn
      ! Outputs
      DOUBLE PRECISION, DIMENSION(0:nQuad,1:nelem), INTENT(OUT) :: quadVals
      DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals
      ! Local Variables

      ! Ansatz value at quad locations for the jth element is just coefficient value
      quadVals = AIn

      edgeVals(0,1:nelem) = AIn(0,:) ! left edge value is just left-most coefficient
      edgeVals(1,1:nelem) = AIn(nNodes,:) ! right edge value is just right-most coefficient

    	! Extend edgeVals periodically
    	edgeVals(:,0) = edgeVals(:,nelem)
    	edgeVals(:,nelem+1) = edgeVals(:,1)

    END SUBROUTINE evalExpansion

    DOUBLE PRECISION FUNCTION dadt(quadVals,flx,u,qWeights,lagDeriv,k,j,nelem,nNodes)
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nelem,nNodes,k,j
      DOUBLE PRECISION, DIMENSION(0:nNodes) :: quadVals,u,qWeights,lagDeriv
      DOUBLE PRECISION, DIMENSION(0:nelem) :: flx

      dadt = SUM( u(:)*quadVals(:)*lagDeriv(:)*qWeights(:) )

      IF( k .eq. 0) THEN
          dadt = dadt + flx(j-1)
      ELSEIF( k .eq. nNodes) THEN
          dadt = dadt - flx(j)
      ENDIF
      dadt = 2D0*dadt/qWeights(k)
    END FUNCTION dadt

    SUBROUTINE computeAverages(avgs,coeffs,quadWeights,nelem,nNodes,nQuadNodes)
      IMPLICIT NONE
      ! Inputs
      INTEGER, INTENT(IN) :: nelem,nNodes,nQuadNodes
      DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs
      DOUBLE PRECISION, DIMENSION(0:nQuadNodes), INTENT(IN) :: quadWeights
      ! Outputs
      DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(OUT) :: avgs
      ! Local variables
      INTEGER :: j

      DO j=1,nelem
        avgs(j) = 0.5D0*SUM(quadWeights(:)*coeffs(:,j))
      ENDDO

    END SUBROUTINE computeAverages

    DOUBLE PRECISION FUNCTION vel_update(t)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: t
      REAL(KIND=8) :: pi
      REAL(KIND=8), parameter :: t_period = 5.d0

      pi = DACOS(-1D0)
      vel_update = DCOS(pi*t/t_period)
    END FUNCTION vel_update
END MODULE splitFunctions
