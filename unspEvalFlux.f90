MODULE unspEvalFlux
  IMPLICIT NONE
  PUBLIC :: evalFluxes
CONTAINS
  SUBROUTINE evalFluxes(Fhat,Ghat,coeffs,quadVals,elemAvg,u,v,qWeights,nex,ney,N,nQuad,dt,dxel,dyel)
    ! ---
    ! TODO: Documentation
    ! ---
    USE testParameters
    USE unspFunctions
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,N,nQuad
    DOUBLE PRECISION, INTENT(IN) :: dt,dxel,dyel
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvg
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nQuad,1:nex,1:ney), INTENT(IN) :: u,v

    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nex,1:ney), INTENT(OUT) :: Fhat
    DOUBLE PRECISION, DIMENSION(0:nQuad,1:nex,0:ney), INTENT(OUT) :: Ghat
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nQuad,1:nex,1:ney), INTENT(OUT) :: quadVals
    DOUBLE PRECISION, DIMENSION(0:N,0:N,1:nex,1:ney), INTENT(INOUT) :: coeffs

    ! Local Variables
    INTEGER :: i,j,k
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,0:nex+1,1:ney) :: edgeValsEW
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,1:nex,0:ney+1) :: edgeValsNS


    ! Get edge and quadrature values then compute unmodified fluxes
    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

    IF(.NOT. doFluxMod) THEN
      !RETURN
    ENDIF

    SELECT CASE(limitingMeth)
      CASE(4)
        ! Modified FCT
        WRITE(*,*) 'WARNING -- This is not fully functional yet!'
        WRITE(*,*) 'Corner nodes are currently being modified twice (for F and G fluxes)'
        WRITE(*,*) 'CURRENT FORMULATION DOES *NOT* GIVE UNIQUE SOLUTION'
        CALL fctPolyMod(coeffs,Fhat,Ghat,elemAvg,u,v,dxel,dyel,dt,nex,ney,N,nQuad)
        STOP
      CASE(5)
        ! Lambda modification
        CALL lambdaPolyMod(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dxel,dyel,dt,nex,ney,N,nQuad)
      CASE DEFAULT
        RETURN
    END SELECT
    ! Conservatively adjust interior nodes
    CALL conservPolyMod(coeffs,edgeValsNS,edgeValsEW,qWeights,nex,ney,N,nQuad)

    ! Update quadurature and edge values using modified polynomial coefficients
    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,n,nex,ney)
  END SUBROUTINE evalFluxes

  SUBROUTINE fctPolyMod(coeffs,Fhat,Ghat,elemAvg,u,v,dx,dy,dt,nex,ney,N,nQuad)
    ! Makes FCT-based correction to fluxes to ensure that element means remain non-negative
    ! Then modifies local polynomials so that FCT fluxes will be consistent
    ! NOTE:: THIS IS NOT CONSERVATIVE! Requires calling conservPolyMod() after to keep conservation
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,N,nQuad
    DOUBLE PRECISION, INTENT(IN) :: dx,dy,dt
    DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvg
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nQuad,1:nex,1:ney), INTENT(IN) :: u,v

    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:N,0:N,1:nex,1:ney), INTENT(INOUT) :: coeffs
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nex,1:ney), INTENT(INOUT) :: Fhat
    DOUBLE PRECISION, DIMENSION(0:nQuad,1:nex,0:ney), INTENT(INOUT) :: Ghat

    ! Local variables
    INTEGER :: i,j,k,m
    DOUBLE PRECISION :: Qij,Pij,eps,mean
    DOUBLE PRECISION, DIMENSION(0:nex+1,0:ney+1) :: fluxRatio

    eps = 1D-10

    DO i=1,nex
      DO j=1,ney
        mean = elemAvg(i,j)
        Qij = mean*dx*dy/dt
        ! Compute net flux out of element (i,j)
        Pij = eps + netFluxOut(i,j,Fhat,Ghat,nex,ney,nQuad)
        fluxRatio(i,j) = MIN(1D0,Qij/Pij)
      ENDDO !j
    ENDDO !i

    fluxRatio(0,:) = fluxRatio(nex,:)
    fluxRatio(nex+1,:) = fluxRatio(1,:)
    fluxRatio(:,0) = fluxRatio(:,ney)
    fluxRatio(:,ney+1) = fluxRatio(:,1)

    ! TODO: This section should modify coefficients so that numerical flux is consistent with
    ! the newly modified polynomials. However, it is currently adjusting the four corner nodes
    ! (which have BOTH a F and G flux) TWICE.
    ! Not only is this overkill in terms of modification, it means that both fluxes aren't
    ! truly consistent with the modified polynomial.
    ! Instead, this may need to adjust both fluxes with the **same** modification factor
    ! and then fix the polynomial so that it is consistent

    ! Adjust horizontal fluxes
    DO i=0,nex
      DO j=1,ney
        DO k=0,N
          IF(Fhat(k,i,j) .ge. 0D0) THEN
            Fhat(k,i,j) = fluxRatio(i,j)*Fhat(k,i,j)
          ELSE
            Fhat(k,i,j) = fluxRatio(i+1,j)*Fhat(k,i,j)
          ENDIF
        ENDDO !k
      ENDDO !j
    ENDDO !i

    ! Multiplicatively adjust upstream DG approximating polynomials so that new fluxes are consistent
    ! Look at inflow boundary on left side (i=0)
    DO j=1,ney
      DO k=0,N
        IF(u(0,k,1,j) .gt. 0D0) THEN
          coeffs(N,k,nex,j) = Fhat(k,i,j)/u(0,k,1,j)
        ELSE IF(u(0,k,1,j) .lt. 0D0) THEN
          coeffs(0,k,1,j) = Fhat(k,i,j)/u(0,k,1,j)
        ELSE
          ! Do nothing
        ENDIF
      ENDDO !k
    ENDDO !j

    ! Look at right edge of remaining elements
    DO i=1,nex
      DO j=1,ney
        DO k=0,N
          IF(u(N,k,i,j) .gt. 0D0) THEN
            coeffs(N,k,i,j) = Fhat(k,i,j)/u(N,k,i,j)
          ELSE IF(u(N,k,i,j) .lt. 0D0) THEN
            coeffs(0,k,i+1,j) = Fhat(k,i,j)/u(N,k,i,j)
          ELSE
            ! Do nothing
          ENDIF
        ENDDO !k
      ENDDO !j
    ENDDO !i

    ! Adjust vertical fluxes
    DO i=1,nex
      DO j=0,ney
        DO k=0,N
          IF(Ghat(k,i,j) .ge. 0D0) THEN
            Ghat(k,i,j) = fluxRatio(i,j)*Ghat(k,i,j)
          ELSE
            Ghat(k,i,j) = fluxRatio(i,j+1)*Ghat(k,i,j)
          ENDIF
        ENDDO !k
      ENDDO !j
    ENDDO !i

    ! Multiplicatively adjust upstream DG approximating polynomials so that new fluxes are consistent
    ! Look at inflow boundary on bottom side (j=0)
    DO i=1,nex
      DO k=0,N
        IF(v(k,0,i,1) .gt. 0D0) THEN
          coeffs(k,N,i,ney) = Ghat(k,i,j)/v(k,0,i,1)
        ELSE IF(v(k,0,i,1) .lt. 0D0) THEN
          coeffs(k,0,i,1) = Ghat(k,i,j)/v(k,0,i,1)
        ELSE
          ! Do nothing if v=0
        ENDIF
      ENDDO !k
    ENDDO !j

    ! Look at top edge of remaining elements
    DO i=1,nex
      DO j=1,ney
        DO k=0,N
          IF(v(k,N,i,j) .gt. 0D0) THEN
            coeffs(k,N,i,j) = Ghat(k,i,j)/v(k,N,i,j)
          ELSE IF(v(k,N,i,j) .lt. 0D0) THEN
            coeffs(k,0,i,j+1) = Ghat(k,i,j)/v(k,N,i,j)
          ELSE
            ! Do nothing
          ENDIF
        ENDDO !k
      ENDDO !j
    ENDDO !i

  END SUBROUTINE fctPolyMod

  SUBROUTINE conservPolyMod(coeffs,edgeValsNS,edgeValsEW,qWeights,nex,ney,N,nQuad)
    ! Makes additive conservative modification to internal (1:nNodes-1) coefficients
    ! corresponding to internal nodal quadrature locations within each element
    ! NOTE:: This only changes the values of INTERNAL points
    ! Inputs -
    !  . edgeValsNS : array of top and bottom edge values of polynomial BEFORE modification
    !  . edgeValsEW : array of left and right edge values of polynomial BEFORE modification
    !  . qWeights : array quadrature weights
    ! Outputs -
    !  . coeffs : array of nodal coefficients after desired edge node modification
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,N,nQuad
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,0:nex+1,1:ney), INTENT(IN) :: edgeValsEW
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,1:nex,0:ney+1), INTENT(IN) :: edgeValsNS

    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:N,0:N,1:nex,1:ney), INTENT(INOUT) :: coeffs
    ! Local variables
    INTEGER :: i,j,l,m
    DOUBLE PRECISION :: massChg,massCor

    ! Correct mass for modified polynomials to maintain conservation
    DO i=1,nex
      DO j=1,ney
        massChg = 0D0
        ! Sum changes over EW faces
        DO l=0,nQuad
          massChg = massChg+(edgeValsEW(0,l,i,j)-coeffs(0,l,i,j))*qWeights(0)*qWeights(l)
          massChg = massChg+(edgeValsEW(1,l,i,j)-coeffs(N,l,i,j))*qWeights(nQuad)*qWeights(l)
        ENDDO !l

        ! To avoid duplication, NS only summed over non-corner edge nodes
        DO l=1,nQuad-1
          massChg = massChg+(edgeValsNS(0,l,i,j)-coeffs(l,0,i,j))*qWeights(0)*qWeights(l)
          massChg = massChg+(edgeValsNS(1,l,i,j)-coeffs(l,N,i,j))*qWeights(nQuad)*qWeights(l)
        ENDDO !l

        ! Compute mass that should be added to each internal node (there are (N-1)^2 of these)
        massCor = massChg/((N-1)*(N-1))
        DO l=1,N-1
          DO m=1,N-1
            coeffs(l,m,i,j) = coeffs(l,m,i,j)+massCor/(qWeights(l)*qWeights(m))
          ENDDO !m
        ENDDO !l
      ENDDO !j
    ENDDO !i
  END SUBROUTINE conservPolyMod

  SUBROUTINE lambdaPolyMod(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dx,dy,dt,nex,ney,N,nQuad)
    !
    ! Lambda polynomial modification for elemenet mean non-negativity
    ! Note: Assumes that incoming polynomial values in coeffs array are all non-negative
    !
    USE testParameters
    USE unspFunctions
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,N,nQuad
    DOUBLE PRECISION, INTENT(IN) :: dx,dy,dt
    DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvg
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nQuad,1:nex,1:ney), INTENT(IN) :: u,v
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:N,0:N,1:nex,1:ney), INTENT(INOUT) :: coeffs
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nex,1:ney), INTENT(INOUT) :: Fhat
    DOUBLE PRECISION, DIMENSION(0:nQuad,1:nex,0:ney), INTENT(INOUT) :: Ghat
    ! Local Variables
    INTEGER :: i,j,k
    DOUBLE PRECISION :: eps,Pij,alpha,avg,lam
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,0:nex+1,1:ney) :: edgeValsEW
    DOUBLE PRECISION, DIMENSION(0:1,0:nQuad,1:nex,0:ney+1) :: edgeValsNS
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nQuad,1:nex,1:ney) :: quadVals

    eps = 1D-10

    ! Corner nodes are examined twice (once for horizontal fluxes and again for vertical fluxes)
    ! The modification to these nodes is determined by the more strict condition from either flux
    ! but are still uniquely determined.
    DO i=1,nex
      DO j=1,ney
        avg = elemAvg(i,j)
        ! Compute net flux out of element (i,j)
        Pij = eps+netFluxOut(i,j,Fhat,Ghat,nex,ney,nQuad)

        ! Make horizontal adjustment
        DO k=0,nQuad
          ! Right edge
          alpha = MAX(0D0,Fhat(k,i,j))/Pij
          lam = lambda(u(nQuad,k,i,j),dt,dx)
          IF(u(nQuad,k,i,j) .gt. 0D0) THEN
            coeffs(N,k,i,j) = MIN(coeffs(N,k,i,j),2D0*alpha*avg/(qWeights(k)*lam))
          ENDIF

          ! Left edge
          alpha = MAX(0D0,-1D0*Fhat(k,i-1,j))/Pij
          lam = lambda(u(0,k,i,j),dt,dx)
          IF(u(0,k,i,j) .lt. 0D0) THEN
            coeffs(0,k,i,j) = MIN(coeffs(0,k,i,j),2D0*alpha*avg/(qWeights(k)*lam))
          ENDIF
        ENDDO !k

        ! Make vertical adjustment
        DO k=0,nQuad
          ! Top edge
          alpha = MAX(0D0,Ghat(k,i,j))/Pij
          lam = lambda(v(k,nQuad,i,j),dt,dy)
          IF(v(k,nQuad,i,j) .gt. 0D0) THEN
            coeffs(k,N,i,j) = MIN(coeffs(k,N,i,j),2D0*alpha*avg/(qWeights(k)*lam))
          ENDIF

          ! Bottom edge
          alpha = MAX(0D0,-1D0*Ghat(k,i,j-1))/Pij
          lam = lambda(v(k,0,i,j),dt,dy)
          IF(v(k,0,i,j) .lt. 0D0) THEN
            coeffs(k,0,i,j) = MIN(coeffs(k,0,i,j),2D0*alpha*avg/(qWeights(k)*lam))
          ENDIF
        ENDDO !k
      ENDDO !j
    ENDDO !i

    ! Update fluxes to reflect modified polynomial values
    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

  END SUBROUTINE lambdaPolyMod

  DOUBLE PRECISION FUNCTION lambda(u,dt,dx)
     IMPLICIT NONE
     ! Inputs
     DOUBLE PRECISION, INTENT(IN) :: u,dt,dx
     ! Local variables
     DOUBLE PRECISION :: eps

     eps = 1D-8
     lambda = MIN(MAX(eps,ABS(u)*dt/dx),1D0)
  END FUNCTION lambda

  DOUBLE PRECISION FUNCTION netFluxOut(i,j,Fhat,Ghat,nex,ney,nQuad)
    !
    ! Computes the net flux OUT of cell (i,j)
    !
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: i,j,nex,ney,nQuad
    DOUBLE PRECISION, DIMENSION(0:nQuad,0:nex,1:ney), INTENT(IN) :: Fhat
    DOUBLE PRECISION, DIMENSION(0:nQuad,1:nex,0:ney), INTENT(IN) :: Ghat
    ! Outputs
    ! Local variables
    INTEGER :: k

    netFluxOut = 0D0
    DO k=0,nQuad
      netFluxOut = netFluxOut + MAX(Fhat(k,i,j),0D0) - MIN(Fhat(k,i-1,j),0D0)
      netFluxOut = netFluxOut + MAX(Ghat(k,i,j),0D0) - MIN(Ghat(k,i,j-1),0D0)
    ENDDO !k
  END FUNCTION netFluxOut

END MODULE unspEvalFlux
