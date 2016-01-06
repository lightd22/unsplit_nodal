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
      RETURN
    ENDIF

    SELECT CASE(limitingMeth)
      CASE(1,4)
      ! FCT
        CALL fctMod(coeffs,Fhat,Ghat,elemAvg,u,v,dxel,dyel,dt,nex,ney,N,nQuad,qWeights)
      CASE(5)
        ! Lambda modification
        CALL lambdaLimiting(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dxel,dyel,dt,nex,ney,N,nQuad)
        ! Update quadurature and edge values using modified polynomial coefficients
        CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
        ! Update fluxes to reflect modified polynomial values
        CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

      CASE(6)
        ! Sliced Lambda modification
        CALL splitLambdaMod(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dxel,dyel,dt,nex,ney,N,nQuad)
        ! Update quadurature and edge values using modified polynomial coefficients
        CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
        ! Update fluxes to reflect modified polynomial values
        CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

      CASE DEFAULT
        RETURN
    END SELECT

  END SUBROUTINE evalFluxes

  SUBROUTINE fctMod(coeffs,Fhat,Ghat,elemAvg,u,v,dx,dy,dt,nex,ney,N,nQuad,qWeights)
    ! Applies Flux Corrected Transport (FCT) modification to fluxes to ensure that element
    ! means remain non-negative.
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

    ! Local variables
    INTEGER :: i,j,k,m
    DOUBLE PRECISION, DIMENSION(0:nex,1:ney) :: Fbar
    DOUBLE PRECISION, DIMENSION(1:nex,0:ney) :: Gbar
    DOUBLE PRECISION :: Qij,Pij,eps,mean
    DOUBLE PRECISION, DIMENSION(0:nex+1,0:ney+1) :: fluxRatio

    eps = 1D-10

    ! Get horizontal mean fluxes through east/west interfaces
    DO j=1,ney
      DO i=0,nex
        Fbar(i,j) = 0.5D0*sum(qweights(:)*Fhat(:,i,j))
      ENDDO !i
    ENDDO !j

    ! Get vertical mean fluxes through north/south interfaces
    DO j=0,ney
      DO i=1,nex
        Gbar(i,j) = 0.5D0*sum(qweights(:)*Ghat(:,i,j))
      ENDDO !i
    ENDDO !j

    DO j=1,ney
      DO i=1,nex
        mean = elemAvg(i,j)
        Qij = mean*dx*dy/dt ! Maximum permissible net flux out of S(i,j)
        ! Compute actual net flux out of element S(i,j)
!        Pij = eps + netFluxOut(i,j,Fhat,Ghat,nex,ney,nQuad)
        Pij = eps + meanNFO(i,j,Fbar,Gbar,dx,dy,nex,ney)
!        DO k=0,N
!          Pij = Pij + qWeights(k)*dy*(max(0D0,Fhat(k,i,j))-min(0D0,Fhat(k,i-1,j)))
!          Pij = Pij + qWeights(k)*dx*(max(0D0,Ghat(k,i,j))-min(0D0,Ghat(k,i,j-1)))
!        ENDDO !N
!        Pij = 0.5D0*Pij
        fluxRatio(i,j) = MIN(1D0,Qij/Pij)
      ENDDO !i
    ENDDO !j

    ! Periodically extend flux ratios to ghost cells
    fluxRatio(0,:) = fluxRatio(nex,:)
    fluxRatio(nex+1,:) = fluxRatio(1,:)
    fluxRatio(:,0) = fluxRatio(:,ney)
    fluxRatio(:,ney+1) = fluxRatio(:,1)

    ! Adjust horizontal fluxes
    DO j=1,ney
      DO i=0,nex
        IF(Fbar(i,j) .ge. 0D0) THEN
          Fhat(:,i,j) = fluxRatio(i,j)*Fhat(:,i,j)
        ELSE
          Fhat(:,i,j) = fluxRatio(i+1,j)*Fhat(:,i,j)
        ENDIF
      ENDDO !i
    ENDDO !j

    ! Adjust vertical fluxes
    DO j=0,ney
      DO i=1,nex
        IF(Gbar(i,j) .ge. 0D0) THEN
          Ghat(:,i,j) = fluxRatio(i,j)*Ghat(:,i,j)
        ELSE
          Ghat(:,i,j) = fluxRatio(i,j+1)*Ghat(:,i,j)
        ENDIF
      ENDDO !i
    ENDDO !j

  END SUBROUTINE fctMod

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
    DOUBLE PRECISION :: massCor,massChg

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

  SUBROUTINE localConserv(coeffs,massChg,i,j,qWeights,nex,ney,N,nQuad)
    ! Makes additive conservative modification to internal (1:nNodes-1) coefficients
    ! corresponding to internal nodal quadrature locations within each element
    ! NOTE:: This only changes the values of INTERNAL points
    ! Inputs -
    !  . massChg : ammount of mass to be added back to element
    !  . i,j : element indicies
    !  . qWeights : array quadrature weights
    ! Outputs -
    !  . coeffs : array of nodal coefficients after edge node modification
    ! Inputs
    INTEGER, INTENT(IN) :: i,j,nex,ney,N,nQuad
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    DOUBLE PRECISION, INTENT(IN) :: massChg
    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:N,0:N,1:nex,1:ney), INTENT(INOUT) :: coeffs
    ! Local variables
    INTEGER :: l,m
    DOUBLE PRECISION :: massCor

    massCor = massChg/((N-1)*(N-1))
    DO m=1,N-1
      DO l=1,N-1
        coeffs(l,m,i,j) = coeffs(l,m,i,j)+massCor/(qWeights(l)*qWeights(m))
      ENDDO !m
    ENDDO !l

  END SUBROUTINE localConserv

  SUBROUTINE lambdaLimiting(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dx,dy,dt,nex,ney,N,nQuad)
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
    DOUBLE PRECISION :: eps,Pij,alpha,avg,lam,dc,dM

    eps = 1D-10
    ! Corner nodes are examined twice (once for horizontal fluxes and again for vertical fluxes)
    ! The modification to these nodes is determined by the more strict condition from either flux
    ! and is still uniquely determined.
    DO j=1,ney
      DO i=1,nex
        dM = 0D0
        avg = elemAvg(i,j)
        ! Compute net flux out of element (i,j)
        Pij = eps+netFluxOut(i,j,Fhat,Ghat,nex,ney,nQuad)

        DO k=0,nQuad
          ! Horizontal flux adjustment
          ! Left edge
          IF(u(0,k,i,j) .lt. 0D0) THEN
            dc = coeffs(0,k,i,j)

            alpha = MAX(0D0,-1D0*Fhat(k,i-1,j))/Pij
            lam = lambda(u(0,k,i,j),dt,dx)
            coeffs(0,k,i,j) = MIN(coeffs(0,k,i,j),2D0*alpha*avg/(qWeights(k)*lam))

            dc = dc - coeffs(0,k,i,j)
            dM = dM+dc*qweights(k)
          ENDIF

          ! Right edge
          IF(u(nQuad,k,i,j) .gt. 0D0) THEN
            dc = coeffs(N,k,i,j)

            alpha = MAX(0D0,Fhat(k,i,j))/Pij
            lam = lambda(u(nQuad,k,i,j),dt,dx)
            coeffs(N,k,i,j) = MIN(coeffs(N,k,i,j),2D0*alpha*avg/(qWeights(k)*lam))

            dc = dc - coeffs(N,k,i,j)
            dM = dM+dc*qweights(k)
          ENDIF

          ! Vertical flux adjustment
          ! Bottom edge
          IF(v(k,0,i,j) .lt. 0D0) THEN
            dc = coeffs(k,0,i,j)

            alpha = MAX(0D0,-1D0*Ghat(k,i,j-1))/Pij
            lam = lambda(v(k,0,i,j),dt,dy)
            coeffs(k,0,i,j) = MIN(coeffs(k,0,i,j),2D0*alpha*avg/(qWeights(k)*lam))

            dc = dc - coeffs(k,0,i,j)
            dM = dM+dc*qweights(k)
          ENDIF

          ! Top edge
          IF(v(k,nQuad,i,j) .gt. 0D0) THEN
            dc = coeffs(k,N,i,j)

            alpha = MAX(0D0,Ghat(k,i,j))/Pij
            lam = lambda(v(k,nQuad,i,j),dt,dy)
            coeffs(k,N,i,j) = MIN(coeffs(k,N,i,j),2D0*alpha*avg/(qWeights(k)*lam))

            dc = dc - coeffs(k,N,i,j)
            dM = dM+dc*qweights(k)
          ENDIF
        ENDDO !k
        dM = dM*qWeights(0)
        CALL localConserv(coeffs,dM,i,j,qWeights,nex,ney,N,nQuad)
      ENDDO !i
    ENDDO !j

    ! Update fluxes to reflect modified polynomial values
!    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
!    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

  END SUBROUTINE lambdaLimiting

  SUBROUTINE splitLambdaMod(coeffs,Fhat,Ghat,elemAvg,u,v,qWeights,dx,dy,dt,nex,ney,N,nQuad)
    !
    ! Second Lambda polynomial modification for element mean non-negativity
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
    DOUBLE PRECISION :: eps,Pij,fluxRatio,avg,lam,dc,dM

    eps = 1D-10

    DO j=1,ney
      DO i=1,nex
        avg = elemAvg(i,j)
        dM = 0D0
        DO k=0,nQuad
          ! Horizontal slice adjustment
!          avg = 0.5D0*SUM(qWeights(:)*coeffs(:,k,i,j))
          Pij = eps+MAX(0D0,Fhat(k,i,j))-MIN(0D0,Fhat(k,i-1,j))
          IF(u(0,k,i,j) .lt. 0D0) THEN
            dc = coeffs(0,k,i,j)

            fluxRatio = MAX(0D0,-1D0*Fhat(k,i-1,j))/Pij
            call coeffAdjust(coeffs(0,k,i,j),u(0,k,i,j),avg,fluxRatio,dt,dx)

            dc = dc - coeffs(0,k,i,j)
            dM = dM+dc*qWeights(k)
          ENDIF

          IF(u(nQuad,k,i,j) .gt. 0D0) THEN
            dc = coeffs(nQuad,k,i,j)

            fluxRatio = MAX(0D0,Fhat(k,i,j))/Pij
            call coeffAdjust(coeffs(nQuad,k,i,j),u(nQuad,k,i,j),avg,fluxRatio,dt,dx)

            dc = dc - coeffs(nQuad,k,i,j)
            dM = dM+dc*qWeights(k)
          ENDIF

          ! Vertical slice adjustment
!          avg = 0.5D0*SUM(qWeights(:)*coeffs(k,:,i,j))
          Pij = eps+MAX(0D0,Ghat(k,i,j))-MIN(0D0,Ghat(k,i,j-1))
          IF(v(k,0,i,j) .lt. 0D0) THEN
            dc = coeffs(k,0,i,j)

            fluxRatio = MAX(0D0,-1D0*Ghat(k,i,j-1))/Pij
            call coeffAdjust(coeffs(k,0,i,j),v(k,0,i,j),avg,fluxRatio,dt,dy)

            dc = dc - coeffs(k,0,i,j)
            dM = dM+dc*qWeights(k)
          ENDIF

          IF(v(k,nQuad,i,j) .gt. 0D0) THEN
            dc = coeffs(k,nQuad,i,j)

            fluxRatio = MAX(0D0,Ghat(k,i,j))/Pij
            call coeffAdjust(coeffs(k,nQuad,i,j),v(k,nQuad,i,j),avg,fluxRatio,dt,dy)

            dc = dc - coeffs(k,nQuad,i,j)
            dM = dM+dc*qWeights(k)
          ENDIF
        ENDDO !k
        dM = dM*qWeights(0)
        CALL localConserv(coeffs,dM,i,j,qWeights,nex,ney,N,nQuad)
      ENDDO !i
    ENDDO !j

    ! Update fluxes to reflect modified polynomial values
!    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuad,N,nex,ney)
!    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuad,N,nex,ney)

  END SUBROUTINE splitLambdaMod

  SUBROUTINE coeffAdjust(coeff,u,avg,fluxRatio,dt,dx)
    ! Inputs
    DOUBLE PRECISION, INTENT(IN) :: u,avg,fluxRatio,dt,dx
    ! Outputs
    DOUBLE PRECISION, INTENT(INOUT) :: coeff
    ! Local Variables
    DOUBLE PRECISION :: lam

    lam = lambda(u,dt,dx)
    coeff = MIN(coeff,fluxRatio*avg/lam)
  END SUBROUTINE coeffAdjust

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

  DOUBLE PRECISION FUNCTION meanNFO(i,j,Fbar,Gbar,dx,dy,nex,ney)
    !
    ! Computes the net AVERAGE flux out of cell (i,j)
    !
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: i,j,nex,ney
    DOUBLE PRECISION, INTENT(IN) :: dx,dy
    DOUBLE PRECISION, DIMENSION(0:nex,1:ney), INTENT(IN) :: Fbar
    DOUBLE PRECISION, DIMENSION(1:nex,0:ney), INTENT(IN) :: Gbar
    ! Outputs
    ! Local Variables
    DOUBLE PRECISION :: meanEW,meanNS
    meanEW = dy*(max(0D0,Fbar(i,j))-min(0D0,Fbar(i-1,j))) ! net average flux through east/west interfaces
    meanNS = dx*(max(0D0,Gbar(i,j))-min(0D0,Gbar(i,j-1))) ! net average flux through north/south interfaces
    meanNFO = meanEW + meanNS
  END FUNCTION meanNFO

END MODULE unspEvalFlux
