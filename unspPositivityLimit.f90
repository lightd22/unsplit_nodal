MODULE unspPositivityLimit
  IMPLICIT NONE
  PUBLIC :: limitMeanPositivity,limitNodePositivity
CONTAINS
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

      ! Local variables
      INTEGER :: i,j,l,beta,p,q
      DOUBLE PRECISION :: avg,valMin,leftTrace,rightTrace,topTrace,botTrace,magicPt
      DOUBLE PRECISION :: eps,theta,traceMin,massChg,weight,dCoeff
      DOUBLE PRECISION, DIMENSION(0:gqOrder) :: magicTmp

      eps = epsilon(1D0)
      SELECT CASE(limitingMeth)
      CASE(2)
        DO j = 1,ney
          DO i = 1,nex
            ! Step 1: Pull element average
            avg = elemAvgs(i,j)

            ! Step 2: Truncate along edges and evaluate magic point
            massChg = 0D0
            magicPt = 0D0

            ! Top and bottom edges
            DO q=0,nOrder,nOrder
              DO p=0,nOrder
                weight = 0.25*gllWeights(p)*gllWeights(0)
                dCoeff = 0.5*(ABS(coeffs(p,q,i,j))-coeffs(p,q,i,j))
                massChg = massChg + weight*dCoeff
                coeffs(p,q,i,j) = coeffs(p,q,i,j)+dCoeff

  !                massChg = massChg - weight*MIN(coeffs(p,q,i,j),0D0)
  !                coeffs(p,q,i,j) = MAX(coeffs(p,q,i,j),0D0)

                ! Track contribution to magic point
                magicPt = magicPt + weight*coeffs(p,q,i,j)
              ENDDO !p
            ENDDO !q

            ! Left and right edges
            DO q = 1,nOrder-1
              DO p = 0,nOrder,nOrder
                weight = 0.25*gllWeights(p)*gllWeights(q)
                dCoeff = 0.5*(ABS(coeffs(p,q,i,j))-coeffs(p,q,i,j))
                massChg = massChg + weight*dCoeff
                coeffs(p,q,i,j) = coeffs(p,q,i,j)+dCoeff

  !                massChg = massChg - weight*MIN(coeffs(p,q,i,j),0D0)
  !                coeffs(p,q,i,j) = MAX(coeffs(p,q,i,j),0D0)

                ! Track contribution to magic point
                magicPt = magicPt + weight*coeffs(p,q,i,j)
              ENDDO !p
            ENDDO !q

            ! Check magic point for positivity
            magicPt = avg - magicPt
            IF(magicPt .lt. 0D0) THEN
              ! Update mass change due to truncating magic pt
              massChg = massChg - magicPt

              ! Truncate interior nodes
              coeffs(1:nOrder-1,1:nOrder-1,i,j) = 0D0
            ENDIF

            ! Rescale remaining coefficents to maintain conservation
            theta = avg/MAX(avg+massChg,TINY(1D0))
            coeffs(:,:,i,j) = theta*coeffs(:,:,i,j)
          ENDDO !j
        ENDDO !i
      CASE(1,3)
        DO j=1,ney
          DO i=1,nex
            ! Step 1: Initialize mean and minimum value
            avg = elemAvgs(i,j)
            valMin = HUGE(1D0)

            ! Step 2: Compute minimum value between traces and "magic point"
            DO beta=0,gqOrder
              ! Look at left, right, top, and bottom trace values for this
              ! Gauss quadrature point (indexed by beta)
              leftTrace = SUM(coeffs(0,:,i,j)*lagGaussVal(:,beta))
              rightTrace = SUM(coeffs(nOrder,:,i,j)*lagGaussVal(:,beta))
              botTrace = SUM(coeffs(:,0,i,j)*lagGaussVal(:,beta))
              topTrace = SUM(coeffs(:,nOrder,i,j)*lagGaussVal(:,beta))

              ! Update minimum
              valMin = MIN(valMin,leftTrace,rightTrace,botTrace,topTrace)

              magicTmp(beta) = gqWeights(beta)*(mu1*(rightTrace+leftTrace)+mu2*(topTrace+botTrace))
            ENDDO !beta
            traceMin = valMin

            ! Compute "magic point" value
            magicPt = avg-0.25D0*zsMinWeight*SUM(magicTmp)
            magicPt = magicPt/(1D0-zsMinWeight)

            ! Update minimum
            !valMin = MIN(valMin,magicPt,MINVAL(coeffs(:,:,i,j)))
            valMin = MIN(valMin,magicPt)

            !valMin = MINVAL(coeffs(:,:,i,j))
            valMin = valMin - eps

            ! Step 3: Compute theta and rescale polynomial
            theta = MIN(abs(avg)/(abs(valMin-avg)),1D0)

            IF(magicPt .lt. traceMin .AND. theta .lt. 1D0) THEN
              magCount = magCount+1
            ENDIF

            ! Rescale reconstructing polynomial for (i,j)th element
            coeffs(:,:,i,j) = theta*(coeffs(:,:,i,j) - avg) + avg
          ENDDO !i
        ENDDO !j
      CASE DEFAULT
        ! Do nothing
      END SELECT

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
    ! Local variables
    INTEGER :: i,j,p,q
    DOUBLE PRECISION :: valMin,eps,theta,avg
    DOUBLE PRECISION :: dCoeff,massChg,weight

    eps = epsilon(1d0)

    SELECT CASE(limitingMeth)
    CASE(2:5)
      ! ================================================================================
      ! Use TMAR limiting at GLL nodes
      ! ================================================================================
      DO j=1,ney
        DO i=1,nex
          avg = elemAvgs(i,j)
          massChg = 0D0
          DO q=0,nOrder
            DO p=0,nOrder
              weight = 0.25*gllWeights(p)*gllWeights(q)
              dCoeff = 0.5*(ABS(coeffs(p,q,i,j))-coeffs(p,q,i,j))
              massChg = massChg + weight*dCoeff
              coeffs(p,q,i,j) = coeffs(p,q,i,j)+dCoeff

  !            coeffs(p,q,i,j) = MAX(coeffs(p,q,i,j),0D0) ! Truncate negative nodes
  !            Mp = Mp+weight*coeffs(p,q,i,j)
            ENDDO !q
          ENDDO !p
  !        theta = avg/MAX(Mp,TINY(1D0)) ! Linear reduction factor
          theta = avg/MAX(avg+massChg,TINY(1D0))
          coeffs(:,:,i,j) = theta*coeffs(:,:,i,j) ! Reduce remaining (positive) nodes by reduction factor
        ENDDO!i
      ENDDO !j
      CASE(1)
        ! ================================================================================
        ! Use ZS (2010)-style linear rescaling
        ! ================================================================================
        DO j=1,ney
          DO i=1,nex
            ! Compute minimum value amongst nodal values (read from coefficients)
            valMin = MINVAL(coeffs(:,:,i,j))
            valMin = valMin - eps
            avg = elemAvgs(i,j)

            ! Compute rescaling factor
            theta = MIN(abs(avg)/(abs(valMin-avg)),1D0)
            ! Rescale approximating polynomial
            coeffs(:,:,i,j) = theta*(coeffs(:,:,i,j) - avg) + avg
          ENDDO !i
        ENDDO !j
      CASE DEFAULT
        ! Do nothing
      END SELECT
  END SUBROUTINE limitNodePositivity

  DOUBLE PRECISION FUNCTION massDiff(coeff,weight)
    IMPLICIT NONE
    ! Inputs
    DOUBLE PRECISION, INTENT(IN) :: coeff,weight
    ! Local Variables

    massDiff = -0.25D0*weight*MIN(coeff,0D0)
  END FUNCTION massDiff
END MODULE unspPositivityLimit
