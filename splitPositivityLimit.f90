MODULE splitPositivityLimit
  IMPLICIT NONE
  PUBLIC :: limitNodePositivity,limitMeanPositivity
CONTAINS
  SUBROUTINE limitNodePositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
    ! Modifies approximating polynomial so that nodal values are non-negative for output
    USE testParameters
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
    DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
    DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: avgVals
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    ! Local Variables
    INTEGER :: j,l
    DOUBLE PRECISION :: avg,theta,valMin,valMax,eps,Mt,Mp

    eps = epsilon(1d0)
    SELECT CASE(limitingMeth)
    CASE(1)
      ! Use linear rescaling for polynomial modification
      DO j=1,nelem
        avg = avgVals(j)
        valMin = MINVAL(qIn(:,j))-eps
        ! -- Compute rescaling factor
        theta = MIN( abs(avg/(valMin-avg)),1D0 )

        ! -- Rescale polynomial
        qIn(:,j) = theta*(qIn(:,j)-avg) + avg
      ENDDO!j
    CASE(2:3)
      ! Use "mass aware truncation" for for polynomial modification
      DO j=1,nelem
        Mp = 0D0
        Mt = 0D0
        avg = avgVals(j)

        DO l=0,nNodes
  !          Mt = Mt + qWeights(l)*qIn(l,j)
            qIn(l,j) = MAX(0D0,qIn(l,j)) ! Zero out negative nodes
            Mp = Mp + qWeights(l)*qIn(l,j)
        ENDDO !l
  !      theta = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
        theta = 2D0*ABS(avg)/MAX(Mp,TINY(1D0))
        qIn(:,j) = theta*qIn(:,j) ! Reduce remaining positive nodes by reduction factor
      ENDDO !j
    CASE DEFAULT
      ! Do not rescale polynomial

    END SELECT !limitingMeth
  END SUBROUTINE limitNodePositivity

  SUBROUTINE limitMeanPositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
    ! Modifies approximating polynomial so that element mean value remains non-negative
    ! according to Zhang and Shu (2010) Thm 2.2
    USE testParameters
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
    DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
    DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: avgVals
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
    ! Local Variables
    INTEGER :: j
    DOUBLE PRECISION :: avg,theta,valMin,eps,qStar

    eps = epsilon(1d0)
    DO j=1,nelem
      avg = avgVals(j)

      SELECT CASE(limitingMeth)
        CASE(1)
          ! Compute 'magic point' value from Zhang and Shu 2011: Survey and New Developments
          qStar = avg-qWeights(0)*qIn(0,j)-qWeights(0)*qIn(nNodes,j)
          qStar = qStar/(1D0-2D0*qWeights(0))

          valMin = MIN(qStar,qIn(0,j),qIn(nNodes,j))-eps
        CASE(4)
          ! Use all points to determine new rescaling
          valMin = MINVAL(qIn(:,j))-eps
        CASE DEFAULT
          ! Do nothing
          RETURN
      END SELECT !limitingMeth

      ! Compute rescaling factor
      theta = MIN(avg/abs(valMin-avg),1D0)

      ! Rescale polynomial
      qIn(:,j) = theta*(qIn(:,j)-avg) + avg
    ENDDO
  END SUBROUTINE limitMeanPositivity
END MODULE splitPositivityLimit
