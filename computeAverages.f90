SUBROUTINE computeAverages(elemAvgs,gllWeights,coeffs,nex,ney,nOrder,nQuad)
  USE testParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nex,ney,nOrder,nQuad
  DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder,1:nex,1:ney), INTENT(IN) :: coeffs
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: gllWeights
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(OUT) :: elemAvgs
  ! Local variables
  INTEGER :: i,j,l
  DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder) :: tmpArray

  DO i=1,nex
    DO j=1,ney
      DO l=0,nQuad
        tmpArray(:,l) = gllWeights(l)*gllWeights(:)*coeffs(:,l,i,j)
      ENDDO!l
      elemAvgs(i,j) = 0.25D0*SUM(tmpArray)
    ENDDO !j
  ENDDO !i

END SUBROUTINE computeAverages
