SUBROUTINE limitMeanPositivity(coeffs,elemAvgs,lagGaussVal,gqWeights,gllWeights,&
                    nex,ney,nOrder,gqOrder,nQuad)
    ! Modifies approximating polynomial so that element mean value remains non-negative
    ! after forward update according to Zhang and Shu (2010) Thm 5.2
    USE testParameters
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,nOrder,gqOrder,nQuad
    DOUBLE PRECISION, DIMENSION(0:gqOrder), INTENT(IN) :: gqWeights
    DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: gllWeights
    DOUBLE PRECISION, DIMENSION(0:norder,0:gqOrder), INTENT(IN) :: lagGaussVal
    DOUBLE PRECISION, DIMENSION(1:nex,1:ney), INTENT(IN) :: elemAvgs
    ! Outputs
    DOUBLE PRECISION, DIMENSION(1:nex,1:ney,0:nOrder,0:nOrder), INTENT(INOUT) :: coeffs
    ! Local variables
    INTEGER :: i,j,l,beta,index
    DOUBLE PRECISION :: avg,valMin,leftTrace,rightTrace,topTrace,botTrace,magicPt
    DOUBLE PRECISION :: eps,theta,traceMin,massChg,weight
    DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder) :: tmpArray
    DOUBLE PRECISION, DIMENSION(0:gqOrder) :: magicTmp

    ! External function interface
    INTERFACE
      DOUBLE PRECISION FUNCTION massDiff(coeff,weight)
        IMPLICIT NONE
        ! Inputs
        DOUBLE PRECISION, INTENT(IN) :: coeff,weight
      END FUNCTION massDiff

    END INTERFACE

    eps = epsilon(1D0)
    SELECT CASE(limitingMeth)
    CASE(1,3)
        DO i=1,nex
          DO j=1,ney
            ! Step 1: Initialize mean and minimum value
            avg = elemAvgs(i,j)
            valMin = HUGE(1D0)

            ! Step 2: Compute minimum value between traces and "magic point"
            DO beta=0,gqOrder
              ! Look at left, right, top, and bottom trace values for this
              ! Gauss quadrature point (indexed by beta)
              leftTrace = SUM(coeffs(i,j,0,:)*lagGaussVal(:,beta))
              rightTrace = SUM(coeffs(i,j,nOrder,:)*lagGaussVal(:,beta))
              botTrace = SUM(coeffs(i,j,:,0)*lagGaussVal(:,beta))
              topTrace = SUM(coeffs(i,j,:,nOrder)*lagGaussVal(:,beta))

              ! Update minimum
              valMin = MIN(valMin,leftTrace,rightTrace,botTrace,topTrace)

              magicTmp(beta) = gqWeights(beta)*(mu1*(rightTrace+leftTrace)+mu2*(topTrace+botTrace))
            ENDDO !beta
            traceMin = valMin

            ! Compute "magic point" value
            magicPt = avg-0.25D0*zsMinWeight*SUM(magicTmp)
            magicPt = magicPt/(1D0-zsMinWeight)

            ! Update minimum
            !valMin = MIN(valMin,magicPt,MINVAL(coeffs(i,j,:,:)))
            valMin = MIN(valMin,magicPt)

            !valMin = MINVAL(coeffs(i,j,:,:))
            valMin = valMin - eps

            ! Step 3: Compute theta and rescale polynomial
            theta = MIN(abs(avg)/(abs(valMin-avg)),1D0)

            IF(magicPt .lt. traceMin .AND. theta .lt. 1D0) THEN
              magCount = magCount+1
            ENDIF

            ! Rescale reconstructing polynomial for (i,j)th element
            coeffs(i,j,:,:) = theta*(coeffs(i,j,:,:) - avg) + avg

          ENDDO !j
        ENDDO !i
      CASE(2)
        ! -- TODO (Devin): Make this more efficient / readable
        DO i = 1,nex
          DO j = 1,ney
            ! Step 1: Pull element average
            avg = elemAvgs(i,j)

            ! Step 2: Truncate along edges and evaluate magic point
            massChg = 0D0
            magicPt = 0D0
            DO beta = 0,nOrder-1

              index = nOrder-beta
              leftTrace = coeffs(i,j,0,beta)
              rightTrace = coeffs(i,j,nOrder,index)
              topTrace = coeffs(i,j,beta,nOrder)
              botTrace = coeffs(i,j,index,0)

              weight = gllWeights(0)*gllWeights(beta)

              ! Truncate along left edge
              massChg = massChg + massDiff(leftTrace,weight)
              coeffs(i,j,0,beta) = MAX(leftTrace,0D0)

              ! Truncate along right edge
              massChg = massChg + massDiff(rightTrace,weight)
              coeffs(i,j,nOrder,index) = MAX(rightTrace,0D0)

              ! Truncate along top edge
              massChg = massChg + massDiff(topTrace,weight)
              coeffs(i,j,beta,nOrder) = MAX(topTrace,0D0)

              ! Truncate along bottom edge
              massChg = massChg + massDiff(botTrace,weight)
              coeffs(i,j,index,0) = MAX(botTrace,0D0)

              ! Track modification to magic point
!              magicPt = magicPt + gllWeights(beta)*(mu1*(rightTrace+leftTrace)+mu2*(topTrace+botTrace))
            ENDDO !beta
!            magicPt = avg - 0.25D0*gllWeights(0)*magicPt
!            IF(magicPt .lt. 0D0) THEN
!              magCount = magCount + 1

              ! Update mass change due to truncating magic point
!              massChg = massChg - magicPt

              ! Truncate internal nodes to zero
!              coeffs(i,j,1:nOrder-1,1:nOrder-1) = 0D0
!            ENDIF

            ! Rescale remaining coefficents to maintain conservation
            theta = avg/MAX(avg+massChg,TINY(1D0))
            coeffs(i,j,:,:) = theta*coeffs(i,j,:,:)

          ENDDO !j
        ENDDO !i
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
  DOUBLE PRECISION, DIMENSION(1:nex,1:ney,0:nOrder,0:nOrder), INTENT(INOUT) :: coeffs
  ! Local variables
  INTEGER :: i,j,p,q
  DOUBLE PRECISION :: valMin,eps,theta,avg
  DOUBLE PRECISION :: Mt,mP
  DOUBLE PRECISION, DIMENSION(0:nOrder,0:nOrder) :: tmpArray

  eps = epsilon(1d0)

  SELECT CASE(limitingMeth)
    CASE(1)
      ! ================================================================================
      ! Use ZS (2010)-style linear rescaling
      ! ================================================================================
      DO i=1,nex
        DO j=1,ney
          ! Compute minimum value amongst nodal values (read from coefficients)
          valMin = MINVAL(coeffs(i,j,:,:))
          valMin = valMin - eps
          avg = elemAvgs(i,j)

          ! Compute rescaling factor
          theta = MIN(abs(avg)/(abs(valMin-avg)),1D0)
          ! Rescale approximating polynomial
          coeffs(i,j,:,:) = theta*(coeffs(i,j,:,:) - avg) + avg
        ENDDO !j
      ENDDO !i
    CASE(2,3)
      ! ================================================================================
      ! Use TMAR limiting at GLL nodes
      ! ================================================================================
      DO i=1,nex
        DO j=1,ney
          Mt = 0D0
          Mp = 0D0
          DO p=0,nOrder
            DO q=0,nOrder
              Mt = Mt+gllWeights(p)*gllWeights(q)*coeffs(i,j,p,q)
              coeffs(i,j,p,q) = MAX(coeffs(i,j,p,q),0D0) ! Truncate negative nodes
              Mp = Mp+gllWeights(p)*gllWeights(q)*coeffs(i,j,p,q)
            ENDDO !q
          ENDDO !p
          theta = MAX(Mt,0D0)/MAX(Mp,TINY(1D0)) ! Linear reduction factor
          coeffs(i,j,:,:) = theta*coeffs(i,j,:,:) ! Reduce remaining (positive) nodes by reduction factor
        ENDDO!j
      ENDDO !i
    CASE DEFAULT
      write(*,*) '***** ERROR *****'
      write(*,*) 'IN limitNodePositivity()... NO SUCH LIMITING METHOD'
      write(*,*) 'Stopping...'
      STOP
    END SELECT
END SUBROUTINE limitNodePositivity

DOUBLE PRECISION FUNCTION massDiff(coeff,weight)
  IMPLICIT NONE
  ! Inputs
  DOUBLE PRECISION, INTENT(IN) :: coeff,weight
  ! Local Variables

  massDiff = -0.25D0*weight*MIN(coeff,0D0)
END FUNCTION massDiff
