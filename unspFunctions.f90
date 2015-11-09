MODULE unspFunctions
  IMPLICIT NONE
  PUBLIC :: numFlux,evalExpansion,dadt,computeAverages,vel_update
CONTAINS
  DOUBLE PRECISION FUNCTION dadt(i,j,l,m,quadVals,Fhat,Ghat,uIn,vIn,gllWeight,lagrangeDeriv,dxel,dyel,nQuadNodes,norder,nex,ney)
  ! Computes dadt for element (i,j) for use in SSPRK3 update step
  ! Does this in 3 steps: dadt = 2*(A+B+C) where
  !	1) A+B = Int(Int( dyel*dPs/dxi*Pt*F+dxel*dPt/deta*Ps*G )) -- Interior contribution to change
  !	2) C = -dxel*Int( Ps*(Ghat(eta=1) - (-1)**t*Ghat(eta=-1)) ) -- Contribution through North/South faces
  !	3) D = -dyel*Int( Pt*(Fhat(xi=1) - (-1)**s*Fhat(xi=-1)) ) -- Contribution through East/West faces
  !
  	IMPLICIT NONE
  	! Inputs
  	INTEGER, INTENT(IN) :: i,j,l,m,nQuadNodes,norder,nex,ney
  	REAL(KIND=8), INTENT(IN) :: dxel,dyel
   	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes), INTENT(IN) :: uIn,vIn,quadVals
  	REAL(KIND=8), DIMENSION(0:norder,0:nQuadNodes), INTENT(IN) :: lagrangeDeriv
  	REAL(KIND=8), DIMENSION(0:nQuadNodes), INTENT(IN) :: gllWeight
  	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nex,1:ney), INTENT(IN) :: Fhat
  	REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex,0:ney), INTENT(IN) :: Ghat

  	! Local variables
  	REAL(KIND=8) :: A,B,C,D

  	! Step 1: Compute A & B
    A = (dyel/gllWeight(l))*SUM( gllWeight(:)*uIn(:,m)*quadVals(:,m)*lagrangeDeriv(l,:) )
    B = (dxel/gllWeight(m))*SUM( gllWeight(:)*vIn(l,:)*quadVals(l,:)*lagrangeDeriv(m,:) )

  	! Step 2: Compute NS Flux Contribution
    C = 0D0
    IF(m .eq. nQuadNodes) THEN
      C = C + Ghat(l,i,j)
    ELSEIF( m .eq. 0) THEN
      C = C - Ghat(l,i,j-1)
    ENDIF

  	C = (-1D0*dxel/gllWeight(m))*C

  	! Step 3: Compute EW Flux Contribution
    D = 0D0
    IF(l .eq. nQuadNodes) THEN
      D = D + Fhat(m,i,j)
    ELSEIF( l .eq. 0) THEN
      D = D - Fhat(m,i-1,j)
    ENDIF

    D = (-1D0*dyel/gllWeight(l))*D

  	dadt = (2D0)*(A+B+C+D)

  END FUNCTION dadt

  SUBROUTINE numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,nQuadNodes,norder,nex,ney)
  	IMPLICIT NONE
  	! Inputs
  	INTEGER, INTENT(IN) :: nQuadNodes,norder,nex,ney
    REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney) :: u, v
  	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,0:nex+1,1:ney), INTENT(IN) :: edgeValsEW
  	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,1:nex,0:ney+1), INTENT(IN) :: edgeValsNS

  	! Outputs
  	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nex,1:ney), INTENT(OUT) :: Fhat
  	REAL(KIND=8), DIMENSION(0:nQuadNodes,1:nex,0:ney), INTENT(OUT) :: Ghat

  	! Local Variables
  	INTEGER :: i,j,k,l,m,s,t
  	REAL(KIND=8), DIMENSION(0:nQuadNodes) :: u_tmp, v_tmp
  	INTEGER :: which_el,which_sign

  	DO j=1,ney
  		DO i=1,nex

  			u_tmp = u(nQuadNodes,:,i,j)
  			v_tmp = v(:,nQuadNodes,i,j)

  			DO k=0,nQuadNodes
  				! Compute numerical flux through North edges (Ghat) (There are nQuadNodes+1 nodes per element)
  				which_sign = 1-(1-INT(SIGN(1D0,v_tmp(k))) )/2
  				which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2
  				Ghat(k,i,j) = v_tmp(k)*edgeValsNS(which_sign,k,i,which_el)

  				! Compute numerical flux through East edges (Fhat) (There are nQuadNodes+1 nodes per element)
  				which_sign = 1-(1-INT(SIGN(1D0,u_tmp(k))) )/2
  				which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
  				Fhat(k,i,j) = u_tmp(k)*edgeValsEW(which_sign,k,which_el,j)
  			ENDDO !k
  		ENDDO !i
  	ENDDO !j

      j = 0
      DO i=1,nex
  			v_tmp = v(:,0,i,j+1)
  			DO k=0,nQuadNodes
    			! Compute numerical flux through North edges (Ghat) (There are nQuadNodes+1 nodes per element)
    			which_sign = 1-(1-INT(SIGN(1D0,v_tmp(k))) )/2
    			which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2
      		Ghat(k,i,j) = v_tmp(k)*edgeValsNS(which_sign,k,i,which_el)
        ENDDO!k
      ENDDO!i

      i=0
      DO j=1,ney
        u_tmp = u(0,:,i+1,j)
        DO k=0,nQuadNodes
    			! Compute numerical flux through East edges (Fhat) (There are nQuadNodes+1 nodes per element)
    			which_sign = 1-(1-INT(SIGN(1D0,u_tmp(k))) )/2
    			which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
    			Fhat(k,i,j) = u_tmp(k)*edgeValsEW(which_sign,k,which_el,j)
  	    ENDDO !k
      ENDDO !j

  END SUBROUTINE numFlux

  SUBROUTINE evalExpansion(quadVals,edgeValsNS,edgeValsEW,coeffs,nQuadNodes,norder,nex,ney)
  ! Evaluate current expansion at interior quadrature locations and quadrature locations along EW and NS edges of each element
  	IMPLICIT NONE
  	! Inputs
  	INTEGER, INTENT(IN) :: nQuadNodes,norder,nex,ney
  	REAL(KIND=8), DIMENSION(0:norder,0:norder,1:nex,1:ney), INTENT(IN) :: coeffs

  	! Outputs
  	REAL(KIND=8), DIMENSION(0:nQuadNodes,0:nQuadNodes,1:nex,1:ney), INTENT(OUT) :: quadVals
  	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,0:nex+1,1:ney), INTENT(OUT) :: edgeValsEW
  	REAL(KIND=8), DIMENSION(0:1,0:nQuadNodes,1:nex,0:ney+1), INTENT(OUT) :: edgeValsNS

  	! Local Variables
  	INTEGER :: i,j

  	DO j=1,nex
  		DO i=1,nex
        ! Evaluate expansion at interior GLL quad points
        quadVals(:,:,i,j) = coeffs(:,:,i,j)

  	    ! Evaluate expansion at NS edges
  			edgeValsNS(1,:,i,j) = coeffs(:,norder,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:),2))
  			edgeValsNS(0,:,i,j) = coeffs(:,0,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:)*tmp2,2))

  			! Evaluate expansion at EW edges
  			edgeValsEW(1,:,i,j) = coeffs(norder,:,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:),1))
  			edgeValsEW(0,:,i,j) = coeffs(0,:,i,j) !SUM(Leg(t,:)*SUM(coeffs(i,j,:,:)*tmp1,1))
  		ENDDO !j
  	ENDDO !i

  	! Extend edge values periodically
    edgeValsEW(:,:,0,:) = edgeValsEW(:,:,nex,:)
    edgeValsEW(:,:,nex+1,:) = edgeValsEW(:,:,1,:)

    edgeValsNS(:,:,:,0) = edgeValsNS(:,:,:,ney)
    edgeValsNS(:,:,:,ney+1) = edgeValsNS(:,:,:,1)

  END SUBROUTINE evalExpansion

  DOUBLE PRECISION FUNCTION vel_update(t)
  		IMPLICIT NONE
  		REAL(KIND=8), INTENT(IN) :: t
      REAL(KIND=8) :: pi
      REAL(KIND=8), parameter :: t_period = 5.d0

      pi = DACOS(-1D0)
      vel_update = DCOS(pi*t/t_period)
  END FUNCTION vel_update

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
END MODULE unspFunctions
