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

SUBROUTINE coeff_update(A,u0,v0,gllNodes,gllWeights,gaussNodes,lagrangeDeriv,time,dt,dxel,dyel,nex,ney,norder,dgorder,&
                        dozshulimit,transient)
    IMPLICIT NONE

	! External functions
	REAL(KIND=8), EXTERNAL :: dadt ! RHS function for evolution ODE for expansion coefficent
	REAL(KIND=8), EXTERNAL :: vel_update ! Velocity update function

    ! Inputs
    INTEGER, INTENT(IN) :: nex,ney,norder,dgorder
    REAL(KIND=8), INTENT(IN) :: dxel,dyel,dt,time
    REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: gllNodes,gllWeights,gaussNodes
    REAL(KIND=8), DIMENSION(0:norder,0:dgorder), INTENT(IN) :: lagrangeDeriv
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(IN) :: u0,v0
    LOGICAL, INTENT(IN) :: dozshulimit,transient

    ! Outputs
    REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(INOUT) :: A

    ! Local Variables
    INTEGER :: i,j,l,m,stage
    REAL(KIND=8) :: coef1
    REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder) :: A1,A2
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder) :: quadVals,u,v
	REAL(KIND=8), DIMENSION(0:nex+1,1:ney,0:1,0:dgorder) :: edgeValsEW
	REAL(KIND=8), DIMENSION(1:nex,0:ney+1,0:1,0:dgorder) :: edgeValsNS
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder) :: Ghat


	! ######################
	! Time step using Strong-Stability Preserving RK3
	! ######################
	coef1 = dt/dxel/dyel
    u = u0
    v = v0

    A1 = A
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
	    CALL evalExpansion(quadVals,edgeValsNS,edgeValsEW,A1,dgorder,norder,nex,ney,dozshulimit)
	    CALL numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,dgorder,norder,nex,ney)

        IF( isnan(max(maxval(abs(Fhat)),maxval(abs(Ghat)) )) )  then
           write(*,*) 'fuck'
            stop
        endif

        ! Forward step of SSPRK3
        	DO i=1,nex
            	DO j=1,ney

        		DO l=0,norder
            		DO m=0,norder
	            		A2(i,j,l,m) = A1(i,j,l,m)+ &
	        				coef1*dadt(i,j,l,m,quadVals(i,j,:,:),Fhat,Ghat,u(i,j,:,:),v(i,j,:,:),gllWeights,&
                                   lagrangeDeriv,dxel,dyel,dgorder,norder,nex,ney)
	            	ENDDO !m
	        	ENDDO!l
            	ENDDO !j
        	ENDDO!i

		SELECT CASE(stage)
		CASE(1)
			A1 = A2
		CASE(2)
			A1 = 0.75D0*A + 0.25D0*A2
		CASE(3)
			A1 = A/3d0 + 2D0*A2/3D0
		END SELECT
    ENDDO !stage
    A = A1

END SUBROUTINE coeff_update

REAL(KIND=8) FUNCTION dadt(i,j,l,m,quadVals,Fhat,Ghat,uIn,vIn,gllWeight,lagrangeDeriv,dxel,dyel,dgorder,norder,nex,ney)
! Computes dadt for element (i,j) for use in SSPRK3 update step
! Does this in 3 steps: dadt = [2/(dxel*dyel)](A+B+C) where
!	1) A = Int(Int( dyel*dPs/dxi*Pt*F+dxel*dPt/deta*Ps*G )) -- Interior contribution to change
!	2) B = -dxel*Int( Ps*(Ghat(eta=1) - (-1)**t*Ghat(eta=-1)) ) -- Contribution through North/South faces
!	3) C = -dyel*Int( Pt*(Fhat(xi=1) - (-1)**s*Fhat(xi=-1)) ) -- Contribution through East/West faces
!
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: i,j,l,m,dgorder,norder,nex,ney
	REAL(KIND=8), INTENT(IN) :: dxel,dyel
 	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: uIn,vIn,quadVals
	REAL(KIND=8), DIMENSION(0:dgorder,0:dgorder), INTENT(IN) :: lagrangeDeriv
	REAL(KIND=8), DIMENSION(0:dgorder), INTENT(IN) :: gllWeight
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(IN) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(IN) :: Ghat

	! Local variables
	REAL(KIND=8) :: A,B,C,D

	! Step 1: Compute A
    A = (dyel/gllWeight(l))*SUM( gllWeight(:)*uIn(:,m)*quadVals(:,m)*lagrangeDeriv(l,:) )
    B = (dxel/gllWeight(m))*SUM( gllWeight(:)*vIn(l,:)*quadVals(l,:)*lagrangeDeriv(m,:) )

	! Step 2: Compute NS Flux Contribution
    C = 0D0
    IF( m .eq. dgorder) THEN
        C = C + Ghat(i,j,l)
    ELSEIF( m .eq. 0) THEN
        C = C - Ghat(i,j-1,l)
    ENDIF
	C = (-1D0*dxel/gllWeight(m))*C

	! Step 3: Compute EW Flux Contribution
    D = 0D0
    IF( l .eq. dgorder) THEN
        D = D + Fhat(i,j,m)
    ELSEIF( l .eq. 0) THEN
        D = D - Fhat(i-1,j,m)
    ENDIF
    D = (-1D0*dyel/gllWeight(l))

	dadt = (2D0)*(A+B+C+D)

END FUNCTION dadt

SUBROUTINE numFlux(Fhat,Ghat,u,v,edgeValsNS,edgeValsEW,dgorder,norder,nex,ney)
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: dgorder,norder,nex,ney
    REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder) :: u, v
	REAL(KIND=8), DIMENSION(0:nex+1,1:ney,0:1,0:dgorder), INTENT(IN) :: edgeValsEW
	REAL(KIND=8), DIMENSION(1:nex,0:ney+1,0:1,0:dgorder), INTENT(IN) :: edgeValsNS

	! Outputs
	REAL(KIND=8), DIMENSION(0:nex,1:ney,0:dgorder), INTENT(OUT) :: Fhat
	REAL(KIND=8), DIMENSION(1:nex,0:ney,0:dgorder), INTENT(OUT) :: Ghat

	! Local Variables
	INTEGER :: i,j,k,l,m,s,t
	REAL(KIND=8), DIMENSION(0:dgorder) :: u_tmp, v_tmp
	INTEGER :: which_el,which_sign

	DO i=1,nex
		DO j=1,ney
	
			u_tmp = u(i,j,dgorder,:) 
			v_tmp = v(i,j,:,dgorder) 

			DO k=0,dgorder	
				! Compute numerical flux through North edges (Ghat) (There are dgorder+1 nodes per element)
				which_sign = 1-(1-INT(SIGN(1D0,v_tmp(k))) )/2
				which_el = j + (1-INT(SIGN(1D0,v_tmp(k))) )/2
				Ghat(i,j,k) = v_tmp(k)*edgeValsNS(i,which_el,which_sign,k)

				! Compute numerical flux through East edges (Fhat) (There are dgorder+1 nodes per element)
				which_sign = 1-(1-INT(SIGN(1D0,u_tmp(k))) )/2
				which_el = i + (1-INT(SIGN(1D0,u_tmp(k))) )/2
				Fhat(i,j,k) = u_tmp(k)*edgeValsEW(which_el,j,which_sign,k)
			ENDDO !k
		ENDDO !j
	ENDDO !i

	! Extend fluxes periodically 
	Ghat(1:nex,0,:) = Ghat(1:nex,ney,:)
	Fhat(0,1:ney,:) = Fhat(nex,1:ney,:)

END SUBROUTINE numFlux

SUBROUTINE evalExpansion(quadVals,edgeValsNS,edgeValsEW,Ain,dgorder,norder,nex,ney,dozshulimit)
! Evaluate current expansion at interior quadrature locations and quadrature locations along EW and NS edges of each element
	IMPLICIT NONE
	! Inputs
	INTEGER, INTENT(IN) :: dgorder,norder,nex,ney
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:norder,0:norder), INTENT(INOUT) :: Ain
    LOGICAL, INTENT(IN) :: dozshulimit

	! Outputs
	REAL(KIND=8), DIMENSION(1:nex,1:ney,0:dgorder,0:dgorder), INTENT(OUT) :: quadVals
	REAL(KIND=8), DIMENSION(0:nex+1,1:ney,0:1,0:dgorder), INTENT(OUT) :: edgeValsEW
	REAL(KIND=8), DIMENSION(1:nex,0:ney+1,0:1,0:dgorder), INTENT(OUT) :: edgeValsNS

	! Local Variables
	INTEGER :: i,j,s,t,l
	REAL(KIND=8), DIMENSION(0:norder) :: arr2
	REAL(KIND=8), DIMENSION(0:norder,0:norder) :: tmp1, tmp2
    REAL(KIND=8) :: theta,elemAvg,valMin

	arr2 = (/ ((-1D0)**i , i=0,norder) /) ! P_l(-1) = (-1)**l

	DO l=0,norder
		tmp1(:,l) = arr2
		tmp2(l,:) = arr2
	ENDDO !l

	DO i=1,nex
		DO j=1,ney
            ! Evaluate expansion at interior quad points
            quadVals(i,j,:,:) = Ain(i,j,:,:)

		    ! Evaluate expansion at NS edges
			edgeValsNS(i,j,1,:) = Ain(i,j,:,norder) !SUM(Leg(t,:)*SUM(Ain(i,j,:,:),2))
			edgeValsNS(i,j,0,:) = Ain(i,j,:,0) !SUM(Leg(t,:)*SUM(Ain(i,j,:,:)*tmp2,2))

			! Evaluate expansion at EW edges
			edgeValsEW(i,j,1,:) = Ain(i,j,norder,:) !SUM(Leg(t,:)*SUM(Ain(i,j,:,:),1))
			edgeValsEW(i,j,0,:) = Ain(i,j,0,:) !SUM(Leg(t,:)*SUM(Ain(i,j,:,:)*tmp1,1))
		ENDDO !j
	ENDDO !i

	IF(dozshulimit) THEN ! Use polynomial modifications from Zhang and Shu (2009)
        write(*,*) 'Warning: Not Implemented yet!' 
        stop
	ENDIF

	! Extend edge values periodically
	edgeValsEW(0,:,:,:) = edgeValsEW(nex,:,:,:)
	edgeValsEW(nex+1,:,:,:) = edgeValsEW(1,:,:,:)

	edgeValsNS(:,0,:,:) = edgeValsNS(:,ney,:,:)
	edgeValsNS(:,ney+1,:,:) = edgeValsNS(:,1,:,:)

END SUBROUTINE evalExpansion

REAL(KIND=8) FUNCTION vel_update(t)
		IMPLICIT NONE
		REAL(KIND=8), INTENT(IN) :: t
	    REAL(KIND=8) :: pi
	    REAL(KIND=8), parameter :: t_period = 5.d0

	    pi = DACOS(-1D0)
	    vel_update = DCOS(pi*t/t_period)
END FUNCTION vel_update
