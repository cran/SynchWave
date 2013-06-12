!   synsq_cwt_squeeze
!   This function is only to be used from synsq_cwt_squeeze()
!
!--------------------------------------------------------------------------------
!    Synchrosqueezing Toolbox
!    Original Authors: Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/),
!                      Hau-Tieng Wu
!--------------------------------------------------------------------------------
!
!    Fortran port by Dongik Jang (dongik.s.jang@gmail.com)

SUBROUTINE synsq_cwt_squeeze(Wx, na, N, w, as, fs, dfs, lfm1, lfM2, out)
    IMPLICIT NONE
    
    INTEGER :: ai, bi, k, N, na
    DOUBLE PRECISION :: lfM2, lfm1, eps
    DOUBLE COMPLEX :: Wx(na, N), out(na, N)
    DOUBLE PRECISION :: w(na, N), Wxr(na, N), Wxi(na, N) 
    DOUBLE PRECISION :: as(na), fs(na), dfs(na), Txr(na, N), Txi(na, N)
    DOUBLE PRECISION :: dfsinv(na), asrtinv(na)
    DOUBLE PRECISION :: Wxbr(na), Wxbi(na), Txbr(na), Txbi(na), wab(na)

    
    Wxr = DBLE(Wx)
    Wxi = AIMAG(Wx)
    eps = 0.00000001

    DO ai = 1, na
        dfsinv(ai) = 1/dfs(ai)
        asrtinv(ai) = 1/dsqrt( as(ai) )
    END DO
    
    
    DO bi = 1, N
        DO ai = 1, na
            Wxbr(ai) = Wxr(ai, bi)
            Wxbi(ai) = Wxi(ai, bi)
            Txbr(ai) = 0.0d0
            Txbi(ai) = 0.0d0
            wab(ai) = w(ai, bi)
        END DO
        
        DO ai = 1, na
            IF( disfinite(wab(ai)) .AND. (wab(ai) .GT. 0) ) THEN
!                k = INT( FLOOR( 0.5d0 + (DBLE(na-1))/(lfM2 - lfm1) * (DLOG(wab(ai))/LOG(2.0) - lfm1) ) ) + 1
                k = INT( FLOOR( 0.5d0 + (DBLE(na-1))/(lfM2 - lfm1) * (log2(wab(ai)) - lfm1) ) ) + 1
            	IF( isfinite(k) .AND. (k .GE. 1) .AND. (k .LE. na) ) THEN
                	Txbr(k) = Txbr(k) + (Wxbr(ai) * asrtinv(ai) * dfsinv(k))
                	Txbi(k) = Txbi(k) + (Wxbi(ai) * asrtinv(ai) * dfsinv(k))
                END IF
			END IF
        END DO
        
        DO ai = 1, na
            Txr(ai, bi) = Txbr(ai)
            Txi(ai, bi) = Txbi(ai)        
        END DO
                    
    END DO
    out = DCMPLX(Txr, Txi)
    
    CONTAINS
    
    LOGICAL FUNCTION isfinite(a)
        IMPLICIT NONE
        INTEGER :: a
        isfinite = (a-a) .EQ. 0
    END FUNCTION isfinite
    
    LOGICAL FUNCTION disfinite(a)
        IMPLICIT NONE
        DOUBLE PRECISION :: a
        disfinite = (a-a) .EQ. 0
    END FUNCTION disfinite
    
    DOUBLE PRECISION FUNCTION log2(a)
        IMPLICIT NONE
        DOUBLE PRECISION :: a
        
        log2 = DLOG(a)/DLOG( 2.0d0)
    
    END FUNCTION log2
    
END SUBROUTINE synsq_cwt_squeeze

!   diff_W
!
!--------------------------------------------------------------------------------
!    Synchrosqueezing Toolbox
!    Original Authors: Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/),
!                      Hau-Tieng Wu
!--------------------------------------------------------------------------------
!
!    Fortran port by Dongik Jang (dongik.s.jang@gmail.com)

SUBROUTINE diff_W (W, na, N, dt, dorder, out)
    IMPLICIT NONE
    
    INTEGER :: ai, bi, k, N, na
    DOUBLE COMPLEX :: W(na, N), out(na, N)
    DOUBLE PRECISION :: dt
    INTEGER :: dorder, iedval
    DOUBLE PRECISION :: Wr(na, N), Wi(na, N), dWr(na, N), dWi(na, N) 
    
    Wr = DBLE(W)
    Wi = AIMAG(W)
    IF (SUM(DABS(Wi)) .EQ. DBLE(0.0) ) THEN
        iedval = 1
    END IF
    
    ! Take derivatives of a certain order, in the column direction
    IF (dorder .EQ. 1) THEN
        ! 1st order Forward Difference 
        DO bi = 1, (N-1)
            DO ai = 1, na
                dWr(ai, bi) = (Wr(ai, bi+1) - Wr(ai, bi))/dt
                IF (iedval .EQ. 1) THEN
                    dWi(ai, bi) = (Wi(ai, bi+1) - Wi(ai, bi))/dt
                END IF
            END DO
        END DO
    ELSE IF (dorder .EQ. 2) THEN
        ! 2nd order Forward Difference
        DO bi = 1, (N-2)
            DO ai = 1, na
                dWr(ai, bi) = (4*Wr(ai, bi+1) - Wr(ai, bi+2) - 3*Wr(ai, bi))/(2*dt)
                IF (iedval .EQ. 1) THEN
                    dWi(ai, bi) = (4*Wi(ai, bi+1) - Wi(ai, bi+2) - 3*Wi(ai, bi))/(2*dt)
                END IF
            END DO
        END DO
    ELSE IF (dorder .EQ. 4) THEN
        ! 4nd order Forward Difference
        DO bi = 3, (N-2)
            DO ai = 1, na
                dWr(ai, bi) = (-Wr(ai, bi+2) + 8*Wr(ai, bi+1) - &
                    8*Wr(ai, bi-1) + Wr(ai, bi-2))/(12*dt)
                IF (iedval .EQ. 1) THEN
                    dWi(ai, bi) = (-Wi(ai, bi+2) + 8*Wi(ai, bi+1) - &
                        8*Wi(ai, bi-1) + Wi(ai, bi-2))/(12*dt)
                END IF  
            END DO
        END DO
    END IF
    out = DCMPLX(dWr, dWi)
    
END SUBROUTINE diff_W

! curve_ext 
! Original Author: Jianfeng Lu (now at NYU CIMS)
! Maintainer: Eugene Brevdo

!--------------------------------------------------------------------------------
!    Synchrosqueezing Toolbox
!    Original Authors: Eugene Brevdo (http://www.math.princeton.edu/~ebrevdo/),
!                      Hau-Tieng Wu
!--------------------------------------------------------------------------------
!
!    Fortran port by Dongik Jang (dongik.s.jang@gmail.com)

SUBROUTINE curve_ext (Tx, na, N, fs, lambda, FreqOut, EnergyOut)
    IMPLICIT NONE
    
    INTEGER ::  N, na
    DOUBLE COMPLEX :: Tx(na, N)
    DOUBLE PRECISION :: Txr(na, N), Txi(na, N), fs(na), sumTxr, sumTxi
    DOUBLE PRECISION :: FreqOut(N), EnergyOut
    DOUBLE PRECISION :: lambda, eps, val, mval
    DOUBLE PRECISION :: Energy(na, N), FVal(na, N), sumval
    INTEGER :: Freq(N)
    INTEGER :: i, j ,k, ind
    
    
    eps = 1e-8
    Txr = DBLE(Tx)
    Txi = AIMAG(Tx)
    
    sumTxr = sum(abs(Txr))
    sumTxi = sum(abs(Txi))
    
!   Calculate in terms of energy - more numerically stable    
    Energy(:, :) = 0.0d0
    sumval = 0.0d0
    DO i = 1, N
        DO j = 1, na
            IF ( sumTxr .GT. 0 ) THEN
                Energy(j, i) = Txr(j, i) * Txr(j, i)
            END IF
                
            IF ( sumTxi .GT. 0 ) THEN
                Energy(j, i) = Energy(j, i) + Txi(j, i) * Txi(j, i)
            END IF
            sumval = sumval + Energy(j, i)
        END DO
    END DO
    
    DO i = 1, N
        DO j = 1, na
            Energy(j, i) = -log(Energy(j, i)/sumval + eps)
        END DO
    END DO
    
!   Simple first order dynamic program: find minimum energy to
!   traverse from frequency bin j at time i to frequency bin k at
!   time i+1, where the energy cost is proportional to
!   lambda*dist(j,k)
    
    DO j = 1, na
        Fval(j, 1) = Energy(j, 1)
    END DO
    
    DO i = 2, N
        DO j = 1, na
            FVal(j, i) = huge(0.0d0)
        END DO
        
        DO j = 1, na
            DO k = 1, na
                val = fs(j)-fs(k)
                val = FVal(k, i-1) + lambda * val * val
                IF ( FVal(j, i) .GT. val ) THEN
                    FVal(j, i) = val
                END IF
            END DO
                
            FVal(j, i) = FVal(j, i) + Energy(j, i)
        END DO
    END DO
    
!   Find the locations for the minimum values of our energy
!   functional at each time location i.  Store in freq.
    
    mval = FVal(1, N)
    DO i = 2, N
        Freq(i) = -1
    END DO
    
    Freq(N) = 0
    DO j = 2, na
        IF ( FVal(j, N) .LT. mval ) THEN
            mval = FVal(j, N)
            Freq(N) = j
        END IF
    END DO
    
    
    DO i = N-1, 1, -1
        IF ( Freq(i+1) .NE. -1 ) THEN
            ind = Freq(i+1)
            val = FVal(ind, i+1) - Energy(ind, i+1)
            DO j = 1, na
                mval = (fs(j) - fs(ind))
                mval = (val - FVal(j, i)) - lambda * mval * mval
                IF ( DABS(mval) .LT. eps) THEN
                    Freq(i) = j
                    EXIT
                END IF            
            END DO       
        END IF
    END DO 
    
    DO i = 1, N
!        FreqOut(i) = DBLE(Freq(i)) + 1.0d0 
        FreqOut(i) = DBLE(Freq(i)) 
    END DO
    
    EnergyOut = Fval(Freq(N), N)
    
!    CONTAINS
!    
!    DOUBLE PRECISION FUNCTION fdist (fs, na, i, j)
!        IMPLICIT NONE
!        INTEGER :: i, j, na
!        DOUBLE PRECISION :: fs(na)
!        
!        fdist = (fs(i)-fs(j)) * (fs(i)-fs(j))    
!        
!    END FUNCTION fdist
END SUBROUTINE curve_ext

!   imdilate
!
!    Author : Dongik Jang (dongik.s.jang@gmail.com)

SUBROUTINE imdilate (A, m, n, SE, ms, ns, out)
    IMPLICIT NONE
    
    INTEGER :: m, n, ms, ns
    INTEGER :: A(m,n), SE(ms, ns), out(m,n)
    INTEGER :: i, j, tmp(ms, ns), d1, d2, xc, yc
    INTEGER :: jid, jed, iid, ied
    
    xc = ceiling( DBLE(ms + 1)/2.0)
    yc = ceiling( DBLE(ns + 1)/2.0)
    DO i = 1, m
        DO j = 1, n
            tmp(:, :) = 0
            iid = max( 1, i - xc + 1 )
            ied = min(m, i - xc + ms)
            jid = max( 1, j - yc + 1 )
            jed = min(n, j - yc + ns)
            
            tmp(max(1, xc-i+iid):min(ms, xc-i+ied), max(1, yc-j+jid):min(ns, yc-j+jed)) = A(iid:ied, jid:jed)
             
            tmp = tmp * SE
            out(i,j) = maxval(tmp)
        END DO
    END DO
END SUBROUTINE imdilate 
