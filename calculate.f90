
subroutine calcu()
    use params
    implicit none

    integer :: i, j, n
    real(16) :: cn, cs, ce, cw, dn, ds, de, dw
    real(16) :: yp, tmult, resor

    ! Chapter 1, assembly of coefficients
    do i = 3, nim1
        do j = 2, njm1
            ! convection coefficients
            cn = 0.5*densit*(V(i  , j+1) + V(i-1, j+1))*SEWU(i)
            cs = 0.5*densit*(V(i  , j  ) + V(i-1, j  ))*SEWU(i)
            ce = 0.5*densit*(U(i+1, j  ) + U(i  , j  ))*SNS(j)
            cw = 0.5*densit*(U(i  , j  ) + U(i-1, j  ))*SNS(j)
            ! diffusion coefficients
            dn = viscos*SEWU(i)/DYNP(j)
            ds = viscos*SEWU(i)/DYPS(j)
            de = viscos*SNS(j)/DXEPU(i)
            dw = viscos*SNS(j)/DXPWU(i)
            ! assemble
            AN(i, j) = AMAX1(ABS(0.5*cn), dn) - 0.5*cn
            AS(i, j) = AMAX1(ABS(0.5*cs), ds) + 0.5*cs
            AE(i, j) = AMAX1(ABS(0.5*ce), de) - 0.5*ce
            AW(i, j) = AMAX1(ABS(0.5*cw), dw) + 0.5*cw
            DU(i, j) = SNS(j)
            SU(i, j) = DU(i, j)*(P(i-1, j) - P(i, j))
            SP(i, j) = 0.0
        end do
    end do

    ! Chapter 2, Problem modifications
    ! Top wall
    yp = YV(nj) - Y(njm1)
    j = njm1
    do i = 3, nim1
        tmult = viscos/yp
        SP(i, j) = SP(i, j) - tmult*SEWU(i)
        SU(i, j) = SU(i, j) + tmult*SEWU(i)*U(i, j+1)
        AN(i, j) = 0.0
    end do
    ! Buttom wall
    yp = Y(2) - YV(2)
    j = 2
    do i = 3, nim1
        tmult = viscos/yp
        SP(i, j) = SP(i, j) - tmult*SEWU(i)
        SU(i, j) = SU(i, j) + tmult*SEWU(i)*U(i, j-1)
        AS(i, j) = 0.0
    end do
    ! Side wall
    do j = 2, njm1
        AE(NIM1, j) = 0.0
        AW(3, j) = 0.0
    end do

    ! Chapter 3, final coef
    resoru = 0.0
    do i = 3, nim1
        do j = 2, njm1
            AP(i, j) = AN(i, j) + AS(i, j) + AE(i, j) + AW(i, j) - SP(i, j)
            DU(i, j) = DU(i, j)/AP(i, j)
            resor    = AN(i, j)*U(i, j+1) + AS(i, j)*U(i, j-1) & 
                     + AE(i, j)*U(i+1, j) + AW(i, j)*U(i-1, j) &
                     - AP(i, j)*U(i, j) + SU(i, j)
            resoru = resoru + abs(resor)
            ! under relaxation
            AP(i, j) = AP(i, j)/urfu
            SU(i, j) = SU(i, j) + (1.-urfu)*AP(i, j)*U(i, j)
            DU(i, j) = DU(i, j)*urfu
        end do
    end do            
    
    ! Chapter 4, Solution of difference equations
    do n = 1, nswpu
        call lisolv(3, 2, ni, nj, U)
    end do      

end subroutine calcu

subroutine calcv()
    use params
    implicit none

    integer :: i, j, n
    real(16) :: cn, cs, ce, cw, dn, ds, de, dw
    real(16) :: xp, tmult, resor

    ! Chapter 1, assembly of coefficients
    do i = 2, nim1
        do j = 3, njm1
            ! convection coefficients
            cn = 0.5*densit*(V(i  , j+1) + V(i  , j  ))*SEW(i)
            cs = 0.5*densit*(V(i  , j  ) + V(i  , j-1))*SEW(i)
            ce = 0.5*densit*(U(i+1, j  ) + U(i+1, j-1))*SNSV(j)
            cw = 0.5*densit*(U(i  , j  ) + U(i  , j-1))*SNSV(j)
            ! diffusion coefficients
            dn = viscos*SEW(i)/DYNPV(j)
            ds = viscos*SEW(i)/DYPSV(j)
            de = viscos*SNSV(j)/DXEP(i)
            dw = viscos*SNSV(j)/DXPW(i)
            ! assemble
            AN(i, j) = AMAX1(ABS(0.5*cn), dn) - 0.5*cn
            AS(i, j) = AMAX1(ABS(0.5*cs), ds) + 0.5*cs
            AE(i, j) = AMAX1(ABS(0.5*ce), de) - 0.5*ce
            AW(i, j) = AMAX1(ABS(0.5*cw), dw) + 0.5*cw
            DV(i, j) = SEW(i)
            SU(i, j) = DV(i, j)*(P(i, j-1) - P(i, j))
            SP(i, j) = 0.0
        end do
    end do

    ! Chapter 2, Problem modifications
    ! Side wall
    xp = X(2) - XU(2)
    i = 2
    do j = 3, njm1
        tmult = viscos/xp
        SP(i, j) = SP(i, j) - tmult*SNSV(j)
        SU(i, j) = SU(i, j) + tmult*SNSV(j)*V(i-1, j)
        AW(i, j) = 0.0
    end do
    ! East wall
    xp = XU(ni) - X(nim1)
    i = nim1
    do j = 3, njm1
        tmult = viscos/xp
        SP(i, j) = SP(i, j) - tmult*SNSV(j)
        SU(i, j) = SU(i, j) + tmult*SNSV(j)*V(i+1, j)
        AE(i, j) = 0.0
    end do
    ! Top wall
    do i = 2, nim1
        AS(i, 3) = 0.0
        AN(i, njm1) = 0.0
    end do

    ! Chapter 3, final coef
    resorv = 0.0
    do i = 2, nim1
        do j = 3, njm1
            AP(i, j) = AN(i, j) + AS(i, j) + AE(i, j) + AW(i, j) - SP(i, j)
            DV(i, j) = DV(i, j)/AP(i, j)
            resor    = AN(i, j)*V(i, j+1) + AS(i, j)*V(i, j-1) & 
                     + AE(i, j)*V(i+1, j) + AW(i, j)*V(i-1, j) &
                     - AP(i, j)*V(i, j) + SU(i, j)
            resorv = resorv + abs(resor)
            ! under relaxation
            AP(i, j) = AP(i, j)/urfv
            SU(i, j) = SU(i, j) + (1.-urfv)*AP(i, j)*V(i, j)
            DV(i, j) = DV(i, j)*urfv
        end do
    end do    
    
    ! Chapter 4, Solution of difference equations
    do n = 1, nswpv
        call lisolv(2, 3, ni, nj, V)
    end do

end subroutine calcv

subroutine calcp()
    use params
    implicit none

    integer :: i, j, n
    real(16) :: cn, cs, ce, cw, smp
    real(16) :: yp, tmult, resor, ppref

    resorm = 0.0

    ! Chapter 1, assembly of coefficients
    do i = 2, nim1
        do j = 2, njm1
            ! calculate coefficients
            AN(i, j) = densit*SEW(i)*DV(i, j+1)
            AS(i, j) = densit*SEW(i)*DV(i, j)
            AE(i, j) = densit*SNS(j)*DU(i+1, j)
            AW(i, j) = densit*SNS(j)*DU(i, j)
            ! source terms
            cn = densit*V(i  , j+1)*SEW(i)
            cs = densit*V(i  , j  )*SEW(i)
            ce = densit*U(i+1, j  )*SNS(j)
            cw = densit*U(i  , j  )*SNS(j)
            smp = cn - cs + ce - cw
            SP(i, j) = 0.0
            SU(i, j) = - smp
            ! comupute sum of absolute mass sources
            resorm = resorm + abs(smp)
        end do
    end do

    ! Chapter 3, final coef
    do i = 2, nim1
        do j = 2, njm1
            AP(i, j) = AN(i, j) + AS(i, j) + AE(i, j) + AW(i, j) - SP(i, j)
        end do
    end do
    
    ! Chapter 4, Solution of difference equations
    do n = 1, nswpp
        call lisolv(2, 2, ni, nj, PP)
    end do

    ! Chapter 5, Correct velocities and pressure
    ! velocities
    do i = 2, nim1
        do j = 2, njm1
            if (i .NE. 2) U(i, j) = U(i, j) + DU(i, j)*(PP(i-1, j) - PP(i, j))
            if (j .NE. 2) V(i, j) = V(i, j) + DV(i, j)*(PP(i, j-1) - PP(i, j))
        end do
    end do
    ! pressures
    ppref = PP(ipref, jpref)
    do i = 2, nim1
        do j = 2, njm1
            P(i, j) = P(i, j) + urfp*(PP(i, j) - ppref)
            PP(i, j) = 0.0
        end do
    end do

end subroutine calcp

subroutine lisolv(istart, jstart, ni, nj, PHI)
    use params, only: nx, ny, AP, AN, AS, AE, AW, SU, SP
    implicit none

    integer, intent(in) :: istart, jstart, ni, nj
    real(16), intent(inout) :: PHI(nx, ny)
    integer :: nim1, njm1, jstm1, i, j, jj
    real(16) :: A(ny), B(ny), C(ny), D(NY), term

    nim1 = ni - 1
    njm1 = nj - 1
    jstm1 = jstart - 1
    A(jstm1) = 0.0
    ! commence W-E sweep
    do i = istart, nim1
        C(jstm1) = PHI(i, jstm1)
        ! commence S-N traverse
        do j = jstart, njm1
            ! assemble TDMA coefficents
            A(j) = AN(i, j)
            B(j) = AS(i, j)
            C(j) = AE(i, j)*PHI(i+1, j) + AW(i, j)*PHI(i-1, j) + SU(i, j)
            D(j) = AP(i, j)
            ! calculate coefficients of recurrence formula
            term = 1./(D(j) - B(j)*A(j-1))
            A(j) = A(j)*term
            C(j) = (C(j) + B(j)*C(j-1))*term
        end do
        ! obtain new PHI's
        do jj = jstart, njm1
            j = nj + jstm1 - jj
            PHI(i, j) = A(j)* PHI(i, j+1) + C(j)
        end do
    end do

end subroutine lisolv