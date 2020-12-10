subroutine grid()
    use params
    implicit none

    integer :: i, j, nig, njg
    real(16) :: dx, dy

    nig = 10
    ni = nig + 2         ! = 12
    nim1 = ni - 1        ! = 11
    dx = xmax/float(nig) ! 1.0/10 = 0.1

    ! X = [-0.05, 0.05, 0.15, ... 1.05]
    X(1) = -0.5*dx
    X(2) = -X(1)
    do i = 3, ni
        X(i) = X(i-1) + dx
    end do

    njg = 10
    nj = njg + 2
    njm1 = nj - 1
    dy = ymax/float(njg)

    ! Y = [-0.05, 0.05, 0.15, ... 1.05]
    Y(1) = -0.5*dy
    Y(2) = -Y(1)
    do j = 3, nj
        Y(j) = Y(j-1) + dy
    end do

    ! print '(12(f6.2))', X
    ! print '(12(f6.2))', Y

end subroutine grid


subroutine init()
    use params
    implicit none
    
    integer :: i, j

    ! Assign dx value
    DXPW(1) = 0.0
    DXEP(ni) = 0.0
    do i = 1, nim1
        DXEP(i) = X(i+1) - X(i)
        DXPW(i+1) = DXEP(i)
    end do

    ! Assign dy value
    DYPS(1) = 0.0
    DYNP(nj) = 0.0
    do j = 1, njm1
        DYNP(j) = Y(j+1) - Y(j)
        DYPS(j+1) = DYNP(j)
    end do
    ! print '(12(f6.2))', DXPW
    ! print '(12(f6.2))', DXEP
    ! print '(12(f6.2))', DYPS
    ! print '(12(f6.2))', DYNP

    ! Averaged value of dx
    SEW(1) = 0.0
    SEW(ni) = 0.0
    do i = 2, nim1
        SEW(i) = 0.5*(DXEP(i) + DXPW(i))
    end do

    ! Averaged value of dy
    SNS(1) = 0.0
    SNS(nj) = 0.0
    do j = 2, njm1
        SNS(j) = 0.5*(DYNP(j) + DYPS(j))
    end do
    ! print '(12(f6.2))', SEW
    ! print '(12(f6.2))', SNS

    ! Posion of midpoint
    XU(1) = 0.0
    do i = 2, ni
        XU(i) = 0.5*(X(i) + X(i-1))
    end do
    ! print '(12(f6.2))', XU

    DXPWU(1) = 0.0
    DXPWU(2) = 0.0
    DXEPU(1) = 0.0
    DXEPU(ni) = 0.0
    do i = 1, nim1
        DXEPU(i) = XU(i+1) - XU(i)
        DXPWU(i+1) = DXEPU(i)
    end do
    ! print '(12(f6.2))', DXEPU
    ! print '(12(f6.2))', DXPWU

    SEWU(1) = 0.0
    do i = 2, ni
        SEWU(i) = (X(i) - X(i-1))
    end do
    ! print '(12(f6.2))', SEWU

    YV(1) = 0.0
    do j=2, nj
        YV(j) = 0.5*(Y(j) + Y(j-1))
    end do
    ! print '(12(f6.2))', YV

    DYPSV(1) = 0.0
    DYPSV(2) = 0.0
    DYNPV(1) = 0.0
    DYNPV(nj) = 0.0
    do j = 2, njm1
        DYNPV(j) = YV(j+1) - YV(j)
        DYPSV(j+1) = DYNPV(j)
    end do
    ! print '(12(f6.2))', DYPSV
    ! print '(12(f6.2))', DYNPV

    SNSV(1) = 0.0
    do j = 2, nj
        SNSV(j) = Y(j) - Y(j-1)
    end do
    ! print '(12(f6.2))', SNSV

    do i = 1, ni
        do j = 1, nj
            U(i, j) = 0.0
            V(i, j) = 0.0
            P(i, j) = 0.0
            PP(i, j) = 0.0
            SU(i, j) = 0.0
            SP(i, j) = 0.0
        end do
    end do

    do i = 1, ni
        U(i, nj) = uwall
    end do

    do i = 1, ni
        do j = 1, nj
            DU(i, j) = 0.0
            DV(i, j) = 0.0
        end do
    end do   

end subroutine init