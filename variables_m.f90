module params
    ! Length of variable array
    integer, parameter :: nx = 12
    integer, parameter :: ny = 12

    ! module all
    ! Unknown
    integer ni, nj, nim1, njm1
    
    ! module geom
    ! Geometry related parameters
    ! DX means dx, EP and PW means positions 
    ! S means What?
    ! NS means North-South?
    ! EW means East-West?
    ! XU is midpoint in X positions
    real(16) :: X(nx), Y(ny), DXEP(nx), DXPW(nx), DYNP(ny), DYPS(ny)
    real(16) :: SNS(ny), SEW(nx), XU(nx), YV(ny)
    
    ! module uvel
    ! X-direction velocity = U
    ! defined at xu (midpoint of x)
    ! DXEPU means dx between E-P, DXPWU means dx between P-W 
    ! S means small?
    ! SEWU is dx at e and u defined in X (NOT XU)
    integer :: nswpu = 2
    real(16) :: urfu = 0.5
    real(16) :: resoru, DXEPU(nx), DXPWU(nx), SEWU(nx)

    ! module vvel
    ! Y-direction velocity = V
    ! nswpv means the number fo iteration?
    integer :: nswpv = 2
    real(16) :: urfv = 0.5
    real(16) :: resorv, DYNPV(ny), DYPSV(ny), SNSV(ny)

    ! module pcor
    ! Unknown
    integer :: nswpp = 3
    integer :: ipref = 2
    integer :: jpref = 2
    real(16) :: urfp = 0.8
    real(16) :: resorm, DU(nx, ny), DV(nx, ny)

    ! module var
    ! Array of the variables which should be solved
    ! PP is what?
    real(16) :: U(nx, ny), V(nx, ny), P(nx, ny), PP(nx, ny)

    ! module fluid
    ! Parameter which is constant in the field
    real(16) :: viscos = 1.0e-3
    real(16) :: densit = 1000.0

    ! module coef
    ! Direction (N-S-E-W) related coefficients
    ! A means TDMA coefficients like something
    ! S means source term
    ! SU related to Pressure gradient
    real(16) :: AP(nx, ny), AN(nx, ny), AS(nx, ny), AE(nx, ny), AW(nx, ny)
    real(16) :: SU(nx, ny), SP(nx, ny)

    ! module cavsiz
    ! Boundary values
    real(16) :: xmax = 1.0
    real(16) :: ymax = 1.0
    real(16) :: uwall = 1.0e-4
end module params
