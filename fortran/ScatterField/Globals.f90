MODULE Globals
    IMPLICIT NONE;
    real*8 a, h, f, theta, cp(2), cs(2), rho(2), lamda(2), mu(2), pi, k
    real*8 Kappa(2), Kappa_(2),omega
    complex*16 ci
    parameter (pi=3.141592653589793d0)
    parameter (ci = (0d0,1d0))
    namelist /media/ cp, cs, rho, theta, a, f, h, omega
    CONTAINS  
        SUBROUTINE InitGlobals
            open(unit=1,file='input.txt',status='old')
            read(1, media);
            close(1);
            mu(1) = rho(1)*cs(1)**2; mu(2) = rho(2)*cs(2)**2; 
            lamda(1) = rho(1)*(cp(1)**2 - 2d0*cs(1)**2); lamda(2) = rho(2)*(cp(2)**2 - 2d0*cs(2)**2);
            !k = 2*pi*f/cp(1);
            !Kappa(1) = 2d0*pi*f/cp(1); Kappa(2) = 2d0*pi*f/cs(1);
!            Kappa_(1) = 2d0*pi*f/cp(2); Kappa_(2) = 2d0*pi*f/cs(2);
             k = omega/cp(1);
            Kappa(1) = omega/cp(1); Kappa(2) = omega/cs(1);
            Kappa_(1) = omega/cp(2); Kappa_(2) = omega/cs(2);
            
            print*, 'Globals have been inited!', '  mu=', mu,  ' lamda = ', lamda, kappa  
        END SUBROUTINE InitGlobals
END  MODULE Globals   