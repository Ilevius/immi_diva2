SUBROUTINE uMinusP
    use Globals;
    use MakeK
    implicit none;
    real*8 psi, psiStep, Rmin, Rstep, Rmax
    real*8, allocatable:: psi2h(:), psih(:)
    real*8 t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    integer dotsNumber, i, j, psiNumber
    real*8, allocatable:: x(:), x2h(:), xh(:)
    real*8, allocatable:: z(:), z2h(:), zh(:)
    real*8, allocatable:: R(:)
    real*8, allocatable:: R2h(:), Rh(:)
    complex*16, allocatable:: field_int(:), field_sp(:), field_sp1(:), field_sp2(:)
    real*8 currentR, currentPsi

    
    namelist /rayInput/ psi, psiNumber, psiStep, Rmin, Rstep, Rmax, t1, t2, t3, t4, tm, tp, eps, step, IntLimit
    open(unit=1,file='input.txt',status='old')
        read(1, rayInput)
    close(1)
     

    dotsNumber = floor((Rmax-Rmin)/Rstep)+1;
    allocate(R(dotsNumber));  allocate(x(dotsNumber)); allocate(z(dotsNumber)); 
    allocate(R2h(dotsNumber)); allocate(x2h(dotsNumber)); allocate(z2h(dotsNumber)); allocate(psi2h(dotsNumber));
    allocate(Rh(dotsNumber)); allocate(xh(dotsNumber)); allocate(zh(dotsNumber)); allocate(psih(dotsNumber));
    allocate(field_int(dotsNumber)); allocate(field_sp(dotsNumber)); allocate(field_sp1(dotsNumber)); allocate(field_sp2(dotsNumber));
    
    open(1001, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\ReDots.txt', FORM='FORMATTED')
    open(1002, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\ReDots1.txt', FORM='FORMATTED')
    open(1, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\dots.txt', FORM='FORMATTED')
    open(2, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_real.txt', FORM='FORMATTED')
    open(3, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_real.txt', FORM='FORMATTED')   
    open(4, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_imag.txt', FORM='FORMATTED')
    open(5, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_imag.txt', FORM='FORMATTED') 
    open(6, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\integral_abs.txt', FORM='FORMATTED')
    open(7, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\asymptotics_abs.txt', FORM='FORMATTED') 
    
    t1 = Kappa(1)*0.1; t2 = t1; t3 = t1; t4 = (Kappa(2)*1.4d0+1d0);
        
        
        
    do j = 1, psiNumber
        psi = psi + psiStep
        do i = 1, dotsNumber
            R(i) = Rmin + (i-1)*Rstep;
            x(i) = R(i)*cosd(psi)
            z(i) = R(i)*sind(psi) -h
              
            R2h(i) = sqrt(x(i)**2+(z(i)+2d0*h)**2)
            psi2h(i) = atan2(z(i)+2d0*h,x(i))
            Rh(i) = sqrt(x(i)**2+(z(i)+h)**2)
            psih(i) = atan2(z(i)+h,x(i))
        enddo 
        
        call bipolarTest
        call dinn5(UminusPInt,t1,t2,t3,t4,tm,tp,eps,step,IntLimit,dotsNumber,field_int)
        field_int = field_int/(2d0*pi);
        
        call uMinusPSp1
        call uMinusPSp2
        
        field_sp = field_sp1 + field_sp2
         
        do i = 1, dotsNumber
            write(1,*) x(i),  z(i), psi
            
            write(2,*) R(i),  real(field_int(i)), psi 
            write(4,*) R(i),  imag(field_int(i)), psi
            write(6,*) R(i),  abs(field_int(i)), psi
            
            write(3,*) R(i),  real(field_sp(i)), Rh(i) 
            write(5,*) R(i),  imag(field_sp(i)), Rh(i)
            write(7,*) R(i),  abs(field_sp(i)), Rh(i)
        enddo
    enddo    
        
        
        
        
        
        
    close(1001); close(1002); close(1); close(2); close(3); close(4); close(5); close(6); close(7);
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    CONTAINS
        
        SUBROUTINE UminusPInt(alfa, s, n)
        implicit none;
        integer n
        complex*16 alfa, s(n), sigma(2)
            do i = 1, dotsNumber
                sigma = MakeSigma(alfa)
                s(i) = 1d0/(delta(alfa)*CramDelta(0, 0, alfa))*(CramDelta(1, 1, alfa)*exp(-sigma(1)*h) + CramDelta(1, 2, alfa)*exp(-sigma(2)*h))*exp(-sigma(1)*(z(i)+h) - ci*alfa*x(i))
                s(i) = s(i) + 1d0/(delta(-alfa)*CramDelta(0, 0, -alfa))*(CramDelta(1, 1, -alfa)*exp(-sigma(1)*h) + CramDelta(1, 2, -alfa)*exp(-sigma(2)*h))*exp(-sigma(1)*(z(i)+h) + ci*alfa*x(i))
            enddo 
        END SUBROUTINE UminusPInt
        
        
        SUBROUTINE UminusPInt1(alfa, s, n)
        implicit none;
        integer n
        complex*16 alfa, s(n), sigma(2)
            do i = 1, dotsNumber
                sigma = MakeSigma(alfa)
                s(i) = 1d0/(delta(alfa)*CramDelta(0, 0, alfa))*CramDelta(1, 1, alfa)*exp(-sigma(1)*h)*exp(-sigma(1)*(z(i)+h) - ci*alfa*x(i))
                s(i) = s(i) + 1d0/(delta(-alfa)*CramDelta(0, 0, -alfa))*CramDelta(1, 1, -alfa)*exp(-sigma(1)*h)*exp(-sigma(1)*(z(i)+h) + ci*alfa*x(i))
            enddo 
        END SUBROUTINE UminusPInt1
        
        
        SUBROUTINE UminusPInt2(alfa, s, n)
        implicit none;
        integer n
        complex*16 alfa, s(n), sigma(2)
            do i = 1, dotsNumber
                sigma = MakeSigma(alfa)
                s(i) = 1d0/(delta(alfa)*CramDelta(0, 0, alfa))*CramDelta(1, 2, alfa)*exp(-sigma(2)*h)*exp(-sigma(1)*(z(i)+h) - ci*alfa*x(i)) + 1d0/(delta(-alfa)*CramDelta(0, 0, -alfa))*CramDelta(1, 2, -alfa)*exp(-sigma(2)*h)*exp(-sigma(1)*(z(i)+h) + ci*alfa*x(i))
            enddo 
        END SUBROUTINE UminusPInt2
        
        
        
        SUBROUTINE uMinusPSp1
        IMPLICIT NONE;
        integer jj
        complex*16 alfa0         
            do jj = 1, dotsNumber
                alfa0 = cmplx(-Kappa(1)*cos(psi2h(jj)))
                field_sp1(jj) = sqrt(Kappa(1)*sin(psi2h(jj))**2/(2d0*pi*R2h(jj)))*1d0/(delta(alfa0)*CramDelta(0, 0, alfa0))*CramDelta(1, 1, alfa0)*exp((ci*R2h(jj)*Kappa(1)-ci*pi/4d0))
            enddo       
        END SUBROUTINE uMinusPSp1
        
                
        SUBROUTINE uMinusPSp2
        IMPLICIT NONE;
        integer jj, p, stPointsNumber 
        real*8 stPoints(10), alfa0, theta, thetaDerDerSign, thetaDerDer
        complex*16 alfa0c
            do jj = 1, dotsNumber
                currentPsi = psih(jj) 
                currentR = Rh(jj)
                call Halfc(uMinusPSp2ThetaDer, -100d0, 100d0, 2d-2, 1d-8, 10, stPoints, stPointsNumber)
                field_sp2(jj) = 0d0
                do p = 1, stPointsNumber
                    alfa0 =   stPoints(p)
                    alfa0c = cmplx(alfa0)
                    theta = uMinusPSp2Theta(alfa0)
                    thetaDerDer = uMinusPSp2ThetaDerDer(alfa0)
                    thetaDerDerSign = abs(uMinusPSp2ThetaDerDer(alfa0))/uMinusPSp2ThetaDerDer(alfa0)
                    field_sp2(jj) =  sqrt(1d0/(2d0*pi*currentR))*sqrt(1d0/abs(thetaDerDer))*1d0/(delta(alfa0c)*CramDelta(0, 0, alfa0c ))*CramDelta(1, 2, alfa0c )*exp( ci*currentR*theta + ci*pi/4d0*thetaDerDerSign )
                enddo
            enddo      
        END SUBROUTINE uMinusPSp2
        
        
        FUNCTION uMinusPSp2Theta(alfa)
            implicit none
            real*8 alfa
            real*8  uMinusPSp2Theta   
                uMinusPSp2Theta = sqrt( kappa(1)**2 - alfa**2 )*sin(currentPsi) + sqrt(kappa(2)**2 - alfa**2)*h/currentR - alfa*cos(currentPsi)
        END FUNCTION uMinusPSp2Theta
        

        FUNCTION uMinusPSp2ThetaDer(alfa)
        implicit none
        complex*16 uMinusPSp2ThetaDer
        real*8 alfa, temp
            temp = -alfa*sin(currentPsi)/(sqrt( kappa(1)**2 - alfa**2 )) - alfa*h/(currentR*sqrt( kappa(2)**2 - alfa**2 )) - cos(currentPsi)
            uMinusPSp2ThetaDer = cmplx(temp)
        END FUNCTION uMinusPSp2ThetaDer
       

        FUNCTION uMinusPSp2ThetaDerDer(alfa)
        implicit none
        real*8 alfa, uMinusPSp2ThetaDerDer
            uMinusPSp2ThetaDerDer = -sin(currentPsi)*kappa(1)**2/((sqrt( kappa(1)**2 - alfa**2 ))**3) - h*kappa(2)**2/(((sqrt( kappa(2)**2 - alfa**2 ))**3)*currentR)
        END FUNCTION uMinusPSp2ThetaDerDer
        
     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        SUBROUTINE  bipolarTest
            implicit none
            integer i

            do i = 1, dotsNumber
                write(1001,*) R2h(i)*cos(psi2h(i)), R2h(i)*sin(psi2h(i))-2d0*h
                write(1002,*) Rh(i)*cos(psih(i)), Rh(i)*sin(psih(i))-h
             enddo   
        END SUBROUTINE  bipolarTest
    
    
END SUBROUTINE uMinusP