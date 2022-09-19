SUBROUTINE UfurTest
    use Globals
    use MakeK
    implicit none
    integer i, alfa_number
    complex*16 alfa_min, alfa_step, alfa
    
    alfa_min = cmplx(-2d0)
    alfa_step = cmplx(4d-2)
    alfa_number = 100
    
    
    print*, 'U0 + U- = U+ at z = -h equation test has been started'
    open(1, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\UscTest.txt', FORM='FORMATTED')
    open(2, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\VscTest.txt', FORM='FORMATTED')
    do i = 1, alfa_number
        alfa = alfa_min + alfa_step*i   
        write(1,*) real(alfa),  abs(U0(alfa,-h) + Uminus(alfa,-h)), abs(Uplus(alfa,-h))
        write(2,*) real(alfa),  abs(V0(alfa,-h) + Vminus(alfa,-h)), abs(Vplus(alfa,-h))
    enddo    
    close(1)  
    close(2)
    
    
    print*, 'U0 field test has been started'
    open(1, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\U0TestGr1.txt', FORM='FORMATTED')
    open(2, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\U0TestGr2.txt', FORM='FORMATTED')
    do i = 1, alfa_number
        alfa = alfa_min + alfa_step*i   
        write(1,*) real(alfa),  abs(  derU0(alfa, 0d0)  ), abs( ci*alfa*V0(alfa, 0d0) )
        write(2,*) real(alfa),  abs(  lamda(1)*(-ci*alfa*U0(alfa, 0d0)) + (lamda(1)+2d0*mu(1))*derV0(alfa, 0d0)  ), 1d0 
    enddo    
    close(1)  
    close(2)
    
    print*, 'T test has been started'
    open(1, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\TTestGr1.txt', FORM='FORMATTED')
    open(2, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\TTestGr2.txt', FORM='FORMATTED')
    do i = 1, alfa_number
        alfa = alfa_min + alfa_step*i   
        write(1,*) real(alfa),  abs(  mu(1)*( derU0(alfa, -h) + derUminus(alfa, -h) - ci*alfa*(V0(alfa, -h) + Vminus(alfa, -h)) )  ), abs( mu(2)*( derUplus(alfa, -h) - ci*alfa*Vplus(alfa, -h) ) )
        write(2,*) real(alfa),  abs(  -ci*alfa*lamda(1)*(U0(alfa, -h) + Uminus(alfa, -h)) + (lamda(1)+2d0*mu(1))*(derV0(alfa, -h)+derVminus(alfa, -h))  ),  abs(  -ci*alfa*lamda(2)*(Uplus(alfa, -h)) + (lamda(2)+2d0*mu(2))*(derVplus(alfa, -h))  ) 
    enddo    
    close(1)  
    close(2)
    
    
                                                                        ! Fields equation for STAR5 and Cramer's rule test
    
    print*, 'Uminus by Cramer s rule test has been started'
    open(1, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\CrUminus.txt', FORM='FORMATTED')
    open(2, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\CrVminus.txt', FORM='FORMATTED')
    open(3, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\CrUplus.txt', FORM='FORMATTED')
    open(4, file='C:\Users\tiama\OneDrive\Рабочий стол\IMMI\DIVA2\data\CrVplus.txt', FORM='FORMATTED')
    do i = 1, alfa_number
        alfa = alfa_min + alfa_step*i   
        write(1,*) real(alfa),  imag(  Uminus(alfa, 0d0)  ), imag( CrUminus(alfa, 0d0) )
        write(2,*) real(alfa),  imag(  Vminus(alfa, -2d0)  ), imag( CrVminus(alfa, -2d0) )
        write(3,*) real(alfa),  imag(  Uplus(alfa, -6d0-h)  ), imag( CrUplus(alfa, -6d0-h) )
        write(4,*) real(alfa),  imag(  Vplus(alfa, -6d0-h)  ), imag( CrVplus(alfa, -6d0-h) )
    enddo    
    close(1)  
    close(2)
    close(3)
    close(4)
    
END SUBROUTINE UfurTest