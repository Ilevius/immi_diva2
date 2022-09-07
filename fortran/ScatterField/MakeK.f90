MODULE MakeK
    use Globals;
    implicit none;
    
    CONTAINS      
            FUNCTION MakeSigma(alfa)
            implicit none;
            complex*16 alfa, MakeSigma(2)

                if (Imag(alfa) == 0d0)  then 
                    if (abs(alfa) < Kappa(1)) then
                        MakeSigma(1) = -ci*sqrt(Kappa(1)**2 - real(alfa)**2)
                    else 
                        MakeSigma(1) = sqrt(real(alfa)**2 - Kappa(1)**2)
                    endif    
                    
                    if (abs(alfa) < Kappa(2)) then
                        MakeSigma(2) = -ci*sqrt(Kappa(2)**2 - real(alfa)**2)
                    else 
                        MakeSigma(2) = sqrt(real(alfa)**2 - Kappa(2)**2)
                    endif  
                    
                else 
                    MakeSigma(1) = sqrt(alfa**2 - cmplx(Kappa(1)**2))
                    MakeSigma(2) = sqrt(alfa**2 - cmplx(Kappa(2)**2))
                endif              
            END FUNCTION MakeSigma
            
        
            FUNCTION MakeSigma_(alfa)
            implicit none;
            complex*16 alfa, MakeSigma_(2)

                if (Imag(alfa) == 0d0)  then 
                    if (abs(alfa) < Kappa_(1)) then
                        MakeSigma_(1) = -ci*sqrt(Kappa_(1)**2 - real(alfa)**2)
                    else 
                        MakeSigma_(1) = sqrt(real(alfa)**2 - Kappa_(1)**2)
                    endif    
                    
                    if (abs(alfa) < Kappa_(2)) then
                        MakeSigma_(2) = -ci*sqrt(Kappa_(2)**2 - real(alfa)**2)
                    else 
                        MakeSigma_(2) = sqrt(real(alfa)**2 - Kappa_(2)**2)
                    endif  
                    
                else 
                    MakeSigma_(1) = sqrt(alfa**2 - cmplx(Kappa_(1)**2))
                    MakeSigma_(2) = sqrt(alfa**2 - cmplx(Kappa_(2)**2))
                endif              
            END FUNCTION MakeSigma_
    
        
        FUNCTION Delta(alfa)
            implicit none;
            complex*16 alfa, Delta, Sigma(2)
                Sigma = MakeSigma(alfa)
                Delta = 2d0*mu(1)*(-(alfa**2-0.5d0*Kappa(2)**2)**2 + alfa**2*Sigma(1)*Sigma(2));    
        END FUNCTION Delta
        
        
        FUNCTION U0(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, U0, Sigma(2)
            real*8 z
                Sigma = MakeSigma(alfa)
                U0 = ci*alfa*(alfa**2 - 0.5d0*(Kappa(2))**2)/Delta(alfa)*exp(sigma(1)*z)-ci*alfa*Sigma(1)*Sigma(2)/Delta(alfa)*exp(sigma(2)*z)    
        END FUNCTION U0
        
        
        FUNCTION V0(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, V0, Sigma(2)
            real*8 z
                Sigma = MakeSigma(alfa);
                V0 = -Sigma(1)*(alfa**2 - 0.5d0*Kappa(2)**2)/Delta(alfa)*exp(sigma(1)*z) + Sigma(1)*alfa**2/Delta(alfa)*exp(sigma(2)*z)
        END FUNCTION V0
        
        
        FUNCTION derU0(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derU0, Sigma(2)
            real*8 z
                Sigma = MakeSigma(alfa)
                derU0 = ci*alfa*(alfa**2 - 0.5d0*(Kappa(2))**2)*sigma(1)/Delta(alfa)*exp(sigma(1)*z)-ci*alfa*Sigma(1)*Sigma(2)**2/Delta(alfa)*exp(sigma(2)*z)    
        END FUNCTION derU0
        
        
        FUNCTION derV0(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derV0, Sigma(2)
            real*8 z
                Sigma = MakeSigma(alfa);
                derV0 = -Sigma(1)**2*(alfa**2 - 0.5d0*Kappa(2)**2)/Delta(alfa)*exp(sigma(1)*z) + Sigma(1)*sigma(2)*alfa**2/Delta(alfa)*exp(sigma(2)*z)
        END FUNCTION derV0
        
        
        FUNCTION matrixA(alfa)
            IMPLICIT NONE;
            complex*16 alfa, matrixA(4,4), Sigma(2), Sigma_(2)
                sigma = MakeSigma(alfa); sigma_ =  MakeSigma_(alfa);   
            
                matrixA(1,1) = (-1d0,0d0);  matrixA(1,2) = sigma(2);  matrixA(1,3) = (1d0,0d0);  matrixA(1,4) = sigma_(2);
                matrixA(2,1) = sigma(1);  matrixA(2,2) =-alfa**2 ;  matrixA(2,3) = sigma_(1);  matrixA(2,4) = alfa**2;  
                matrixA(3,1) = mu(1)*(ci*alfa*sigma(1)-sigma(1));  matrixA(3,2) = mu(1)*(sigma(2)**2 - ci*alfa**3);  
                matrixA(3,3) = mu(2)*(ci*alfa*sigma_(1)-sigma_(1));  matrixA(3,4) = -mu(2)*(sigma_(2)**2 - ci*alfa**3); 
                matrixA(4,1) = -ci*alfa*lamda(1) + (lamda(1) + 2d0*mu(1))*sigma(1)**2;  matrixA(4,2) = sigma(2)*alfa*(ci*lamda(1) - (lamda(1)+2d0*mu(1))*alfa);  
                matrixA(4,3) = ci*alfa*lamda(2) - (lamda(2) + 2d0*mu(2))*sigma_(1)**2;  matrixA(4,4) = sigma_(2)*alfa*(ci*lamda(2) - (lamda(2)+2d0*mu(2))*alfa);  
        END FUNCTION matrixA
        
        
        FUNCTION matrixB(alfa)
            IMPLICIT NONE;
            complex*16 alfa, matrixB(4,1), Sigma(2), Sigma_(2)
                sigma = MakeSigma(alfa); sigma_ =  MakeSigma_(alfa);
                
                matrixB(1,1) = U0(alfa, -h)
                matrixB(2,1) = V0(alfa, -h)
                matrixB(3,1) = ci*alfa*mu(1)*V0(alfa, -h) - mu(1)*derU0(alfa, -h)
                matrixB(4,1) = ci*lamda(1)*alfa*U0(alfa, -h) - (lamda(1)+2*mu(1))*derV0(alfa, -h)
        END FUNCTION matrixB                  
        
        
        FUNCTION Uminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, Uminus, Sigma(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma = MakeSigma(alfa);
                Uminus = t(1,1)*exp(-sigma(1)*(z+h)) - t(2,1)*sigma(2)*exp(-sigma(2)*(z+h))
        END FUNCTION Uminus
        
        
        FUNCTION Vminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, Vminus, Sigma(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma = MakeSigma(alfa);
                Vminus = -t(1,1)*sigma(1)*exp(-sigma(1)*(z+h)) + t(2,1)*alfa**2*exp(-sigma(2)*(z+h))
        END FUNCTION Vminus
        
        
        FUNCTION Uplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, Uplus, Sigma_(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma_ = MakeSigma_(alfa);
                Uplus = t(3,1)*exp(sigma_(1)*(z+h)) + t(4,1)*sigma_(2)*exp(sigma_(2)*(z+h))
        END FUNCTION Uplus
        
        
        FUNCTION Vplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, Vplus, Sigma_(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma_ = MakeSigma_(alfa);
                Vplus = t(3,1)*sigma_(1)*exp(sigma_(1)*(z+h)) + t(4,1)*alfa**2*exp(sigma_(2)*(z+h))
        END FUNCTION Vplus
        
   !                                                                                              Usc D E R I V A T I V E S
        FUNCTION derUminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derUminus, Sigma(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma = MakeSigma(alfa);
                derUminus = -sigma(1)*t(1,1)*exp(-sigma(1)*(z+h)) + sigma(2)*t(2,1)*sigma(2)*exp(-sigma(2)*(z+h))
        END FUNCTION derUminus
        
        
        FUNCTION derVminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derVminus, Sigma(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma = MakeSigma(alfa);
                derVminus = t(1,1)*sigma(1)**2*exp(-sigma(1)*(z+h)) - sigma(2)*t(2,1)*alfa**2*exp(-sigma(2)*(z+h))
        END FUNCTION derVminus
        
        
        FUNCTION derUplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derUplus, Sigma_(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma_ = MakeSigma_(alfa);
                derUplus = sigma_(1)*t(3,1)*exp(sigma_(1)*(z+h)) + sigma_(2)*t(4,1)*sigma_(2)*exp(sigma_(2)*(z+h))
        END FUNCTION derUplus
        
        
        FUNCTION derVplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, derVplus, Sigma_(2), t(4,1)
            real*8 z
            complex*16 A(4,4), C(4,4), R(4,7)         
                A = matrixA(alfa)
                t = matrixB(alfa)
                call STAR5(A,t,C,R,4,4,1,2)
                Sigma_ = MakeSigma_(alfa);
                derVplus = t(3,1)*sigma_(1)**2*exp(sigma_(1)*(z+h)) + sigma_(2)*t(4,1)*alfa**2*exp(sigma_(2)*(z+h))
        END FUNCTION derVplus
    
        
        
                        !                                       Cramer's rule
        
        FUNCTION partB(i, alfa)
            IMPLICIT NONE;
            integer i
            complex*16 alfa, partB(4,1), Sigma(2), Sigma_(2)
                sigma = MakeSigma(alfa); sigma_ =  MakeSigma_(alfa);     
                if (i == 1) then     
                    partB(1,1) = ci*alfa*(alfa**2 - 0.5d0*kappa(2)**2)
                    partB(2,1) = -(alfa**2 - 0.5d0*kappa(2)**2)*Sigma(1)
                    partB(3,1) = -2d0*ci*alfa*mu(1)*sigma(1)*(alfa**2 - 0.5d0*kappa(2)**2)
                    partB(4,1) = ((lamda(1)+2d0*mu(1))*sigma(1)**2 - lamda(1)*alfa**2)*(alfa**2 - 0.5d0*kappa(2)**2)
                else 
                    partB(1,1) = -ci*alfa*sigma(1)*sigma(2)
                    partB(2,1) = alfa**2*sigma(1)
                    partB(3,1) = ci*mu(1)*sigma(1)*alfa**3 + ci*alfa*mu(1)*sigma(1)*sigma(2)**2
                    partB(4,1) = -2d0*mu(1)*sigma(1)*sigma(2)*alfa**2
                endif
        END FUNCTION partB  
        
        
        FUNCTION CramDelta(i, j, alfa)
        implicit none
        complex*16 CramDelta, alfa, matrix(4,4), column(4,1)
        complex*16 C(4,4), S(4), SM(4)
        integer i, j     
            matrix = matrixA(alfa)         
            if (i>0 .AND. j>0) then 
                column = partB(j, alfa)
                matrix(:,i) = column(:,1)
            endif          
            call DSTAR(matrix,C,S,SM,CramDelta,4,4,2)      
        END FUNCTION CramDelta
        
        
        FUNCTION CrUminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, CrUminus, Sigma(2), t1, t2
            real*8 z
                Sigma = MakeSigma(alfa);
                t1 = CramDelta(1, 1, alfa)*exp(-sigma(1)*h) + CramDelta(1, 2, alfa)*exp(-sigma(2)*h)
                t2 = CramDelta(2, 1, alfa)*exp(-sigma(1)*h) + CramDelta(2, 2, alfa)*exp(-sigma(2)*h)
                CrUminus = t1*exp(-sigma(1)*(z+h)) - t2*sigma(2)*exp(-sigma(2)*(z+h))
                CrUminus = CrUminus/(delta(alfa)*CramDelta(0, 0, alfa))
        END FUNCTION CrUminus
        
        FUNCTION CrVminus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, CrVminus, Sigma(2), t1, t2
            real*8 z
                Sigma = MakeSigma(alfa);
                t1 = CramDelta(1, 1, alfa)*exp(-sigma(1)*h) + CramDelta(1, 2, alfa)*exp(-sigma(2)*h)
                t2 = CramDelta(2, 1, alfa)*exp(-sigma(1)*h) + CramDelta(2, 2, alfa)*exp(-sigma(2)*h)
                CrVminus = -sigma(1)*t1*exp(-sigma(1)*(z+h)) + t2*alfa**2*exp(-sigma(2)*(z+h))
                CrVminus = CrVminus/(delta(alfa)*CramDelta(0, 0, alfa))
        END FUNCTION CrVminus
        
        
        FUNCTION CrUplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, CrUplus, Sigma(2), Sigma_(2), t3, t4
            real*8 z
                Sigma = MakeSigma(alfa); Sigma_ = MakeSigma_(alfa);
                t3 = CramDelta(3, 1, alfa)*exp(-sigma(1)*h) + CramDelta(3, 2, alfa)*exp(-sigma(2)*h)
                t4 = CramDelta(4, 1, alfa)*exp(-sigma(1)*h) + CramDelta(4, 2, alfa)*exp(-sigma(2)*h)
                CrUplus = t3*exp(sigma_(1)*(z+h)) + t4*sigma_(2)*exp(sigma_(2)*(z+h))
                CrUplus = CrUplus/(delta(alfa)*CramDelta(0, 0, alfa))
        END FUNCTION CrUplus
        
        
        FUNCTION CrVplus(alfa, z)
            IMPLICIT NONE;
            complex*16 alfa, CrVplus, Sigma(2), Sigma_(2), t3, t4
            real*8 z
                Sigma = MakeSigma(alfa); Sigma_ = MakeSigma_(alfa);
                t3 = CramDelta(3, 1, alfa)*exp(-sigma(1)*h) + CramDelta(3, 2, alfa)*exp(-sigma(2)*h)
                t4 = CramDelta(4, 1, alfa)*exp(-sigma(1)*h) + CramDelta(4, 2, alfa)*exp(-sigma(2)*h)
                CrVplus = t3*sigma_(1)*exp(sigma_(1)*(z+h)) + t4*alfa**2*exp(sigma_(2)*(z+h))
                CrVplus = CrVplus/(delta(alfa)*CramDelta(0, 0, alfa))
        END FUNCTION CrVplus
        
        
        !                                                           Stationary point calculations
        


        
        
END MODULE MakeK    