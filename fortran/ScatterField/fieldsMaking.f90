SUBROUTINE fieldsMaking
    use Globals
    use MakeK
    implicit none
        print*, 'Field makers has been started'
        !call star5test
    
    
    
    CONTAINS
        SUBROUTINE star5test       
        implicit none
        complex*16 A(2,2), B(2,1), C(2,2), R(2,4)
        
        
        
            A(1,1) = (1d0,0d0); A(1,2) = (1d0,0d0);
            A(2,1) = (2d0,0d0); A(2,2) = (-1d0,0d0);
            
            B(1,1) = (8,0d0); B(2,1) = (-11d0,0d0);
        
            call STAR5(A,B,C,R,2,2,1,2)
            
            print*, B
            
            pause
        
            !    STAR5(A,B,C,R,ND,N,M,kort) is a solver for a complex linear algebraic system Ax=b 
            !       by the double-orthogonalization method 
            ! N - system size; M - number of the right-hand parts;
            ! Nd - dimension (number of the lines) declared in fact; 
            ! A(Nd,N) - N*N complex matrix of the system; 
            ! B(Nd,M) - N*M complex matrix of the right-hand parts;
            ! C(Nd,N),R(Nd,M+3) - complex work arrays
            ! KORT - number of the ortogonalization cycles 
            !(2 - 3 are usually quite enough)
        

        
        
        
        END SUBROUTINE star5test
    
    
    
END SUBROUTINE fieldsMaking