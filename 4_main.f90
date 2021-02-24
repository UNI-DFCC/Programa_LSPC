
! Main Program - Solucion Numerica para 'r' puntos y 'm' modos

program main
    use constants
    implicit none
    integer :: i,j
    integer :: r,n,m,k,pa
    real(dp) :: radio,h,error,w,paso,rd,th,p_th,temperatura
    complex(dp),dimension(:,:),allocatable :: coef
    real(dp),dimension(:,:),allocatable :: A
    real(dp),dimension(:),allocatable :: sol_aux1,sol_aux2,B
    

    write(*,*) '----------SOLUCION DE LA ECUACION DE CALOR----------'
    write(*,*) ' '
    write(*,*) 'Ingrese el numero de puntos en la parte radial:'
    read(*,*) r
    write(*,*) 'Ingrese el numero maximo de modos para resolver:'
    read(*,*) m
    
    radio=1.0 ; error=(10.0**(-6)) ; w=1.5 ; pa=100 
    p_th=2.0*pi/real(pa-1)
    
    n = r+1
    k = m+1
    h = radio/real(n)

    allocate(A(n*2,n*2))
    allocate(B(n*2))
    allocate(coef(k,n))
    allocate(sol_aux1(n*2))
    allocate(sol_aux2(pa))

    
    call m_zeros(A,n*2,n*2) ! Se llena de ceros la matriz
    call v_zeros(B,n*2)   ! Se llena de ceros el vector
    call m_zeros_complex(coef,k,n)
    !call m_print(A,n*2,n*2) ! Se imprime la matriz
    !call v_print(B,n*2)   ! Se imprime el vector
    !call m_print_complex(coef,k,n)
    
    do i=0,m,1
        !print'(A,3/,A)', ' ',' '
        !write(*,*) '---------------------------------------------------'
        !write(*,*) 'Para k=',i
        !write(*,*) ''
        call matrix_problem(A,B,n,i,h)
        !write(*,*) 'Matriz A'
        !call m_print(A,n*2,n*2)
        !write(*,*) 'Vector B'
        !call v_print(B,n*2)
        call sor(A,B,n*2,error,w,.true.,sol_aux1)
        !write(*,*) 'solucion'
        !call v_print(sol_aux1,n*2)
        do j=1,n,1
            coef(i+1,j)=complex(sol_aux1(j),sol_aux1(j+n))
        end do
    end do
    !print'(A,1/,A)', ' ',' '
    !write(*,*) 'matriz de coeficientes'
    !call m_print_complex(coef,k,n)
    
    open(20,file='file1.dat',status='replace')
    open(30,file='file2.dat',status='replace')
    open(40,file='file3.dat',status='replace')

    write(40,*) '#  Matriz de temperaturas T[r,th] ---> T(th,r)'
    do i=1,n,1
        rd= h*real(i-1)
        write(20,*) rd
        do j=0,pa-1,1
            th=p_th*real(j)
            if (i==1) then
                write(30,*) th
            end if
            sol_aux2(j+1) = temperatura(coef,k,n,i,th)
        end do
        write(40,*) sol_aux2
    end do
    close(20)
    close(30)
    close(40)
    
end program main



!--------------------------------------------------------------------
!--------------------------------------------------------------------

