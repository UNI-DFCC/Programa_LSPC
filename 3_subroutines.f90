! REAL SUBROUTINES --------------------------------------------------

subroutine regresion_lineal(vx,vy,n,a0,a1,r)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: suma,sum_v,sum_vu,sum_v2,media,Sr,St
    real(dp) :: funcion_regresion,sr_valor,st_valor
    real(dp) :: a0,a1,r,den,num
    real(dp),dimension(n),intent(in) :: vx,vy
    real(dp),dimension(n) :: y

    call v_zeros(y,n)
    
    num=(real(n)*sum_vu(vx,vy))-(sum_v(vx)*sum_v(vy))
    den=(real(n)*sum_v2(vx))-(sum_v(vx)**2)
    a1 = num/den
    a0 = media(vy)-(a1*media(vx))
    do i=1,n,1
        y(i) = funcion_regresion(a0,a1,vx(i))
    end do
    sr_valor = Sr(vy,y)
    st_valor = St(vy)
    r = ((st_valor-sr_valor)/st_valor)**0.5

    return
end subroutine

!--------------------------------------------------------------------

subroutine matrix_problem(A,B,n,k,h)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n,k
    real(dp),intent(in) :: h
    real(dp),dimension(n*2,n*2),intent(inout) :: A
    real(dp),dimension(n*2),intent(inout) :: B
    real(dp) :: alpha,aux1,aux2,aux3
    complex(dp) :: aux4,contorno_rieman
    
    alpha=0.0
    if (mod(k,4)==0) then
        alpha=4.0
    end if
    
    A(1,1) = 4.0
    A(1+n,1+n) = 4.0
    A(1,2) = (-1.0)*alpha
    A(1+n,2+n) = (-1.0)*alpha
    do i=2,n,1
        aux1 = 1.0/(h**2)
        aux2 = 1.0/(2.0*(real(i-1)*h)*h)
        aux3 = (real(k)**2)/((real(i-1)*h)**2)
        A(i,i-1) = (aux1-aux2)
        A(i+n,i+n-1) = (aux1-aux2)
        A(i,i) = ((-2.0*aux1)-aux3)
        A(i+n,i+n) = ((-2.0*aux1)-aux3)
        if (i/=n) then
            A(i,i+1) = (aux1+aux2)
            A(i+n,i+n+1) = (aux1+aux2)
        end if
        if (i==n) then
            aux4 = contorno_rieman(k)
            B(n) = (-1.0)*(aux1+aux2)*real(aux4)
            B(2*n) = (-1.0)*(aux1+aux2)*aimag(aux4)
        end if
    end do

    return
end subroutine

!--------------------------------------------------------------------

subroutine sor(A,B,n,error,w,logic,sol)
    use constants
    implicit none
    integer,intent(in) :: n
    integer :: i,j,ite
    real(dp),intent(in) :: error,  w
    real(dp) :: delta,cuad,suma1,suma2,suma3
    real(dp),dimension(n,n),intent(in) :: A
    real(dp),dimension(n),intent(in) :: B
    real(dp),dimension(n),intent(inout) :: sol
    real(dp),dimension(n) :: aux,aux1,dif2
    logical,intent(in) :: logic
    ite=0 ; delta=100.0

    call v_zeros(sol,n)
    call v_zeros(aux,n)
    call v_zeros(aux1,n)
    call v_zeros(dif2,n)

    do while (delta>=error)
        aux1=sol
        do i=1,n,1
            suma1 = 0.0
            do j=1,i-1,1
                suma1 = suma1 + (A(i,j)*sol(j))
            end do
            suma2 = 0.0
            do j=i+1,n,1
                suma2 = suma2 + (A(i,j)*sol(j))
            end do
            suma3 = 0.0
            do j=1,i-1,1
                suma3 = suma3 + (A(i,j)*aux(j))
            end do
            aux(i)=(B(i)-((1.0-w)*suma1)-suma2-(w*suma3))/A(i,i)
        end do
        sol=aux
        cuad = 0.0
        dif2=(sol-aux1)
        do i=1,n,1
            cuad = cuad + (dif2(i)**2)
        end do
        delta=(cuad**(0.5))
    end do

    return
end subroutine

!--------------------------------------------------------------------

subroutine m_zeros(A,n,m)
    use constants
    implicit none
    integer,intent(in) :: n,m
    real(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = 0.0
        end do
    end do
    
    return
end subroutine

!---------------------------------------------------------------------

subroutine v_zeros(A,n)
    use constants
    implicit none
    integer,intent(in) :: n
    real(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = 0.0
    end do
    
    return
end subroutine

!---------------------------------------------------------------------

subroutine m_print(A,n,m)
    use constants
    implicit none
    integer,intent(in) :: n,m
    real(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j
    
    write(*,*) '---------------------------------------------------'
    do i=1,n,1
        write(*,*) A(i,:)
    end do
    write(*,*) '---------------------------------------------------'

    return
end subroutine

!---------------------------------------------------------------------

subroutine v_print(A,n)
    use constants
    implicit none
    integer,intent(in) :: n
    real(dp),dimension(n),intent(inout) :: A

    write(*,*) '---------------------------------------------------'
    write(*,*) A
    write(*,*) '---------------------------------------------------'

    return
end subroutine




!---------------------------------------------------------------------
!---------------------------------------------------------------------

! COMPLEX SUBROUTINES -----------------------------------------------

subroutine m_zeros_complex(A,n,m)
    use constants
    implicit none
    integer,intent(in) :: n,m
    complex(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = (0.0,0.0)
        end do
    end do
    
    return
end subroutine

!---------------------------------------------------------------------

subroutine v_zeros_complex(A,n)
    use constants
    implicit none
    integer,intent(in) :: n
    complex(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = (0.0,0.0)
    end do
    
    return
end subroutine

!---------------------------------------------------------------------

subroutine m_print_complex(A,n,m)
    use constants
    implicit none
    integer,intent(in) :: n,m
    complex(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j
    
    write(*,*) '---------------------------------------------------'
    do i=1,n,1
        write(*,*) A(i,:)
    end do
    write(*,*) '---------------------------------------------------'

    return
end subroutine

!---------------------------------------------------------------------

subroutine v_print_complex(A,n)
    use constants
    implicit none
    integer,intent(in) :: n
    complex(dp),dimension(n),intent(inout) :: A

    write(*,*) '---------------------------------------------------'
    write(*,*) A
    write(*,*) '---------------------------------------------------'

    return
end subroutine

!---------------------------------------------------------------------
