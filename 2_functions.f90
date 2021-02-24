! REAL FUNCTIONS -------------------------------------------------

function sol_analitica(r,th) result(aux)
    use constants
    implicit none
    integer :: i
    real(dp),intent(in) :: r,th
    real(dp) :: k,aux,aux1,aux2,aux3,aux4,suma

    aux1 = (-2.0)*pi
    aux2 = r*((-1.0)*pi*cos(th))+((3.0+(8.0*(pi**2)))*sin(th)/6.0)
    suma=0.0
    do i=2,100,1
        k = real(i)
        aux3 = (4.0*pi)*cos(k*th)/((k**2)-1.0)
        aux4 = (8.0*k)*sin(k*th)/((1.0-(k**2))**2)
        suma = suma + ((r**i)*(aux3+aux4))
    end do
    aux = aux1+aux2+suma

    return
end function
!####################################################################
function funcion_contorno(x) result(aux)
    use constants
    implicit none
    real(dp),intent(in) :: x
    real(dp) :: aux

    aux = (1.0+(x**2))*sin(x)
    !aux = sin(x)

    return
end function
!####################################################################
function funcion_regresion(a0,a1,x) result(aux)
    use constants
    implicit none
    real(dp),intent(in) :: a0,a1,x
    real(dp) :: aux

    aux = a0+(a1*x)

    return
end function
!####################################################################
function sum_vu(vx,vy,n) result(suma)
    use constants
    implicit none
    integer,intent(in) :: n
    integer :: i
    real(dp) :: suma
    real(dp),dimension(n),intent(in) :: vx,vy

    suma=0.0
    do i=1,n,1
        suma = suma + (vx(i)*vy(i))
    end do

    return
end function
!####################################################################
function sum_v(v,n) result(suma)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: suma
    real(dp),dimension(n),intent(in) :: v
    
    suma=0.0
    do i=1,n,1
        suma = suma + (v(i))
    end do

    return
end function
!####################################################################
function media(v,n) result(aux)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: aux,sum_v
    real(dp),dimension(n),intent(in) :: v
    
    aux = sum_v(v,n) / real(n)
    
    return
end function
!####################################################################
function sum_v2(v,n) result(suma)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: suma
    real(dp),dimension(n),intent(in) :: v
    
    suma=0.0
    do i=1,n,1
        suma = suma + (v(i)**2)
    end do

    return
end function
!####################################################################
function St(v,n) result(suma)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: suma,media,aux
    real(dp),dimension(n),intent(in) :: v
    
    suma=0.0
    aux=media(v,n)
    do i=1,n,1
        suma = suma + ((v(i)-aux)**2)
    end do

    return
end function
!####################################################################
function Sr(v,u,n) result(suma)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: n
    real(dp) :: suma
    real(dp),dimension(n),intent(in) :: v,u
    
    suma=0.0
    do i=1,n,1
        suma = suma + ((v(i)-u(i))**2)
    end do

    return
end function
!####################################################################



!####################################################################
!--------------------------------------------------------------------
! COMPLEX FUNCTIONS -------------------------------------------------

function contorno_rieman(k) result(suma)
    use constants
    implicit none
    integer,intent(in) :: k
    integer :: n,j
    real(dp) :: x,fun,ex_re,ex_im,funcion_contorno
    complex(dp) :: ex,suma

    n=1000
    suma=(0.0,0.0)
    do j=0,n,1
        x = 2.0*pi*real(j)/real(n)
        ex_re=cos(real(k)*x)
        ex_im=(-1.0)*sin(real(k)*x)
        ex= complex(ex_re,ex_im)
        fun = funcion_contorno(x)
        suma = suma + (ex*(fun/real(n)))
    end do

    return
end function
!####################################################################

function temperatura(coef,k,n,p,th) result(t)
    use constants
    implicit none
    integer :: i
    integer,intent(in) :: k,n,p
    real(dp),intent(in) :: th
    real(dp) :: t,cz,sz
    complex(dp),dimension(k,n),intent(in) :: coef
    complex(dp) :: aux,z,C
    
    aux=coef(1,p)
    do i=2,k,1
        C = coef(i,p)
        cz = cos(real(i-1)*th)
        sz = sin(real(i-1)*th)
        z = complex(cz,sz)
        aux=aux+(C*z)+(conjg(C)*conjg(z))
    end do
    t= real(aux)

    return
end function

!####################################################################


