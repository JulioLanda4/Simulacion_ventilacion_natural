PROGRAM Flujo_en_un_canal_obs_sol
IMPLICIT NONE

! Declaración de variables
CHARACTER*50 name, itchar
INTEGER nx, ny, i, j, max_iter, itmax, it, maxiter_S, ei, ej, counter
REAL*4 x0, xl, y0, yl, dx, dy, dv, Se, So, Pe, gamma, Sn, Ss, Sw, tolerance, tolerance_S, residual, time, dt
REAL*4 ue,uw,vn,vs, Div, Re, dbkx, dbky 
REAL*4, ALLOCATABLE :: T(:), aP(:,:), aN(:,:), aS(:,:), aE(:,:), aW(:,:), sP(:,:)
REAL*4, ALLOCATABLE :: x(:), xc(:), y(:), yc(:), phi(:,:), phia(:,:), uc(:,:), vc(:,:)
REAL*4, ALLOCATABLE :: u1(:,:), v1(:,:), P(:,:), PP(:,:), de(:,:),dn(:,:), u(:,:), v(:,:)

REAL*4, ALLOCATABLE :: x0_obs, y0_obs, xl_obs, yl_obs, u_obs, v_obs
INTEGER, ALLOCATABLE :: mark_cells(:,:)

!u1, v1 tiempo presente y u, v tiempo anterior


! Inicializar variables
nx = 100; ny = 20
x0 = 0.0; xl = 20.0
y0 = -2.5; yl = 2.5
Re=100.0

time=200
dt=0.001
itmax=int(time/dt)+1

max_iter = 1000
maxiter_s=100
tolerance = 1e-05
tolerance_S = 1e-5

dx = (xl - x0)/float(nx)
dy = (yl - y0)/float(ny)

!las caras que tiene cuando es en una dimension es cuando ponemos igual a uno
Se = dy; Sw = dy
Sn = dx; Ss = dx

dv = dx*dy

gamma=1.0/Re



! Alojamiento de la memoria dinámica: Para la matriz
ALLOCATE(aP(nx,ny), aE(nx,ny), aW(nx,ny), aN(nx,ny), aS(nx,ny), sP(nx,ny), uc(0:nx+1,0:ny+1), vc(0:nx+1,0:ny+1))
ALLOCATE(x(0:nx), xc(0:nx+1), y(0:ny), yc(0:ny+1), p(0:nx+1, 0:ny+1), pp(0:nx+1, 0:ny+1))
ALLOCATE(u(0:nx, 0:ny+1), u1(0:nx, 0:ny+1), v(0:nx+1, 0:ny), v1(0:nx+1, 0:ny), de(0:nx+1, 0:ny+1), dn(0:nx+1, 0:ny+1))
ALLOCATE(mark_cells(0:nx+1,0:ny+1))

! Creación de la malla
Call MESH_1D(nx,x0,xl,x,xc)
Call MESH_1D(ny,y0,yl,y,yc)

! Condiciones de frontera
u=1.0; v=0.0; P=0.0; Pp=0.0
de= 0.0; dn=0.0


! Definiendo dimensiones del obstaculo
x0_obs=4.0; xl_obs=5.0; y0_obs=-0.5; yl_obs=0.5
mark_cells = 0

u(0,:)=1.0 !flujo tapon en la cara oeste

u(:,0) = 0.0
u(:,ny+1) = 0.0


u1=u; v1=v 

OPEN(3, FILE='Obstaculo.txt',STATUS='REPLACE')
! Ciclo para marcar el obstaculo
DO i=1, nx
    DO j=1, ny
        if((xc(i) .gt. x0_obs .and. xc(i) .lt. xl_obs) .and. (yc(j) .gt. y0_obs .and. yc(j) .lt. yl_obs)) then
            mark_cells(i,j)=1; u(i,j)=0.0; v(i,j)=0.0
            WRITE(3,*)xc(i),yc(j)
        end if
    END DO
END DO
CLOSE(3)


DO i=1, nx
    DO j=1, ny
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .eq. 1) then
            u(i,j) = 0.0;
        end if
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .eq. 1) then
            v(i,j) = 0.0;
        end if
    END DO
END DO

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!Creacion de arcivho de animacion 
OPEN(2,FILE='anim.gnp',STATUS='REPLACE')
WRITE(2,*)'set xrange[0:20]; set yrange[-2.5:2.5]; set view map; set size ratio 0.25'
!WRITE(2,*)"set xlabel 'x'; set ylabel 'y' rotate by 0"


!////////////////////////////////////////////////////
!AQUI VAMOS 02/05/2023
!///////////////////////////////////////////////////

!Iniciar ciclo temporal 
do it=1, itmax
Div= 1.0
counter=1

!INICIA CICLO SIMPLEC
DO WHILE ((counter .LT. maxiter_s) .and. (Div .GT. tolerance))


! Algoritmo de solución
! Cálculo de los coeficientes

! Ecuacion de la componente U
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0

ei=ubound(u,1)-1 
ej=ubound(u,2)-1

do i=1, ei
    do j=1, ej
        ue=0.5*(u(i,j)+u(i+1,j))
        uw=0.5*(u(i,j)+u(i-1,j))
        vn=0.5*(v(i,j)+v(i+1,j))
        vs=0.5*(v(i,j-1)+v(i+1,j-1))

        ! aE(i,j) = gamma*Se/dx -MIN(ue*se,0.0)
        ! aW(i,j) = gamma*Sw/dx +MAX(uw*sw,0.0)
        ! aN(i,j) = gamma*Sn/dy -MIN(vn*sn,0.0)
        ! aS(i,j) = gamma*Ss/dy +MAX(vs*sn,0.0)

        aE(i,j) = gamma*Se/dx - 0.5*ue*se !MIN(ue*se,0.0)
        aW(i,j) = gamma*Sw/dx + 0.5*uw*sw !+MAX(uw*sw,0.0)
        aN(i,j) = gamma*Sn/dy - 0.5*vn*sn !-MIN(vn*sn,0.0)
        aS(i,j) = gamma*Ss/dy + 0.5*vs*ss !+MAX(vs*sn,0.0)

        aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j) + dv/dt
        sP(i,j) = u(i,j) * dv/dt - (P(i+1,j)-p(i,j)) * dv/dx

        ! Correciones de frontera para el obstaculo //////////////////////////////////////////////////////////////////
        ! Estas correciones son respecto al flujo y no al obstaculo, es decir, las fronteras estan invertidas
        ! La malla esta desplazada medio volumen de control a la derecha

        !   Para el vector U -------------------------------------------------------------------------------

        ! Este
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+2,j) .eq. 1) then
        sp(i,j)=sp(i,j)+ae(i,j)*u1(i+1,j)
        ae(i,j)=0.0
        end if

        ! Oeste
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i-1,j) .eq. 1) then
        sp(i,j)=sp(i,j)+aw(i,j)*u1(i-1,j)
        aw(i,j)=0.0
        end if

        ! Norte
        if(mark_cells(i,j) .eq. 0 .and. (mark_cells(i+1,j+1) .eq. 1 .or. mark_cells(i,j+1) .eq. 1)) then
        ap(i,j)=ap(i,j)+an(i,j)
        sp(i,j)=sp(i,j)+2.0*an(i,j)*u1(i,j+1)
        an(i,j)=0.0
        end if

        ! Sur
        if(mark_cells(i,j) .eq. 0 .and. (mark_cells(i+1,j-1) .eq. 1 .or. mark_cells(i,j-1) .eq. 1)) then
        ap(i,j)=ap(i,j)+as(i,j)
        sp(i,j)=sp(i,j)+2.0*as(i,j)*u1(i,j-1)
        as(i,j)=0.0
        end if

        ! Obstaculo
        if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .eq. 1)) then
        aw(i,j)=0.0; ae(i,j)=0.0; as(i,j)=0.0; an(i,j)=0.0; u1(i,j)=0.0
        ap(i,j)=1.0 ! también se puede hacer ap(i,j)=1.0 con sp(i,j)=u1(i,j)     puede tomar cualquier valor diferente de cero
        sp(i,j)=0.0   ! sp(i,j) = 0.0
        end if
        ! ---------------------------------------------------------------------------------------------

        ! ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    end do
end do



! Corrección por condiciones de frontera
dbkx=0.0; dbky=1.0



! Oeste
aP(1, 1:ej) = aP(1,1:ej) + dbkx * aW(1,1:ej)
sP(1, 1:ej) = sP(1,1:ej) + (1.0+dbkx)*aW(1,1:ej)*u(0,1:ej)
aW(1, 1:ej) = 0.0

! Este
aP(ei,1:ej) = aP(ei,1:ej) - dbkx * aE(ei,1:ej)
sP(ei,1:ej) = sP(ei,1:ej)  !+ (1.0+dbkx)* aE(nx,:)*u(ei+1,1:ej)
aE(ei,1:ej) = 0.0

! Norte
aP(1:ei,ej) = aP(1:ei,ej)+ dbky*aN(1:ei,ej)
sP(1:ei,ej) = sP(1:ei,ej)+ (1.0+dbky)* aN(1:ei, ej)*u(1:ei, ej+1)
aN(1:ei,ej)=0.0

! Sur
aP(1:ei,1) = aP(1:ei,1)+ dbky* aS(1:ei,1)
sP(1:ei,1) = sP(1:ei,1)+ (1.0+dbky) *aS(1:ei,1)*u(1:ei,0)
aN(1:ei,1)=0.0

de(1:ei, 1:ej)=Se/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))



!Método de solución
Call Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance, residual)
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0


!Componente v
ei=ubound(v,1)-1 
ej=ubound(v,2)-1

do i=1, ei
    do j=1, ej
        ue=0.5*(u(i,j)+u(i,j+1))
        uw=0.5*(u(i-1,j)+u(i-1,j+1))
        vn=0.5*(v(i,j)+v(i,j+1))
        vs=0.5*(v(i,j)+v(i,j-1))

        ! aE(i,j) = gamma*Se/dx -MIN(ue*se,0.0)
        ! aW(i,j) = gamma*Sw/dx +MAX(uw*sw,0.0)
        ! aN(i,j) = gamma*Sn/dy -MIN(vn*sn,0.0)
        ! aS(i,j) = gamma*Ss/dy +MAX(vs*sn,0.0)

        aE(i,j) = gamma*Se/dx - 0.5*ue*se !MIN(ue*se,0.0)
        aW(i,j) = gamma*Sw/dx + 0.5*uw*sw !+MAX(uw*sw,0.0)
        aN(i,j) = gamma*Sn/dy - 0.5*vn*sn !-MIN(vn*sn,0.0)
        aS(i,j) = gamma*Ss/dy + 0.5*vs*ss !+MAX(vs*sn,0.0)

        aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j) + dv/dt
        sP(i,j) = v(i,j) * dv/dt - (P(i,j+1)-p(i,j)) * dv/dy

        ! Correcion de las fronteras del obstaculo en la componente v

        ! Este
        if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j+1) .eq. 1 .or. mark_cells(i+1,j) .eq. 1)) then
        ap(i,j)=ap(i,j)+ae(i,j)
        sp(i,j)=sp(i,j)+2.0*ae(i,j)*v1(i+1,j)
        ae(i,j)=0.0
        end if

        ! Oeste
        if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i-1,j+1) .eq. 1 .or. mark_cells(i-1,j) .eq. 1)) then
        ap(i,j)=ap(i,j)+aw(i,j)
        sp(i,j)=sp(i,j)+2.0*aw(i,j)*v1(i-1,j)
        aw(i,j)=0.0
        end if!! !

        ! Norte
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+2) .eq. 1) then
        sp(i,j)=sp(i,j)+an(i,j)*v1(i,j+1)
        an(i,j)=0.0
        end if

        ! Sur
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j-1) .eq. 1) then
        sp(i,j)=sp(i,j)+as(i,j)*v1(i,j-1)
        as(i,j)=0.0
        end if

        ! Obstaculo
        if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .eq. 1)) then
        aw(i,j)=0.0; ae(i,j)=0.0; as(i,j)=0.0; an(i,j)=0.0; v1(i,j)=0.0
        ap(i,j)=1.0
        sp(i,j)=0.0
        end if

    end do
end do

!write(*,*) aP
!stop 

! Corrección por condiciones de frontera
dbkx=1.0; dbky=0.0

! Oeste
aP(1, 1:ej) = aP(1,1:ej) + dbkx * aW(1,1:ej)
sP(1, 1:ej) = sP(1,1:ej) + (1.0+dbkx)*aW(1,1:ej)*v(0,1:ej)
aW(1, 1:ej) = 0.0

! Este
aP(ei,1:ej) = aP(ei,1:ej) - dbkx * aE(ei,1:ej)
sP(ei,1:ej) = sP(ei,1:ej) !+ (1.0+dbkx)* aE(nx,:)*v(ei+1,1:ej)
aE(ei,1:ej) = 0.0

! Norte
aP(1:ei,ej) = aP(1:ei,ej)+ dbky*aN(1:ei,ej)
sP(1:ei,ej) = sP(1:ei,ej)+ (1.0+dbky)* aN(1:ei, ej)*v(1:ei, ej+1)
aN(1:ei,ej)=0.0

! Sur
aP(1:ei,1) = aP(1:ei,1)+ dbky* aS(1:ei,1)
sP(1:ei,1) = sP(1:ei,1)+ (1.0+dbky) *aS(1:ei,1)*v(1:ei,0)
aN(1:ei,1)=0.0

dn(1:ei, 1:ej)=Sn/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))



!Método de solución
Call Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance, residual)
! P=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; sP=0.0

aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0; Pp=0.0

! Correcciones de la presion
ei=ubound(Pp,1)-1
ej=ubound(Pp,2)-1

do i=1, ei
    do j=1, ej
        ue=u1(i,j)
        uw=u1(i-1,j)
        vn=v1(i,j)
        vs=v1(i,j-1)

        aE(i,j) = de(i,j)*Se
        aW(i,j) = de(i-1,j)*Sw
        aN(i,j) = dn(i,j)*sN
        aS(i,j) = dn(i,j-1)*Ss

        aP(i,j) = aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j)
        sP(i,j) = -ue*Se+uw*Sw-vn*Sn+vs*Ss !divergencia de u discretizada
    end do
end do

Call Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance, residual)
!

!PASO 4 CORRECCION DE LAS PRESIONES
P=P+Pp

!PASO 5, CORREGIR LAS VELOCIDADES
!COMPONRENTE U 
DO i=1, nx-1
    DO j=1, ny
        if ((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j) .ne. 1)) then
            u1(i,j)=u1(i,j)+de(i,j)*(Pp(i,j)-Pp(i+1,j))
        end if
END DO
END DO

!COMPONRENTE V
DO i=1, nx
    DO j=1, ny-1
        if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .ne. 1)then
            v1(i,j)=v1(i,j)+dn(i,j)*(Pp(i,j)-Pp(i,j+1))
        end if
END DO
END DO

!::::PASO 5.5:::::
! !Vel en frontera Newman(este)
u1(nx,:)=u1(nx-1,:)
v1(nx+1,:) = v1(nx,:)

!PASO 6. ITERAR HASTA LLEGAR A LA SOLUCION . VERIFICAR LA DIVERGENCIA DE U
Div=MAXVAL(ABS(Sp(1:ei,1:ej)))

counter=counter+1


!Para que nos imprima los resultados cada 100 iteraciones

! if (MOD(it,100) .EQ. 0) then 
! call WriteScalarField2D('temp',it,phi,xc,yc,nx,ny) 
! WRITE(itchar,'(i6)')it !i6 es el numero de caracteres que tendremos, itchar es la variable que creara  caracteres
! itchar=ADJUSTL(itchar) !Cadena de caracteres a la izquierda y no queden espacios vacios 
! name='temp'//itchar(1:LEN_TRIM(itchar))//'.txt'
! WRITE(2,*)'sp "'// name(1:LEN_TRIM(name))//'" w pm3d'
! WRITE(2,*)'pause 0.1' !ir viendo la evolucion temporal del problema

end do !Termina ciclo SIMPLEC

u=u1; v=v1
WRITE (*,*)it, counter, Div 

! Interpolar velocidades
CALL interpolateToNodesVs(vc,v,nx,ny)
CALL interpolateToNodesUs(uc,u,nx,ny)

! Para la animacion
if (mod(it,1000) .eq. 0) then



    ! Escribir campo de valocidades
    CALL WriteVectorField('vel',it,uc,vc,xc,yc,nx,ny)

    write(itchar, '(i6)')it
    itchar = adjustl(itchar)


    name='vel'//itchar(1:LEN_TRIM(itchar))//'.txt'
    WRITE(2,*)'sp "'// name(1:LEN_TRIM(name))//'" w vec lc -1'
    WRITE(2,*)'pause 0.1' !ir viendo la evolucion temporal del problema

!     itchar = itchar(1:len_trim(itchar))
!     write(3,*)"p 'vel"//itchar(1:len_trim(itchar))//".txt' u 1:2:(0.2*$3):(0.2*$4) ev 1:3 w vec lc -1"
!     write(3,*)'pause 0.05'
end if

! if (MOD(it,100) .EQ. 0) then
! call WriteScalarField2D('vel',it,phi,xc,yc,nx,ny)
! WRITE(itchar,'(i6)')it !i6 es el numero de caracteres que tendremos, itchar es la variable que creara  caracteres
! itchar=ADJUSTL(itchar) !Cadena de caracteres a la izquierda y no queden espacios vacios
! name='temp'//itchar(1:LEN_TRIM(itchar))//'.txt'
! WRITE(2,*)'sp "'// name(1:LEN_TRIM(name))//'" w pm3d'
! WRITE(2,*)'pause 0.1' !ir viendo la evolucion temporal del problema
! end if



END DO !Terminamos el ciclo temporal 
close(2)


CALL WriteVectorField('vel',0,uc,vc,xc,yc,nx,ny)
end program 

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER nx,i
REAL*4 dx,x0,xl
REAL*4 x(0:nx),xc(0:nx+1)

dx=1.0/FLOAT(nx)

DO i=0,nx
  x(i)=x0+FLOAT(i)*(xl-x0)*dx
END DO

xc(0)=x(0); xc(nx+1)=x(nx)

DO i=1,nx
  xc(i)=0.5*(x(i)+x(i-1))
END DO

END SUBROUTINE


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine TDMA(x,a,b,c,d,n)

implicit none
 integer n,k
 real a(n),b(n),c(n),d(n),x(n),m

 do k=2,N
  m=a(k)/b(k-1)
  b(k)=b(k)-m*c(k-1)
  d(k)=d(k)-m*d(k-1)
 end do

 x(n)=d(n)/b(n)

 do k=n-1,1,-1
  x(k)=(d(k)-c(k)*x(k+1))/b(k)
 end do

end subroutine

!***************************************************************

subroutine lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ei),b(ei),c(ei),d(ei)
bi=1; bj=1
do j=bj,ej
do i=bi,ei
    a(i)=-aW(i,j)
    b(i)=aP(i,j)
    c(i)=-aE(i,j)
    d(i)=sp(i,j) + aN(i,j) * phi(i,j+1) + aS(i,j) * phi(i,j-1)
end do
 call TDMA(phi(bi:ei,j), a, b, c ,d ,ei)
end do

end subroutine

!****************************************************************

subroutine lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ej),b(ej),c(ej),d(ej)
bi=1; bj=1

do i=bi,ei
do j=bj,ej
    a(j)=-aS(i,j)
    b(j)=aP(i,j)
    c(j)=-aN(i,j)
    d(j)=sp(i,j) + aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) 
end do
 call TDMA(phi(i,bj:ej), a, b, c ,d ,ej)
end do

end subroutine

!****************************
!****************************

real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 acum(ei,ej),residual,NINV
bi=1; bj=1
acum=0
NINV = 1.0 / dfloat(ei*ej)
do i=bi,ei
do j=bj,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function

!************************************************************************************************************

subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer bi,ei,bj,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance
    
    interface
        real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
        implicit none
        integer bi,ei,bj,ej,i,j,nx,ny
        real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
        end function
    end interface

count_iter=0;  residual=1.0

do while((count_iter <= max_iter).and.(residual > tolerance)) 
call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do

end subroutine


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
Subroutine WriteScalarField2D(Name,kx,T,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
    do j=0,ny+1
    do i=0,nx+1
    write(10,*)xc(i),yc(j),T(i,j)
    end do
    write(10,*)''
    end do
close(10)
End Subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine interpolateToNodesUs(uc,us,nx,ny)

    implicit none
    integer nx,ny
    real:: us(0:nx,0:ny+1),uc(0:nx+1,0:ny+1)
    integer bi,ei,bj,ej,i,j
    bi=1; bj=1;ej=ny; ei=nx
    
    ! Internal points
    do i=bi,ei
    do j=bj-1,ej+1
        uc(I,J) = ( us(I-1,J) + us(I,J) ) * 0.5
    end do 
    end do
    
    uc(bi-1,:)=us(0,:)
    uc(ei+1,:)=us(nx,:)
    
    end subroutine
    
    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    subroutine interpolateToNodesVs(vc,vs,nx,ny)
    implicit none
    integer nx,ny
    real:: vs(0:nx+1,0:ny),vc(0:nx+1,0:ny+1)
    integer bi,ei,bj,ej,i,j
    bi=1; bj=1;ej=ny; ei=nx
    
    do i=bi-1,ei+1
    do j=bj,ej
        vc(I,J) = ( vs(I,J) + vs(I,J-1) ) * 0.5
    end do 
    end do
    
    vc(:,bj-1)=vs(:,0)
    vc(:,ej+1)=vs(:,ny)
    end subroutine
    
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
    integer i,j,nx,ny,kx
    real*4 uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
    character*(*)Name
    character*50 txt,Filename
    write(txt,'(i6)')kx
    txt=ADJUSTL(txt)
    Filename=name//txt(1:len_trim(txt))//".txt"
    
    open(11,file=Filename(1:len_trim(Filename)))
        do i=0,nx+1
        do j=0,ny+1
        write(11,*)xc(i),yc(j),uc(i,j),vc(i,j)
        end do
        write(11,*)''
        end do
    close(11)
    End Subroutine
