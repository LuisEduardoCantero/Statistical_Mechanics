!!
!!		    Autor:			      Fecha:
!! 
!!	Miguel Ángel Sandoval Puentes		   27-ENERO-2018
!! 
!!	Como parte del trabajo para la obtención del grado de maestria en ciencias físicas. 
!! Se usó como base un programa proporcionado por el Dr. José Miguel Méndez Alcaraz, asesor.

Module Var_Prog
	Implicit None
!!
!!	Archivos de salida
!!
	character*20 din	!! Datos del programa
	character*20 gout,eout,sout	!! funciones g(r)
!!
!!	Arreglos
!!
	Real*8, Allocatable, Dimension (:) :: x, y, z	!! Posiciones
	Real*8, Allocatable, Dimension (:,:) :: cfgx, cfgy, cfgz	!! Configuraciones
	Real*8, Allocatable, Dimension (:) :: g,s,q			!! funciones g(r),s(q)

!!
!!	Variables de la simulación
!!
   
	Real*8 :: pi	!! Pi=3.141569...
	Real*8 :: rho,phi	!! densidad reducida, fracción de llenado 
	Real*8 :: bl,rc,qmax,dd,ad	!! ancho de caja, radio de corte, q máximo, distancia promedio, parámetro de escalamiento
	Real*8 :: dr,dq		!! Incrementos en g(r) y s(q)
	Real*8 :: x0, y0, z0	!! Posiciones base
	Real*8 :: xx, yy, zz	!! Posiciones nuevas

	Real*8 :: escR,escS,escD,escE	!! Escalas
	Integer*4 :: nct,ncp,ncet,ncep	!! ciclos e incrementos en ciclos
	Integer*4 :: iseed,nss,nt,nr,nq	!! semilla y # de elementos de arreglos
	Integer*4 :: np,nprom

	Real*8 :: drmax			!! Desplazamiento máximo para MC
	Real*8 :: rza			!! Razón de aceptanción de configuraciones
	Real*8 :: V,V0,Vn,dV,Vp		!! Energía de la configuración
	Integer*4 :: ncfa		!! Número de confiduraciones aceptadas
	Integer*4 :: nintr		!! Número de interacciones
	Integer*4 :: idr		!! # ciclos para actualizar drmax
	End Module Var_Prog

	Module Var_Pot
	Implicit None
	Real*8 :: rmin,rmin2,rmax	!! Contactos de potencial
	Real*8 :: eps1,eps2,H,H2,diam	!! Intensidades de pozo y barrera, alcances de potencial, diámetro

	End Module Var_Pot

	Program DB_FSS
	Use Var_Prog
	Use Var_Pot
	Implicit None
	Real*8 :: ini,fin,m1,m2	!! tiempo de ejecución
	Real*8 :: densSo	!! Densidad numérica
	Real*8 :: ran2,uij,uu,rij	!! Funciones
	Integer*4 :: i,j,k
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Inicio de tiempo de simulación
!!
	call cpu_time(ini)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Variables de entrada	
!!
	pi=4.d0*datan(1.d0)
	phi=5.0d-2		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! factor de llenado Luis
	diam=1.0d-6
	densSo=6*phi/(pi*(diam**3.d0))			!! densidad total 
	dd=(1.d0/densSo)**(1.d0/3)		!! distancia prometio general
	ad=diam			!! Adimensionalización
	NSS=5
	iseed=-50
	rho=1.d0
	bl=7.d0!!!!!!!!!!!!!!!!!!!!!!!!!!VeranoUG
	rc=bl/2.d0
	np=idint(bl**3)
	escR=dd/ad
	escS=ad/dd
	escD=(ad/dd)**2.d0
	escE=1.d0
	nct=200		!!!!!!! MOVER PARA PRECTICAR 200000
 	ncp=nct+10	!!!!!!! MOVER PARA PRACTICAR 1000000
	ncet=10
	ncep=10
	qmax=25.d0
	dr=0.01d0
	dq=0.05d0
	nr=idint(rc/dr)
	nq=idint(qmax/dq)
	nt=(ncp-nct)/ncep
	idr=ncp/4
	drmax=0.1*(diam/dd)
	ncfa=0
	nintr=0
	nprom=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Dimensiones
!!
	Allocate(x(np),y(np),z(np))
	Allocate(cfgx(nt,np),cfgy(nt,np),cfgz(nt,np))
	Allocate(g(nr))
	Allocate(q(nq),s(nq))
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Variables de potencial
!!	
	H=0.25d0*diam/dd!!!!!!!!!!!!!VeranoUG
	H2=0.03d0*diam/dd!!!!!!!!!!!!!!VeranoUG
	rmin=diam/dd+H

	eps1=-1.73091d0
	eps2=0.0

	rmin2=diam/dd+H+H2
	rmax=rmin2+0.01d0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Vamos a verificar los potenciales
!!
open(1,file='pot.dat',status='unknown')
Do i=1,10*nr
	rij=i*0.1d0*dr
	uu=uij(rij)
	write(1,*)EscR*rij,uu
enddo
 close(1)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!		Ficheros de salida
!! 
	din='datos.dat' 
	gout='g.dat'
	eout='e.dat'
	sout='s.dat'
	open(10,file=gout,status='unknown')
	open(13,file=eout,status='unknown') 
	open(1111,file=sout,status='unknown')
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Esta sección guarda los datos usados para generar cada paquete de resultados
!!
	open(14,file=din,status='unknown')
	write(14,*)'Semilla=',iseed
	write(14,*)'Número de partículas: np=',np
	write(14,*)'Pasos en g(r): nr=',nr
	write(14,*)'Configuraciones guardadas: nt=',nt
	write(14,*)'Diametro: diam=',diam
	write(14,*)'Distanca pormedio: dd=',dd
	write(14,*)'Factor de llenado: phi=',phi
	write(14,*)'Número de configuraciones: ncp=',ncp
	write(14,*)'Configuraciones de termalizacion: nct=',nct
	write(14,*)'Cada cuando analizar configuración: ncep=',ncep

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Inicio de simulación
!!
	open(1,file='conf_fin.dat',status='unknown')
	open(2,file='configuraciones.dat',status='unknown')
	call iniconfigC()	!! Configuración inicial
	call energia()		!! Energía de configuración inicial
	Do i=1,ncp		!! Configuraciones
		Do j=1,Np	!! Partículas
			x0=x(j)
			y0=y(j)	
			z0=z(j)
			call energiaU0(j)	!! Energía de confguración inicial
			xx=x0+(2.d0*ran2(iseed)-1.d0)*drmax
			yy=y0+(2.d0*ran2(iseed)-1.d0)*drmax	!! Movimiento aleatorio
			zz=z0+(2.d0*ran2(iseed)-1.d0)*drmax
			xx=xx-bl*dnint(xx/bl)
			yy=yy-bl*dnint(yy/bl)	!! Convencion de imagen mínima
			zz=zz-bl*dnint(zz/bl)
			call energiaU(j)	!! Energía de configuración final
			dV=Vn-V0
			if(dV.lt.100.d0)then
			if(dV.le.0.d0.or.((Dexp(-dV)).gt.ran2(iseed)))then	!!Condición para aceptar configuración
				V=V+dV
				x(j)=xx
				y(j)=yy
				z(j)=zz
				ncfa=ncfa+1	!! Número de configuraciones aceptadas
			endif
			endif
			nintr=nintr+1		!! Número de interacciones
			Vp=V/dfloat(Np)	!! Energía por partícula
		if(i.eq.ncp)then
			write(1,*)x(j),y(j),z(j)
		endif
		End do
		if(mod(i,idr).eq.0)then
		rza=ncfa/dfloat(Np*idr)
		if(rza.gt.0.5d0)then
		drmax=drmax*1.05d0
		else
		drmax=drmax*0.95d0
		endif
		ncfa=0
		endif
		if((mod(i,ncep).eq.0.d0).and.(i.gt.nct))then
		nprom=nprom+1
		Do k=1,np
			cfgx(nprom,k)=x(k)
			cfgy(nprom,k)=y(k)
			cfgz(nprom,k)=z(k)
!			write(2,*)i,nprom,x(k),y(k),z(k)
		end do
		end if
	write(13,*)i,Vp
	enddo
	close(1)
	close(13)
	call cpu_time(m2)
	write(14,*)'Configuraciones',m2-ini
!!
	!!	Promedios
!!
!!  g(r)
	m1=m2
	call gr()
	call cpu_time(m2)
!	write(*,*)'funciones g(r)',m2-m1
	write(14,*)'funciones g(r)',m2-m1
!! 
!!  S(q)
!!
	m1=m2
!	call sq()
	call cpu_time(m2)
	write(14,*)'Factores de estructura S(q)',m2-m1
!	write(*,*)'Factores de estructura S(q)',m2-m1
	call cpu_time(fin)
!	write(*,*)'tiempo de ejecución',fin-ini
	write(14,*)'tiempo de ejecución',fin-ini
	close(14)
	stop
100   format(2f12.6)
200   format(3f12.6) 
	End Program DB_FSS
!!
!!============================================================================================================================
!!
!!	Potencial PHS
!!
REAL*8 FUNCTION uij(rij)
Use Var_Prog
Use Var_Pot
Implicit None
Real*8, Intent(in) :: rij
	if(rij.le.diam/dd)then
	uij=1.0d10
	elseif(rij.gt.diam/dd.and.rij.le.rmin)then
		uij=eps1
	elseif(rij.gt.rmin.and.rij.le.rmin2)then
		uij=eps2
	elseif(rij.gt.rmin2)then
		uij=0.0d0
	endif
RETURN
end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina genera una la configuración inicial aleatoria.
!!
	subroutine iniconfigR()
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Real*8 :: ran2
	Real*8 :: dist
	Real*8 :: xR,yR,zR	!! separación entre partículas
	Integer*4 :: i,j!,k
	Do i=1,np
2		x(i)=(ran2(iseed)-0.5d0)*bl
		y(i)=(ran2(iseed)-0.5d0)*bl
		z(i)=(ran2(iseed)-0.5d0)*bl
		do j=1,i-1
			xR=x(i)-x(j)
			yR=y(i)-y(j)
			zR=z(i)-z(j)
			dist=dsqrt(xR**2+yR**2+zR**2)
			if(dist.le.diam)then
			GO TO 2
			end if
		end do
	end do
	return 
	end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina genera una la configuración inicial cúbica.
!!
	subroutine iniconfigC()
	USE Var_Prog
	Implicit None
	Real*8 :: bl2
	Integer*4 :: i,j!,k
	bl2=bl/2.d0
	x(1)=-bl2+0.5d0
	y(1)=-bl2+0.5d0
	z(1)=-bl2+0.5d0
	do i=1,np-1
		x(i+1)=x(i)+1.d0
		y(i+1)=y(i)
		z(i+1)=z(i)
		if (x(i+1).gt.bl2) then
			x(i+1)=-bl2+0.5d0
			y(i+1)=y(i+1)+1.d0
			if(y(i+1).gt.bl2)then
				y(i+1)=-bl2+0.5d0
				z(i+1)=z(i+1)+1.d0
			end if
		endif
	enddo 
	return 
	end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina calcula la energía de la configuración, tomando en cuanta las interacciones dentro del radio de corte.
!!
	subroutine energia()
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Real*8 :: rij,rij2,xij,yij,zij
	Real*8 :: xi,yi,zi
	Real*8 :: uuij,uij
	Integer*4 :: i,j!,k
	
	V=0.d0
	do i=1,np-1
		xi=x(i)
		yi=y(i)
		zi=z(i)
		do j=i+1,np
			xij=x(j)-xi
			yij=y(j)-yi
			zij=z(j)-zi
			xij=xij-bl*dnint(xij/bl)
			yij=yij-bl*dnint(yij/bl)
			zij=zij-bl*dnint(zij/bl)
			rij2=xij*xij+yij*yij+zij*zij
			rij=dsqrt(rij2)
			if (rij.le.rc) then
			uuij=uij(rij)
			V=V+uuij
			endif
		enddo      
	enddo	
	return 
	end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina calcula la energía sobre una partícula, tomando en cuanta las interacciones dentro del radio de corte.
!!
	subroutine energiaU0(i)
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Integer*4, Intent(in) :: i 
	Real*8 :: rij,rij2,xij,yij,zij
	Real*8 :: uuij,uij
	Integer*4 :: j
	
	V0=0.d0
	do j=1,np
		if(j.NE.i)then
		xij=x(j)-x0
		yij=y(j)-y0
		zij=z(j)-z0
		xij=xij-bl*dnint(xij/bl)
		yij=yij-bl*dnint(yij/bl)
		zij=zij-bl*dnint(zij/bl)
		rij2=xij*xij+yij*yij+zij*zij
		rij=dsqrt(rij2)
		if (rij.le.rc) then
		uuij=uij(rij)
		V0=V0+uuij
		endif
		endif
	enddo      
	return 
	end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina calcula la energía sobre una partícula, tomando en cuanta las interacciones dentro del radio de corte.
!!
	subroutine energiaU(i)
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Integer*4, Intent(in) :: i 
	Real*8 :: rij,rij2,xij,yij,zij
	Real*8 :: uuij,uij
	Integer*4 :: j
	
	Vn=0.d0
	do j=1,np
		if(j.NE.i)then
		xij=x(j)-xx
		yij=y(j)-yy
		zij=z(j)-zz
		xij=xij-bl*dnint(xij/bl)
		yij=yij-bl*dnint(yij/bl)
		zij=zij-bl*dnint(zij/bl)
		rij2=xij*xij+yij*yij+zij*zij
		rij=dsqrt(rij2)
		if (rij.le.rc) then
		uuij=uij(rij)
		Vn=Vn+uuij
		endif
		endif
	enddo      
	return 
	end
!!
!!==============================================================================
!!
!!	Funciones de distribución radial g(r)
!!
	subroutine gr()
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Integer*4 :: i,j,k,m
	Real*8 :: da,rij,rij2,xij,yij,zij,aux
	Real*8, Dimension (nr) :: r
	do i=1,nr
		g(i)=0.d0
	enddo
	do k=1,nprom-1 !!calculo de desplazamiento
		do i=1,np-1
			do j=i+1,np
				xij=cfgx(k,j)-cfgx(k,i)
				yij=cfgy(k,j)-cfgy(k,i)
				zij=cfgz(k,j)-cfgz(k,i)
				xij=xij-bl*dnint(xij/bl)
				yij=yij-bl*dnint(yij/bl)
				zij=zij-bl*dnint(zij/bl)
				rij2=xij*xij+yij*yij+zij*zij
				rij=dsqrt(rij2)
				if (rij.le.rc) then
							m=idint(rij/dr)+1
							g(m)=g(m)+2.d0	
				endif
			enddo      
		enddo	   
	enddo
	do i=1,nr
		r(i)=dr/2.d0+dr*(i-1)
		da=4.d0*pi*dr*r(i)*r(i)
	!!		Escalamiento
		r(i)=escR*r(i)
		g(i)=g(i)/(da*np*nprom)
		write(10,100)r(i),g(i)
	enddo
100   format(2f12.6)
200   format(3f12.6) 
	return
	end
!!
!!==============================================================================
!!
!!	Factores de estructura S(q)
!!
	subroutine sq()
	USE Var_Prog
	USE Var_Pot
	Implicit None
	Integer*4 :: i,j,k,ia,ib,ic,ii,jj,kk
	Real*8 :: parti,suma
	Real*8 :: xaux,yaux,zaux,raux,raux2
	Real*8 :: aux
	Real*8, Dimension (nq) :: sigs
	Real*8 :: qx,qy,qz,qq
	Real*8, Dimension (nq) :: ssc,sss
	do i=1,nq
		q(i)=i*dq
		s(i)=0.d0
		sigs(i)=0.d0
	enddo
	do k=1,nprom-1 !!calculo de desplazamiento
		do i=1,nq
			parti=0.d0 
			suma=0.d0
			do j=1,3*NSS**3-2
				ssc(j)=0.d0
				sss(j)=0.d0
			enddo  
			do j=1,np
				xaux=cfgx(k,j)-bl*dnint(cfgx(k,j)/bl)
				yaux=cfgy(k,j)-bl*dnint(cfgy(k,j)/bl)
				zaux=cfgz(k,j)-bl*dnint(cfgz(k,j)/bl)
				raux2=xaux*xaux+yaux*yaux+zaux*zaux
				raux=dsqrt(raux2)
				if (raux.le.rc) then
					parti=parti+1.d0
						Do ia=0,NSS-1
						Do ib=0,NSS-1
						Do ic=0,NSS-1
						if (ia.eq.0.and.ib.eq.0.and.ic.eq.0)then
					ssc(1)=0.d0
					sss(1)=0.d0
						else
		qx=dfloat(ia)
		qy=dfloat(ib)
		qz=dfloat(ic)
		qq=q(i)/dsqrt(qx*qx+qy*qy+qz*qz)
		ii=ia*(NSS**2)+ib*NSS+ic+1
		jj=NSS**3-1+ii
		kk=2*(NSS**3-1)+ii
		ssc(ii)=ssc(ii)+dcos(qq*(qx*xaux+qy*yaux+qz*zaux))
		ssc(jj)=ssc(jj)+dcos(qq*(qx*xaux+qy*yaux-qz*zaux))
		ssc(kk)=ssc(kk)+dcos(qq*(qx*xaux-qy*yaux+qz*zaux))
		sss(ii)=sss(ii)+dsin(qq*(qx*xaux+qy*yaux+qz*zaux))
		sss(jj)=sss(jj)+dsin(qq*(qx*xaux+qy*yaux-qz*zaux))
		sss(kk)=sss(kk)+dsin(qq*(qx*xaux-qy*yaux+qz*zaux))
						endif
						enddo
						enddo
						enddo				
				endif
			enddo
			do j=1,3*NSS**3-2
				suma=suma+ssc(j)*ssc(j)+sss(j)*sss(j)
			enddo
			aux=suma/(parti*(3.d0*NSS**3-2.d0))
			s(i)=s(i)+aux
			sigs(i)=sigs(i)+aux*aux
		enddo
	End do
	do i=1,nq
!! 	Escalamiento
		q(i)=escS*q(i)
		s(i)=s(i)/nprom
		sigs(i)=sigs(i)/nprom
		sigs(i)=sigs(i)-s(i)**2
		sigs(i)=dsqrt(sigs(i)/nprom)
		write(1111,200)q(i),s(i),sigs(i)
	enddo
	return
100   format(2f12.6)
200   format(3f12.6)    
	return
	end
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	ran2
!!
!!	Esta subrutina genera números aleatorios uniformes entre 0 y 1
!!
	REAL*8 FUNCTION ran2(idum)
	implicit integer*4(i-n),real*8(a-h,o-z)
	PARAMETER (im1=2147483563,im2=2147483399,am=1.d0/im1)
	PARAMETER (imm1=im1-1,ia1=40014,ia2=40692)
	PARAMETER (iq1=53668,iq2=52774,ir1=12211,ir2=3791)
	PARAMETER (ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7)
	PARAMETER (rnmx=1.d0-eps)
	dimension iv(ntab)
	SAVE iv,idum2,iy
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then	
		idum=max(-idum,1)	
		idum2=idum
		do j=ntab+8,1,-1	
			k=idum/iq1
			idum=ia1*(idum-k*iq1)-k*ir1
			if (idum.lt.0) idum=idum+im1
			if (j.le.NTAB) iv(j)=idum
		enddo
		iy=iv(1)
	endif
	k=idum/iq1
	idum=ia1*(idum-k*iq1)-k*ir1
	if (idum.lt.0) idum=idum+im1
	k=idum2/iq2
	idum2=ia2*(idum2-k*IQ2)-k*ir2
	if (idum2.lt.0) idum2=idum2+im2
	j=1+iy/ndiv
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+imm1
	ran2=dmin1(am*iy,rnmx)
	return
	END
