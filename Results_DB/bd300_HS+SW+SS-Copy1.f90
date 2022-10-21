!! 	Este programa se usa para calcular la función g(r), el s(q), la función Fself y el <(dr)^2> 
!! para un sistema monodisperso 3D de particulas en un potencial tipo SW+SS mediante dinámica browniana. 
!! Todas las longitudes están escaladas con la distancia promedio entre partículas d=n^-1/3.
!!
!! ===============================================================================
!!
!!		    Autor:			      Fecha:
!! 
!!	Miguel Ángel Sandoval Puentes		   27-ENERO-2018
!! 
!!	Como parte del trabajo para la obtención del grado de maestria en ciencias físicas. 
!! Se usó como base un programa proporcionado por el Dr. José Miguel Méndez Alcaraz, asesor.
!! ===============================================================================
!	PROGRAM Dinamica_Browniana
!! 
!!Declaración de variables
!!
	implicit integer*4(i-n),real*8(a-h,o-z)
	character*20 gout,dout,eout,din
	parameter(mp=343,mt=20000,mr=350)
	dimension x(mp),y(mp),z(mp),r(mr),t(mt),t1(mt),d(mt)
	dimension cfgx(mt,mp),cfgy(mt,mp),cfgz(mt,mp) !!mt Numero de pasos en tiempo
	common/bp/bl,rc,np
	common/cp/dt,nct,ncet
	common/strp/dr,nr,nss
	common/strf/g(mr)
	common/f/fx(mp),fy(mp),fz(mp)
	common/ran/iseed
	common/par/diam,dd
	common/sp/eps,eps1,eps2,H,H2,rho,eex,sige
	common/pot/rlmb,FactFP,dk,dlt,nk
	common/rms/rmin,rmin2,rmin3,rmin4,rmax
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Variables de entrada
!!
	pi=4.d0*datan(1.d0)
	phi=0.4d0	!! factor de llenado
	rho=1.d0
	NSS=3
	bl=7.d0
	rc=bl/2.d0
	np=idint(bl**3)
	diam=1.0d-6
	dd=(pi*diam*diam*diam/(6*phi))**(1.d0/3)
	dtt=3.81d-6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	dtp=1.0d-6
	nct=2000000 
	ncp=20000000
	ncet=1000 
	ncep=1000
	dr=0.01d0
	iseed=50
	nr=idint(rc/dr)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Parámetros del potencial
!!
	rlmb=50.d0
	eps=1.0d0/1.474d0
	factFP=rlmb*eps*(rlmb/(rlmb-1.0d0))**(rlmb-1.0d0)
	rmin=rlmb*diam/((rlmb-1.0d0)*dd)
	dk=10000.d0
	A=1.6654
	B=3.6975
	C=1.3414
	al=A/(eps2**C+B)
	d0=dsqrt(-log(al)/dk)

	eps1=-1.25d0
	eps2=0.1!!!!5.0d0
	H=0.25d0*diam/dd
	H2=0.03d0*diam/dd!+d0
	nk=(dd/diam)*300
	dlt=pi/nk
	ad=diam			!! Adimensionalización
	escR=dd/ad

	rmin2=diam/dd+H-0.394916451232002*dlt
	rmin3=diam/dd+H+(1.d0-0.394916451232002)*dlt
	rmin4=diam/dd+H+H2-d0
	rmax=dsqrt(100.d0/dk)+diam/dd+H+H2
	Write(*,*)escR*rmin,diam/dd
	Write(*,*)escR*rmin2,escR*rmin3,rmin4,rmax
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Escalamiento
!!

	escS=ad/dd
	escD=(dd/ad)**2.d0
	escE=1.d0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Vamos a verificar los potenciales
!!
	open(1,file='pot300.dat',status='unknown')
	Do i=1,100*nr
		rij=i*0.01d0*dr
		u11=uij(rij,aux)
		write(1,*)escR*rij,u11,aux
	enddo
	close(1)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Archivos de salida
!!
	din='datos2.dat'
	gout='g2.dat'
	dout='d2.dat'
	eout='e2.dat'
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
	do i=1,nr
		g(i)=0.d0
	enddo
	eex=0.d0
	sige=0.d0
	nprom=0 
	dt=dtt     
	open(10,file=gout,status='unknown')
	open(12,file=dout,status='unknown')
	open(13,file=eout,status='unknown') 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	Esta sección guarda los datos usados para generar cada paquete de resultados
!!    
	open(14,file=din,status='unknown') 
	write(14,*)'diametro=',diam
	write(14,*)'factor de llenado=',phi
	write(14,*)'semilla=',iseed
	close(14)

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
	!!!!!call iniconfig(x,y,z)
	open(1320,file='conf_fin.dat', status='unknown')
    
	do i=1,mp
		READ(1320,*) x(i),y(i),z(i)
	enddo
	do i=1, mp
	    !write(*,*)x(i),y(i),z(i)
	enddo
        
	do i=1,nct
		call thermal(x,y,z,i)
		call newconfig(x,y,z)
	enddo
	dt=dtp
	do i=1,ncp
		decimal=dfloat(i)/dfloat(ncep)-i/ncep
		if (decimal.lt.1.d-6) then
			call promedio(x,y,z,i)
			nprom=nprom+1
!			write(*,*)nprom
			t(nprom)=dt*ncep*(nprom-1)
			do j=1,np
				cfgx(nprom,j)=x(j)
				cfgy(nprom,j)=y(j)
				cfgz(nprom,j)=z(j)
			enddo
		else
			call thermal(x,y,z,0)
		endif
		call newconfig(x,y,z)
	enddo
	open(1,file='conf_fin2.dat',status='unknown')
	Do i=1,np

		write(1,*)x(i),y(i),z(i)
	end do
	eex=eex/nprom
	sige=sige/nprom
	sige=sige-eex**2
	sige=dsqrt(sige/nprom)
	write(10,*)'excess energy/kTN=',eex
	write(10,*)'sigma(excess energy/kTN)=',sige
	do i=1,nr
		r(i)=dr/2.d0+dr*(i-1)
		da=4.d0*pi*dr*r(i)*r(i)
!!		Escalamiento
		r(i)=EscR*r(i)
		g(i)=g(i)/(da*np*nprom)
		write(10,100)r(i),g(i)
	enddo
	open(15,file='Fself2.dat',status='unknown')
	qqs=2.0d0*pi*EscR
	do i=1,nprom-1 !!calculo de desplazamiento
		Fs=0.d0
		Fs2=0.d0
		dif2=0.d0
		sigd=0.d0
		do j=1,nprom-i
			do k=1,np
				difx=cfgx(j+i,k)-cfgx(j,k)
				dify=cfgy(j+i,k)-cfgy(j,k)
				difz=cfgz(j+i,k)-cfgz(j,k)
				Do ia=0,NSS-1
				Do ib=0,NSS-1
				Do icc=0,NSS-1
			if (ia.ne.0.or.ib.ne.0.or.icc.ne.0)then
			qqx=dfloat(ia)
			qqy=dfloat(ib)
			qqz=dfloat(icc)
			qq=qqs/dsqrt(qqx*qqx+qqy*qqy+qqz*qqz)
			Fs=Fs+Dcos(qq*(qqx*difx+qqy*dify+qqz*difz))
			Fs=Fs+Dcos(qq*(qqx*difx+qqy*dify-qqz*difz))
			Fs=Fs+Dcos(qq*(qqx*difx-qqy*dify+qqz*difz))
			endif
				enddo
				enddo
				enddo
			aux=difx*difx+dify*dify+difz*difz
			Fs2=Fs2+dsin(qqs*dsqrt(aux))/(qqs*dsqrt(aux))
			dif2=dif2+aux
			sigd=sigd+aux*aux
			enddo

		enddo
		d(i)=dif2
!		d(i)=dif2/(4.d0*t(i+1))
		d(i)=d(i)/(np*(nprom-i))
!		sigd=sigd/(16.d0*t(i+1)**2)
		sigd=sigd/(np*(nprom-i))
		sigd=sigd-d(i)**2
		sigd=dsqrt(sigd/(np*(nprom-i))) 
!! 	Escalamiento
		t1(i+1)=EscD*t(i+1) 
	write(15,*)EscD*t(i+1),Fs/(Np*(nprom-i)*(3.d0*(NSS**3-1.d0))),Fs2/(Np*(nprom-i))       
		write(12,200)t1(i+1),EscD*d(i),sigd     
	enddo
	close(15)
	close(10)
	close(11)
	close(12)
	close(13)
	stop
100   format(2f12.6)
200   format(3f12.6)      
	end  
!!
!!============================================================================================================================
!!
!!	Potencial SW+SS
!!
	REAL*8 FUNCTION uij(rij,aux)
	implicit integer*4(i-n),real*8(a-h,o-z)
	common/par/diam,dd
	common/sp/eps,eps1,eps2,H,H2,rho,eex,sige
	common/pot/rlmb,FactFP,dk,dlt,nk
	common/rms/rmin,rmin2,rmin3,rmin4,rmax
	if(rij.le.rmin)then
	uij=factFP*((diam/(rij*dd))**rlmb-(diam/(rij*dd))**(rlmb-1.0))+eps+eps1
	fact1=(rlmb*dd/diam)*(diam/(rij*dd))**(rlmb+1.0d0)
	fact2=((rlmb-1.0d0)*dd/diam)*(diam/(rij*dd))**rlmb
	aux=factFP*(fact1-fact2)
	aux=aux/rij
	elseif(rij.gt.rmin.and.rij.le.rmin2)then
		uij=eps1
		aux=0.d0
	elseif(rij.gt.rmin2.and.rij.le.rmin3)then
	uij=(eps1+eps2)/2-((abs(eps1)+abs(eps2))/2)*Dcos(nk*(rij-rmin2))
	aux=-(abs(eps1)+abs(eps2))*nk*Dsin(nk*(rij-rmin2))/(2.0d0*rij)
	elseif(rij.gt.rmin3.and.rij.le.rmin4)then
		uij=eps2
		aux=0.d0
	elseif(rij.gt.rmin4)then
		if(dk*(rij-rmin4)**2.le.100.d0)then
		uij=eps2*dexp(-dk*(rij-rmin4)**2)
		aux=2.d0*dk*(rij-rmin4)*uij
		aux=aux/rij
		else
		uij=0.d0
		aux=0.d0
		endif
	endif
		RETURN
	end
!!
!!============================================================================================================================
!!
!!	Esta subrrutina genera una la configuración inicial del sistema, ordena las moleculas en una red.
!!
	subroutine iniconfig(x,y,z)
	implicit integer*4(i-n),real*8(a-h,o-z)
	parameter(mp=343)
	dimension x(mp),y(mp),z(mp)
	common/bp/bl,rc,np
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
!! Esta subrrutina calcula la fuerza total ejercida sobre cada partícula por el resto
!!     
	subroutine thermal(x,y,z,ic)
	implicit integer*4(i-n),real*8(a-h,o-z)
	parameter(mp=343)
	dimension x(mp),y(mp),z(mp)
	common/sp/eps,eps1,eps2,H,H2,rho,eex,sige
	common/bp/bl,rc,np
	common/cp/dt,nct,ncet
	common/f/fx(mp),fy(mp),fz(mp)
	common/par/diam,dd
	pi=4.d0*datan(1.d0)
	eexi=0.d0
	do i=1,np
		fx(i)=0.d0
		fy(i)=0.d0
		fz(i)=0.d0
	enddo
	do i=1,np-1
		do j=i+1,np
			xij=x(j)-x(i)
			yij=y(j)-y(i)
			zij=z(j)-z(i)
			xij=xij-bl*dnint(xij/bl)
			yij=yij-bl*dnint(yij/bl)
			zij=zij-bl*dnint(zij/bl)
			rij2=xij*xij+yij*yij+zij*zij
			rij=dsqrt(rij2)
			if (rij.le.rc) then
				uuij=uij(rij,aux)
				eexi=eexi+uuij
				fxij=aux*xij
				fyij=aux*yij
				fzij=aux*zij
				fx(i)=fx(i)-fxij
				fx(j)=fx(j)+fxij
				fy(i)=fy(i)-fyij
				fy(j)=fy(j)+fyij
				fz(i)=fz(i)-fzij
				fz(j)=fz(j)+fzij
			endif
		enddo      
	enddo
	if (ic.eq.0) return
	decimal=dfloat(ic)/dfloat(ncet)-ic/ncet
	if (decimal.lt.1.d-6) then
		write(13,100)ic,eexi
	endif
	return
100   format(i12,f12.4)      
	end 
!!
!!==============================================================================
!!
!! Esta subrrutina calcula las nuevas pocisiones de las partículas debido a la fuerza ejercida sobre ellas y al termino estocástico.
!!           
	subroutine newconfig(x,y,z)
	implicit integer*4(i-n),real*8(a-h,o-z)
	parameter(mp=343)
	dimension x(mp),y(mp),z(mp)
	common/bp/bl,rc,np
	common/cp/dt,nct,ncet
	common/f/fx(mp),fy(mp),fz(mp)
	common/ran/iseed
	sigma=dsqrt(2.d0*dt)
	do i=1,np
		dx=sigma*gasdev(iseed)
		dy=sigma*gasdev(iseed)
		dz=sigma*gasdev(iseed)
		x(i)=x(i)+fx(i)*dt+dx
		y(i)=y(i)+fy(i)*dt+dy
		z(i)=z(i)+fz(i)*dt+dz
	enddo
	return 
	end      
!!
!!==============================================================================
!!
!! Esta subrrutina tiene una función similar a Thermal, pero se usa para cuando el sistema ya termializó, es decir, cuando ya es de interés estudiar el sistema en si y por tanto también una sección para guardar los datos relacionados con la funcion g(r)
!!      
	subroutine promedio(x,y,z,ic)
	implicit integer*4(i-n),real*8(a-h,o-z)
	parameter(mp=343,mr=350)
	dimension x(mp),y(mp),z(mp)
	common/sp/eps,eps1,eps2,H,H2,rho,eex,sige
	common/bp/bl,rc,np
	common/cp/dt,nct,ncet
	common/strp/dr,nr,nss
	common/strf/g(mr)
	common/f/fx(mp),fy(mp),fz(mp)
	common/par/diam,dd
	pi=4.d0*datan(1.d0)
	eexi=0.d0
	do i=1,np
		fx(i)=0.d0
		fy(i)=0.d0
		fz(i)=0.d0
	enddo
	do i=1,np-1
		do j=i+1,np
			xij=x(j)-x(i)
			yij=y(j)-y(i)
			zij=z(j)-z(i)
			xij=xij-bl*dnint(xij/bl)
			yij=yij-bl*dnint(yij/bl)
			zij=zij-bl*dnint(zij/bl)
			rij2=xij*xij+yij*yij+zij*zij
			rij=dsqrt(rij2)
			if (rij.le.rc) then
				uuij=uij(rij,aux)
				eexi=eexi+uuij
				fxij=aux*xij
				fyij=aux*yij
				fzij=aux*zij
				fx(i)=fx(i)-fxij
				fx(j)=fx(j)+fxij
				fy(i)=fy(i)-fyij
				fy(j)=fy(j)+fyij
				fz(i)=fz(i)-fzij
				fz(j)=fz(j)+fzij
				m=idint(rij/dr)+1
				g(m)=g(m)+2.d0
			endif
		enddo      
	enddo
	eex=eex+eexi
	sige=sige+eexi**2
	write(13,100)nct+ic,eexi
	return
100   format(i12,f12.4)            
	end 
!!
!!================================================================================================================
!!
!!	Esta subrutina genera números aleatorios con distribución gaussiana a partir de números generados de forma uniforme
!!
	function gasdev(idum)
	implicit integer*4(i-n),real*8(a-h,o-z)
	save iset,gset
	data iset/0/
	if (iset.eq.0) then
10       v1=2.d0*ran1(idum)-1.d0
		v2=2.d0*ran1(idum)-1.d0
		rsq=v1*v1+v2*v2
		if (rsq.ge.1.d0.or.rsq.eq.0.d0) goto 10
		fac=dsqrt(-2.d0*dlog(rsq)/rsq)
		gset=v1*fac
		gasdev=v2*fac
		iset=1
	else
		gasdev=gset
		iset=0
	endif
	return
	end
!!
!!================================================================================================================
!!
!!	Genera numeros aleatorios con distribución uniforme
!!             
	function ran1(idum)
	implicit integer*4(i-n),real*8(a-h,o-z)
	parameter(ia=16807,im=2147483647,am=1.d0/im,iq=127773)
	parameter(ir=2836,ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7)
	parameter(rnmx=1.d0-eps)
	dimension iv(ntab)
	save iv,iy
	data iv /ntab*0/, iy /0/
	if (idum.le.0.or.iy.eq.0) then
		idum=max0(-idum,1)
		do j=ntab+8,1,-1
			k=idum/iq
			idum=ia*(idum-k*iq)-ir*k
			if (idum.lt.0) idum=idum+im
			if (j.le.ntab) iv(j)=idum
		enddo
		iy=iv(1)
	endif
	k=idum/iq
	idum=ia*(idum-k*iq)-ir*k
	if (idum.lt.0) idum=idum+im
	j=1+iy/ndiv
	iy=iv(j)
	iv(j)=idum
	ran1=dmin1(am*iy,rnmx)
	return
	end

