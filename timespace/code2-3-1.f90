# Wave field extrapolation program
implicit undefined (a-z)
complex cd(48),ce(48),cf(48),q(48),aa,a,b,c,cshift
real p(96,48,12),phase,pi2,dx,dz,v,z0,x0,dt,dw,lambda,w,wov,x
integer ix,nx,iz,nz,iw,nw,it,nt
open(3,file='plot30',status='new',access='direct',form='unformatted',recl=l)

nt=12; nx=48; nz=96; dx=2; dz=1; pi2=2.*3.141592
v=1; lambda=nz*dz/4; dw=v*pi2/lambda; dt=pi2/(nt*dw); nw=2

do iz=l,nz; do ix=l,nx; do it=l,nt { p(iz,ix,it) = 0. )
do iw = 1,nw { # superimpose nw frequencies
  w = iw*dw; wov = w/v # frequency / velocity
  xO = nx*dx/3; z0 = nz*dz/3
  do ix = 1,nx { # initial conditions for a
    x = ix*dx-x0; # collapsing spherical wave
    phase = -wov*sqrt(z0**2+x**2)
    q(ix) = cexp(cmplx(O.,phase))}

  aa = dz/(4.*(0.,-l.)*wov*dx**2) # tridiagonal matrix coefficients
  a = -aa; b = 1.+2.*aa; c = -aa
  do iz = 1,nz { # extrapolation in depth
    do ix = 2,nx-1 # diffraction term
      cd(ix) = aa*q(ix+ 1) + (1.-2.*aa)*q(ix) + aa*q(ix-l)
    cd(1) = 0.; cd(nx) = 0.
    call ctris(nx,-a,a,b,c,-~,cd,~,ce,cf)
      # "ctrisn solves complex tridiagonal equations
      # i.e. "rtrisn with complex variables
    cshift = cexp(cmplx(0.,wov*dz))
    do ix = 1,nx # shifting term
      q(ix) = q(ix) * cshift
    do it=l,nt { # evolution in time
      cshift = cexp(cmplx(O.,-w*it*dt))
      do ix = 1,nx
        p(iz,ix,it) = p(iz,ix,it)+q(ix)*cshift
      }}}
write(3,rec=l) (((p(iz,ix,it),iz=l ,nz),ix=l,nx),it=l,nt)
stop; end
