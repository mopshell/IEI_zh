# Time Domain 15-degree Diffraction Movie
# Star: w=p(t  ,z) y=p(t  ,z+1)
# Star: u=p(t+1,z) v=p(t+1,z+1)
real p(36,96),u(36),w(36),v(36),y(36),e(36),f(36),d(36),2(96),alfa,beta
integer ix,nx,iz,nz,it,nt,kbyte
nx = 36; nz = 96; n t = 96; kbyte=1
alfa = .125 # v*dz*dt/(8*dx*dx)
beta = .140 # accurate x derivative parameter; simplest case b=O.
open(3,file='plot40',status='new',access='direct',form='unformatted',recl=1)
do iz=1,nz; do ix=1,nx; p(ix,iz) = 0. # clear space
do iz=nz/5,nz,nz/4 # Set up initial model
  do it=1,15 # of 4 band limited
      do ix=1,4 # "point" scatterers.
        p(ix,it+iz) = (5.-ix)*(8-it)*exp(-.1*(it-8)**2)
apb = alfa+beta; amb = alfa-beta # tridiagonal coefficients
diag = 1.+2.*amb; offdi = -amb
do iz=nz,2,-2 { # Climb up in steps of 2 z-levels
  do i=1,nz; z(i)=O.; z(iz)=1. # Pointer t o current z-level
    write(3,rec=kbyte) (z(i),i=1,nz),((p(ix,i),i=1,nz),ix=1,nx)
    kbyte = kbyte + nx*nz*4 + nz*4
      do ix=1,nx
        { u(ix) = p(ix,iz-1); v(ix) = u(ix) )
      do it=iz,nt {
        do ix=1,nx #update the differencing star
          { w(ix) = u(ix); y(ix) = v(ix); v(ix) = p(ix,it) )
        dd = (1.-apb)*(v(1)+w(1))+apb*(v(2)+~(2))
        d(1) = dd-diag*y(1)-offdi*(y(1)+y(2))
        do ix=2,nx-1 {
          dd = (1.-2.*apb)*(v(ix)+w(ix))
          dd = dd + apb*(v(ix-1)+w(ix-1)+v(ix+1)+w(ix+1))
          d(ix) = dd-diag*y(ix)-offdi*(y(ix-1)+y(ix+1)) }
        dd = (1.-apb)*(v(nx)+w(nx))+apb*(v(nx-1)+w(nx-1))
        d(nx) = dd-diag*y(nx)-offdi*(y(nx)+y(nx-1))
        call rtris(nx,diag+offdi,offdi,diag,offdi,diag+offdi,d,u,e,f)
        do ix=1,nx
          p(ix,it) = u(ix)
        }
    }
do i=1,nz; z(i)=0.; z(1)=1.
write(3,rec=kbyte) (z(i),i=1,nz),((p(ix,i),i=1,nz),ix=1,nx)
stop; end
