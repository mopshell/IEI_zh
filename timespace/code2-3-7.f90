# Migration in the (omega,x,z)-domain
real q(48,64),pi2,alpha,dt,dtau,dw
complex cq(48,64),cd(48),ce(48),cf(48),aa,a,b,c,cshift
integer ix,nx,iz,nz,iw,nw,it,nt
open(4,file='plot36', status='new', access='direct',form= 'unformatted', recl=l)

nt = 64; nz = nt; nx = 48; pi2=2.*3.141592
dt=1.; dtau=1.; dw=pi2/(dt*nt); nw=nt/2;
alpha = .25       # alpha = v*v*dtau/(4*dx*dx)
do iz=l,nz; do ix=l,nx; { q(ix,iz) = 0.; cq(ix,iz)=O. )
do it=nt/3,nt,nt/4
  do ix=1,4 # Broadened impulse source
    { cq(iu,it) = ( 5 . 4 ~ ) ; cq(ix,it+l) = (5.-ix) )
call rowcc(nx,nt,cq,+l.,+l.) # F.T. over time.
do iz = 1,nz { # iz and iw loops interchangeable
do iw = 2,nw { # iz and iw loops interchangeable
  aa =- alpha /( (O.,-l.)*(iw-l)*dw )
  a = -aa; b = 1.+2.+aa; c = -aa
  do ix = 2,nx-1
    cd(ix) = aa*cq(ix+l,iw) + (1.-2.*aa)*cq(ix,iw)+ aa*cq(ix-1,iw)
  cd(1) = 0.; cd(nx) = 0.
  call ctris(nx,-a,a7b,c,-c,cd,cq(l,iw),ce,cf)
  cshift = cexp(cmplx(O.,-(iw-l)*dw*dtau))
  do ix=l,nx
    cq(ix,iw) = cq(ix,iw) * cshift
  do ix = 1,nx
    q(ix,iz) = q(ix,iz)+cq(ix,iw) # q(t=O) = C Q(w)
  }}

write(4,rec=l) ((q(ix,iz),iz=l,nz),ix=l,nx)
stop; end
