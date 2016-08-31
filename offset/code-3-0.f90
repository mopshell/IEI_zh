# Synthetic marine data tape movie generation
integer kbyte,it,nt,ih,nh,is,ns,iz,nz,itO,iy
real p(512),b(512),refl(25116),z(25),geol(25),random
open(3,file="plot" ,status='new',acce='direct',form='unformatted',recl=1)
nt = 512; nh = 48; ns = 10; nz = 25;kbyte = 1
do iz=1,nz # Reflector depth
  z(iz) = nt*random() # random() is on the interval (O.,1.)
do iz=1,nz # Reflector strength with depth.
  geol(iz) = 2.*random()-1.
do is = 1,ns # Give texture to the Geology
  do iz = 1,nz
    refl(iz,is) = (1.+random())*geol(iz)
do it = 1,nt # Prepare a wavelet
  b(it) = exp(-it*.08)*sin(.5*it-.5)
do is = ns,1,-1 { # Shots. Run backwards.
  do ih = 1,nh { # down cable h = (g-s)/2
    iy = (is-1)+(ih-1) # y = midpoint
    iy = 1+ (iy-ns*(iy/ns)) # periodic with midpoint
    do it = 1,nt
      p(it) = 0.
      do iz = 1,nz { # Add in a hyperbola for each layer
        it0 = sqrt( z(iz)**2 + 100.*(ih-1)**2 )
        do it = 1,nt-it0 { # Add in the wavelet
          p(it+it0)= p(it+itO) + refl(iz,iy)*b(it)
          }
        }
      write(3,rec=kbyte) (p(it),it=1,nt); kbyte = kbyte+nt*4
      }
  }
stop; end
