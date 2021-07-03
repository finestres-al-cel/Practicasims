c      Generation of initial conditions file (ic.dat)
c      Homogeneous sphere.

       INCLUDE '../nbody1.h' !Parameters file
       real*8 x1,y1,z1,x,y,z,ran
       integer n,i,j,seed
       dimension x(1000000),y(1000000),z(1000000),ran(1000000)
       i=0
       n=np
       seed=86456

c       do j=1,3*n
         call RANDOM_NUMBER(ran)
c         ran(j)=rand(seed)!int(3.145654*j**2/pi**3)+int(j/pi*(pi+j))) !Generation of random numbers for positions and
c velocities
c       do j=1,3*n
c         write(*,*)ran(j)
c       enddo

c      Generate an uniform distribution of positions in a sphere with R=1

       do j=0,3*n,3
          x1=ran(j)/0.5-1.
          y1=ran(j+1)/0.5-1.
          z1=ran(j+2)/0.5-1.

         sumt=x1*x1+y1*y1+z1*z1 !Radius
       
        if (sumt.le.1.) then
          i=i+1
          x(i)=x1
          y(i)=y1
          z(i)=z1
        endif

       enddo

       open(1,file='../DATA/ic.dat',status='unknown')
       do j=1,i
         write(1,101)x(j),y(j),z(j),0.,0.,0.,0.,0.,0.,1.
       enddo

       close(1)
 100   format(i10,2f18.10,i10,3f18.10)
       
 101   format(10f18.10)


       end
    
       
       

