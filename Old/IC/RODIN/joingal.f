      implicit none
      real*8 pos,vel,mass,no,pos1,vel1,mass1
      integer i,j,k,n
      dimension pos(1000000,3),vel(1000000,3),mass(1000000)
      dimension pos1(1000000,3),vel1(1000000,3),mass1(1000000)

      open(1,file='ic.dat',status='unknown')
      n=0
      do i=1,10000000
      read(1,*,end=10)(pos(i,j),j=1,3),(vel(i,j),j=1,3),no,no,no,mass(i)
      n=n+1
      enddo
 
  10  continue
      close(1)

      do i=1,n
       do j=1,3
       pos1(i,j)=pos(i,j)+15.0d0
       enddo
       vel1(i,1)=vel(i,1)-6.0d0
       vel1(i,2)=vel(i,2)-0.0d0
       vel1(i,3)=vel(i,3)-0.0d0
       vel(i,1)=vel(i,1)+2.0d0
       vel(i,2)=vel(i,2)+0.0d0
       vel(i,3)=vel(i,3)+0.0d0
      enddo
       
      open(2,file='ic_new.dat',status='unknown')
       do i=1,n
         write(2,100)(pos(i,j),j=1,3),(vel(i,j),j=1,3),no,no,no,mass(i)
       enddo
       do i=1,n
       write(2,100)(pos1(i,j),j=1,3),(vel1(i,j),j=1,3),no,no,no,mass(i)
       enddo
      close(2)
 100  format(10(f20.8,1x))
      end
      
