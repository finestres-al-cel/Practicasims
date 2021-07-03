c ---------Practica simulacions Joves i Ciencia 2015----------------
c ---------------Dr. Santi Roca-Fàbrega-----------------------------
c --------------Universitat de Barcelona----------------------------
c ------------------sroca@am.ub.es----------------------------------
c-------------------------------------------------------------------- 

c                    Run simple N-body simulations
c                    using leap-frog integration scheme
c
c                    Logic: - read simulation parameters from nbody1.h 
c                           - read initial conditions (coordinates, velocities, 
c                             masses) from file 'ic.dat' 
c                           - run for total time RunTime (defined in nbody1.h)
c                           - save the data for continuation of the run
c                           in particles.dat 
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h' !Parameters file
      REAL*8    Ekin,Epot,Etot,yc(3)
      integer  npi
   
      COMMON /Particles/ npi

                          ! open files
      Open(1,file='DATA/ic.dat',status='unknown') ! Initial conditions file
      Open(2,file=path//'trajectories.dat',status='unknown') !trajectories control
      Open(3,file=path//'energies.dat',position='append') ! Energy
control
      Open(8,file=path//'analysis.dat',position='append') ! analysis file

      time=time0
      PrintTime=0.
      istep=istep0 
      npi=0
      write(*,*) ' Time=',time
      write(*,*) ' Time step          =',dt,' Step=',iStep
      write(*,*) ' Force softening    =',epsilon
      write(*,*) ' Run Time           =',RunTime

c    Read data from initial conditions file 'ic.dat'

         Do i=1,Np
            read(1,*,end=10) (Coords(k,i),k=1,10)
            npi=npi+1
         EndDo 


  10     continue

         write(*,*)'Real number of particles =', npi

         E0 =-0.88451

                          ! set parameters
      Nsteps     = INT(RunTime/dt)  ! number of steps for this run
      PrintTime  = time +dPrintTime ! next output time
      Ntraj      = 5                ! save this number of trajectories
      iSnapshot  = 1                ! current snapshot
      ienergy    = 50              ! save energies this often
      indx       = 0

c     Main simulation loop
 
      Do i=1,Nsteps            

c     Compute accelerations at each position in a subroutine
 
         Call GetAccelerations 

c     Move particles to new positions (t+dt) using accelerations

         Call MoveParticles

c     Actualize time and step number for the next computations
         time = time + dt
         iStep= iStep+ 1

c     Write particle trajectories for orbital analysis

         Call WriteTrajectories(Ntraj,time)

c     Compute kinetic and potential energy

         Call Energies(Ekin,Epot)

c     Compute total energy to check energy conservation

            Etot = Ekin+Epot

c     Write energy control file, also window output

            write(*,20) time,iStep,Ekin,Epot,(Etot/E0-1)*100.,Ekin/Etot
            write(3,20) time,iStep,Ekin,Epot,Etot,Ekin/Etot

c     Make first analysis

            Call Analyze(indx,yc,time)

c     Write positions, velocities and masses in a file

         If(time.ge.PrintTime)Then
            Call WriteSnapshot(yc,time)
            PrintTime  = PrintTime + dPrintTime !actualization for new 
c print time
         EndIf 
      EndDo 

c     Save final positions and velocities

      Call SaveMoment


 20   format(f8.3,i8,6g13.5)

      Stop
      End
c-------------------------------------------------------------------- 
c               Find Acceleations for all particles
      SUBROUTINE GetAccelerations
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'
      REAL*8   acc,rr

      COMMON /particles/ npi

      Do i=1,Npi         ! set acceleration counters to zero
         Do k=7,9
            Coords(k,i) =0.
         EndDo 
      EndDo 
       Do i=1,Npi-1       ! sum contributions pair-wise
         Do j=i+1,Npi
            rr = (sqrt(
     &           (Coords(1,i)-Coords(1,j))**2+
     &           (Coords(2,i)-Coords(2,j))**2+
     &           (Coords(3,i)-Coords(3,j))**2+ eps2 ) )**3
            Do k=1,3
               acc           =  (Coords(k,i)-Coords(k,j))/rr
               Coords(k+6,i) =   Coords(k+6,i) -Coords(10,j)*acc
               Coords(k+6,j) =   Coords(k+6,j) +Coords(10,i)*acc
            EndDo 
         EndDo 
      EndDo 
      Return
      End
c-------------------------------------------------------------------- 
c     Equation of motion pos=posini+vel*dt          
      SUBROUTINE MoveParticles
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'

      COMMON /particles/ npi

      Do i=1,Npi         
         Do k=1,3
            Coords(k+3,i) = Coords(k+3,i) + Coords(k+6,i)*dt
            Coords(k  ,i) = Coords(k  ,i) + Coords(k+3,i)*dt

         EndDo 
      EndDo 
      Return
      End
c-------------------------------------------------------------------- 
c               write Ntraj trajectories into file 2
      SUBROUTINE WriteTrajectories(Ntraj,time)
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'

      write(2,10)time,iStep,((Coords(i,j),i=1,3),j=1,Ntraj)


 10   format(g12.5,i8,20(2x,3f9.4))
      Return
      End
