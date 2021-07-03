      REAL*8           time
      REAL*8           Coords,dt,epsilon,time0,eps2,Rmaxinit
      REAL*8           RunTime,dPrinttime,amasstotal0
      REAL*8           rs,Rmininit
      INTEGER          np,istep0
      INTEGER          istep
      character*7       path

      PARAMETER(path='./DATA/') ! location to place output files

      PARAMETER (Nsnapshots =   100)   ! maximum number of snapshots

      PARAMETER (Np         =9000)! Max number of particles(max 150000)

      PARAMETER (time0       =  0.)  !Initial time
      PARAMETER (dt        = 0.003)  !Time step
      PARAMETER (epsilon   = 0.003)  !Force softening
      PARAMETER (eps2      = epsilon*epsilon)

      PARAMETER (RunTime    = 0.1)  !final time
      PARAMETER (istep0     =  0   )
      PARAMETER (dprinttime =  0.003)  !time when to print

      PARAMETER (pi         =  3.14159265)   
      PARAMETER (pi4        =  4.*pi)   
      CHARACTER*80     FileNAMES(Nsnapshots)
      COMMON /MAINDATA/Coords(10,Np)
      DATA FileNAMES/'Snap.1.dat','Snap.2.dat','Snap.3.dat',
     &               'Snap.4.dat','Snap.5.dat','Snap.6.dat',
     &               'Snap.7.dat','Snap.8.dat','Snap.9.dat',
     &               'Snap.10.dat',
     &               'Snap.11.dat','Snap.12.dat','Snap.13.dat',
     &               'Snap.14.dat','Snap.15.dat','Snap.16.dat',
     &               'Snap.17.dat','Snap.18.dat','Snap.19.dat',
     &               'Snap.20.dat',
     &               'Snap.21.dat','Snap.22.dat','Snap.23.dat',
     &               'Snap.24.dat','Snap.25.dat','Snap.26.dat',
     &               'Snap.27.dat','Snap.28.dat','Snap.29.dat',
     &               'Snap.30.dat',
     &               'Snap.31.dat','Snap.32.dat','Snap.33.dat',
     &               'Snap.34.dat','Snap.35.dat','Snap.36.dat',
     &               'Snap.37.dat','Snap.38.dat','Snap.39.dat',
     &               'Snap.40.dat',
     &               'Snap.41.dat','Snap.42.dat','Snap.43.dat',
     &               'Snap.44.dat','Snap.45.dat','Snap.46.dat',
     &               'Snap.47.dat','Snap.48.dat','Snap.49.dat',
     &               'Snap.50.dat',
     &               'Snap.51.dat','Snap.52.dat','Snap.53.dat',
     &               'Snap.54.dat','Snap.55.dat','Snap.56.dat',
     &               'Snap.57.dat','Snap.58.dat','Snap.59.dat',
     &               'Snap.60.dat',
     &               'Snap.61.dat','Snap.62.dat','Snap.63.dat',
     &               'Snap.64.dat','Snap.65.dat','Snap.66.dat',
     &               'Snap.67.dat','Snap.68.dat','Snap.69.dat',
     &               'Snap.70.dat',
     &               'Snap.71.dat','Snap.72.dat','Snap.73.dat',
     &               'Snap.74.dat','Snap.75.dat','Snap.76.dat',
     &               'Snap.77.dat','Snap.78.dat','Snap.79.dat',
     &               'Snap.80.dat',
     &               'Snap.81.dat','Snap.82.dat','Snap.83.dat',
     &               'Snap.84.dat','Snap.85.dat','Snap.86.dat',
     &               'Snap.87.dat','Snap.88.dat','Snap.89.dat',
     &               'Snap.90.dat',
     &               'Snap.91.dat','Snap.92.dat','Snap.93.dat',
     &               'Snap.94.dat','Snap.95.dat','Snap.96.dat',
     &               'Snap.97.dat','Snap.98.dat','Snap.99.dat',
     &               'Snap.100.dat'/

c     Initial conditions parameters (for equilibriu sphere)
      
      PARAMETER (Rmaxinit= 10.0)
      PARAMETER (Rmininit= 0.001)
      PARAMETER (aMasstotal0= 1.)
      PARAMETER (rs = 0.25 ) 
