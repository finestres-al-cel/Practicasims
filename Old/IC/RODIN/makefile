
#opt =  -byteswapio -O3  -mcmodel=medium -Mlarge_arrays
opt =  -O3 #-mcmodel=medium -i-dynamic 
#opt = -O3 -xW -openmp -ipo -no-prec-div -i-dynamic -mcmodel=medium -convert big_endian
#opt = -O3 -openmp -i-dynamic -mcmodel=medium -convert big_endian

lib =      
#fort = ifort                    
fort = gfortran               
#fort = ifc
load = $(fort) $(opt)  $(olist) $(lib) -o Rodin1.x  

comp = $(name).o : $(name).f $(dep);          \
                    $(fort) $(opt) -c  $(name).f $(lib) 


Rodin2 : PMparameters.h  Rodin2.h Rodin2.o; \
          $(fort) $(opt) Rodin2.o -o Rodin2 

Rodin2.o : Rodin2.f PMparameters.h  Rodin2.h; \
            $(fort) $(opt) -c Rodin2.f 




