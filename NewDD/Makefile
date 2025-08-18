#CC = mpicc
#OPT = -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DNMEG=16000 -DWGROUPSIZE=5  -DUSE_MPI
CC = mpicc
OPT = -O3 -DNENER=0 -DNPRE=8 -DNMEG=20000 -DWGROUPSIZE=5  -DUSE_MPI -DQUADHILBERT -DREAD_SINK -DNCHEM=9 -DNDUST=4
#OPT = -O3 -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DNMEG=20000 -DWGROUPSIZE=5  -DUSE_MPI
#CC = icc
#OPT = -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DNMEG=20000 -DWGROUPSIZE=8  


RANLIB = ranlib


CFLAGS = $(OPT) 
CCLAGS = $(OPT)
#CAMB Makefile
#Edit for your compiler
#Note there are many ifc versions, some of which behave oddly

#Intel , -openmp toggles mutli-processor:
#note version 10.0 gives wrong result for lensed when compiled with -openmp [fixed in 10.1]
#F90C     = ifort
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0

# Intel 9 on IA-64 (eg. COSMOS)
# (do "module load icomp90" before compiling)
#F90C = ifort
#FFLAGS = -openmp -fpp2 -w -O3 -ip -mP2OPT_hlo_prefetch=F

#Intel ifc, add -openmp for multi-processor (some have bugs):
#F90C     = ifc
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -quiet -fpp2
#some systems can can also add e.g. -tpp7 -xW

#G95 compiler
#F90C   = g95
#FFLAGS = -O2

#Gfortran compiler: if pre v4.3 add -D__GFORTRAN__
#F90C     = gfc
#FFLAGS =  -O2 

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#F90C     = f90
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1

#Absoft ProFortran, single processor:
#F90C     = f95
#FFLAGS = -O2 -cpu:athlon -s -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_"

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3

#PGF90
#F90C = pgf90
#FFLAGS = -O2 -DESCAPEBACKSLASH

#Sun V880
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O2 -xarch=native64 -openmp -ftrap=%none
#try removing -openmp if get bus errors. -03, -04 etc are dodgy.

#IBM XL Fortran, multi-processor (run gmake)
#F90C     = xlf90_r
#FFLAGS  = -DESCAPEBACKSLASH -DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3


LIBS       = libmyram.a

#Shouldn't need to change anything else...


#TIMEROBJ       = ost.o ost.ver1.o ost.ver2.o Model4TreeForce.o fof.o nnost.o
#TIMEROBJ       = rd_amr.o  rd_info.o  rd_part.o header.o utils.o Memory2.o find_leaf_gas.o rd_grav.o rd_hydro.o rd_sink.o
TIMEROBJ       = rd_amr.o  rd_info.o  rd_part.o header.o utils.o Memory2.o find_leaf_gas.o rd_hydro.o rd_sink.o


#TEST	= test.exe
#TESTOBJ	= test.o

NEWDD	= newdd.exe
NEWDDOBJ	= newdd.o

INCLUDES = -I./ -I../

default: 
	$(MAKE) $(LIBS)
	$(RANLIB) $(LIBS)


all: 
	$(MAKE) $(LIBS)
#	$(MAKE) $(TEST)
	$(MAKE) $(NEWDD)

	$(MAKE) $(DMBASE)

new: 
	$(MAKE) clean  
	$(MAKE) $(LIBS)

#test:
#	-rm -f test.exe
#	$(MAKE) $(TESTOBJ)
#	$(MAKE) $(TEST)

newdd:
	-rm -f newdd.exe
	$(MAKE) $(NEWDDOBJ)
	$(MAKE) $(NEWDD)


$(LIBS): $(TIMEROBJ)
	ar -rcv $@ $?
	ranlib $@


OSTLIB = -L. -lmyram  -L/home/kjhan/local/lib -lm

%.c: 
	$(CC) $(CFLAGS) -c $*.c


clean:
	-rm -f *.o *.a *.d core  *.exe


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c 

#rules to build binary from source


.c.o : 
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.c.exe :
#	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -o $*.exe $(TESTOBJ) $(OSTLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -o $*.exe $(NEWDDOBJ) $(OSTLIB)


