#FFLAGS=-O3 -fopenmp
# FFLAGS=-fbounds-check -fbacktrace -g
FFLAGS=-mcmodel=medium -O3 -fopenmp

F90SRCS =           \
bndset.f90    		\
cal_adv.f90			\
cal_coef.f90		\
cal_dt.f90			\
cal_rs.f90			\
cal_sh.f90  		\
cal_thetas.f90		\
cal_vel.f90			\
cal_xs_ys.f90 		\
sor2.f90			\
main.f90            \

SRCS  =  $(F90SRCS)

.SUFFIXES: .o .f .f90
# FCOBJS = $(FCSRCS:.f=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  =  $(F90OBJS)

.f90.o:
	gfortran -c $(FFLAGS) $<

.f.o:
	gfortran -c $(FFLAGS) $<

a.out: $(OBJS)
	gfortran $(OBJS) $(FFLAGS) -o a.out

clean:
	rm $(OBJS)
	rm a.out
