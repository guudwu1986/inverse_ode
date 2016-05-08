all: rkf45_test.out

PATH_ROOT=$(INVERSE_ODE)
PATH_RKF45=$(PATH_ROOT)/solver/rkf45
PATH_TEST=$(PATH_ROOT)/test/rkf45

COMPILER=gfortran -fdefault-real-8
STATIC_LIB_COMPILER=$(COMPILER) -c

rkf45_test.out: rkf45.o machine.o rkf45_test.o
	$(COMPILER) rkf45.o machine.o rkf45_test.o \
	  -o rkf45_test.out

rkf45.o: $(PATH_RKF45)/rkf45.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/rkf45.f \
	  -o rkf45.o
machine.o: $(PATH_RKF45)/machine.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/machine.f \
	  -o machine.o
rkf45_test.o: $(PATH_TEST)/rkf45_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/rkf45_test.f90 \
	  -o rkf45_test.o

clear:
	rm *.o *.out
