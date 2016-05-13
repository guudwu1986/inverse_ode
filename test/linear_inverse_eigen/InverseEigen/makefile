all: inverse_eigen_test.out

PATH_ROOT=$(INVERSE_ODE)
PATH_NEWUOA=$(PATH_ROOT)/solver/newuoa_dp
PATH_RKF45=$(PATH_ROOT)/solver/rkf45
PATH_LINEAR_INVERSE=$(PATH_ROOT)/solver/linear_inverse
PATH_TEST=$(PATH_ROOT)/test/linear_inverse_eigen/InverseEigen

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

inverse_eigen_test.out: \
  bigden.o biglag.o newuoa.o newuob.o trsapp.o update.o \
  rkf45.o machine.o \
  linear_inverse_eigen_mod.o \
  inverse_eigen_test.o
	$(COMPILER) \
	  bigden.o biglag.o newuoa.o newuob.o \
	  trsapp.o update.o \
	  rkf45.o machine.o \
	  linear_inverse_eigen_mod.o \
	  inverse_eigen_test.o \
	  -llapack \
	  -o inverse_eigen_test.out

bigden.o: $(PATH_NEWUOA)/bigden.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/bigden.f \
	  -o bigden.o
biglag.o: $(PATH_NEWUOA)/biglag.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/biglag.f \
	  -o biglag.o
newuoa.o: $(PATH_NEWUOA)/newuoa.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/newuoa.f \
	  -o newuoa.o
newuob.o: $(PATH_NEWUOA)/newuob.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/newuob.f \
	  -o newuob.o
trsapp.o: $(PATH_NEWUOA)/trsapp.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/trsapp.f \
	  -o trsapp.o
update.o: $(PATH_NEWUOA)/update.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/update.f \
	  -o update.o
rkf45.o: $(PATH_RKF45)/rkf45.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/rkf45.f \
	  -o rkf45.o
machine.o: $(PATH_RKF45)/machine.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/machine.f \
	  -o machine.o
linear_inverse_eigen_mod.o: \
  $(PATH_LINEAR_INVERSE)/linear_inverse_eigen_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE)/linear_inverse_eigen_mod.f90 \
	  -o linear_inverse_eigen_mod.o
inverse_eigen_test.o: $(PATH_TEST)/inverse_eigen_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/inverse_eigen_test.f90 \
	  -o inverse_eigen_test.o

clear:
	rm *.o *.out *.mod
