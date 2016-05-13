all: construct_test.out

PATH_ROOT=$(INVERSE_ODE)
PATH_NEWUOA=$(PATH_ROOT)/solver/newuoa_dp
PATH_TEST=$(PATH_ROOT)/test/linear_inverse_eigen/ConstructODE
PATH_LINEAR_INVERSE_EIGEN=$(PATH_ROOT)/solver/linear_inverse

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

construct_test.out: \
  bigden.o biglag.o newuoa.o newuob.o trsapp.o update.o \
  linear_inverse_eigen_mod.o \
  construct_test.o
	$(COMPILER) \
	  bigden.o biglag.o newuoa.o newuob.o \
	  trsapp.o update.o \
	  linear_inverse_eigen_mod.o \
	  construct_test.o \
	  -llapack \
	  -o construct_test.out

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
linear_inverse_eigen_mod.o: \
  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90 \
	  -o linear_inverse_eigen_mod.o
construct_test.o: $(PATH_TEST)/construct_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/construct_test.f90 \
	  -o construct_test.o

clear:
	rm *.o *.out *.mod
