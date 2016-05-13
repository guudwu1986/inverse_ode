all: newuoa_objective_test.out

PATH_ROOT=$(INVERSE_ODE)
PATH_TEST=$(PATH_ROOT)/test/linear_inverse_eigen/TestResidueNewuoaObjective
PATH_LINEAR_INVERSE_EIGEN=$(PATH_ROOT)/solver/linear_inverse

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

newuoa_objective_test.out: \
  linear_inverse_eigen_mod.o \
  newuoa_objective_test.o
	$(COMPILER) \
	  linear_inverse_eigen_mod.o \
	  newuoa_objective_test.o \
	  -llapack \
	  -o newuoa_objective_test.out

linear_inverse_eigen_mod.o: \
  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90 \
	  -o linear_inverse_eigen_mod.o
newuoa_objective_test.o: $(PATH_TEST)/newuoa_objective_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/newuoa_objective_test.f90 \
	  -o newuoa_objective_test.o

clear:
	rm *.o *.out *.mod
