all: residue_test.out

PATH_ROOT=$(INVERSE_ODE)
PATH_TEST=$(PATH_ROOT)/test/linear_inverse_eigen/ResidualSumOfSquares
PATH_LINEAR_INVERSE_EIGEN=$(PATH_ROOT)/solver/linear_inverse

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

residue_test.out: \
  linear_inverse_eigen_mod.o \
  residue_test.o
	$(COMPILER) \
	  linear_inverse_eigen_mod.o \
	  residue_test.o \
	  -llapack \
	  -o residue_test.out

linear_inverse_eigen_mod.o: \
  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE_EIGEN)/linear_inverse_eigen_mod.f90 \
	  -o linear_inverse_eigen_mod.o
residue_test.o: $(PATH_TEST)/residue_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/residue_test.f90 \
	  -o residue_test.o

clear:
	rm *.o *.out *.mod
