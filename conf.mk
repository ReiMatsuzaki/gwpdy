# -- Directories --
SRC=${PROJ_ROOT}/src
BUILD=${PROJ_ROOT}/build
VPATH=${BUILD}:${SRC}

# -- Common --
ifeq (${FC},gfortran)
	FF=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace -g -cpp -ffree-line-length-512 -fopenmp
	LDFLAGS=-llapack -lblas
endif
INCLUDE=-I${BUILD} -I${SRC}
MODS0=const sys math timer strutil

# -- utility function
mod2obj = $(addprefix ${BUILD}/, $(addsuffix .o, $(1)))

# -- compile --
%.x:
	${FC} ${FF} $^ -o $@ -cpp ${LIBS} ${LDFLAGS} ${INCLUDE}

${BUILD}/%.o: ${SRC}/%.f90
	@if [ ! -d ${BUILD} ]; \
	   then mkdir -p ${BUILD}; \
	fi
	cd ${BUILD}; ${FC} ${FF} ${INCLUDE} -c $< -o $@

# -- objects --
OBJ_MONO_DY=$(call mod2obj, ${MODS0} spline mole_fit dy_mono)
OBJ_BRANCH_DY=$(call mod2obj, ${MODS0} gwp spline mole_fit dy_branch)
OBJ_BRANCH2_DY=$(call mod2obj, ${MODS0} pwgto spline mole_fit dy_branch2)
