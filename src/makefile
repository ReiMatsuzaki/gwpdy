include ../local.mk

# ==== Common ====
ifeq (${FC},gfortran)
	FF=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace -g -cpp -ffree-line-length-512 -fopenmp
	LDFLAGS=-llapack -lblas
endif
MODS0=const math

# -- Directories --
SRC=${PROJ_ROOT}/src
BUILD=${PROJ_ROOT}/build
TESTS=${PROJ_ROOT}/tests

# -- compile --
${BUILD}/%.x:
	@if [ ! -d ${BUILD} ]; \
	   then mkdir -p ${BUILD}; \
	    fi
	${FC} ${FF} $^ -o $@ -cpp ${LIBS} ${LDFLAGS}


${BUILD}/%.o: ${SRC}/%.f90
	@if [ ! -d ${BUILD} ]; \
	   then mkdir -p ${BUILD}; \
	    fi
	${FC} ${FF} ${INCLUDE} -c $< -o $@


# -- Clean --
clean:
	rm -f ${BUILD}/*.o
	rm -f ${BUILD}/*.mod
	rm -f ${BUILD}/*.x

# -- utility function
mod2obj = $(addprefix ${BUILD}/, $(addsuffix .o, $(1)))

# ==== utest(gwp) ====
${BUILD}/utest_gwp.x: $(call mod2obj, ${MODS0} gwp utest_gwp)
check_gwp: ${BUILD}/utest_gwp.x
	$<

${BUILD}/utest_dy.x: $(call mod2obj, ${MODS0} gwp dy_mono utest_dy)
check_dy: ${BUILD}/utest_dy.x
	$<
