include ../../local.mk
include ../../conf.mk

TARGET=run.x

${TARGET}: ${OBJ_MONO_DY} ${TARGET:%.x=%.f90}
all: ${TARGET}
	./$< > out.log
	python plot.py
	imgcat fig/prob.pdf
