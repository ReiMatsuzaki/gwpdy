include ../../../local.mk
include ../../../conf.mk

TARGET=run.x

${TARGET}: ${OBJ_MONO_DY} run.f90
all: ${TARGET}
	python pot.py
	./$< > run.log
	python plot.py
	imgcat fig/prob.pdf

check: ${TARGET}
	./$<
