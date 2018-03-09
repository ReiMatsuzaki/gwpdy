include ../../local.mk
include ../../conf.mk

TARGET=mono_harm1d.x

${TARGET}: ${OBJ_MONO_DY} ${TARGET:%.x=%.f90}
all: ${TARGET}
	./$< > out.log
	python plot.py
	imgcat fig/r.pdf

