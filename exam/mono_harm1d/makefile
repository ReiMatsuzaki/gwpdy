include ../../local.mk
include ../../conf.mk
VPATH=${BUILD}:${SRC}

TARGET=mono_harm1d.x

${TARGET}: $(call mod2obj, ${MODS0} gwp dy_mono) mono_harm1d.f90
all: ${TARGET}
	./$< > out.log
	python plot.py
	imgcat fig/r.pdf

