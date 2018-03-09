include ../../local.mk
include ../../conf.mk

TARGET=run.x

${TARGET}: ${OBJ_BRANCH_DY} run.f90
all: ${TARGET}
	python pot.py
	./$< > run.log
	python plot.py
	imgcat fig/prob.pdf

run: ${TARGET}
	./$<

clean:
	rm ${TARGET}
	rm -rf fig
	rm -rf out
	rm *.csv
	rm *.log
