include ../local.mk
include ../conf.mk

# -- Clean --
clean:
	rm -f ${BUILD}/*.o
	rm -f ${BUILD}/*.mod
	rm -f *.mod
	rm -f ${BUILD}/*.x

# ==== exe =====
# -- utest(strutil) --
${BUILD}/utest_strutil.x: $(call mod2obj, ${MODS0} utest utest_strutil)
check_strutil: ${BUILD}/utest_strutil.x
	$<

# -- utest(gwp) --
${BUILD}/utest_gwp.x: $(call mod2obj, ${MODS0} gwp utest_gwp)
check_gwp: ${BUILD}/utest_gwp.x
	$<

# -- utest(dy) --
${BUILD}/utest_dy.x: $(call mod2obj, ${MODS0} gwp dy_mono utest utest_dy)
check_dy: ${BUILD}/utest_dy.x
	$<

# -- example(Harmonic oscillator) --


# -- main --
${BUILD}/gwpdy.x: ${BUILD}/gwpdy.a
	${FF} ${FF} $^ -o $@ -cpp
check_main: ${BUILD}/gwpdy.x 
	$< -nf 1 -ne 1


