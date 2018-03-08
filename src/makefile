include ../local.mk
include ../conf.mk

# -- Clean --
clean:
	rm -f ${BUILD}/*.o
	rm -f ${BUILD}/*.mod
	rm -f *.mod
	rm -f ${BUILD}/*.x

# -- quick lunch  --
check_%: ${BUILD}/utest_%.x
	$<

# ==== exe =====
${BUILD}/utest_strutil.x: $(call mod2obj, ${MODS0} utest utest_strutil)
${BUILD}/utest_spline.x: $(call mod2obj, ${MODS0} spline utest utest_spline)
${BUILD}/utest_gwp.x: $(call mod2obj, ${MODS0} gwp utest_gwp)
${BUILD}/utest_dy.x: $(call mod2obj, ${MODS0} gwp dy_mono utest utest_dy)

# -- main --
${BUILD}/gwpdy.x: ${BUILD}/gwpdy.a
	${FF} ${FF} $^ -o $@ -cpp
check_main: ${BUILD}/gwpdy.x 
	$< -nf 1 -ne 1


