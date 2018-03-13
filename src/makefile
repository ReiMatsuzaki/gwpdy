include ../local.mk
include ../conf.mk

# -- Clean --
clean:
	rm -f ${BUILD}/*.o
	rm -f ${BUILD}/*.mod
	rm -f ${BUILD}/*.x

# -- quick lunch  --
check_%: ${BUILD}/utest_%.x
	$<

# ==== exe =====
${BUILD}/utest_strutil.x: $(call mod2obj, ${MODS0} utest utest_strutil)
${BUILD}/utest_spline.x: $(call mod2obj, ${MODS0} spline utest utest_spline)
${BUILD}/utest_gwp.x: $(call mod2obj, ${MODS0} gwp utest utest_gwp)
${BUILD}/utest_dy.x: $(call mod2obj, ${MODS0} spline mole_fit gwp dy_mono dy_branch utest utest_dy)

check: check_strutil check_spline check_gwp check_dy


