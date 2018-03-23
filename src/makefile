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
${BUILD}/utest_math.x: $(call mod2obj, ${MODS0} utest utest_math)
${BUILD}/utest_strutil.x: $(call mod2obj, ${MODS0} utest utest_strutil)
${BUILD}/utest_spline.x: $(call mod2obj, ${MODS0} spline utest utest_spline)
${BUILD}/utest_gwp.x: $(call mod2obj, ${MODS0} gwp utest utest_gwp)
${BUILD}/utest_pwgto.x: $(call mod2obj, ${MODS0} gwp pwgto utest utest_pwgto)
${BUILD}/utest_dy.x: $(call mod2obj, ${MODS0} spline mole_fit gwp dy_mono dy_branch utest utest_dy)
${BUILD}/utest_dyb2.x: $(call mod2obj, ${MODS0} spline mole_fit pwgto dy_branch2 utest utest_dyb2)
${BUILD}/utest_cl.x: $(call mod2obj, ${MODS0}  cl_psanb utest utest_cl)

check: check_strutil check_spline check_gwp check_dy


