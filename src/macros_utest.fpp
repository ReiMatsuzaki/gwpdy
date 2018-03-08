#define EXPECT_NEAR_D(a,b,eps,ierr) call utest_near_d(a,b,eps,__FILE__,__LINE__,ierr)
#define EXPECT_NEAR_C(a,b,eps,ierr) call utest_near_c(a,b,eps,__FILE__,__LINE__,ierr)
#define EXPECT_EQ_D(a,b,ierr) call utest_near_d(a,b,1.0d-15,__FILE__,__LINE__,ierr)
#define EXPECT_EQ_C(a,b,ierr) call utest_near_c(a,b,1.0d-15,__FILE__,__LINE__,ierr)
#define EXPECT_EQ_I(a,b,ierr) call utest_eq_i(a,b,__FILE__,__LINE__,ierr)
#define EXPECT_EQ_S(a,b,ierr) call utest_eq_s(a,b,__FILE__,__LINE__,ierr)


