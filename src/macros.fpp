
#define MSG_ERR(msg) write(0, '(A, ":", I0, ": ", A)') __FILE__, __LINE__, msg
#define CHK_ERR(ierr) if(ierr .ne. 0) MSG_ERR(""); if(ierr .ne. 0) return

