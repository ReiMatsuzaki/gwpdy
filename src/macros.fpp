
#define begin_err(msg) write(0, '(A, ":", I0, ": ", A)') __FILE__, __LINE__, msg
#define check_err(ierr) if(ierr .ne. 0) begin_err(""); if(ierr .ne. 0) return

