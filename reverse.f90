subroutine reverse(A,AT,N)
implicit none 
integer n,m,info,k,i,j
integer ipiv(n)
double precision A(n,n)
double precision work(n)
double precision AT(n,n)
AT=A
call dgetrf( n, n,AT , n,ipiv,info)
call dgetri( n, AT, n,ipiv, work, n,info)
end

subroutine MKL_AL_AR(A,AR,n)
include 'mkl.fi' 
use   lapack95
!use mkl95_lapack
!!use mymkl
!use mkl95_precision

	implicit none
	integer n,ilo,ihi,info
	double precision A(n,n),AL(n,n),AR(n,n),A1(n,n),AL1(n,n),AR1(n,n),abnrm
double precision scalev(n),wr(n),wi(n),work(n),iwork(n),rconde(n),rcondv(n)

A1=CMPLX(A)
!call geevx(a, w[,vl] [,vr] [,balanc] [,ilo] [,ihi] [,scale] [,abnrm] [,rconde] [,rcondv] [,info])
!call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
!call dgeev('V', 'V', 4, A1, 4,wr,wi, AL1,4, AR1, 4,work, 16, info)   !,ilo,ihi,scalev,abnrm, rconde,rcondv, info)
!call geev( A1, wr, AL,AR,info)   !,ilo,ihi,scalev,abnrm, rconde,rcondv, info)
call	 dgeevx('P', 'V', 'V', 'V', 4, A1, 4, wr, wi, Al, 4, AR, 4, & 
			ilo, ihi, scalev, abnrm, rconde, rcondv, work, 40,iwork, info)
write(*,*) info
end



