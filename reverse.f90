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
