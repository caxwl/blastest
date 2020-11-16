module ISO_C_BINDING
  ! Dummy module providing sizes of datatypes
  integer, parameter :: C_INT=4
  integer, parameter :: C_DOUBLE=8
end module ISO_C_BINDING 


! Computes dot product using BLAS routine
function dot_blas(n,x1,x2) result(d)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: n
  real(C_DOUBLE), intent(in), dimension(n) :: x1, x2
  real(C_DOUBLE) :: d, ddot
  integer, parameter :: inc=1

  d=ddot(n,x1,inc,x2,inc)
end function dot_blas

! Computes dot product using for loop
function dot_for(n,x1,x2) result(d)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: n
  real(C_DOUBLE), intent(in), dimension(n) :: x1, x2
  real(C_DOUBLE) :: d
  integer(C_INT) :: i

  ! Compute dot product with for loop
  d=0.0
  do i=1,n
     d=d+x1(i)*x2(i)
  end do
end function dot_for

! Computes dot product using intrinsic Fortran function
function dot_intrinsic(n,x1,x2) result(d)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: n
  real(C_DOUBLE), intent(in), dimension(n) :: x1, x2
  real(C_DOUBLE) :: d

  d=dot_product(x1,x2)
end function dot_intrinsic

! Computes matrix-vector product y=A^T x using BLAS routine
subroutine matvec_blas(m,n,mat,vec,res)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n
  real(C_DOUBLE), intent(in), dimension(m,n) :: mat
  real(C_DOUBLE), intent(in), dimension(n) :: vec
  real(C_DOUBLE), intent(out), dimension(m) :: res
  integer, parameter:: inc=1

  ! SLOW!
!  call dgemv('n',m,n,1.0,mat,m,vec,inc,0.0,res,inc)
  ! OK
  call dgemv('t',m,n,1.0,mat,m,vec,inc,0.0,res,inc)
end subroutine matvec_blas

! Computes matrix-vector product using for loop
subroutine matvec_for(m,n,mat,vec,res)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n
  real(C_DOUBLE), intent(in), dimension(m,n) :: mat
  real(C_DOUBLE), intent(in), dimension(n) :: vec
  real(C_DOUBLE), intent(out), dimension(m) :: res
  integer, parameter:: inc=1
  integer(C_INT) :: i, j

  ! Initialize result vector
  res=0.0

  ! Perform product, remember correct memory order
  do i=1,n
     do j=1,m
        res(i)=res(i)+mat(j,i)*vec(j)
     end do
  end do
end subroutine matvec_for

! Computes matrix-vector product using intrinsic function
subroutine matvec_intrinsic(m,n,mat,vec,res)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n
  real(C_DOUBLE), intent(in), dimension(m,n) :: mat
  real(C_DOUBLE), intent(in), dimension(n) :: vec
  real(C_DOUBLE), intent(out), dimension(m) :: res

  ! Initialize result vector
  res=matmul(transpose(mat),vec)
end subroutine matvec_intrinsic


! Computes matrix-matrix product C=A B using BLAS routine
subroutine matmat_blas(m,n,k,A,B,C)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n, k
  real(C_DOUBLE), intent(in), dimension(m,k) :: A
  real(C_DOUBLE), intent(in), dimension(k,n) :: B
  real(C_DOUBLE), intent(out), dimension(m,n) :: C
  integer, parameter:: inc=1

  ! 42.6 sec for 1000x1000
  call dgemm('n','n',m,n,k,1.0,a,m,b,k,0.0,c,m)

  ! 42.1 sec
!  call dgemm('n','t',m,n,k,1.0,a,m,b,k,0.0,c,m)

  ! 42.4 sec
!  call dgemm('t','n',m,n,k,1.0,a,m,b,k,0.0,c,m)

  ! 42.8 sec
!  call dgemm('t','t',m,n,k,1.0,a,m,b,k,0.0,c,m)
end subroutine matmat_blas

! Computes matrix-matrix product using for loop
subroutine matmat_for(m,n,k,A,B,C)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n, k
  real(C_DOUBLE), intent(in), dimension(m,k) :: A
  real(C_DOUBLE), intent(in), dimension(k,n) :: B
  real(C_DOUBLE), intent(out), dimension(m,n) :: C
  integer, parameter:: inc=1

  integer(C_INT) :: im, in, ik

  c=0.0
  do im=1,m
     do in=1,n
        do ik=1,k
           c(im,in)=c(im,in)+a(im,ik)*b(ik,in)
        end do
     end do
  end do
end subroutine matmat_for

! Computes matrix-matrix product using intrinsic function
subroutine matmat_intrinsic(m,n,k,A,B,C)
  use ISO_C_BINDING
  implicit none
  integer(C_INT), intent(in) :: m, n, k
  real(C_DOUBLE), intent(in), dimension(m,k) :: A
  real(C_DOUBLE), intent(in), dimension(k,n) :: B
  real(C_DOUBLE), intent(out), dimension(m,n) :: C
  integer, parameter:: inc=1

  c=matmul(a,b)
end subroutine matmat_intrinsic
