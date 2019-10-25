module matrix
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

contains
  ! invert matrix
  subroutine invert(a,ainv,numNodes) !switch to QR decomp
    implicit none

    integer :: p,i,j,l,l2,numNodes
    real(dp) :: factor
    real(dp), dimension (:,:), intent (in) :: a
    real(dp), dimension (:,:), intent (out) :: ainv
    real(dp), dimension (numNodes,numNodes*2) :: aa
    aa(1:,1:) = 0
    l = numNodes
    l2 = 2*numNodes

  ! set up (n)x(2n) matrix to invert
    do i = 1,l,1
      do j = 1,l,1
        aa(i,j) = a(i,j)
      end do
        aa(i,i+numNodes) = 1
    end do

  ! p is pivot point
    do p = 1,l,1
      ! make pivot row start with 1
      do j = 1,l2,1
        if(j/=p)then
          aa(p,j) = aa(p,j)/aa(p,p)
        end if
      end do !
      aa(p,p)= aa(p,p)/aa(p,p)
      ! make zeros in pivot column
      do i = 1,l,1
        if(i/=p)then
          factor = -aa(i,p)/aa(p,p)
          do j = 1,l2,1
            aa(i,j) = aa(i,j)+(factor*aa(p,j))
          end do !j
        end if
      end do !i
    end do !p
  ! separate out n x n inverse matrix
  do i = 1,l,1
    do j = 1,l,1
      ainv(i,j) = aa(i,j+l)
    end do
  end do
  end subroutine invert

  function det(A,n) ! Through LU decomposition (dont find L because not needed for det)
    implicit none
    integer :: p,i,j,n
    real(dp), dimension (n,n), intent (in) :: A
    real(dp), dimension (n,n) :: U
    real(dp) :: det,factor
    U = A

    ! this prevents errors from having points at zero
    do j = 1,n,1
      U(3,j) = U(3,j)+U(2,j)
      U(3,j) = U(3,j)+U(1,j)
      U(1,j) = U(1,j)+U(3,j)
      U(2,j) = U(2,j)+U(3,j)
    end do !j

    ! form upper triangle matrix
    do p = 1,n,1
      do i = p+1,n,1
        factor = -U(i,p)/U(p,p)
        do j = 1,n,1
          U(i,j) = U(i,j)+(factor*U(p,j))
        end do !j
      end do !i
    end do !p
    ! prod diagonal
    det = 1
    do i = 1,n,1
      det = det * U(i,i)
    end do !i
  end function det

  subroutine printArray(A)
    implicit none

    integer :: i
    real(dp), dimension (:,:), intent (in) :: A
    integer, dimension(2) :: limits
    limits = shape(A)
    ! print each row 1 at a time
    do i = 1,limits(1),1
      print *,A(i,:)
    end do
    print *,

  end subroutine printArray

end module matrix

module finiteElement
  use matrix
  implicit none

contains
  subroutine assembleMatrix(x,n,numNodes,numEls,a,Q)
    implicit none

    real(dp), intent (in) :: a  !thermal conductivity
    real(dp), intent (in) :: Q
    integer, intent (in) :: numNodes, numEls
    integer :: el,i,j
    integer, dimension(:,:), intent (in) :: n
    real(dp), dimension(:,:), intent (in) :: x
    real(dp), dimension(numNodes,numNodes) :: gMatrix
    real(dp), dimension(3,3) :: p = 1
    real(dp), dimension(3,3) :: element = 1
    real(dp) :: determinate

    do el = 1,numEls,1
      ! Define P matrix of Interest
      p(1,1:2) = x(n(el,1),:)
      p(2,1:2) = x(n(el,2),:)
      p(3,1:2) = x(n(el,3),:)

      ! Assemble Element matrix
      element(1,1) = (p(2,2)-p(3,2))*(p(2,2)-p(3,2)) + (p(3,1)-p(2,1))*(p(3,1)-p(2,1))
      element(1,2) = (p(2,2)-p(3,2))*(p(3,2)-p(1,2)) + (p(3,1)-p(2,1))*(p(1,1)-p(3,1))
      element(1,3) = (p(2,2)-p(3,2))*(p(1,2)-p(2,2)) + (p(3,1)-p(2,1))*(p(2,1)-p(1,1))
      element(2,1) = (p(3,2)-p(1,2))*(p(2,2)-p(3,2)) + (p(1,1)-p(3,1))*(p(3,1)-p(2,1))
      element(2,2) = (p(3,2)-p(1,2))*(p(3,2)-p(1,2)) + (p(1,1)-p(3,1))*(p(1,1)-p(3,1))
      element(2,3) = (p(3,2)-p(1,2))*(p(1,2)-p(2,2)) + (p(1,1)-p(3,1))*(p(2,1)-p(1,1))
      element(3,1) = (p(1,2)-p(2,2))*(p(2,2)-p(3,2)) + (p(2,1)-p(1,1))*(p(3,1)-p(2,1))
      element(3,2) = (p(1,2)-p(2,2))*(p(3,2)-p(1,2)) + (p(2,1)-p(1,1))*(p(1,1)-p(3,1))
      element(3,3) = (p(1,2)-p(2,2))*(p(1,2)-p(2,2)) + (p(2,1)-p(1,1))*(p(2,1)-p(1,1))
      element = -a*element/(2* det(p,3));

      ! build Global matrix
      do i = 1,3,1
        do j = 1,3,1
          gMatrix(n(el,i),n(el,j)) = gMatrix(n(el,i),n(el,j)) + element(i,j)
        end do !j
      end do !i
    end do !el

  end subroutine assembleMatrix



end module finiteElement


program main
  use matrix
  use finiteElement
  implicit none

  real(dp), parameter :: a = 1. !thermal conductivity
  real(dp), parameter :: Q = 1. ! known constant
  integer, dimension(3,3) :: n
  real(dp), dimension(4,2) :: x
  integer :: numNodes, numEls
  integer, dimension(2) :: temp

  n(1,:) = (/1,2,3/)
  n(2,:) = (/2,4,3/)
  n(3,:) = (/1,3,4/)
  x(:,1) = (/1.315,0.,1.315,2.63/)
  x(:,2) = (/2.27,0.,0.76,0./)
  temp = shape(n)
  numEls = temp(1)
  temp = shape(x)
  numNodes = temp(1)
  call assembleMatrix(x,n,numNodes,numEls,a,Q)
  ! integer ::numNodes = 3
  ! integer, dimension(2)::i
  ! real(dp), dimension (3,3) :: a
  ! real(dp), dimension (3,3) :: ainv
  ! a(1,:) = (/1, -1, 3/)
  ! a(2,:) = (/2, 1, 2/)
  ! a(3,:) = (/-2, -2, 1/)
  ! i  = shape(a)
  ! numNodes = i(1)
  ! ! print *, numNodes
  !
  ! call invert(a,ainv,numNodes)
  ! call printArray(a)

end program main
