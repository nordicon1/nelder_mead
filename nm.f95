! http://www.jasoncantarella.com/downloads/NelderMeadProof.pdf !
!-------------------------------------------
!             simple Nelder-Mead           !
!-------------------------------------------


!-------------------------------------------
! Module : commons
!-------------------------------------------
module commons

  implicit none
contains
  !---------------------------------------
  ! Name : get_random
  !---------------------------------------
  function get_random(MN, MX) result(rand)
    implicit none
    real(8),intent(in) :: MN
    real(8),intent(in) :: MX
    real(8) rand
    real(8),save:: area
    area = abs(MX-MN)

    call random_seed() 
    call random_number(rand)
    rand = (rand * area) - (area / 2.0d0)
  end function get_random

  !---------------------------------------
  ! Name : linespace
  !---------------------------------------
  function linspace(MN, MX, STEP) result(y)
    implicit none
    ! Argument
    real(8), intent(in) :: MN
    real(8), intent(in) :: MX
    integer, intent(in) :: STEP
    ! Output
    real(8) :: y(STEP)
    ! Working memory
    integer i
    real(8),save :: width_pix
    real(8),save :: center_pix
    real(8),save :: org_pix_x   ! x-origin of the pixel
                                ! (at left-edge of pixel)

    width_pix = (MX - MN) / STEP
    center_pix = width_pix / 2.0d0

    do i = 1, size(y)
      org_pix_x = mn + width_pix * (i - 1)
      y(i) = org_pix_x + center_pix   ! x-position
    end do

    return
  end function linspace 
end module commons


program main
    use commons
    implicit none
    ! input
    real(8),allocatable :: x0(:)
    real(8),allocatable :: x1(:)
    ! output
    real(8),allocatable :: y(:,:)
    !
    ! Constant value
    !
    !  Function search area
    integer, parameter :: MAX_ITER = 10
    integer,parameter :: select_obj = 3  
    logical,parameter :: RANDOM_ON = .false.
    real(8),parameter :: BND_X0_MIN = 0.0d0 ! Search area
    real(8),parameter :: BND_X0_MAX = 5.0d0 ! search area
    real(8),parameter :: BND_X1_MIN = 0.0d0 ! search area
    real(8),parameter :: BND_X1_MAX = 5.0d0 ! search area
    integer,parameter :: NUM_PIX_X0  = 255 
    integer,parameter :: NUM_PIX_X1  = 255
    
    !logical,parameter :: ALLOW_CROSS_BOUNDARY = .true.
    logical,parameter :: ALLOW_CROSS_BOUNDARY = .false.
    character(8),parameter :: GOAL = 'MINIMIZE' 
    !character(8),parameter :: GOAL = 'MAXIMIZE' 
    
    ! Working variable
    integer,save :: i,j
    real(8),save :: rand_value
    real(8),save :: model(2)
    real(8),save :: V(3,2)    ! Three Vertex  
    
    !-----------------------------------------------
    ! Check Goal Setting
    !-----------------------------------------------
    If (GOAL == 'MINIMIZE') Then
      write(*,*) 'Goal: ', 'MINIMIZE'

    Else If (GOAL == 'MAXIMIZE') Then
      write(*,*) 'Goal: ', 'MAXIMIZE'
    
    Else
      write(*,*) 'Error: Goal=', GOAL, ' is invalid'
      return
    
    End if
    
    
    !-----------------------------------------------
    ! Make Model
    !-----------------------------------------------
    
    x0 = linspace(BND_X0_MIN, BND_X0_MAX, NUM_PIX_X0)
    x1 = linspace(BND_X1_MIN, BND_X1_MAX, NUM_PIX_X1)
    
    allocate(y(size(x0),size(x1)))
    
    open(10,file='f0.csv')
    
    Do i = 1, size(x0)
      Do j = 1, size(x1)
        model(1) = x0(i)
        model(2) = x1(j)
        y(i,j) = f(model)
        write(10,*) x0(i), x1(j), y(i,j)
        !write(*,*)  x0(i), x1(j), y(i,j)
      End Do
      ! brank line for GNUPLOT 3D plot
      write(10,*) ''
    End do
    
    close(10)

    deallocate(y)
    
    !-----------------------------------------------
    ! Optimization
    !-----------------------------------------------
    call NelderMead(          &  ! Melder-Mead
      BND_X0_MIN, BND_X0_MAX, &  !  X0:Min,Max
      BND_X1_MIN, BND_X1_MAX  &  !  X1:Min,Max
    )                            !


contains
  !--------------------------------------------
  ! Objective
  !--------------------------------------------
  function objective_1(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = x0 ** 2 + x1 ** 2
  end function objective_1
  

  function objective_2(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = 2.0d0*(x0**2) -1.05d0*(x0**4)+((x0**6)/6.d0)+(x0*x1)+(x1**2)
  end function objective_2


  !----------------------------------------------------------------
  ! Nelder-Mead BenchMark
  !  0 < x0,x1 < 5
  !----------------------------------------------------------------
  function objective_3(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = (x0**2) - (4.0d0 * x0) + (x1**2) - x1 - (x0*x1)
  end function objective_3


  !----------------------------------------------------------------
  ! sphere
  !  -2 < x0,x1 < 2
  !----------------------------------------------------------------
  function objective_4(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = (x0 ** 2) + (x1 ** 2)
  end function objective_4

  
  !----------------------------------------------------------------
  ! LED
  !  0 < x0,x1 < 1
  !----------------------------------------------------------------
  function objective_5(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = (x0-(3.0d0/4.0d0))**2 + (x1/100d0)
  end function objective_5


  !----------------------------------------------------------------
  ! Goldstein-price
  !  -3 < x,y < 3
  !----------------------------------------------------------------
  function objective_6(x, y) result(u)
    implicit none
    real(8), intent(in) :: x   ! input
    real(8), intent(in) :: y   ! input
    real(8) u                   ! output
    u = (   ((1.0d0+((x+y+1.0d0))**2))  &
          * (                           &
               19.0d0-(14.0d0*x)        &
              +((3.0d0*x)**2)           &
              -(14.0d0*y)               &
              +(6.0d0*x*y)              &
              +((3.0d0*y)**2)           &
            )                           &
          * (30.0d0                     &
              + ((2.0d0*x-3.0d0*y)**2)  &
              * (                       &
                   (18.0d0-(32.0d0*x))  &
                  +(12.0d0*x**2)        &
                  +(48.d0*y)            &
                  -(36.0d0*x*y)         &
                  +(27.0d0*y**2)        &
                )                       &
            )                           &
        )
  end function objective_6

  !----------------------------------------------------------------
  !  0 < x0,x1 < 1
  !----------------------------------------------------------------
  function objective_7(x0, x1) result(y)
    implicit none
    real(8), intent(in) :: x0   ! input
    real(8), intent(in) :: x1   ! input
    real(8) y                   ! output
    y = (x0*x1)/(x0+x1)
  end function objective_7

  function f(x) result(y)
    implicit none
    real(8), intent(in) :: x(2)   ! input
    real(8) y                     ! output
    real(8),save :: x0
    real(8),save :: x1
    real(8), parameter :: oo = huge(1.0)
    x0 = x(1)
    x1 = x(2)
   

    If (select_obj == 1) Then
      y = objective_1(x0, x1)

    Else If (select_obj == 2) Then
      y = objective_2(x0, x1)
    
    Else If (select_obj == 3) Then
      y = objective_3(x0, x1)
    
    Else If (select_obj == 4) Then
      y = objective_4(x0, x1)
    
    Else If (select_obj == 5) Then
      y = objective_5(x0, x1)
    
    Else If (select_obj == 6) Then
      y = objective_6(x0, x1)

    Else If (select_obj == 7) Then
      y = objective_7(x0, x1)

    End If
    
    
    If (.not.ALLOW_CROSS_BOUNDARY) Then
      If ((BND_X0_MIN > X0 .or. BND_X0_MAX < X0) .or. &
          (BND_X1_MIN > X1 .or. BND_X1_MAX < X1)) Then
        If (GOAL == 'MINIMIZE') Then
          y = oo
        Else if (GOAL == 'MAXIMIZE') Then
          y = -1.0d0 * oo
        End If
      End If
    End If
  end function
  

  subroutine NelderMead(&
    bnd_x0_min, &
    bnd_x0_max, &
    bnd_x1_min, &
    bnd_x1_max  &
    )
    implicit none
    ! Arguments
    real(8), intent(in) :: bnd_x0_min ! bound x min
    real(8), intent(in) :: bnd_x0_max ! bound x max
    real(8), intent(in) :: bnd_x1_min ! bound y min
    real(8), intent(in) :: bnd_x1_max ! bound y max
    ! Parameter
    integer, parameter :: x = 1     ! index 
    integer, parameter :: y = 2     ! index
    ! Working Memory
    integer, save :: k,i,n        ! Loop variable
    real(8),save :: M(2)      ! Midpoint
    real(8),save :: B(2)      ! Best vertex
    real(8),save :: G(2)      ! Good Vertex (next to best)
    real(8),save :: W(2)      ! Worst vertex
    real(8),save :: R(2)      ! Refrection vertex
    real(8),save :: E(2)      ! Extend vertex
    real(8),save :: C(2)      ! Contranction 
    real(8),save :: C1(2)     ! Contranction 
    real(8),save :: C2(2)     ! Contranction
    real(8),save :: S(2)      ! 
    real(8),save :: V(3,2)    ! Three Vertex  
    real(8),save :: TMP(2)    ! 

    !-------------------------------------------------
    ! Initial Triangle BGW
    !  -Let f (x, y) be the function that is to be minimized. 
    !   To start, we are given three vertices
    !    of a triangle: Vk = (xk , yk ), k = 1, 2, 3. 
    !   The function f (x, y) is then evaluated at each
    !    of the three points: zk = f (xk , yk ) for k = 1, 2, 3. 
    !   The subscripts are then reordered
    !    so that z1 ≤ z2 ≤ z3. We use the notation
    !-------------------------------------------------

    If (RANDOM_ON) Then
      V(1,x) = get_random(bnd_x0_min, bnd_x0_max)
      V(1,y) = get_random(bnd_x1_min, bnd_x1_max)
      V(2,x) = get_random(bnd_x0_min, bnd_x0_max)
      V(2,y) = get_random(bnd_x1_min, bnd_x1_max)
      V(3,x) = get_random(bnd_x0_min, bnd_x0_max)
      V(3,y) = get_random(bnd_x1_min, bnd_x1_max)
    Else
      V(1,x) = 0.0d0
      V(1,y) = 0.0d0
      V(2,x) = 1.2d0
      V(2,y) = 0.0d0
      V(3,x) = 0.0d0
      V(3,y) = 0.8d0
    End If

    Do i=1,3
      Do j = 1, 3
        If ( f(V(i,:)) < f(V(j,:)) ) Then
          TMP = V(i,:)
          V(i,:) = V(j,:)
          V(j,:) = TMP
        End If
      End Do
    End Do
    
    B = V(1,:)
    G = V(2,:)
    W = V(3,:)
    
    write(*,*) 'Initial Triangle BGW'
    write(*,*) 'B(x,y)',B(x),B(y), 'value:',f(B)
    write(*,*) 'G(x,y)',G(x),G(y), 'value:',f(G)
    write(*,*) 'W(x,y)',W(x),W(y), 'value:',f(W)
    write(*,*) ''

    !-------------------------------------------------
    ! Midpoint of the Good Side
    !  -The construction process uses the midpoint of 
    !    the line segment joining B and G . It is
    !    found by averaging the coordinates
    !-------------------------------------------------
    M = (B + G) / 2.0d0
    write(*,*) 'Midpoint of the Good Side'
    write(*,*) 'M(x,y)',M(x),M(y)
    write(*,*) ''

    !-------------------------------------------------
    ! Reflection Using the Point R
    !-------------------------------------------------
    R = (2.0d0 * M) - W
    write(*,*) 'Reflection Using the Point R'
    write(*,*) 'R(x,y)',R(x),R(y)
    write(*,*) ''

    !-------------------------------------------------
    ! Expansion Using the Point E
    !-------------------------------------------------
    If ( f(R) < f(W) ) Then
      E = (2.0d0 * R) - M
      write(*,*) 'Expanson using the Point E'
      write(*,*) 'E(x,y)',E(x),E(y)
      write(*,*) ''
    End If

    !-------------------------------------------------
    ! Contraction Using the Point C
    !  -If the function value at R and W are the same,
    !    another point must be tested.
    !   Perhaps the function is smaller at M, 
    !    but we cannot replace W with M because we must
    !    have a triangle.
    !   Consider the two midpoints C1 and C2 of the line
    !    segments WM and MR, respectively.
    !   The point with the smaller functin value is called C,
    !    and the new triangle is BGC.
    !   Note. The choice between C1 and C2 might seem inappropriate
    !    for the two-dimensinal case, but it is important in higher
    !    dimentions.
    !-------------------------------------------------

    !-------------------------------------------------
    ! Shrink toward B
    !  -If the function value at C is not less than at W,
    !    the points G and W must be shrunk toward B.
    !   The point G is replaces with M, and W is replacesd
    !    with S, switch is the midpoint of the line segment
    !    joining B with W
    !-------------------------------------------------

    open(10, file='vertex.csv')
    
    write(10,*) V(1,:), f(V(1,:))
    write(10,*) V(2,:), f(V(2,:))
    write(10,*) V(3,:), f(V(3,:))
    
    

    Do k=2, MAX_ITER


      M = (B + G) / 2.0d0
      R = (2.0d0 * M) - W
      If ( f(R) < f(G) ) Then
        ! case(i) {either reflect or extend}
        If ( f(B) < f(R) ) Then
          ! Refrlect
          !  replace W with R
          W = R
        Else
          ! Extend
          !  Compute E and f(E)
          E = (2.0d0 * R) - M
          If ( f(E) < f(B) ) Then
            ! replace W with E
            W = E
          Else
            ! relpace W with R
            W = R
          End If
        End If
      Else
        ! case(ii) {either contract or shrink}
        If ( f(R) < f(W) ) Then
          ! replace W with R
          W = R
        End If
        ! Compute C1 = (W + M)/2 or 
        !         C2 = (M + R)/2 and
        !         f(C)
        C1 = (W + M)/2.0d0
        C2 = (M + R)/2.0d0
        If ( f(C1) < f(C2) ) Then
          C = C1
        Else
          C = C2
        End If

        If ( f(C) < f(W) ) Then
          ! replace W with C
          W = C
        Else
          ! Compute S and f(S)
          S = (B + W)/2.0d0
          ! replace W with S
          W = S
          ! replace G with M
          G = M
        End If
      End If

     
      !
      ! Evaluate Three Vertex,
      ! Set B, G, E
      !

      V(1,:) = B
      V(2,:) = G
      V(3,:) = W
      If (GOAL == 'MINIMIZE') Then
        Do i=1,3
          Do j = 1, 3
            If ( f(V(i,:)) < f(V(j,:)) ) Then
              TMP = V(i,:)
              V(i,:) = V(j,:)
              V(j,:) = TMP
            End If
          End Do
        End Do
      
      Else If (GOAL == 'MAXIMIZE') Then
         Do i=1,3
          Do j = 1, 3
            If ( f(V(i,:)) > f(V(j,:)) ) Then
              TMP = V(i,:)
              V(i,:) = V(j,:)
              V(j,:) = TMP
            End If
          End Do
        End Do
      
      End if
      B = V(1,:)
      G = V(2,:)
      W = V(3,:)

      ! Print result 
      write(*,*) '---------------------------&
                  ---------------------------'
      write(*,*) 'times: ', k

      write(*,*) 'Best point'
      write(*,*) 'B(x,y)',B(x),B(y), f(B)
      write(*,*) ''
      
      write(*,*) 'Good point'
      write(*,*) 'G(x,y)',G(x),G(y), f(G)
      write(*,*) ''

      write(*,*) 'Worst point'
      write(*,*) 'W(x,y)',W(x),W(y), f(W)
      write(*,*) ''


      write(*,*) 'Midpoint of the Good Side'
      write(*,*) 'M(x,y)',M(x),M(y)
      write(*,*) ''

      ! Reflection Using the Point R
      write(*,*) 'Reflection Using the Point R'
      write(*,*) 'R(x,y)',R(x),R(y)
      write(*,*) ''

      ! Expansion Using the Point E
      write(*,*) 'Expanson using the Point E'
      write(*,*) 'E(x,y)',E(x),E(y)
      write(*,*) ''

      ! output to csv
      write(10,*) B, f(B)
      write(10,*) G, f(G)
      write(10,*) W, f(W)

    Enddo

    close(10)
  end subroutine NelderMead
end program main
