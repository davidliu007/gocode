SUBROUTINE GP_tree_calculation_FULL_CODE(n_levels,n_nodes,node_type,node_values,L_echo)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    write out a file that is a fortran program containing all of the information
!    needed to calculate the solutions to an individual output code (tree) from GPCODE
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      12/14/2020  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

! Data Dictionary
INTEGER (KIND=i4b), PARAMETER :: null_node=-9999

INTEGER (KIND=i4b), INTENT(IN) :: n_levels ! total number of levels
INTEGER (KIND=i4b), INTENT(IN) :: n_nodes  ! total number of nodes

INTEGER (KIND=i4b), INTENT(IN) :: Node_Type(n_nodes)

REAL (KIND=r8b), INTENT(INOUT) :: Node_Values(n_nodes)

REAL (KIND=r8b) :: left_node_value,right_node_value,cff

INTEGER (KIND=i4b) :: i_function_node,i_node_left,i_node_right,i_node

LOGICAL, INTENT(IN) :: L_echo

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! calculate the tree function/value
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IF (L_echo) THEN
  DO i_node=1,n_nodes
    IF (NODE_Type(i_node) .ne. null_node) THEN
      WRITE(*,*) 'R: ',i_node,Node_Type(i_node),Node_Values(i_node)
    END IF
  END DO
END IF

DO i_function_node=(2**(n_levels-1))-1,1,-1  ! process all nodes from 1 level up from n_levels

  i_node_left=i_function_node*2    ! sets the 'left terminal node's index'
  i_node_right=i_node_left+1  ! sets the 'right terminal node's index'

  left_node_value=Node_Values(i_node_left)
  right_node_value=Node_Values(i_node_right)

  IF (Node_Type(i_function_node) .gt. 0) THEN ! it is a math function node, do the math

    SELECT CASE(Node_Type(i_function_node))

    CASE(1) ! LHS + RHS
      Node_Values(i_function_node)=left_node_value+right_node_value
    CASE(2) ! LHS - RHS
      Node_Values(i_function_node)=left_node_value-right_node_value
    CASE(3) ! LHS * RHS
      IF (ISNAN(left_node_value) .or. ISNAN(right_node_value)) THEN                                                                           
        Node_Values(i_function_node)=0.0D+00
      ELSE
        Node_Values(i_function_node)=left_node_value*right_node_value
      END IF
    CASE(4) ! protected: LHS/RHS
      IF (ABS(right_node_value) > 0.0D+00) THEN
        Node_Values(i_function_node)=left_node_value/right_node_value
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(5) ! Ivlev Grazing Function: (1.0D+00-EXP(-ABS(LHS*RHS)))
      cff=ABS(left_node_value*right_node_value)
      IF (cff < 100.0D+00) THEN
        Node_Values(i_function_node)=1.0D+00-EXP(-1.0D+00*cff)
      ELSE
        Node_Values(i_function_node)=1.0D+00
      END IF
    CASE(6) ! 'Michealis-Menton (ABS(RHS) / (ABS(LHS) + ABS(RHS)))'
      cff=ABS(left_node_value)+ABS(right_node_value)
      IF (cff .gt. 0.0D+00) THEN
        Node_Values(i_function_node)=ABS(right_node_value)/cff
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(7) ! Mayzaud-Poulet Grazing Function: ABS(LHS*RHS)*(1.0D+00 - EXP(-ABS(LHS*RHS)))'
      cff=ABS(left_node_value*right_node_value)
      IF (cff < 100.0D+00) THEN
        Node_Values(i_function_node)=cff*(1.0D+00-EXP(-1.0D+00*cff))
      ELSE
        Node_Values(i_function_node)=cff
      END IF
    CASE(8) ! 'LHS**RHS
      IF (ISNAN (left_node_value) .or. ISNAN (right_node_value)) THEN
        Node_Values(i_function_node)=0.0D+00
      ELSE
        IF (ABS(left_node_value) <= 1.0D-99) THEN
          Node_Values(i_function_node)=0.0D+00
        ELSE
          IF (ABS(right_node_value) <= 1.0D-99) THEN
            Node_Values(i_function_node)=1.0D+00
          ELSE
            IF (ABS(left_node_value-right_node_value) <= 1.0D-99 ) THEN ! try to eliminate a**a functions
              Node_Values(i_function_node)=0.0D+00
            ELSE
              Node_Values(i_function_node)=ABS(left_node_value)**right_node_value
              Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
              Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
            END IF
          END IF
        END IF
      END IF
    CASE(9) ! 'EXP(-ABS(LHS*RHS)'
      IF (ISNAN(left_node_value) .or. ISNAN (right_node_value)) THEN
        Node_Values(i_function_node)=0.0D+00
      ELSE
        cff=ABS(left_node_value*right_node_value)
        IF (cff < 100.0D+00) THEN
          Node_Values(i_function_node)=EXP(-1.0D+00*cff)
        ELSE
          Node_Values(i_function_node)=0.0D+00
        END IF 
      END IF 
    CASE(10) !  MIN(LHS,RHS)
      Node_Values(i_function_node)=MIN(left_node_value,right_node_value)
    CASE(11) !  MAX(LHS,RHS)
      Node_Values(i_function_node)=MAX(left_node_value,right_node_value)
    CASE(12) !   if a /= 0.0D+00, THEN b, ELSE 0.0D+00
      IF (left_node_value .ne. 0.0D+00) THEN
        Node_Values(i_function_node)=right_node_value
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(13) !   if LHS >  RHS    , THEN 1.0D+00, ELSE 0.0D+00
      IF (left_node_value .gt. right_node_value) THEN
        Node_Values(i_function_node)=1.0D+00
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(14) !   if LHS >=  RHS    , THEN 1.0D+00, ELSE 0.0D+00
      IF (left_node_value .ge. right_node_value) THEN
        Node_Values(i_function_node)=1.0D+00
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(15) !   if LHS < RHS    , THEN 1.0D+00, ELSE 0.0D+00
      IF (left_node_value .lt. right_node_value) THEN
        Node_Values(i_function_node)=1.0D+00
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(16) !   if LHS <=  RHS    , THEN 1.0D+00, ELSE 0.0D+00
      IF (left_node_value .le. right_node_value) THEN
        Node_Values(i_function_node)=1.0D+00
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(17) ! EXP(LHS)
      Node_Values(i_function_node)=EXP(left_node_value)
      Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
      Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
    CASE(18) ! EXP(RHS)
      Node_Values(i_function_node)=EXP(right_node_value)
      Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
      Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
    CASE(19) ! EXP(-LHS)
      Node_Values(i_function_node)=EXP(-1.0D+00*left_node_value)
      Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
      Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
    CASE(20) ! EXP(-RHS)
      Node_Values(i_function_node)=EXP(-1.0D+00*right_node_value)
      Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
      Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
    CASE(21) ! LOG_rhs(ABS(LHS)) ! protected routine with ABS values
      IF (left_node_value .gt. 0.0D+00 .and. right_node_value .gt. 0.0D+00) THEN
        Node_Values(i_function_node)=LOG10(left_node_value)/LOG10(right_node_value)
      ELSE
        Node_Values(i_function_node)=0.0D+00
      END IF
    CASE(22) ! LHS
      IF (ISNAN(left_node_value)) THEN
        Node_Values(i_function_node)=0.0D+00
      ELSE
        Node_Values(i_function_node)=left_node_value
      END IF
    CASE(23) ! LHS**2
      IF (ISNAN(left_node_value)) THEN
        Node_Values(i_function_node)=0.0D+00
      ELSE
        IF (ABS(left_node_value) <= 1.0D-99) THEN
          Node_Values(i_function_node)=0.0D+00
        ELSE
          Node_Values(i_function_node)=ABS(left_node_value)**2
          Node_Values(i_function_node)=MIN(Node_Values(i_function_node),1.0D+19)
          Node_Values(i_function_node)=MAX(Node_Values(i_function_node),1.0D-19)
        END IF
      END IF
    CASE DEFAULT
      WRITE(*,'(//A//)') 'Wrong case number chosen in tree evaluations'
      STOP
    END SELECT

    IF (Node_Type(i_function_node) .gt. 0 .and. L_echo) THEN
      WRITE(*,*)
      WRITE(*,*) 'i_function_node: ',i_function_node
      WRITE(*,*) 'i_node_left, i_node_right: ',i_node_left,i_node_right
      WRITE(*,*) 'Node_Type(i_function_node): ',Node_Type(i_function_node)
      WRITE(*,*) 'left_node_value and right_node_value: ', &
        left_node_value,right_node_value
      WRITE(*,*) 'Node_Values(i_function_node): ',Node_Values(i_function_node)
    END IF

  END IF

END DO

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE GP_tree_calculation_FULL_CODE
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890