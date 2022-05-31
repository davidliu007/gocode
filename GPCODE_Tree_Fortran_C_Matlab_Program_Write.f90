PROGRAM GPCODE_Tree_Fortran_C_Matlab_Program_Write
!  Purpose:
!    1) to read in the GPCODE output trees and parameters from either .dot or tree_string files
!    2) read in the IOCCG chlorophyll a and Rrs test data and calculate the OC4 chla estimate
!    3) to process the GPCODE trees and calculate the chlorophyll estimates
!    4) to write out the Fortran/C/Matlab codes for specific tree strings/parameter sets
!    5) to test the Fotran/C/Matlab codes
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      01/13/2021  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! several notes regarding the GPCODE assumptions
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! -9999 is the default value for NULL nodes
! 1 to n_node_functions for the math function calculations
! 0 is the Node_Type for a parameter
! n_levels=11                 ! number of tree levels
! n_node_functions=23         ! number of function cases
! n_variables=6               ! number of input variables
! -1 seawifs_rrs411
! -2 seawifs_rrs443
! -3 seawifs_rrs490
! -4 seawifs_rrs510
! -5 seawifs_rrs555
! -6 seawifs_rrs670

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

INTEGER (KIND=i4b), PARAMETER :: n_levels=11                 ! number of tree levels
INTEGER (KIND=i4b), PARAMETER :: n_nodes=(2**n_levels)-1     ! total number of nodes
INTEGER (KIND=i4b), PARAMETER :: n_node_functions=23         ! number of function cases
INTEGER (KIND=i4b), PARAMETER :: n_Dot_File_Cases=15         ! number of GPCODE cases
INTEGER (KIND=i4b), PARAMETER :: n_Tree_String_File_Cases=12 ! number of GPCODE cases
INTEGER (KIND=i4b), PARAMETER :: n_IOCCG_data=500            ! number of IOCCG training data points
INTEGER (KIND=i4b), PARAMETER :: n_SeaBASS_data=1016         ! number of SeaBASS chla/Rrs data points
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999        ! tree node number for "NULL" nodes

REAL (KIND=r8b), PARAMETER :: zero=0.0D+00
REAL (KIND=r8b), PARAMETER :: one=1.0D+00

INTEGER (KIND=i4b) :: Node_Type(n_nodes)

REAL (KIND=r8b) :: Node_Parameters(n_nodes)
REAL (KIND=r8b) :: Node_Values(n_nodes)

REAL (KIND=r8b) :: cff,xval,R4S,v1,v2,v3,v4,v5,v6
REAL (KIND=r8b) :: a_0,a_1,a_2,a_3,a_4

REAL (KIND=r8b), ALLOCATABLE :: chlorophyll(:)
REAL (KIND=r8b), ALLOCATABLE :: OC4_Answer(:)
REAL (KIND=r8b), ALLOCATABLE :: GPCODE_Answer(:)
REAL (KIND=r8b), ALLOCATABLE :: all_variables(:,:)
REAL (KIND=r8b), ALLOCATABLE :: variables(:)

REAL (KIND=r8b) :: OC4_MAPE,OC4_SSE,chla_SSE0,OC4_chla,chla
REAL (KIND=r8b) :: GPCODE_MAPE,GPCODE_SSE
REAL (KIND=r8b) :: seawifs_rrs411,seawifs_rrs443,seawifs_rrs490
REAL (KIND=r8b) :: seawifs_rrs510,seawifs_rrs555,seawifs_rrs670

INTEGER (KIND=i4b) :: i_node,i_variable,n_parameters,i_base,n_data,n_variables
INTEGER (KIND=i4b) :: status,i_case,i_data,n_File_cases

CHARACTER (LEN=200) :: CMD
CHARACTER (LEN=200) :: File_Name,Training_Data_File_Name
CHARACTER (LEN=200) :: Code_Root_File_Name,Code_File_Name
CHARACTER (LEN=100) :: Case_Dot_File_Names(n_Dot_File_Cases)
CHARACTER (LEN=100) :: Case_Tree_String_File_Names(n_Tree_String_File_Cases)

LOGICAL :: L_Fortran_Program_Write
LOGICAL :: L_C_Program_Write
LOGICAL :: L_Matlab_Program_Write
LOGICAL :: L_Fortran_Program_Run
LOGICAL :: L_C_Program_Run
LOGICAL :: L_Matlab_Program_Run
LOGICAL :: L_echo
LOGICAL :: L_Read_in_Dot_Tree
LOGICAL :: L_IOCCG_data
LOGICAL :: L_Seawifs,L_MODIS,L_VIIRS

!  'MERIS.dot           ', &  ! 1
!  'VIIRS.dot           ', &  ! 15 YES
!  'gw_Log08032020.dot  ', &  ! 16
!  'OLCI.dot            ', &  ! 3
DATA Case_Dot_File_Names / &
  'VIIRS.dot           ', &  ! 15 YES
  'SeaWifs2.dot        ', &  ! 5 YES
  'MODI_Aqua.dot       ', &  ! 2 YES
  'SeaWiFS_tree.dot    ', &  ! 4
  'SeaWifs_BS.dot      ', &  ! 6
  'SeaWifs_BS0363.dot  ', &  ! 7
  'SeaWifs_BS0364.dot  ', &  ! 8
  'SeaWifs_BS0416.dot  ', &  ! 9
  'SeaWifs_BS0794.dot  ', &  ! 10
  'SeaWifs_BS1124.dot  ', &  ! 11
  'SeaWifs_BS_A0405.dot', &  ! 12
  'SeaWifs_BS_B0663.dot', &  ! 13
  'SeaWifs_BS_C0612.dot', &  ! 14
  'seaWifs_1151.dot    ', &  ! 17
  'seaWifs_BS0393.dot  '/    ! 18

DATA Case_Tree_String_File_Names / &
  'Case_01_GPCODE_Best_Tree.txt', &
  'Case_02_GPCODE_Best_Tree.txt', &
  'Case_03_GPCODE_Best_Tree.txt', &
  'Case_04_GPCODE_Best_Tree.txt', &
  'Case_05_GPCODE_Best_Tree.txt', &
  'Case_06_GPCODE_Best_Tree.txt', &
  'Case_07_GPCODE_Best_Tree.txt', &
  'Case_08_GPCODE_Best_Tree.txt', &
  'Case_09_GPCODE_Best_Tree.txt', &
  'Case_10_GPCODE_Best_Tree.txt', &
  'Case_11_GPCODE_Best_Tree.txt', &
  'Case_12_GPCODE_Best_Tree.txt'/

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the program write and execution logicals
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

L_Fortran_Program_Write=.true.
L_Fortran_Program_Run=.false.

L_C_Program_Write=.true.
L_C_Program_Run=.false.

L_Matlab_Program_Write=.false.
L_Matlab_Program_Run=.false.

L_Read_in_Dot_Tree=.false.  ! read in the Tree Strings not the .dot files else read in the Tree strings

IF (L_Read_in_Dot_Tree) THEN
  n_File_Cases=n_Dot_File_Cases
ELSE
  n_File_Cases=n_Tree_String_File_Cases
END IF

L_IOCCG_data=.true.

IF (L_IOCCG_data) THEN
  n_data=n_IOCCG_data
ELSE
  n_data=n_SeaBASS_data
END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the 'read_GPCODE_tree' subroutine to read in the tree node_type and node_parameters
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

DO i_case=12,12  ! ,3  ! n_File_Cases

  L_Seawifs=.false.
  L_MODIS=.false.
  L_VIIRS=.false.
  IF (i_case .eq. 1) THEN
    L_VIIRS=.true.
    n_variables=5
  ELSE IF (i_case .eq. 2) THEN
    L_Seawifs=.true.
    n_variables=6
  ELSE IF (i_case .eq. 3) THEN
    L_MODIS=.true.
    n_variables=6
  ELSE IF (i_case .eq. 12) THEN
    L_Seawifs=.true.
    n_variables=6
  END IF

  IF (ALLOCATED(chlorophyll)) DEALLOCATE(chlorophyll)
  ALLOCATE (chlorophyll(n_data))
  chlorophyll=zero

  IF (ALLOCATED(all_variables)) DEALLOCATE(all_variables)
  ALLOCATE (all_variables(n_data,n_variables))
  all_variables=zero

  IF (ALLOCATED(variables)) DEALLOCATE(variables)
  ALLOCATE (variables(n_variables))
  variables=zero

  IF (ALLOCATED(OC4_Answer)) DEALLOCATE(OC4_Answer)
  ALLOCATE (OC4_Answer(n_data))
  OC4_Answer=zero

  IF (ALLOCATED(GPCODE_Answer)) DEALLOCATE(GPCODE_Answer)
  ALLOCATE (GPCODE_Answer(n_data))
  GPCODE_Answer=zero

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! get the training data set
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! set the metrics for the OC4 Algorithm to zero
  OC4_MAPE=0.0D+00
  OC4_SSE=0.0D+00
  chla_SSE0=0.0D+00

  IF (L_IOCCG_data) THEN ! for the IOCCG data
    Training_Data_File_Name='IOCCG_SeaWifs.txt'
  ELSE
    Training_Data_File_Name='SeaBASS_chla_Rrs.txt'
  END IF
  OPEN(UNIT=10,FILE=Training_Data_File_Name,STATUS='old',FORM='formatted', &
    ACTION='read',IOSTAT=status)
  WRITE(*,*) Training_Data_File_Name
  WRITE(*,*) Training_Data_File_Name(1:LEN_TRIM(Training_Data_File_Name))

  IF (status /= 0) THEN
    WRITE(*,*) 'Error opening the training data file for the first time'
    WRITE(*,1010) Training_Data_File_Name(1:LEN_TRIM(Training_Data_File_Name)),status
    STOP
  END IF

  IF (L_IOCCG_data) THEN ! for the IOCCG data
    READ(10,*) ! skip the first line
  ELSE ! for the SeaBASS data
    READ(10,*) ! skip the first line
    READ(10,*) ! skip the first line
    READ(10,*) ! skip the first line
  END IF

  DO i_data=1,n_data

    READ(10,*) chla,v1,v2,v3,v4,v5,v6

    seawifs_rrs411=v1
    seawifs_rrs443=v2
    seawifs_rrs490=v3
    seawifs_rrs510=v4
    seawifs_rrs555=v5
    seawifs_rrs670=v6

    R4S=MAX(seawifs_rrs443,seawifs_rrs490,seawifs_rrs510)/seawifs_rrs555

    xval=LOG10(R4S)

!   the latest 4th-order polynomial coefficients for the OC-X algorithm
    a_0=0.3272D+00
    a_1=-2.9940D+00
    a_2=2.7218D+00
    a_3=-1.2259D+00
    a_4=-0.5683D+00

    cff=a_0+(a_1*xval)+(a_2*(xval**2))+(a_3*(xval**3))+(a_4*(xval**4))

    OC4_chla=10.0**cff
    OC4_Answer(i_data)=OC4_chla
    OC4_MAPE=OC4_MAPE+ABS((OC4_chla-chla)/chla)
    OC4_SSE=OC4_SSE+((OC4_chla-chla)**2)
    chla_SSE0=chla_SSE0+(chla**2)

  END DO

  OC4_MAPE=100.0D+00*OC4_MAPE/REAL(n_data)

  REWIND(10)
  CLOSE(10)

  WRITE(*,*) 'OC4 SSE: ',OC4_SSE
  WRITE(*,*) 'OC4 MAPE: ',OC4_MAPE
  WRITE(*,*) 'chla_SSE0: ',chla_SSE0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! get the training data set
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  IF (L_IOCCG_data) THEN ! for the IOCCG data
    IF (L_Seawifs) THEN
      Training_Data_File_Name='SeaWifs.txt'
      Training_Data_File_Name='SeaWiFS_IOCCG_Training_data.csv'
    ELSE IF (L_MODIS) THEN
      Training_Data_File_Name='MODIS_IOCCG_Training_data.csv'
    ELSE IF (L_VIIRS) THEN
      Training_Data_File_Name='VIIRS_IOCCG_Training_data.csv'
    END IF
  ELSE
    Training_Data_File_Name='SeaBASS_chla_Rrs.txt'
  END IF
  OPEN(UNIT=10,FILE=Training_Data_File_Name,STATUS='old',FORM='formatted', &
    ACTION='read',IOSTAT=status)

  IF (status /= 0) THEN
    WRITE(*,*) 'Error opening the training data file for the second time'
    WRITE(*,*) Training_Data_File_Name(1:LEN_TRIM(Training_Data_File_Name)),status
    STOP
  END IF

  READ(10,*) ! skip the first line

  DO i_data=1,n_data
    READ(10,*) chlorophyll(i_data),all_variables(i_data,1:n_variables)
  END DO

  REWIND(10)
  CLOSE(10)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  File_Name=''
  Code_Root_File_Name=''
  IF (L_Read_in_Dot_Tree) THEN
    File_Name=Case_Dot_File_Names(i_case)
  ELSE
    File_Name=Case_Tree_String_File_Names(i_case)
  END IF
  Code_Root_File_Name=File_Name(1:LEN_TRIM(File_Name)-4)  ! strip off the '.dot' or '.txt' endings
  WRITE(*,*) Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))

! set the Node_Type and Node_Parameter strings to zero
  Node_Type=0
  Node_Parameters=0.0

  status=0
  IF (L_Read_in_Dot_Tree) THEN
    CALL read_GPCODE_tree_dot_file(Code_Root_File_Name,n_nodes, &
      Node_Type,Node_Parameters,n_parameters,status)
  ELSE
    CALL read_GPCODE_tree_string_file(Code_Root_File_Name,n_nodes, &
      Node_Type,Node_Parameters,n_parameters,status)
  END IF

  IF (status .ne. 0) THEN
    IF (L_Read_in_Dot_Tree) THEN
      WRITE(*,*) 'i_status .ne. 0: error in read_GPCODE_tree_dot_file subroutine'
    ELSE
      WRITE(*,*) 'i_status .ne. 0: error in read_GPCODE_tree_tree_string_file subroutine'
    END IF
    STOP
  END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Calculate the tree values using the test data set chlorophyll a and Rrs(lambda)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  GPCODE_MAPE=0.0D+00
  GPCODE_SSE=0.0D+00

  DO i_data=1,n_data

    DO i_variable=1,n_variables
      variables(i_variable)=all_variables(i_data,i_variable)
    END DO

    Node_Values=0.0D+00
    DO i_node=1,n_nodes
      IF (Node_Type(i_node) .ne. Null_Tree_Node) THEN
        IF (Node_Type(i_node) .eq. 0) THEN
          Node_Values(i_node)=Node_Parameters(i_node)
        ELSE IF (Node_Type(i_node) .lt. 0) THEN

!!!!!!!!!!!!!!!!!!!!!!
!            needed to correct the GPCODE indexing error/issue [-2 <==>-7] to [1 <==> 6]
          i_variable=(-1*Node_Type(i_node))   !  - 1
!!!!!!!!!!!!!!!!!!!!!!

          IF (i_variable .lt. 1 .or. i_variable .gt. 6) THEN
            WRITE(*,*) 'Node_Values(i_node) index: i_node is out-of-range (1-6)'
            STOP
          END IF

          Node_Values(i_node)=variables(i_variable)
        END IF
      END IF
    END DO

    L_echo=.false.
    CALL GP_tree_calculation_FULL_CODE(n_levels,n_nodes,Node_Type,Node_Values,L_echo)

!   I have tried to run this both ways, there is no difference
!off        GPCODE_Answer(i_data)=Node_Values(1)
    GPCODE_Answer(i_data)=ABS(Node_Values(1))

    GPCODE_MAPE=GPCODE_MAPE+ &
      ABS((GPCODE_Answer(i_data)-chlorophyll(i_data))/chlorophyll(i_data))
    GPCODE_SSE=GPCODE_SSE+((GPCODE_Answer(i_data)-chlorophyll(i_data))**2)

    WRITE(*,*) i_data,GPCODE_Answer(i_data),OC4_Answer(i_data),chlorophyll(i_data)

    IF (L_IOCCG_data) THEN
      i_base=0
    ELSE
      i_base=200
    END IF
    IF (L_Read_in_Dot_Tree) THEN
      WRITE(i_base+100+30+i_case,*) i_data,GPCODE_Answer(i_data),OC4_Answer(i_data), &
        chlorophyll(i_data)
    ELSE
      WRITE(i_base+200+30+i_case,*) i_data,GPCODE_Answer(i_data),OC4_Answer(i_data), &
        chlorophyll(i_data)
    END IF

  END DO
  GPCODE_MAPE=100.0D+00*GPCODE_MAPE/REAL(n_data) ! convert to units of [%]

  WRITE(*,*) 'GPCODE SSE: ',GPCODE_SSE
  WRITE(*,*) 'GPCODE MAPE: ',GPCODE_MAPE

  IF (L_Read_in_Dot_Tree) THEN
    WRITE(i_base+100+30+i_case,*) n_parameters
  ELSE
    WRITE(i_base+200+30+i_case,*) n_parameters
  END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Write out the fortran/C/Matlab codes for chosen tree or .dot files
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Code_Root_File_Name=File_Name(1:LEN_TRIM(File_Name)-4)//'_code'
  IF (L_Fortran_Program_Write) THEN
    WRITE(*,*) 'going to the Fortran code writing routine'
    CALL Write_GPCODE_Tree_Test_Fortran_code(Code_Root_File_Name,n_levels,Node_Type, &
      Node_Parameters,n_data,n_variables,status)
      IF (status /= 0) THEN
         L_Fortran_Program_Run=.false.
      END IF
  END IF

  IF (L_C_Program_Write) THEN
    WRITE(*,*) 'going to the C code writing routine'
    CALL Write_GPCODE_Tree_Test_C_code(Code_Root_File_Name,n_levels,Node_Type, &
      Node_Parameters,n_data,n_variables,status)
      IF (status /= 0) THEN
         L_C_Program_Run=.false.
      END IF
  END IF

  IF (L_Matlab_Program_Write) THEN
    WRITE(*,*) 'going to the Matlab code writing routine'
    CALL Write_GPCODE_Tree_Test_Matlab_code(Code_Root_File_Name,n_levels,Node_Type, &
      Node_Parameters,n_data,n_variables,status)
      IF (status /= 0) THEN
         L_Matlab_Program_Run=.false.
      END IF
  END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! test the programs
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  IF (L_Fortran_Program_Run) THEN
    CMD='ncargf77 -o '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//' '//&
       Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.f90'
    CALL system(CMD)
    CMD='./'//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
    CALL system(CMD)
  END IF

  IF (L_C_Program_Run) THEN
    CMD='cc -o '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//' '//&
      Code_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.c'
    CALL system(CMD)
    CMD='./'//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
    CALL system(CMD)
  END IF

  IF (ALLOCATED(chlorophyll)) DEALLOCATE(chlorophyll)
  IF (ALLOCATED(all_variables)) DEALLOCATE(all_variables)
  IF (ALLOCATED(variables)) DEALLOCATE(variables)
  IF (ALLOCATED(OC4_Answer)) DEALLOCATE(OC4_Answer)
  IF (ALLOCATED(GPCODE_Answer)) DEALLOCATE(GPCODE_Answer)

END DO

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! STOP END
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
STOP
END PROGRAM GPCODE_Tree_Fortran_C_Matlab_Program_Write
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
SUBROUTINE read_GPCODE_tree_dot_file(Code_Root_File_Name,n_nodes,Node_Type, &
  Node_Parameters,n_parameters,status)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    program to read in GPCODE tree files in order to convert algorithms to usable application code
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      01/11/2021  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

! Data Dictionary: declare KIND settings
INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

! Data Dictionary: declare constants
CHARACTER (LEN=200), INTENT(IN) :: Code_Root_File_Name
INTEGER (KIND=i4b), INTENT(IN) :: n_nodes
INTEGER (KIND=i4b), INTENT(OUT) :: Node_Type(n_nodes)
REAL (KIND=r8b), INTENT(OUT) :: Node_Parameters(n_nodes)
INTEGER (KIND=i4b), INTENT(OUT) :: n_parameters
INTEGER (KIND=i4b), INTENT(OUT) :: status

INTEGER (KIND=i4b), PARAMETER :: n_math_funcs=23
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999        ! tree node number for "NULL" nodes

INTEGER (KIND=i4b) :: i_ch,i_node,j_node
INTEGER (KIND=i4b) :: i_first_r_bracket,i_first_l_bracket
INTEGER (KIND=i4b) :: i_second_r_bracket,i_second_l_bracket
INTEGER (KIND=i4b) :: i_r_paren,i_l_paren,i_l_quote,i_r_quote
INTEGER (KIND=i4b) :: icff,i_read

INTEGER (KIND=i4b) :: Node_Variable

INTEGER (KIND=i4b) :: i_math_func,i_length

CHARACTER (LEN=1) :: ch_1
CHARACTER (LEN=50) :: ch_50
CHARACTER (LEN=200) :: File_Name
CHARACTER (LEN=200) :: Dot_File_Name
CHARACTER (LEN=7) :: Function_Node_Name(n_math_funcs),ch7

LOGICAL :: first_r_bracket,first_l_bracket,second_r_bracket,second_l_bracket
LOGICAL :: r_paren,l_paren,r_quote,l_quote

DATA Function_Node_Name /'+      ','-      ','*      ','/      ','IGF    ', &
                         'MMT    ','MPGF   ','pow    ','exp    ','min    ', &
                         'max    ','if     ','>      ','>=     ','<      ', &
                         '<=     ','expLP  ','expRP  ','expLM  ','expRM  ', &
                         'LOG_rhs','Mult_1 ','Square '/

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! zero out the node arrays
Node_Type=Null_Tree_Node
Node_Parameters=0.0D+00
n_parameters=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Dot_File_Name=Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.dot'
WRITE(*,*) Dot_File_Name
WRITE(*,*) 'Opening GPCODE Tree File: ',Dot_File_Name
OPEN(UNIT=10,FILE=Dot_File_Name,STATUS='old',FORM='formatted',ACTION='read',IOSTAT=status)
WRITE(*,*) 'status: ',status,' Dot_File_Name: ',Dot_File_Name

IF (status /= 0) THEN
  WRITE(*,1010) file_name(1:LEN_TRIM(Dot_File_Name)),status
  STOP
END IF

! skip the first two lines
READ(10,*)  ! skip the first line
READ(10,*)  ! skip the second line

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! read in the tree
DO j_node=1,n_nodes

  READ(10,1,END=300) ch_50  ! read the tree string
  1 FORMAT(a50)
  IF (j_node .gt. 1) READ(10,*,END=300) ! skip the marker lines between the tree script

! reset the character string read markers and logicals
  first_r_bracket=.false.
  first_l_bracket=.false.
  second_r_bracket=.false.
  second_l_bracket=.false.
  r_paren=.false.
  l_paren=.false.
  r_quote=.false.
  l_quote=.false.
  i_first_l_bracket=0
  i_first_r_bracket=0
  i_second_l_bracket=0
  i_second_r_bracket=0
  i_l_paren=0
  i_r_paren=0
  i_l_quote=0
  i_r_quote=0

  DO i_ch=1,LEN_TRIM(ch_50)

    READ(ch_50(i_ch:i_ch),'(a1)') ch_1

    IF (ch_1 .eq. '[' .and. .not. first_l_bracket) THEN
      first_l_bracket=.true.
      i_first_l_bracket=i_ch
    END IF
    IF (ch_1 .eq. '[' .and. first_l_bracket) THEN
      second_l_bracket=.true.
      i_second_l_bracket=i_ch
    END IF
    IF (ch_1 .eq. ']' .and. .not. first_r_bracket) THEN
      first_r_bracket=.true.
      i_first_r_bracket=i_ch
    END IF
    IF (ch_1 .eq. ']' .and. first_r_bracket) THEN
      second_r_bracket=.true.
      i_second_r_bracket=i_ch
    END IF

    IF (ch_1 .eq. '(' .and. .not. l_paren) THEN
      l_paren=.true.
      i_l_paren=i_ch
    END IF
    IF (ch_1 .eq. ')' .and. l_paren) THEN
      r_paren=.true.
      i_r_paren=i_ch
    END IF
  
    IF (ch_1 .eq. '"' .and. .not. l_quote) THEN
      l_quote=.true.
      i_l_quote=i_ch
    END IF
    IF (ch_1 .eq. '"' .and. l_quote) THEN
      r_quote=.true.
      i_r_quote=i_ch
    END IF
          
  END DO

! read in the node number
  icff=i_first_r_bracket-(i_second_l_bracket+1)
  IF (icff .eq. 1) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i1)') i_node
  ELSE IF (icff .eq. 2) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i2)') i_node
  ELSE IF (icff .eq. 3) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i3)') i_node
  ELSE IF (icff .eq. 4) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i4)') i_node
  ELSE IF (icff .eq. 5) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i5)') i_node
  ELSE IF (icff .eq. 6) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i6)') i_node
  ELSE IF (icff .eq. 7) THEN
    READ(ch_50(i_second_l_bracket+1:i_first_r_bracket-1),'(i7)') i_node
  END IF

! check for a function
  IF (.not. r_paren .and. .not. l_paren) THEN ! the line is a function node
    icff=(i_r_quote-2)-(i_first_r_bracket+1)
    i_read=i_first_r_bracket+2

    DO i_math_func=1,n_Math_Funcs
      i_length=LEN_TRIM(Function_Node_Name(i_math_func))
      ch7=Function_Node_Name(i_math_func)
      IF (ch_50(i_read:i_read+icff) .eq. ch7(1:i_length)) THEN
        Node_Type(i_node)=i_math_func
      END IF
    END DO

  ELSE IF (r_paren .and. l_paren) THEN ! the line is a terminal node
!   determine if it is a (P) "Parameter" or (V) "Variable"
    READ(ch_50(i_l_paren+1:i_l_paren+1),'(a1)') ch_1
    IF (ch_1 .eq. 'V') THEN
      READ(ch_50(i_r_paren+6:i_r_paren+6),'(i1)') Node_Variable
      Node_Type(i_node)=-1*Node_Variable ! negative Node_Types means there are variables at the node terminal
    ELSE IF (ch_1 .eq. 'P') THEN
      n_parameters=n_parameters+1
      READ(ch_50(i_r_paren+3:i_r_quote-1),'(e11.5)') Node_Parameters(i_node) 
      Node_Type(i_node)=0  ! zero means that there is a parameter value on the node site
    END IF
  END IF

END DO  ! loop for reading in the .dot file string code
300 CLOSE(10)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE read_GPCODE_tree_dot_file
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
SUBROUTINE read_GPCODE_tree_string_file(Tree_Root_File_Name,n_nodes,Node_Type, &
  Node_Parameters,n_parameters,status)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    program to read in GPCODE tree files in order to convert algorithms to usable application code
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      01/11/2021  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

! Data Dictionary: declare KIND settings
INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

! Data Dictionary: declare constants
INTEGER (KIND=i4b), INTENT(IN) :: n_nodes
!off INTEGER (KIND=i4b), PARAMETER :: n_math_funcs=23
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999        ! tree node number for "NULL" nodes

REAL (KIND=r8b), INTENT(OUT) :: Node_Parameters(n_nodes)
INTEGER (KIND=i4b), INTENT(OUT) :: Node_Type(n_nodes)
INTEGER (KIND=i4b), INTENT(OUT) :: n_parameters
INTEGER (KIND=i4b), INTENT(OUT) :: status

INTEGER (KIND=i4b) :: i,j,k,i_node,j_node,i_parameter

CHARACTER (LEN=5) :: ch_5
CHARACTER (LEN=5) :: ch_10
CHARACTER (LEN=200), INTENT(IN) :: Tree_Root_File_Name
CHARACTER (LEN=200) :: Tree_String_File_Name

LOGICAL :: L_read_tree_nodes

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! zero out the node arrays
Node_Type=Null_Tree_Node
Node_Parameters=0.0D+00
n_parameters=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Tree_String_File_Name=Tree_Root_File_Name(1:LEN_TRIM(Tree_Root_File_Name))//'.txt'
WRITE(*,*) Tree_String_File_Name
WRITE(*,*) 'Opening GPCODE Tree File: ',Tree_String_File_Name
OPEN(UNIT=10,FILE=Tree_String_File_Name,STATUS='old',FORM='formatted', &
  ACTION='read',IOSTAT=status)
WRITE(*,*) 'status: ',status,' Tree_String_File_Name: ',Tree_String_File_Name

IF (status /= 0) THEN
  WRITE(*,*) 'Error in opening the Tree String File'
  WRITE(*,1010) Tree_String_File_Name(1:LEN_TRIM(Tree_String_File_Name)),status
  STOP
END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! read in the tree nodes
L_read_tree_nodes=.true.
n_parameters=0
DO i_node=1,n_nodes
  IF (L_read_tree_nodes) THEN
    READ(10,'(1x,a5)') ch_5
    IF (ch_5 .ne. 'Nodes' .and. i_node .eq. 1) THEN
      WRITE(*,*) 'Error in tree string file: ', &
        Tree_String_File_Name(1:LEN_TRIM(Tree_String_File_Name))
    END IF
    IF (ch_5 .eq. 'Nodes') THEN ! read in the node
      BACKSPACE(10)
      READ(10,*) ch_5,i,j,k,j_node,Node_Type(j_node)
      WRITE(*,*) ch_5,i,j,k,j_node,Node_Type(j_node)

! fix the index for the variables [-2 <==> -7] to [-1 <==> -6] 
  IF (Node_Type(j_node) .lt. 0) Node_Type(j_node)=Node_Type(j_node)+1 

    ELSE
      L_read_tree_nodes=.false.
      BACKSPACE(10)
      READ(10,*) ch_10,n_parameters
    END IF
  END IF
END DO
DO i_parameter=1,n_parameters
  READ(10,*) ch_10,i,j,k,i_node,Node_Parameters(i_node)
  WRITE(*,*) i_parameter,n_parameters,ch_10,i,j,k,i_node,Node_Parameters(i_node)
END DO

CLOSE(10)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! RETURN END
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE read_GPCODE_tree_string_file
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
SUBROUTINE Write_GPCODE_Tree_Test_Fortran_code(Code_Root_File_Name,n_levels,Node_Type, &
  Node_Parameters,n_data,n_variables,status)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    program to read in GPCODE tree files in order to convert algorithms to usable application code
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      01/11/2021  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

! Data Dictionary: declare KIND settings
INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999        ! tree node number for "NULL" nodes

! Data Dictionary: declare constants
INTEGER (KIND=i4b), INTENT(IN) :: n_levels
INTEGER (KIND=i4b), INTENT(IN) :: n_variables
INTEGER (KIND=i4b), INTENT(IN) :: n_data

INTEGER (KIND=i4b) :: n_nodes    ! total number of nodes
INTEGER (KIND=i4b) :: i_node

REAL (KIND=r8b), INTENT(IN) :: Node_Parameters((2**n_levels)-1)
INTEGER (KIND=i4b), INTENT(IN) :: Node_Type((2**n_levels)-1)
INTEGER (KIND=i4b), INTENT(OUT) :: status

CHARACTER (LEN=200), INTENT(IN) :: Code_Root_File_Name
CHARACTER (LEN=200) :: Tree_Code_File_Name

n_nodes=(2**n_levels)-1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! open up a new Fortran code file and write out the tree string code in Fortran-90
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Tree_Code_File_Name=Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.f90'
WRITE(*,*) 'Fortran 90 Tree_Code_File_Name: ',Tree_Code_File_Name

WRITE(*,*) 'Opening GPCODE Tree Fortran Code File: ',Tree_Code_File_Name
OPEN(UNIT=10,FILE=Tree_Code_File_Name,STATUS='new',FORM='formatted', &
  ACTION='write',IOSTAT=status)
WRITE(*,*) 'status: ',status,' Tree_Code_File_Name: ',Tree_Code_File_Name

IF (status /= 0) THEN
  WRITE(*,*) 'Error in opening the Tree Fortran Code File'
  WRITE(*,1010) Tree_Code_File_Name(1:LEN_TRIM(Tree_Code_File_Name)),status
  STOP
END IF

WRITE(10,*) 'PROGRAM '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '! program to test the tree format to program test'
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) ! add a blank line
WRITE(10,*) 'IMPLICIT NONE'
WRITE(10,*) ! add a blank line
WRITE(10,*) 'INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value'
WRITE(10,*) 'INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value'
WRITE(10,*) ! add a blank line
WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: n_levels=',n_levels, &
  ' ! number of tree levels'
WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: n_data=',n_data, &
  ' ! number of IOCCG data points'
WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: n_nodes=(2**n_levels)-1', &
  ' ! number of tree nodes'
!off WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: n_node_functions=', &
!off  n_node_functions,' ! number of function cases'
WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: n_variables=',n_variables, &
  ' ! number of input variables'
WRITE(10,*) 'INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999 '//&
  ' ! tree node number for "NULL" nodes'
WRITE(10,*) ! add a blank line

WRITE(10,*) ! add a blank line
WRITE(10,*) 'INTEGER (KIND=i4b) :: Node_Type(n_nodes)'
WRITE(10,*) 'REAL (KIND=r8b) :: Node_Parameters(n_nodes)'
WRITE(10,*) 'REAL (KIND=r8b) :: Node_Values(n_nodes)'
WRITE(10,*) 'REAL (KIND=r8b) :: variables(n_variables)'
WRITE(10,*) ! add a blank line
WRITE(10,*) 'REAL (KIND=r8b) :: cff,chlorophyll(n_data)'
WRITE(10,*) 'REAL (KIND=r8b) :: OC4_SSE,OC4_MAPE,chla_SSE0'
WRITE(10,*) 'REAL (KIND=r8b) :: GPCODE_SSE,GPCODE_MAPE'
WRITE(10,*) 'REAL (KIND=r8b) :: all_variables(n_data,6)'
WRITE(10,*) 'REAL (KIND=r8b) :: OC4_Answer(n_data)'
WRITE(10,*) 'REAL (KIND=r8b) :: a_0,a_1,a_2,a_3,a_4,R4S,xval'
WRITE(10,*) 'REAL (KIND=r8b) :: seawifs_rrs411,seawifs_rrs443,seawifs_rrs490'
WRITE(10,*) 'REAL (KIND=r8b) :: seawifs_rrs510,seawifs_rrs555,seawifs_rrs670'
WRITE(10,*) ! add a blank line'

WRITE(10,*) 'INTEGER (KIND=i4b) :: i_variable,i_data,i_node'
WRITE(10,*) 'INTEGER (KIND=i4b) :: status'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'CHARACTER (LEN=200) :: Training_Data_File_Name'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'Logical :: L_echo ! logical to control the data dump'
WRITE(10,*) ! add a blank line

WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '! get the training data set'
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) ! add a blank line

WRITE(10,*) '! set the metrics for the OC4 Algorithm to zero'
WRITE(10,*) 'OC4_MAPE=0.0D+00'
WRITE(10,*) 'OC4_SSE=0.0D+00'
WRITE(10,*) 'chla_SSE0=0.0D+00'
WRITE(10,*) ! add a blank line

WRITE(10,*) "Training_Data_File_Name='IOCCG_SeaWifs.txt'"
WRITE(10,*) "OPEN(UNIT=10,FILE=Training_Data_File_Name,STATUS='old',FORM='formatted', &"
WRITE(10,*) "  ACTION='read',IOSTAT=status)"
WRITE(10,*) ! add a blank line

WRITE(10,*) 'IF (status /= 0) THEN'
WRITE(10,*) '  WRITE(*,1010) &'
WRITE(10,*) '    Training_Data_File_Name(1:LEN_TRIM(Training_Data_File_Name)),status'
WRITE(10,*) '  STOP'
WRITE(10,*) 'END IF'
WRITE(10,*) ! add a blank line

WRITE(10,*) '! set the most recent OC-X 4th-order polynomial coefficients'
WRITE(10,*) '  a_0=0.3272D+00'
WRITE(10,*) '  a_1=-2.9940D+00'
WRITE(10,*) '  a_2=2.7218D+00'
WRITE(10,*) '  a_3=-1.2259D+00'
WRITE(10,*) '  a_4=-0.5683D+00'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'READ(10,*) ! skip the first line'
WRITE(10,*) 'DO i_data=1,n_data'
WRITE(10,*) ! add a blank line

WRITE(10,*) '  READ(10,*) chlorophyll(i_data),all_variables(i_data,1), &'
WRITE(10,*) '    all_variables(i_data,2),all_variables(i_data,3), &'
WRITE(10,*) '    all_variables(i_data,4),all_variables(i_data,5), &'
WRITE(10,*) '    all_variables(i_data,6)'
WRITE(10,*) ! add a blank line

WRITE(10,*) '  seawifs_rrs411=all_variables(i_data,1)'
WRITE(10,*) '  seawifs_rrs443=all_variables(i_data,2)'
WRITE(10,*) '  seawifs_rrs490=all_variables(i_data,3)'
WRITE(10,*) '  seawifs_rrs510=all_variables(i_data,4)'
WRITE(10,*) '  seawifs_rrs555=all_variables(i_data,5)'
WRITE(10,*) '  seawifs_rrs670=all_variables(i_data,6)'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  R4S=MAX(seawifs_rrs443,seawifs_rrs490,seawifs_rrs510)/seawifs_rrs555'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  xval=LOG10(R4S)'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  cff=a_0+(a_1*xval)+(a_2*(xval**2))+(a_3*(xval**3))+(a_4*(xval**4))'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  OC4_Answer(i_data)=10.0**cff'
WRITE(10,*) '  OC4_MAPE=OC4_MAPE+ &'
WRITE(10,*) '    ABS((OC4_Answer(i_data)-chlorophyll(i_data))/chlorophyll(i_data))'
WRITE(10,*) '  OC4_SSE=OC4_SSE+((OC4_Answer(i_data)-chlorophyll(i_data))**2)'
WRITE(10,*) '  chla_SSE0=chla_SSE0+(chlorophyll(i_data)**2)'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'END DO'
WRITE(10,*) 'CLOSE(10)'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'OC4_MAPE=100.0D+00*OC4_MAPE/REAL(n_data)'
WRITE(10,*) ! add a blank line

WRITE(10,*) "WRITE(*,*) 'OC4 SSE: ',OC4_SSE"
WRITE(10,*) "WRITE(*,*) 'OC4 MAPE: ',OC4_MAPE"
WRITE(10,*) "WRITE(*,*) 'chla_SSE0: ',chla_SSE0"
WRITE(10,*) ! add a blank line

WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '! set the Node_Type and Node_Parameters strings'
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

WRITE(10,*) ! add a blank line
WRITE(10,*) '  Node_Type=Null_Tree_Node   ! ID number for NULL nodes'
WRITE(10,*) '  Node_Parameters=0.0D+00'
WRITE(10,*) '  Node_Values=0.0D+00'
WRITE(10,*) ! add a blank line
DO i_node=1,n_nodes
  IF (Node_Type(i_node) .ne. Null_Tree_Node) THEN
    IF (Node_Type(i_node) .eq. 0) THEN
      WRITE(10,*) '  Node_Type(',i_node,')=',Node_Type(i_node)
      WRITE(10,*) '  Node_Parameters(',i_node,')=',Node_Parameters(i_node)
    ELSE IF (Node_Type(i_node) .lt. 0 .or. Node_Type(i_node) .gt. 0) THEN
      WRITE(10,*) '  Node_Type(',i_node,')=',Node_Type(i_node)
    END IF
  END IF
END DO
WRITE(10,*) ! put in a space

WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '! process the GPCODE tree using the IOOCG Rrs data'
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

WRITE(10,*) ! put in a space

WRITE(10,*) 'GPCODE_MAPE=0.0D+00'
WRITE(10,*) 'GPCODE_SSE=0.0D+00'
WRITE(10,*) ! put in a space

WRITE(10,*) 'READ(10,*) ! skip the first line'
WRITE(10,*) 'DO i_data=1,n_data'
WRITE(10,*) ! put in a space

WRITE(10,*) '  variables(1)=all_variables(i_data,1)'
WRITE(10,*) '  variables(2)=all_variables(i_data,2)'
WRITE(10,*) '  variables(3)=all_variables(i_data,3)'
WRITE(10,*) '  variables(4)=all_variables(i_data,4)'
WRITE(10,*) '  variables(5)=all_variables(i_data,5)'
WRITE(10,*) '  variables(6)=all_variables(i_data,6)'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  Node_Values=0.0D+00'
WRITE(10,*) '! fill in the Node_Values tree with the variables and parameters'
WRITE(10,*) '  DO i_node=1,n_nodes'
WRITE(10,*) '    IF (Node_Type(i_node) .ne. Null_Tree_Node) THEN'
WRITE(10,*) '      IF (Node_Type(i_node) .eq. 0) THEN'
WRITE(10,*) '        Node_Values(i_node)=Node_Parameters(i_node)'
WRITE(10,*) '      ELSE IF (Node_Type(i_node) .lt. 0) THEN'
WRITE(10,*) '        i_variable=-1*Node_Type(i_node)'
WRITE(10,*) '        Node_Values(i_node)=variables(i_variable)'
WRITE(10,*) '      END IF'
WRITE(10,*) '    END IF'
WRITE(10,*) '  END DO'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  L_echo=.false.'
WRITE(10,*) '  CALL GP_tree_calculation_FULL_CODE(n_levels,n_nodes, &'
WRITE(10,*) '    node_type,node_values,L_echo)'
WRITE(10,*) !  add a blank line

WRITE(10,*) "  WRITE(*,*) 'GPCODE Fortran-code estimate: ',Node_Values(1)"
WRITE(10,*) '  WRITE(*,*) i_data,chlorophyll(i_data),Node_Values(1)'
WRITE(10,*) '  WRITE(33,*) i_data,chlorophyll(i_data),Node_Values(1)'
WRITE(10,*) !  add a blank line

WRITE(10,*) '  GPCODE_MAPE=GPCODE_MAPE+ &'
WRITE(10,*) '    ABS((Node_Values(1)-chlorophyll(i_data))/chlorophyll(i_data))'
WRITE(10,*) '  GPCODE_SSE=GPCODE_SSE+((Node_Values(1)-chlorophyll(i_data))**2)'
WRITE(10,*) !  add a blank line

WRITE(10,*) 'END DO'
WRITE(10,*) ! add a blank line

WRITE(10,*) 'GPCODE_MAPE=100.0D+00*GPCODE_MAPE/REAL(n_data)'
WRITE(10,*) ! add a blank line

WRITE(10,*) "WRITE(*,*) 'GPCODE SSE: ',GPCODE_SSE"
WRITE(10,*) "WRITE(*,*) 'GPCODE MAPE: ',GPCODE_MAPE"
WRITE(10,*) ! add a blank line

WRITE(10,*) "1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)"
WRITE(10,*) ! add a blank line

WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '! STOP END'
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) 'STOP'
WRITE(10,*) 'END PROGRAM '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
WRITE(10,*) '!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(10,*) '!23456789012345678901234567890123456789012345678901234'// &
  '567890123456789012345678901234567890'
!off WRITE(10,*) "INCLUDE './GP_tree_calculation_FULL_CODE.f90"
WRITE(10,*) "INCLUDE './GPCODE_analysis_code/GP_tree_calculation_FULL_CODE.f90'"

CLOSE(10)

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! RETURN END
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE Write_GPCODE_Tree_Test_Fortran_code
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
SUBROUTINE Write_GPCODE_Tree_Test_C_code(Code_Root_File_Name,n_levels,Node_Type, &
  Node_Parameters,n_data,n_variables,status)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    program to read in GPCODE tree files in order to convert algorithms to usable application code
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      01/14/2021  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

! Data Dictionary: declare KIND settings
INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

! Data Dictionary: declare constants
INTEGER (KIND=i4b), INTENT(IN) :: n_levels
INTEGER (KIND=i4b), INTENT(IN) :: n_data
INTEGER (KIND=i4b), INTENT(IN) :: n_variables
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999

INTEGER (KIND=i4b) :: n_nodes  ! total number of nodes
INTEGER (KIND=i4b) :: i_node

REAL (KIND=r8b), INTENT(IN) :: Node_Parameters((2**n_levels)-1)
INTEGER (KIND=i4b), INTENT(IN) :: Node_Type((2**n_levels)-1)
INTEGER (KIND=i4b), INTENT(OUT) :: status

CHARACTER (LEN=200), INTENT(IN) :: Code_Root_File_Name
CHARACTER (LEN=200) :: Tree_Code_File_Name

n_nodes=(2**n_levels)-1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! create the c code
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Tree_Code_File_Name=Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.c'
WRITE(*,*) Tree_Code_File_Name
WRITE(*,*) 'Opening GPCODE Tree C Code File: ',Tree_Code_File_Name
OPEN(UNIT=20,FILE=Tree_Code_File_Name,STATUS='new',FORM='formatted', &
  ACTION='write',IOSTAT=status)
WRITE(*,*) 'status: ',status,' Tree_Code_File_Name: ',Tree_Code_File_Name

IF (status /= 0) THEN
  WRITE(*,*) 'Error in opening the Tree C Code File'
  WRITE(*,1010) Tree_Code_File_Name(1:LEN_TRIM(Tree_Code_File_Name)),status
  STOP
END IF

WRITE(20,*) '/* C VERSION OF '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'*/'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
          'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* program to test the tree format to test the C program */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*)  ! add a blank line

WRITE(20,*) '#include <stdio.h>'
WRITE(20,*) '#include <stdlib.h>'
WRITE(20,*) '#include <math.h>'
WRITE(20,*)  ! add a blank line

WRITE(20,*) '#define n_levels ',n_levels,' /* number of tree levels */'
WRITE(20,*) '#define n_nodes ',(2**n_levels)-1,' /* number of nodes */'
WRITE(20,*) '#define n_variables ',n_variables,' /* number of input variables */'
WRITE(20,*) '#define n_data ',n_data,' /* number of input variables */'
WRITE(20,*) '#define Null_Tree_Node ',Null_Tree_Node,' /* Null Tree Node default value */'
WRITE(20,*)  ! add a blank line

WRITE(20,*) 'int Node_Type[n_nodes+1];'
WRITE(20,*) 'float Node_Parameters[n_nodes+1];'
WRITE(20,*) 'float Node_Values[n_nodes+1];'
WRITE(20,*) 'float all_variables[n_data][n_variables+1];'
WRITE(20,*) 'float a_0,a_1,a_2,a_3,a_4;'
WRITE(20,*) 'float seawifs_rrs411, seawifs_rrs443, seawifs_rrs490;'
WRITE(20,*) 'float seawifs_rrs510, seawifs_rrs555, seawifs_rrs670;'
WRITE(20,*) 'float var1,var2,var3,var4,var5,var6,R4S,dff;'
WRITE(20,*) 'float right_node_value,left_node_value;'
WRITE(20,*) 'float chla[n_data+1],cff,chla_SSE0;'
WRITE(20,*) 'float OC4_SSE,OC4_MAPE;'
WRITE(20,*) 'float OC4_Answer[n_data+1];'
WRITE(20,*) 'float GPCODE_SSE,GPCODE_MAPE;'
WRITE(20,*) 'float GPCODE_Answer[n_data+1];'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'FILE *file_input_pointer;'
WRITE(20,*) 'FILE *file_output_pointer;'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'int i_function,i_function_start,i_function_stop;'
WRITE(20,*) 'int i_level,i_node_left,i_node_right,i_start;'
WRITE(20,*) 'int i_node,i_variable,i_data;'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'char file_input_name[50]="IOCCG_SeaWifs.txt";'
WRITE(20,*) 'char file_output_name[50]="GPCODE_c_test_output.txt";'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* Main code */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'int main()'
WRITE(20,*) '{'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* NOTE: Read in the IOCCG chlorophyll and Rrs training data */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* set the most recent OC-X 4th-order polynomial coefficients */'
WRITE(20,*) 'a_0=0.3272;'
WRITE(20,*) 'a_1=-2.9940;'
WRITE(20,*) 'a_2=2.7218;'
WRITE(20,*) 'a_3=-1.2259;'
WRITE(20,*) 'a_4=-0.5683;'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'file_input_pointer = fopen(file_input_name,"r");'
WRITE(20,*) 'file_output_pointer = fopen(file_output_name,"w");'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'fprintf(file_output_pointer,"IOCCG vs GPCODE chla: \r\n");'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'if (file_input_pointer == NULL)'
WRITE(20,*) ' {'
WRITE(20,*) '   printf("Error !");'
WRITE(20,*) '   exit(1);'
WRITE(20,*) ' }'
WRITE(20,*) ! add a blank line

WRITE(20,*) ' fscanf(file_input_pointer,"%*[^\n]\n"); /* skip the first line */'
WRITE(20,*) ! add a blank line

WRITE(20,*) ' for (i_data=1; i_data<=n_data; i_data++){'
!off WRITE(20,*) ' for (i_data=1; i_data<=n_data; i_data++){'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  fscanf(file_input_pointer,"%f %f %f %f %f %f %f",'//&
  '&chla[i_data],&var1,&var2,&var3,&var4,&var5,&var6);'
WRITE(20,*) '  printf("%f, %f, %f, %f, %f, %f, %f \n", '//&
  'chla[i_data],var1,var2,var3,var4,var5,var6);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  all_variables[i_data][1]=var1;'
WRITE(20,*) '  all_variables[i_data][2]=var2;'
WRITE(20,*) '  all_variables[i_data][3]=var3;'
WRITE(20,*) '  all_variables[i_data][4]=var4;'
WRITE(20,*) '  all_variables[i_data][5]=var5;'
WRITE(20,*) '  all_variables[i_data][6]=var6;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  seawifs_rrs411=all_variables[i_data][1];'
WRITE(20,*) '  seawifs_rrs443=all_variables[i_data][2];'
WRITE(20,*) '  seawifs_rrs490=all_variables[i_data][3];'
WRITE(20,*) '  seawifs_rrs510=all_variables[i_data][4];'
WRITE(20,*) '  seawifs_rrs555=all_variables[i_data][5];'
WRITE(20,*) '  seawifs_rrs670=all_variables[i_data][6];'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  cff=seawifs_rrs510;'
WRITE(20,*) '  if (seawifs_rrs443 > cff) {'
WRITE(20,*) '    cff=seawifs_rrs443;'
WRITE(20,*) '  }'
WRITE(20,*) '  if (seawifs_rrs490 > cff) {'
WRITE(20,*) '    cff=seawifs_rrs490;'
WRITE(20,*) '  }'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  R4S=cff/seawifs_rrs555;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  printf("%.9f , %.9f \r\n ", R4S,cff);'
WRITE(20,*) '  printf("%.9f , %.9f \r\n ", R4S,cff);'
WRITE(20,*) '  printf("%.9f , %.9f \r\n ", R4S,cff);'
WRITE(20,*) '  printf("%.9f , %.9f \r\n ", R4S,cff);'

WRITE(20,*) '  dff=log10(R4S);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  cff=a_0+(a_1*dff)+(a_2*(pow(dff,2)))+(a_3*(pow(dff,3)))+(a_4*(pow(dff,4)));'
WRITE(20,*) '  printf("%.9f , %.9f \r\n ", cff,dff);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  OC4_Answer[i_data]=pow(10.0,cff);  /* 10.0**cff */'
WRITE(20,*) '  cff=(OC4_Answer[i_data]-chla[i_data])/chla[i_data];'
WRITE(20,*) '  OC4_MAPE=OC4_MAPE+fabs(cff);'
WRITE(20,*) '  cff=(OC4_Answer[i_data]-chla[i_data]);'
WRITE(20,*) '  OC4_SSE=OC4_SSE+pow(cff,2);'
WRITE(20,*) '  chla_SSE0=chla_SSE0+pow(chla[i_data],2);'

WRITE(20,*) '  printf("%.9f , %.9f \r\n ", chla[i_data],OC4_Answer[i_data]);'
WRITE(20,*)    ! add a blank line

WRITE(20,*) '  }' ! end of the data read loop  
WRITE(20,*)    ! add a blank line

WRITE(20,*) '  OC4_MAPE=100.0*OC4_MAPE/n_data;'
WRITE(20,*)    ! add a blank line

WRITE(20,*) ' fclose(file_input_pointer);'
WRITE(20,*) ' fclose(file_output_pointer);'
WRITE(20,*)   ! add a blank line

WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* Initialize the Node_Values array */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*)   ! add a blank line

WRITE(20,*) '/* set the defaults for nodes */'
WRITE(20,*) 'for(i_node=1; i_node<=n_nodes+1;i_node++){Node_Type[i_node]=Null_Tree_Node;}'
WRITE(20,*) 'for(i_node=1; i_node<=n_nodes+1; i_node++){Node_Parameters[i_node]=0.0;}'
WRITE(20,*) 'for(i_node=1; i_node<=n_nodes+1; i_node++){Node_Values[i_node]=0.0;}'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* write out the active node_type and parameters */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) ! add a blank line

DO i_node=1,n_nodes
  IF (Node_Type(i_node) .ne. Null_Tree_Node) THEN
    IF (Node_Type(i_node) .eq. 0) THEN
      WRITE(20,*) 'Node_Type[',i_node,']=',Node_Type(i_node),';'
      WRITE(20,*) 'Node_Parameters[',i_node,']=',Node_Parameters(i_node),';'
    ELSE IF (Node_Type(i_node) .gt. 0) THEN
      WRITE(20,*) 'Node_Type[',i_node,']=',Node_Type(i_node),';'
    ELSE IF (Node_Type(i_node) .lt. 0) THEN
      WRITE(20,*) 'Node_Type[',i_node,']=',Node_Type(i_node),';'
    END IF
  END IF
END DO
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* loop through the data set and calculate the GPCODE tree values */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) ! add a blank line

WRITE(20,*) 'GPCODE_SSE=0.0;'
WRITE(20,*) 'GPCODE_MAPE=0.0;'
WRITE(20,*) ! add a blank line

WRITE(20,*) ' for (i_data=1; i_data<=n_data; i_data++){'
!off WRITE(20,*) ' for (i_data=1; i_data<=n_data; i_data++){'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* fill in the tree_values array */'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* fill in the Node_Values tree with the variables and parameters */'
WRITE(20,*) 'for (i_node=1; i_node<=n_nodes; i_node++){'
WRITE(20,*) '  Node_Values[i_node]=0.0;'
WRITE(20,*) '  if (Node_Type[i_node] > Null_Tree_Node) {'
WRITE(20,*) '    if (Node_Type[i_node] == 0) {'
WRITE(20,*) '      Node_Values[i_node]=Node_Parameters[i_node];'
WRITE(20,*) '    }'
WRITE(20,*) '    if (Node_Type[i_node] < 0) {'
WRITE(20,*) '      i_variable=-1*Node_Type[i_node];'
WRITE(20,*) '      Node_Values[i_node]=all_variables[i_data][i_variable];'
WRITE(20,*) '    }'
WRITE(20,*) '  }'
WRITE(20,*) '}'
WRITE(20,*) ! add a blank line

WRITE(20,*) '/* move up the tree structure levels from (2**(n_levels-1))-1 to 1 */'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  i_start=pow(2,(n_levels-1))-1;'
WRITE(20,*) '  for (i_function=i_start; i_function>=1; i_function--) {'
WRITE(20,*) '    i_node_left=i_function*2; /* set left node index */'
WRITE(20,*) '    i_node_right=i_node_left+1; /* set right node index */'
WRITE(20,*) !    add a blank line
WRITE(20,*) '    left_node_value=Node_Values[i_node_left]; /* set left node value */'
WRITE(20,*) '    right_node_value=Node_Values[i_node_right]; /* set right node value */'
WRITE(20,*) !    add a blank line

WRITE(20,*) '    if (Node_Type[i_function] > 0) {'
WRITE(20,*) !      add a blank line

WRITE(20,*) '      switch(Node_Type[i_function]) {'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 1: /* [LHS + RHS] */'
WRITE(20,*) '          Node_Values[i_function]=left_node_value+right_node_value;'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 2: /* [LHS - RHS] */'
WRITE(20,*) '          Node_Values[i_function]=left_node_value-right_node_value;'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 3: /* [LHS * RHS] */'
WRITE(20,*) '          Node_Values[i_function]=left_node_value*right_node_value;'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 4: /* protected [LHS / RHS] */'
WRITE(20,*) '          if (fabsf(right_node_value) > 0.0E+00) {'
WRITE(20,*) '            Node_Values[i_function]=left_node_value/right_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 5: /* Ivlev Grazing Function: '// &
  '[(1.0E+00-exp(-fabsf(LHS*RHS)))] */'
WRITE(20,*) '          cff=fabsf(left_node_value*right_node_value);'
WRITE(20,*) '          if (cff < 100.0E+00) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00-exp(-1.0E+00*cff);'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 6: /* Michealis-Menton '// &
  '[(fabsf(RHS) / (fabsf(LHS) + fabsf(RHS)))] */'
WRITE(20,*) '          cff=fabsf(left_node_value)+fabsf(right_node_value);'
WRITE(20,*) '          if (cff > 0.0E+00) {'
WRITE(20,*) '            Node_Values[i_function]=fabsf(right_node_value)/cff;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 7: /* Mayzaud-Poulet Grazing Function: '//&
  ' [fabsf(LHS*RHS)*(1.0E+00 - exp(-fabsf(LHS*RHS)))] */'
WRITE(20,*) '          cff=fabsf(left_node_value*right_node_value);'
WRITE(20,*) '          if (cff < 100.0E+00) {'
WRITE(20,*) '            Node_Values[i_function]=cff*(1.0E+00-exp(-1.0E+00*cff));'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=cff;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 8: /* LHS**RHS */'
WRITE(20,*) '          if (fabsf(left_node_value) <= 1.0E-99) {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            if (fabsf(right_node_value) <= 1.0E-99) {'
WRITE(20,*) '              Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '            }'
WRITE(20,*) '            else {'
WRITE(20,*) '              if (fabsf(left_node_value-right_node_value) <= 1.0E-99) {'
WRITE(20,*) '                Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '              }'
WRITE(20,*) '              else {'
WRITE(20,*) '                Node_Values[i_function]='// &
  'pow(fabsf(left_node_value),right_node_value);'
WRITE(20,*) '                if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '                  Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '                }'
WRITE(20,*) '                if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '                Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '                }'
WRITE(20,*) '              }'
WRITE(20,*) '            }'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 9: /* exp(-fabsf(LHS*RHS) */'
WRITE(20,*) '          if (isnan(left_node_value) || isnan(right_node_value)) {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            cff=fabsf(left_node_value*right_node_value);'
WRITE(20,*) '            if (cff < 100.0E+00 ) {'
WRITE(20,*) '              Node_Values[i_function]=exp(-1.0E+00*cff);'
WRITE(20,*) '            }'
WRITE(20,*) '            else {'
WRITE(20,*) '              Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '            }'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 10: /* min(LHS,RHS) */'
WRITE(20,*) '          if (left_node_value < right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=left_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=right_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 11: /* max(LHS,RHS) */'
WRITE(20,*) '          if (left_node_value > right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=left_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=right_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 12: /* if LHS /= 0.0E+00, then RHS, else 0.0E+00 */'
WRITE(20,*) '          if (left_node_value != 0.0E+00) {'
WRITE(20,*) '            Node_Values[i_function]=right_node_value;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 13: /* if LHS > RHS, then 1.0E+00, else 0.0E+00 */'
WRITE(20,*) '          if (left_node_value > right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 14: /* if LHS >= RHS, then 1.0E+00, else 0.0E+00 */'
WRITE(20,*) '          if (left_node_value >= right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 15: /* if LHS < RHS, then 1.0E+00, else 0.0E+00 */'
WRITE(20,*) '          if (left_node_value < right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 16: /* if LHS <= RHS, then 1.0E+00, else 0.0E+00 */'
WRITE(20,*) '          if (left_node_value <= right_node_value) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 17: /* exp(LHS) */'
WRITE(20,*) '          Node_Values[i_function]=exp(left_node_value);'
WRITE(20,*) '          if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '          }'
WRITE(20,*) '          if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 18: /* exp(RHS) */'
WRITE(20,*) '          Node_Values[i_function]=exp(right_node_value);'
WRITE(20,*) '          if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '          }'
WRITE(20,*) '          if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 19: /* exp(-LHS) */'
WRITE(20,*) '          Node_Values[i_function]=exp(-1.0E+00*left_node_value);'
WRITE(20,*) '          if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '          }'
WRITE(20,*) '          if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 20: /* exp(-RHS) */'
WRITE(20,*) '          Node_Values[i_function]=exp(-1.0E+00*right_node_value);'
WRITE(20,*) '          if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '          }'
WRITE(20,*) '          if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '            Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 21: /* log_|RHS|(|LHS|)) */'
WRITE(20,*) '          cff=log10(fabsf(right_node_value));'
WRITE(20,*) '          Node_Values[i_function]=log10(fabsf(left_node_value))/cff;'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 22: /* LHS */'
WRITE(20,*) '          Node_Values[i_function]=left_node_value;'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        case 23: /* LHS**2 */'
WRITE(20,*) '          if (fabsf(left_node_value) <= 1.0E-99) {'
WRITE(20,*) '            Node_Values[i_function]=0.0E+00;'
WRITE(20,*) '          }'
WRITE(20,*) '          else {'
WRITE(20,*) '            Node_Values[i_function]=pow(left_node_value,2);'
WRITE(20,*) '            if (Node_Values[i_function] > 1.0E+19) {'
WRITE(20,*) '              Node_Values[i_function]=1.0E+19;'
WRITE(20,*) '            }'
WRITE(20,*) '            if (Node_Values[i_function] < 1.0E-19) {'
WRITE(20,*) '              Node_Values[i_function]=1.0E-19;'
WRITE(20,*) '            }'
WRITE(20,*) '          }'
WRITE(20,*) '          break;'
WRITE(20,*) ! add a blank line

WRITE(20,*) '        default:'
WRITE(20,*) '          printf("default case\r\n");'
WRITE(20,*) '          printf("Wrong case number in tree evaluations\r\n");'
WRITE(20,*) '          printf("Bad case number\r\n");'
WRITE(20,*) '          exit(0);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '      } /* end switch */'
WRITE(20,*) !      add a blank line
WRITE(20,*) ! add a blank line

WRITE(20,*) '  } /* end if */'
WRITE(20,*) ! add a blank line

WRITE(20,*) '} /* next for tree node fun*/'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  printf("GPCODE C-code tree estimate: %12.9f\r\n", Node_Values[1]);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  printf("Corrected C-code GPCODE tree estimate: %12.9f , '// &
  '%d\r\n", Node_Values[1],i_data);'
WRITE(20,*) '  fprintf(file_output_pointer,"%d, %12.9f, %12.9f\r\n", '//&
  'i_data,chla[i_data],Node_Values[1]);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  GPCODE_Answer[i_data]=Node_Values[1];'
WRITE(20,*) '  cff=fabs((GPCODE_Answer[i_data]-chla[i_data])/chla[i_data]);'
WRITE(20,*) '  GPCODE_MAPE=GPCODE_MAPE+cff;'
WRITE(20,*) '  GPCODE_SSE=GPCODE_SSE+(pow((GPCODE_Answer[i_data]-chla[i_data]),2));'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  }'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  GPCODE_MAPE=100.0*GPCODE_MAPE/n_data;'
WRITE(20,*)    ! add a blank line

WRITE(20,*) '  printf("C-code     i_data: %d\r\n",i_data);'
WRITE(20,*) '  printf("           IOCCG chla: %12.9f\r\n",chla[i_data]);'
WRITE(20,*) '  printf("           OC4 estimate: %12.9f\r\n",OC4_Answer[i_data]);'
WRITE(20,*) '  printf("           GPCODE tree estimate: %12.9f\r\n",Node_Values[1]);'
WRITE(20,*) ! add a blank line

WRITE(20,*) '  printf("%.9f , %.9f \r\n ", GPCODE_SSE,GPCODE_MAPE);'
WRITE(20,*)    ! add a blank line

WRITE(20,*) '  return 0;'
WRITE(20,*) '} /*close main */'
WRITE(20,*) ! add a blank line
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* end C program */'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*)  ! add a blank line

WRITE(20,*) '/* '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'*/'
WRITE(20,*) '/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
  'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */'
WRITE(20,*) '/* 123456789012345678901234567890123456789012345678901234'// &
  '567890123456789012345678901234567890 */'

CLOSE(20)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! RETURN END
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE Write_GPCODE_Tree_Test_C_code
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
SUBROUTINE Write_GPCODE_Tree_Test_Matlab_code(Code_Root_File_Name,n_levels,Node_Type, &
  Node_Parameters,n_data,n_variables,status)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  Purpose:
!    program to read in GPCODE tree files in order to convert algorithms to usable application code
!
!  Record of revisions:
!      Date        Programmer         Description of change
!      ====        ==========         =====================
!      12/14/2015  J.R. Moisan        Original code
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

IMPLICIT NONE

! Data Dictionary: declare KIND settings
INTEGER (KIND=4), PARAMETER :: i4b = 4 ! compiler dependent value
INTEGER (KIND=4), PARAMETER :: r8b = 8 ! compiler dependent value

! Data Dictionary: declare constants
INTEGER (KIND=i4b), INTENT(IN) :: n_levels
INTEGER (KIND=i4b), INTENT(IN) :: n_data
INTEGER (KIND=i4b), INTENT(IN) :: n_variables
INTEGER (KIND=i4b), PARAMETER :: Null_Tree_Node=-9999

INTEGER (KIND=i4b) :: i_node
INTEGER (KIND=i4b) :: n_nodes     ! total number of nodes

INTEGER (KIND=i4b), INTENT(IN) :: Node_Type((2**n_levels)-1)
REAL (KIND=r8b), INTENT(IN) :: Node_Parameters((2**n_levels)-1)

!off REAL (KIND=r8b) :: Variables(n_variables)

INTEGER (KIND=i4b), INTENT(OUT) :: status

CHARACTER (LEN=200), INTENT(IN) :: Code_Root_File_Name
CHARACTER (LEN=200) :: Tree_Code_File_Name

n_nodes=(2**n_levels)-1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Tree_Code_File_Name=Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.m'
WRITE(*,*) Tree_Code_File_Name
WRITE(*,*) 'Opening GPCODE Tree Matlab Code File: ',Tree_Code_File_Name
OPEN(UNIT=30,FILE=Tree_Code_File_Name,STATUS='new',FORM='formatted', &
  ACTION='write',IOSTAT=status)
WRITE(*,*) 'status: ',status,' Tree_Code_File_Name: ',Tree_Code_File_Name

IF (status /= 0) THEN
  WRITE(*,*) 'Error in opening the Tree Matlab Code File'
  WRITE(*,1010) Tree_Code_File_Name(1:LEN_TRIM(Tree_Code_File_Name)),status
  STOP
END IF

WRITE(30,*) '% Matlab VERSION OF '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) '% program to test the tree format to program test'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*)  ! add a blank line
WRITE(30,*) '% integer n_levels;'
WRITE(30,*) '% integer n_nodes;'
WRITE(30,*) '% integer n_node_functions;'
WRITE(30,*) '% integer n_variables;'
WRITE(30,*) '% integer Node_Type(1:n_nodes);'
WRITE(30,*) '% real Tree_Parameters(1:n_nodes);'
WRITE(30,*) '% real Tree_Values(1:n_nodes);'
WRITE(30,*) '% real chla(1:n_data);'
WRITE(30,*) '% real OC4_Answer(1:n_data);'
WRITE(30,*) '% real GP_Answer(1:n_data);'
WRITE(30,*) '% real OC4_SSE,OC4_MAPE,GP_SSE,GP_MAPE,chla_SSE0;'
WRITE(30,*) '% real all_variables(1:n_data,1:n_variables);'
WRITE(30,*)  ! add a blank line

WRITE(30,*) '% real right_node_value,left_node_value;'
WRITE(30,*) '% real cff;'
WRITE(30,*) '% real v0,v1,v2,v3,v4,v5,v6;'
WRITE(30,*) '% real a_0,a_1,a_2,a_3,a_4;'
WRITE(30,*)  ! add a blank line

WRITE(30,*) '% real seawifs_rrs411,seawifs_rrs443,seawifs_rrs490;'
WRITE(30,*) '% real seawifs_rrs510,seawifs_rrs555,seawifs_rrs670;'
WRITE(30,*)  ! add a blank line

WRITE(30,*) '% integer i_function,i_start;'
WRITE(30,*) '% integer i_node_left,i_node_right;'
WRITE(30,*) '% integer i_node,i_variable;'
WRITE(30,*)  ! add a blank line

WRITE(30,*) 'n_levels =',n_levels,'; % number of tree levels'
WRITE(30,*) 'n_nodes =',(2**n_levels)-1,'; % number of nodes'
WRITE(30,*) 'Null_Tree_Node = -9999; % number of nodes '
WRITE(30,*) 'n_variables =',n_variables,'; % number of input variables'
WRITE(30,*) 'n_data =',n_data,'; % number of input training data'
WRITE(30,*)  ! add a blank line

WRITE(30,*) 'char file_input_name[100]="SeaWifs.txt";'
WRITE(30,*) 'char file_output_name[100]="'// &
  Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))//'.txt";'
WRITE(30,*)  ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) '% NOTE: Read in the IOCCG chlorophyll and Rrs training data'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% open the input data and output results files'
WRITE(30,*) "input_file_pointer=fopen(input_file_name,'r');"
WRITE(30,*) "if input_file_pointer < 0, error('Cannot open input file'); end"
WRITE(30,*) ! add a blank line

WRITE(30,*) "output_file_pointer=fopen(output_file_name,'w');"
WRITE(30,*) "if output_file_pointer < 0, error('Cannot open output file'); end"
WRITE(30,*) ! add a blank line

WRITE(30,*) 'formatSpec="%f %f %f %f %f %f %f \n";'
WRITE(30,*) '% read in the IOCCG data set'
WRITE(30,*) ! add a blank line

!off WRITE(30,*) 'fscanf(input_file_pointer,format_Spec,all_variables)'
!off WRITE(30,*) "all_variables=all_variables'"
WRITE(30,*) 'tline=fgets(input_file_pointer); % skip the first line in the SeaWifs.txt file'
WRITE(30,*) 'for i_data=1;n_data'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  fscanf(input_file_pointer,formatSpec,v0,v1,v2,v3,v4,v5,v6;'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  chla(i_data)=v0;'
WRITE(30,*) '  all_variables(i_data,1)=v1;'
WRITE(30,*) '  all_variables(i_data,2)=v2;'
WRITE(30,*) '  all_variables(i_data,3)=v3;'
WRITE(30,*) '  all_variables(i_data,4)=v4;'
WRITE(30,*) '  all_variables(i_data,5)=v5;'
WRITE(30,*) '  all_variables(i_data,6)=v6;'
WRITE(30,*) ! add a blank line

WRITE(30,*) 'end % end of reading in the IOCCG data file'
WRITE(30,*) 'fclose(input_file_pointer); % close the input IOCCG data file'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) '% calculate the OC4 chlorophyl values, SSE and MAPE'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
            'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% set the most recent OC-X 4th-order polynomial coefficients'
WRITE(30,*) 'a_0=0.3272;'
WRITE(30,*) 'a_1=-2.9940;'
WRITE(30,*) 'a_2=2.7218;'
WRITE(30,*) 'a_3=-1.2259;'
WRITE(30,*) 'a_4=-0.5683;'
WRITE(30,*) ! add a blank line

WRITE(30,*) 'for i_data=1;n_data'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  seawifs_rrs411=all_variables[i_data][1];'
WRITE(30,*) '  seawifs_rrs443=all_variables[i_data][2];'
WRITE(30,*) '  seawifs_rrs490=all_variables[i_data][3];'
WRITE(30,*) '  seawifs_rrs510=all_variables[i_data][4];'
WRITE(30,*) '  seawifs_rrs555=all_variables[i_data][5];'
WRITE(30,*) '  seawifs_rrs670=all_variables[i_data][6];'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  cff=seawifs_rrs510;'
WRITE(30,*) '  if (seawifs_rrs443 > cff) {'
WRITE(30,*) '    cff=seawifs_rrs443;'
WRITE(30,*) '  }'
WRITE(30,*) '  if (seawifs_rrs490 > cff) {'
WRITE(30,*) '    cff=seawifs_rrs490;'
WRITE(30,*) '  }'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  R4S=cff/seawifs_rrs555;'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  dff=log10(R4S);'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  cff=a_0+(a_1*dff)+(a_2*(pow(dff,2)))+(a_3*(pow(dff,3)))+(a_4*(pow(dff,4)));'
WRITE(30,*) '  printf("%.9f , %.9f \r\n ", cff,dff);'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  OC4_Answer[i_data]=pow(10.0,cff);  /* 10.0**cff */'
WRITE(30,*) '  cff=(OC4_Answer(i_data)-chla(i_data))/chla(i_data);'
WRITE(30,*) '  OC4_MAPE=OC4_MAPE+fabs(cff);'
WRITE(30,*) '  cff=OC4_Answer(i_data)-chla(i_data);'
WRITE(30,*) '  OC4_SSE=OC4_SSE+pow(cff,2);'
WRITE(30,*) '  chla_SSE0=chla_SSE0+pow(chla(i_data),2);'
WRITE(30,*) !  add a blank line

WRITE(30,*) '  printf("%.9f , %.9f \r\n ", chla(i_data),OC4_Answer(i_data));'
WRITE(30,*)    ! add a blank line

WRITE(30,*) 'end % next for'
WRITE(30,*) ! add a blank line

WRITE(30,*) 'OC4_MAPE=100.0*OC4_MAPE/n_data;'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  printf("%.9f , %.9f \r\n ", OC4_SSE,OC4_MAPE);'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) '% Initialize the Node_Values array'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) ! add a blank line
WRITE(30,*) 'Node_Type(1:n_nodes)=Null_Tree_Node; % default ID for NULL nodes'
WRITE(30,*) 'Node_Parameters(1:n_nodes)=0.0;'
WRITE(30,*) 'Node_Values(1:n_nodes)=0.0;'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx '
WRITE(30,*) '% set the Tree_Type and Tree_Parameters'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) ! add a blank line

DO i_node=1,n_nodes
  IF (Node_Type(i_node) .ne. Null_Tree_Node) THEN
    IF (Node_Type(i_node) .eq. 0) THEN
      WRITE(30,*) 'Node_Type(',i_node,')=',Node_Type(i_node),';'
      WRITE(30,*) 'Node_Parameters(',i_node,')=',Node_Parameters(i_node),';'
    ELSE ! IF (Node_Type(i_node) .lt. 0 .or. Node_Type(i_node) .gt. 0) THEN
      WRITE(30,*) 'Node_Type(',i_node,')=',Node_Type(i_node),';'
    END IF
  END IF
END DO
WRITE(30,*) ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx '
WRITE(30,*) '% Loop through the IOCCG data set and compute the GPCODE chla values'
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'// &
              'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(30,*) ! add a blank line

WRITE(20,*) 'fprintf(output_file_pointer,"IOCCG vs GPCODE chla: \r\n");'
WRITE(20,*) ! add a blank line


WRITE(30,*) 'for i_data=1;n_data'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  % fill in the Node_Values tree with the variables and parameters'
WRITE(30,*) '  for i_node=1:n_nodes'
WRITE(30,*) '    if (Node_Type(i_node) == 0)'
WRITE(30,*) '      Node_Values(i_node)=Node_Parameters(i_node);'
WRITE(30,*) '    end'
WRITE(30,*) '    if (Node_Type(i_node) < 0 && Node_Type(i_node) > Null_Tree_Node)'
WRITE(30,*) '      i_variable=-1*Node_Type(i_node);'
WRITE(30,*) '      Node_Values(i_node)=variables(i_variable);'
WRITE(30,*) '    end'
WRITE(30,*) '  end % next for'
WRITE(30,*)

WRITE(30,*) '  % move up the tree structure levels from n_level-1 to 1'
WRITE(20,*) '  i_start=((n_levels-1)^2)-1;'
WRITE(20,*) '  for i_function=i_start;-1:1'
WRITE(30,*) !    add a blank line

WRITE(30,*) '    i_node_left=i_function*2; % set left node index '
WRITE(30,*) '    i_node_right=i_node_left+1; % set right node index'
WRITE(30,*) !    add a blank line

WRITE(30,*) '    left_node_value=Node_Values(i_node_left); % set left node value'
WRITE(30,*) '    right_node_value=Node_Values(i_node_right); % set right node value'
WRITE(30,*) !    add a blank line

WRITE(30,*) '    if (Node_Type(i_function) > 0) % the node is a function'
WRITE(30,*) !      add a blank line

WRITE(30,*) '      switch Node_Type(i_function)'
WRITE(30,*) !      add a blank line

WRITE(30,*) '      case 1 %: [LHS + RHS]'
WRITE(30,*) '        Node_Values(i_function)=left_node_value+right_node_value;'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 2 % : [LHS - RHS]'
WRITE(30,*) '        Node_Values(i_function)=left_node_value-right_node_value;'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 3 %: [LHS * RHS]'
WRITE(30,*) '        if (isnan(left_node_value) || isnan(right_node_value))'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=left_node_value*right_node_value;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 4 % : protected [LHS / RHS]'
WRITE(30,*) '        if (abs(right_node_value) > 0.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=left_node_value/right_node_value;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 5 % : Ivlev Grazing Function: '// &
                     '[(1.0E+00-exp(-abs(LHS*RHS)))]'
WRITE(30,*) '        cff=abs(left_node_value*right_node_value);'
WRITE(30,*) '        if (cff < 100.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00-exp(-1.0E+00*cff);'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 6 % : Michealis-Menton [(abs(RHS) / (abs(LHS) + abs(RHS)))]'
WRITE(30,*) '        cff=abs(left_node_value)+abs(right_node_value);'
WRITE(30,*) '        if (cff > 0.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=abs(right_node_value)/cff;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 7 % : Mayzaud-Poulet Grazing Function: '// &
                     '[abs(LHS*RHS)*(1.0-exp(-abs(LHS*RHS)))]'
WRITE(30,*) '        cff=abs(left_node_value*right_node_value);'
WRITE(30,*) '        if (cff < 100.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=cff*(1.0E+00-exp(-1.0E+00*cff));'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=cff;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 8 % : LHS**RHS'
WRITE(30,*) '        if (isnan(left_node_value) || isnan(right_node_value))'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          if (abs(left_node_value) <= 1.0E-99)'
WRITE(30,*) '            Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '          else'
WRITE(30,*) '            if (abs(right_node_value) <= 1.0E-99)'
WRITE(30,*) '              Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '            else'
WRITE(30,*) '              if (abs(left_node_value-right_node_value) <= 1.0E-99)'
WRITE(30,*) '                Node_Values[i_function]=1.0E+00;'
WRITE(30,*) '              else'
WRITE(30,*) '                Node_Values(i_function)=left_node_value^right_node_value;'
WRITE(30,*) '                if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '                  Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '                end'
WRITE(30,*) '                if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '                  Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '                end'
WRITE(30,*) '              end'
WRITE(30,*) '            end'
WRITE(30,*) '          end'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 9 % : exp(-abs(LHS*RHS)'
WRITE(30,*) '        if (isnan(left_node_value) || isnan(right_node_value))'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          cff=abs(left_node_value*right_node_value);'
WRITE(30,*) '          if (cff < 100.0E+00 )'
WRITE(30,*) '            Node_Values(i_function)=exp(-1.0E+00*cff);'
WRITE(30,*) '          else'
WRITE(30,*) '            Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '          end'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 10 % : min(LHS,RHS)'
WRITE(30,*) '        if (left_node_value < right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=left_node_value;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=right_node_value;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '        case 11 % : max(LHS,RHS)'
WRITE(30,*) '        if (left_node_value > right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=left_node_value;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=right_node_value;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 12 % : if LHS /= 0.0E+00, then RHS, else 0.0E+00'
WRITE(30,*) '        if (left_node_value ~= 0.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=right_node_value;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 13 % : if LHS > RHS, then 1.0E+00, else 0.0E+00'
WRITE(30,*) '        if (left_node_value > right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 14 % : if LHS >= RHS, then 1.0E+00, else 0.0E+00'
WRITE(30,*) '        if (left_node_value >= right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 15 % : if LHS < RHS, then 1.0E+00, else 0.0E+00'
WRITE(30,*) '        if (left_node_value < right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 16 % : if LHS <= RHS, then 1.0E+00, else 0.0E+00'
WRITE(30,*) '        if (left_node_value <= right_node_value)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 17 % : exp(LHS)'
WRITE(30,*) '        Node_Values(i_function)=exp(left_node_value);'
WRITE(30,*) '        if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '        end'
WRITE(30,*) '        if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 18 % : exp(RHS)'
WRITE(30,*) '        Node_Values(i_function)=exp(right_node_value);'
WRITE(30,*) '        if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '        end'
WRITE(30,*) '        if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 19 % : exp(-LHS)'
WRITE(30,*) '        Node_Values(i_function)=exp(-1.0E+00*left_node_value);'
WRITE(30,*) '        if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '        end'
WRITE(30,*) '        if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 20 % : exp(-RHS)'
WRITE(30,*) '        Node_Values(i_function)=exp(-1.0E+00*right_node_value);'
WRITE(30,*) '        if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '        end'
WRITE(30,*) '        if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '          Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 21 % : log_|RHS|(|LHS|))'
WRITE(30,*) '        if (left_hand_node > 0.0E+00 & right_hand_node > 0.0E+00)'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          cff=log10(abs(right_node_value));'
WRITE(30,*) '          Node_Values(i_function)=log10(abs(left_node_value))/cff;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 22 % : LHS'
WRITE(30,*) '        if (isnan(left_node_value))'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          Node_Values(i_function)=left_node_value;'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      case 23 % : LHS**2'
WRITE(30,*) '        if (isnan(left_node_value))'
WRITE(30,*) '          Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '        else'
WRITE(30,*) '          if (abs(left_node_value) <= 1.0E-99)'
WRITE(30,*) '            Node_Values(i_function)=0.0E+00;'
WRITE(30,*) '          else'
WRITE(30,*) '            Node_Values(i_function)=left_node_value^2;'
WRITE(30,*) '            if (Node_Values(i_function) > 1.0E+19)'
WRITE(30,*) '              Node_Values(i_function)=1.0E+19;'
WRITE(30,*) '            end'
WRITE(30,*) '            if (Node_Values(i_function) < 1.0E-19)'
WRITE(30,*) '              Node_Values(i_function)=1.0E-19;'
WRITE(30,*) '            end'
WRITE(30,*) '          end'
WRITE(30,*) '        end'
WRITE(30,*) !        add a blank line

WRITE(30,*) '      otherwise % default/error detection:'
WRITE(30,*) "        disp('default case')"
WRITE(30,*) "        disp('Wrong case number in tree evaluations')"
WRITE(30,*) "        disp('Bad case number')"
WRITE(30,*) '      end'
WRITE(30,*) !      add a blank line

WRITE(30,*) '    end % end if'
WRITE(30,*) !    add a blank line

WRITE(30,*) '  end % next for in tree node calculations'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  % GPCODE Matlab-code answer is:'
WRITE(30,*) '  GP_Answer(i_data)=Node_Values(1);'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  cff=(GP_Answer(i_data)-chla(i_data))/chla(i_data);'
WRITE(30,*) '  GP_MAPE=GP_MAPE+fabs(cff);'
WRITE(30,*) '  cff=GP_Answer(i_data)-chla(i_data);'
WRITE(30,*) '  GP_SSE=GP_SSE+pow(cff,2);'
WRITE(30,*) ! add a blank line

WRITE(30,*) '  v1=chla(i_data);'
WRITE(30,*) '  v2=OC4_Answer(i_data);'
WRITE(30,*) '  v3=GP_Answer(i_data);'
WRITE(30,*) '  fprintf(output_file_pointer,"%.9f, %.9f, %.9f \r\n ", v1, v2, v3;'
WRITE(30,*) !  add a blank line

WRITE(30,*) 'end % next for in tree node calculations'
WRITE(30,*) ! add a blank line

WRITE(30,*) 'fclose(input_file_pointer) % close the input IOCCG data file'
WRITE(30,*) ! add a blank line

WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
WRITE(20,*) '% end Matlab program: '//Code_Root_File_Name(1:LEN_TRIM(Code_Root_File_Name))
WRITE(30,*) '% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

CLOSE(30)

1010 FORMAT(' ',' Error opening file: ',A120,' IOSTAT = ',I6)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! RETURN END
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
RETURN
END SUBROUTINE Write_GPCODE_Tree_Test_Matlab_code
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
INCLUDE '/Users/jmoisan/NEPAC/GPCODE_analysis_code/GP_tree_calculation_FULL_CODE.f90'
