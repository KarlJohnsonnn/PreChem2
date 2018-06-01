!===================================================================================!
!                                                                                   !
!                     Module for reading the chemical system                        !
!                 file .sys and building the coeficient matrices                    !
!                                                                                   !
!===================================================================================!
MODULE Chemsys_Mod
  !
  USE Kind_Mod
  USE String_Mod
  USE LexicalStringSort
  USE hashtbl
  USE mo_unirnk
  USE mo_control
  USE mo_reac
  USE InputTool_Mod
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PARAMETER ::  maxLENinActDuct=9
  ! 
  TYPE Duct_T
    CHARACTER(LenName) :: Species=''
    CHARACTER(LenType) :: Type
    REAL(dp)           :: Koeff
    INTEGER            :: iSpecies=0
  END TYPE Duct_T

  TYPE Special_T
    INTEGER                         :: nVariables = 0
    INTEGER,            ALLOCATABLE :: iVariables(:)
    CHARACTER(LenName), ALLOCATABLE :: cVariables(:)
    CHARACTER(LenLine)              :: Formula = ''
    LOGICAL                         :: Temp = .FALSE.
  END TYPE Special_T
  !
  ! LIST FORM
  TYPE Reaction_T
    CHARACTER(LenType)        :: Type, TypeConstant
    CHARACTER(LenName)        :: Comment
    CHARACTER(LenLine)        :: Line1, Line2, Line3, Line4
    CHARACTER(LenName)        :: Factor
    TYPE(Duct_T),     POINTER :: Educt(:)=>NULL(), Product(:)=>NULL()
    REAL(dp),     ALLOCATABLE :: Constants(:)
    TYPE(Duct_T),     POINTER :: InActEduct(:)=>NULL(), InActProduct(:)=>NULL()
    TYPE(Special_T)           :: Special
    INTEGER                   :: nInActEd=0, nInActPro=0
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Reaction_T
  !
  ! ARRAY FORM
  TYPE ReactionStruct_T
    CHARACTER(LenType)              :: Type,  TypeConstant
    CHARACTER(LenLine)              :: Line1='' , Line2='' , Line3='', Line4=''
    LOGICAL                         :: bR = .FALSE. , brX = .FALSE. 
    CHARACTER(LenName)              :: Factor = ''
    CHARACTER(LenName)              :: Comment = ''
    CHARACTER(2)                    :: direction = ''
    REAL(dp)                        :: SumAqCoef     
    TYPE(Special_T)                 :: Special
    TYPE(Duct_T)  ,     ALLOCATABLE :: Educt(:), Product(:)
    REAL(dp),           ALLOCATABLE :: Constants(:)
    REAL(dp),           ALLOCATABLE :: LowConst(:), HighConst(:), TroeConst(:) ! combustion press dep reactions
    REAL(dp),           ALLOCATABLE :: InActEduct(:), InActProduct(:)
    INTEGER                         :: nInActEd = 0, nInActPro = 0, nActEd = 0, nActPro = 0
    INTEGER                         :: nConst = 0
    INTEGER                         :: HenrySpc = 0
    LOGICAL                         :: TB = .FALSE. , TBextra=.FALSE.
    INTEGER,            ALLOCATABLE :: TBidx(:)
    CHARACTER(LenName), ALLOCATABLE :: TBspc(:)
    REAL(dp),           ALLOCATABLE :: TBalpha(:)
    CHARACTER(LenName), ALLOCATABLE :: InActEductSpc(:), InActProductSpc(:)
  END TYPE ReactionStruct_T
  !
  !
  TYPE ListReaction_T
    TYPE(Reaction_T), POINTER :: Start=>NULL()
    TYPE(Reaction_T), POINTER :: End=>NULL()
    INTEGER :: LenList=0
  END TYPE ListReaction_T
  !
  TYPE Species_T
    CHARACTER(LenName) :: Species=''
    LOGICAL            :: isHenry=.FALSE.
    REAL(dp)     :: Hf=0.0d0, Gf=0.0d0, Cp=0.0d0
  END TYPE Species_T

  

  
  TYPE(Reaction_T), POINTER   :: System
  TYPE(ListReaction_T), SAVE  :: ListRGas, ListRHenry, ListRAqua,        &
  &                              ListRDiss, ListRSolid, ListRPartic,     &
  &                              ListRMicro
  !
  TYPE(hash_tbl_sll)          :: ListAqua, ListGas, ListSolid,           &
  &                              ListPartic, ListNonReac, ListAtoms
  TYPE(hash_tbl_sll)          :: ListFamilies
  !
  TYPE(Species_T), ALLOCATABLE, TARGET :: ListAqua2(:), ListGas2(:),     &
  &                               ListSolid2(:), ListPartic2(:), &
  &                               ListNonReac2(:)
  INTEGER :: InputUnit=10
  INTEGER, PARAMETER :: MaxEduct=10
  INTEGER, PARAMETER :: MaxProduct=10
  !
  CHARACTER(33), PARAMETER :: SetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZapsc[]()=+*'
  CHARACTER(14), PARAMETER :: SetConstants='ABCDEFGKINMOR/'
  CHARACTER(12), PARAMETER :: SetExponent='0123456789+-'

  TYPE Element_T
    CHARACTER(5) :: Element=''
  END TYPE Element_T

  TYPE(Element_T) :: Elements(11)=(/       &
  &                    Element_T('(')      &
  &                    ,Element_T(')')     &
  &                    ,Element_T('exp')   &
  &                    ,Element_T('+')     &
  &                    ,Element_T('-')     &
  &                    ,Element_T('*')     &
  &                    ,Element_T('/')     &
  &                    ,Element_T('**')    &
  &                    ,Element_T('abs')   &
  &                    ,Element_T('sqrt')  &
  &                    ,Element_T('log')   /)
 
  
  INTEGER :: nsr                        ! # activ species + all Reactions

  INTEGER :: UnitGas=0
  INTEGER :: UnitAqua=0
  !
  CHARACTER(20) :: Filename
  CHARACTER(20) :: IniName
  !
  REAL(dp), PARAMETER :: RGas=8.3145d0
  REAL(dp), PARAMETER :: TRef=280.0d0 !298.15d0
  !
  TYPE(Reaction_T), POINTER :: Current
  TYPE(ReactionStruct_T), ALLOCATABLE :: ReactionSystem(:)
  TYPE(ListReaction_T), ALLOCATABLE :: CompleteReactionList(:)
  !
  !
  REAL(dp), ALLOCATABLE :: Emis(:)          & ! emission values
  !&                            , InitValAct(:)    & ! initial values activ spc
  &                      , InitValInAct(:)    ! initial values inactiv spc
  !
  !
  CHARACTER(LenName), ALLOCATABLE :: RO2spcG(:) , RO2spcA(:)
  INTEGER, ALLOCATABLE :: RO2idxG(:) , RO2idxA(:)
  !
  !
  !REAL(dp) :: aH2O
  !
  REAL(dp), ALLOCATABLE :: sumBAT(:)         ! sum_j=1,n_s (b_ij-a_ij),  i el. N_R

  INTEGER :: fNumber = 0
  !
  CONTAINS
  ! ------------------------------------
  ! -----------SUBROUTINEN--------------
  ! ------------------------------------
  !
  SUBROUTINE SortReactionList(ReacStructOut,ReacStructIn)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructIn(:)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructOut(:)
    !
    INTEGER :: i
    CHARACTER(20), ALLOCATABLE :: ReacTypeSorted(:)
    INTEGER, ALLOCATABLE :: iReacTypeSorted(:)
    !
    ! sort the reaction list --> TypeConstant
    ALLOCATE(ReacTypeSorted(neq))
    ALLOCATE(iReacTypeSorted(neq))  
    DO i=1,neq
      ReacTypeSorted(i)=ReacStructIn(i)%Type
    END DO
    CALL StringSort(ReacTypeSorted,iReacTypeSorted)
    ALLOCATE(ReacStructOut(neq))
    !
    DO i=1,SIZE(ReacStructIN)
      ReacStructOut(i)=ReacStructIn(iReacTypeSorted(i))
    END DO
    DEALLOCATE(ReacStructIn)
    DEALLOCATE(iReacTypeSorted)
  END SUBROUTINE SortReactionList
  !
  !
  SUBROUTINE ReadSpecies(Out)
    LOGICAL :: Out
    !
    CHARACTER(LenName) :: Species
    CHARACTER(LenType) :: Type
    INTEGER :: Pos

    READ(InputUnit,'(a100)',END=1) Species

    DO
      Pos = SCAN( Species , "'" )
      IF ( Pos > 0 ) THEN
        Species(Pos:) = Species(Pos+1:)
      ELSE
        EXIT
      END IF
    END DO
    IF ( Species /= '' ) CALL InsertSpecies(Species,Type)

    Out = .FALSE.
    GO TO 999

  1 CONTINUE
   
    Out = .TRUE.
999 CONTINUE
  END SUBROUTINE ReadSpecies

  !
  SUBROUTINE ReadReaction(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: CLASS
    CHARACTER(40) :: TypeR
    INTEGER :: idxFAC
    
    iLine = 0
    DO

      READ( InputUnit , '(A400)' , IOSTAT=is ) LocLine
      idxFAC = INDEX(LocLine,'$')

      IF ( idxFAC > 0 ) THEN
        SELECT CASE (TRIM(LocLine(idxFAC:)))
          CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','$O2O2','$aH2O','$+M','$(+M)','$RO2','$RO2aq')
            IF ( TRIM(LocLine(idxFAC:)) == '$H2'    ) nr_FAC_H2    = nr_FAC_H2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2N2'  ) nr_FAC_O2N2  = nr_FAC_O2N2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$M'     ) nr_FAC_M     = nr_FAC_M     + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2'    ) nr_FAC_O2    = nr_FAC_O2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$N2'    ) nr_FAC_N2    = nr_FAC_N2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$H2O'   ) nr_FAC_H2O   = nr_FAC_H2O   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2O2'  ) nr_FAC_O2O2  = nr_FAC_O2O2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2'   ) nr_FAC_RO2   = nr_FAC_RO2   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2aq' ) nr_FAC_RO2aq = nr_FAC_RO2aq + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$aH2O'  ) nr_FAC_aH2O  = nr_FAC_aH2O  + 2
            nr_FACTOR = nr_FACTOR + 1
          CASE DEFAULT
            WRITE(*,*) '  Unknown FACTOR:  ', TRIM(LocLine(idxFAC:)), '   at Line:  ???'
        END SELECT
      END IF

      IF ( ABS(is) > 0 ) EXIT
      
      ! if no comment or blank line then
      IF ( ADJUSTL(LocLine(1:1)) /= '#'       .AND.  &
      &    ADJUSTL(LocLine(1:7)) /= 'COMMENT' .AND.  &
      &    LEN(TRIM(LocLine)) > 0 ) THEN
        iLine = iLine + 1
        Line(iLine) = LocLine
        IF ( iLine == 4 ) EXIT
      END IF

    END DO

    IF ( iLine >= 3 ) THEN
      Pos = SCAN(Line(1),'#')
      IF ( Pos > 0 ) Line(1) = Line(1)(:Pos-1)
      
      ! if new reaction line starts go one line back
      IF ( INDEX(Line(4),'CLASS') > 0 ) BACKSPACE(InputUnit)
      
      ! read the reaction TYPE
      PosColon = Index(Line(1),':')
      CLASS    = ADJUSTL(Line(1)(PosColon+1:))
      
      ! count the number of each reaction type
      SELECT CASE (CLASS)
        
        CASE ('GAS')      ! gaseous phase reactions

          nr_gas = nr_gas + 1
          CALL InsertReaction( ListRGas , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3','PHOTAB','PHOTABC','PHOTMCM')
              nr_G_photo = nr_G_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
              IF ( TypeR == 'PHOTO'   ) nr_PHOTOkpp  = nr_PHOTOkpp  + 1
              IF ( TypeR == 'PHOTO2'  ) nr_PHOTO2kpp = nr_PHOTO2kpp + 1
              IF ( TypeR == 'PHOTO3'  ) nr_PHOTO3kpp = nr_PHOTO3kpp + 1
            CASE ('CONST')
              nr_G_const = nr_G_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','TEMP1','TEMP2','TEMP3','TEMP4')
              nr_G_temp  = nr_G_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('TROE','TROEF','TROEQ','TROEQF','TROEXP','TROEMCM')
              nr_G_troe = nr_G_troe + 1
              IF ( TypeR == 'TROE'    ) nr_TROE    = nr_TROE    + 1
              IF ( TypeR == 'TROEF'   ) nr_TROEf   = nr_TROEf   + 1
              IF ( TypeR == 'TROEQ'   ) nr_TROEq   = nr_TROEq   + 1
              IF ( TypeR == 'TROEQF'  ) nr_TROEqf  = nr_TROEqf  + 1
              IF ( TypeR == 'TROEXP'  ) nr_TROExp  = nr_TROExp  + 1
              IF ( TypeR == 'TROEMCM' ) nr_TROEmcm = nr_TROEmcm + 1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM',  &
            &     'SPEC2MCM','SPEC3MCM','SPEC4MCM','SPEC5MCM', &
            &     'SPEC6MCM','SPEC7MCM','SPEC8MCM','SPEC9MCM'  )
              nr_G_spec = nr_G_spec + 1
              IF ( TypeR == 'SPEC1' ) nr_SPEC1 = nr_SPEC1 + 1
              IF ( TypeR == 'SPEC2' ) nr_SPEC2 = nr_SPEC2 + 1
              IF ( TypeR == 'SPEC3' ) nr_SPEC3 = nr_SPEC3 + 1
              IF ( TypeR == 'SPEC4' ) nr_SPEC4 = nr_SPEC4 + 1
              IF ( TypeR == 'SPEC1MCM' ) nr_SPEC1mcm = nr_SPEC1mcm + 1
              IF ( TypeR == 'SPEC2MCM' ) nr_SPEC2mcm = nr_SPEC2mcm + 1
              IF ( TypeR == 'SPEC3MCM' ) nr_SPEC3mcm = nr_SPEC3mcm + 1
              IF ( TypeR == 'SPEC4MCM' ) nr_SPEC4mcm = nr_SPEC4mcm + 1
              IF ( TypeR == 'SPEC5MCM' ) nr_SPEC5mcm = nr_SPEC5mcm + 1
              IF ( TypeR == 'SPEC6MCM' ) nr_SPEC6mcm = nr_SPEC6mcm + 1
              IF ( TypeR == 'SPEC7MCM' ) nr_SPEC7mcm = nr_SPEC7mcm + 1
              IF ( TypeR == 'SPEC8MCM' ) nr_SPEC8mcm = nr_SPEC8mcm + 1
              IF ( TypeR == 'SPEC9MCM' ) nr_SPEC9mcm = nr_SPEC9mcm + 1
            CASE ('S4H2O')
              nr_S4H2O = nr_S4H2O + 1
            CASE ('T1H2O')
              nr_T1H2O = nr_T1H2O + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_G_special = nr_G_special + 1
            CASE ('HOM1')
              nr_HOM1 = nr_HOM1 + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown gaseous reaction: ', TypeR
          END SELECT

        CASE ('HENRY')        ! phase transfer pseudo-reactions 

          nr_henry = nr_henry + 1
          CALL InsertReaction( ListRHenry , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('TEMP','TEMP1','TEMP2','TEMP3')
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
            CASE ('CONST')
              nr_CONST = nr_CONST + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_H_special = nr_H_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown phase transfer reaction: ', TypeR
          END SELECT

        CASE ('AQUA')         ! aquatic phase reactions 

          nr_aqua = nr_aqua + 1
          CALL InsertReaction( ListRAqua , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTAB','PHOTABC','PHOTMCM')
              nr_A_photo = nr_A_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
            CASE ('CONST')
              nr_A_const = nr_A_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','Temp1''TEMP2','TEMP3','TEMP4')
              nr_A_temp  = nr_A_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('ASPEC1','ASPEC2','ASPEC3','ASPEC4')
              nr_A_spec  = nr_A_spec + 1
              IF ( TypeR == 'ASPEC1' ) nr_ASPEC1 = nr_ASPEC1 + 1
              IF ( TypeR == 'ASPEC2' ) nr_ASPEC2 = nr_ASPEC2 + 1
              IF ( TypeR == 'ASPEC3' ) nr_ASPEC3 = nr_ASPEC3 + 1
              IF ( TypeR == 'ASPEC4' ) nr_ASPEC4 = nr_ASPEC4 + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_A_special = nr_A_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown aqueous reaction: ', TypeR
          END SELECT

        CASE ('DISS')        ! fast aquatic phase equil. reactions 

          nr_diss = nr_diss + 1
          CALL InsertReaction( ListRDiss , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('DCONST','DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5','MESKHIDZE')
              IF ( TypeR == 'DCONST'    ) nr_DCONST = nr_DCONST + 1
              IF ( TypeR == 'DTEMP'     ) nr_DTEMP  = nr_DTEMP  + 1
              IF ( TypeR == 'DTEMP2'    ) nr_DTEMP2 = nr_DTEMP2 + 1
              IF ( TypeR == 'DTEMP3'    ) nr_DTEMP3 = nr_DTEMP3 + 1
              IF ( TypeR == 'DTEMP4'    ) nr_DTEMP4 = nr_DTEMP4 + 1
              IF ( TypeR == 'DTEMP5'    ) nr_DTEMP5 = nr_DTEMP5 + 1
              IF ( TypeR == 'MESKHIDZE' ) nr_Meskhidze = nr_Meskhidze + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_D_special = nr_D_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown dissociation reaction: ', TypeR
          END SELECT

        CASE ('SOLID')

          nr_solid = nr_solid + 1
          CALL InsertReaction( ListRSolid , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('EQUI')
              nr_S_equi = nr_S_equi + 1
            CASE ('DTEMP3')
              nr_S_temp = nr_S_temp + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_S_special = nr_S_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown solid reaction: ', TypeR
          END SELECT

        CASE ('PARTI')

          nr_parti = nr_parti + 1
          CALL InsertReaction( ListRPartic , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_P_special = nr_P_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown particle reaction: ', TypeR
          END SELECT

        CASE ('MICROPHYS')

          nr_micphys = nr_micphys + 1
          CALL InsertReaction( ListRMicro , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_M_special = nr_M_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown microphysic reaction: ', TypeR
          END SELECT

        CASE DEFAULT
          WRITE(*,*) '  Unknown reaction CLASS: ', CLASS
          STOP
      END SELECT

      Out = .FALSE.
    ELSE
      Out = .TRUE.
    END IF

  END SUBROUTINE ReadReaction
  

  
  SUBROUTINE PrintSpecies(ListName,Unit,phs)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: Unit
    INTEGER, ALLOCATABLE :: spclist(:)
    CHARACTER(1), OPTIONAL :: phs
    !
    INTEGER :: i

    DO i=1,SIZE(ListName)
      WRITE(Unit,*) "'"//TRIM(ListName(i)%Species)//"'"
    END DO
  END SUBROUTINE PrintSpecies
  !
  !
  SUBROUTINE SpcIdx(ListName,idx)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: idx
    !
    INTEGER :: i
    !
    DO i=1,SIZE(ListName)
      IF (i==idx) THEN
        WRITE(*,*) "'"//TRIM(ListName(idx)%Species)//"'"
      END IF
    END DO
  END SUBROUTINE SpcIdx
  !
  !
  SUBROUTINE PrintHeadSpecies(Filename,Unit)
    INTEGER :: Unit
    CHARACTER(*) :: Filename
    !
    CHARACTER(8) :: Date
    CHARACTER(10) :: Time
    INTEGER(8) :: Value(8)
    !
    CALL DATE_AND_TIME(Date,Time,VALUES=Value)
    !
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ' ========  0-dim Simulation of chemical mechanisms  ========'
    WRITE(Unit,*) ' ========     Output -  Chemical Reaction Data      ========'
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' Created:             ', Date(7:8),'.',Date(5:6),'.',Date(1:4)
    WRITE(Unit,*) ' Chemical Mechanism:  ', TRIM(ADJUSTL(FileName))
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================     Units         ======================='
    WRITE(Unit,*) ''
    IF (UnitGas==0) THEN
      WRITE(Unit,*) ' Gas Phase Units:     molec/cm3'
    ELSE
      WRITE(Unit,*) ' Gas Phase Units:     mol/m3'
    END IF
    IF (UnitAqua==0) THEN
      WRITE(Unit,*) ' Aqueous Phase Units: mol/l'
    END IF
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================    Numbers        ======================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ns_GAS &
                 +ns_AQUA &
                 +ns_PARTI &
                 +ns_KAT &
                 +ns_SOLID,  '     Number of Species'
    WRITE(Unit,*) ns_GAS,    '     No. of gaseous species'
    WRITE(Unit,*) ns_AQUA,   '     No. of aqueous species'
    WRITE(Unit,*) ns_PARTI, '     No. of particular species'
    WRITE(Unit,*) ns_SOLID,  '     No. of solid   species'
    WRITE(Unit,*) ns_KAT,'     Number of Non-reactive Species '
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================   Species Names   ======================='
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadSpecies
  !
  !
  SUBROUTINE PrintFinalReactions(Unit)
    INTEGER :: Unit
    !
    WRITE(Unit,*) ''
    WRITE(Unit,*) ''
    WRITE(Unit,*) '========================================================='
    WRITE(Unit,*) '========              End  TAPE2                 ========'
    WRITE(Unit,*) '========     M3TRAS:  Chemical Reaction Data     ========'
    WRITE(Unit,*) '========================================================='
  END SUBROUTINE PrintFinalReactions
  !
  !
  SUBROUTINE PrintHeadReactions(Unit)
    INTEGER :: Unit
  
    nr    = nr_gas  + 2*nr_henry + nr_aqua   &
    &     + 2*nr_diss + nr_solid + nr_parti &
    &     + nr_micphys

    nr_D_Temp = nr_DTEMP  + nr_DTEMP2 &
    &         + nr_DTEMP3 + nr_DTEMP5
    
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ================   Description of Reactions   =============='
    WRITE(Unit,*) ''
    WRITE(Unit,*) nr,          '        NREAK  : Number of Reactions'
    WRITE(Unit,*) nr_gas,      '        NGAS   : Gas phase reactions'
    WRITE(Unit,*) nr_G_photo,  '           Gaseous PHOTO - type reactions'
    WRITE(Unit,*) nr_G_const,  '           Gaseous CONST - type reactions'
    WRITE(Unit,*) nr_G_temp,   '           Gaseous TEMP - type reactions'
    WRITE(Unit,*) nr_SimpTB,   '           Gaseous Simple three-body - type reactions'
    WRITE(Unit,*) nr_G_lind,   '           Gaseous Lindemann - type reactions'
    WRITE(Unit,*) nr_G_troe,   '           Gaseous TROE - type reactions'
    WRITE(Unit,*) nr_G_spec,   '           Gaseous SPEC - type reactions'
    WRITE(Unit,*) nr_G_special,'           Gaseous SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_henry,    '        NHENRY : Henry Equilib. reactions'
    WRITE(Unit,*) nr_diss,     '        NDISS  : Dissociation reactions'
    WRITE(Unit,*) nr_DCONST,   '           Aqueous DCONST - type reactions'
    WRITE(Unit,*) nr_D_TEMP,   '           Aqueous DTEMP - type reactions'
    WRITE(Unit,*) nr_D_special,'           Aqueous SPECIAL formula reactions'
    WRITE(Unit,*) nr_aqua,     '        NAQUA  : Aquatic Equilib. reactions'
    WRITE(Unit,*) nr_A_photo,  '           Aqueous PHOTO - type reactions'
    WRITE(Unit,*) nr_A_const,  '           Aqueous CONST - type reactions'
    WRITE(Unit,*) nr_A_temp,   '           Aqueous TEMP - type reactions'
    WRITE(Unit,*) nr_A_spec,   '           Aqueous SPEC - type reactions'
    WRITE(Unit,*) nr_A_special,'           Aqueous SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_parti,   '        NPARTI : Particulare reactions   '
    WRITE(Unit,*) nr_P_special,'           Partic SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_solid,    '        NSOLID : Solid Equilib. reactions'
    WRITE(Unit,*) nr_S_temp,   '           Solid DTEMP3 - type reactions'
    WRITE(Unit,*) nr_S_equi,   '           Solid EQUI - type reactions'
    WRITE(Unit,*) nr_S_spec,   '           Solid SPEC - type reactions'
    WRITE(Unit,*) nr_S_special,'           Solid SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_micphys,  '        NMICRO : Microphysical reactions'
    WRITE(Unit,*) nr_M_special,'           Micro SPECIAL formula - type reactions'
    WRITE(Unit,*)
    WRITE(Unit,*) ' ======================  Reactions   ========================'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadReactions
 

  FUNCTION RemoveSpaces(String) RESULT(StringOut)
    ! replaces multiple spaces in string by one space
    CHARACTER(*), INTENT(IN) :: String
    CHARACTER(LEN(String))   :: StringOUT

    INTEGER :: i

    StringOut = TRIM(String)

    DO
      i = INDEX(TRIM(StringOut),'  ')
      IF ( i == 0 ) EXIT
      StringOut = TRIM(ADJUSTL(StringOut(:i-1))) &
      &           //' '//                        &
      &           TRIM(ADJUSTL(StringOut(i+1:)))
    END DO

  END FUNCTION RemoveSpaces

  !
  SUBROUTINE Print_ChemFile(RS,File,Unit,CK)
    ! IN:
    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
    CHARACTER(*) :: File
    INTEGER      :: Unit
    LOGICAL      :: CK
    ! TEMP:
    INTEGER      :: io_stat
    INTEGER      :: i,j,m,iR
    INTEGER      :: nEduct,nProd
    TYPE(Duct_T) :: ActiveEduct(30)
    TYPE(Duct_T) :: ActiveProduct(30)
    !
    INTEGER      :: nnzA, nnzB
    !
    INTEGER,  ALLOCATABLE :: tmpIdx(:)
    REAL(dp), ALLOCATABLE :: tmpVal(:)
    INTEGER,  ALLOCATABLE :: permutation(:)
    INTEGER :: newLen

    
     
    !-----------------------------------------------------------------------
    ! --- Build the reaction system
    !-----------------------------------------------------------------------
  
    CALL AllListsToArray( RS            &
    &                   , ListRGas    , ListRHenry  &
    &                   , ListRAqua   , ListRDiss   &
    &                   , ListRSolid  , ListRPartic &
    &                   , ListRMicro  )
  
    

    OPEN(unit=Unit, file=File, status='replace', action='write', access='sequential', iostat=io_stat)
    IF ( io_stat /= 0 ) WRITE(*,*) '  ERROR creating chem-file :: ',io_stat
    REWIND(ChemUnit)
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies ( File , Unit )
    
    IF ( ns_GAS   > 0 ) CALL PrintSpecies( ListGas2     , Unit )
    IF ( ns_AQUA  > 0 ) CALL PrintSpecies( ListAqua2    , Unit )
    IF ( ns_SOLID > 0 ) CALL PrintSpecies( ListSolid2   , Unit )
    IF ( ns_PARTI > 0 ) CALL PrintSpecies( ListPartic2  , Unit )
    IF ( ns_KAT   > 0 ) CALL PrintSpecies( ListNonReac2 , Unit )
     
    CALL PrintHeadReactions( Unit )


    DO iR=1,neq
      ! count activ educts in reaction iR
      nEduct = 0
      !print*, 'DEBUG::chemsys    sizeRSe,p= ',iR,SIZE(RS(iR)%Educt),SIZE(RS(iR)%Product)
      !print*, 'DEBUG::chemsys    reaktion = ',TRIM(RS(iR)%Line1)
      DO i=1,SIZE(RS(iR)%Educt)
        SELECT CASE(RS(iR)%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nEduct = nEduct + 1
            ActiveEduct(nEduct) = RS(iR)%Educt(i)
            !print*, 'debug::chemssys   ActiveEduct(nEduct)=',ActiveEduct(nEduct)
        END SELECT
      END DO
      ! count activ products in reaction iR
      nProd = 0
      DO i=1,SIZE(RS(iR)%Product)
        SELECT CASE(RS(iR)%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic','Micro','GAS')
            nProd = nProd + 1
            ActiveProduct(nProd) = RS(iR)%Product(i)
            !print*, 'debug::chemssys   ActiveProduct(nProd)=',ActiveProduct(nProd)
        END SELECT
      END DO
   
      !iR = iR + 1
      WRITE(Unit,*)
      WRITE(Unit,'(A,I6,A)') '#-----------', iR ,'. Reaction ----------- '
    
      WRITE(Unit,*) TRIM(RS(iR)%Type)//'   '//TRIM(RS(iR)%TypeConstant)
     
      WRITE(Unit,*) SIZE(RS(iR)%Educt), SIZE(RS(iR)%Product), nEduct, nProd
     

      ! ----------------------------------------------------
      ! SpeziesIndx Edukt=> 1:#Edukt von Reaktion iR
      ! SpeziesIndx Produkt=> #Edukt+1:#Edukt+#Produkt von Reaktion iR
      ! #aktiver Stoffe der Reaktion
      WRITE(Unit,*) (PositionSpeciesAll(RS(iR)%Educt(i)%Species),  &
      &             i=1,SIZE(RS(iR)%Educt)),                       &
      &             (PositionSpeciesAll(RS(iR)%Product(i)%Species),&
      &             i=1,SIZE(RS(iR)%Product)),                     &
      &             nEduct+nProd
      ! 
      !----------------------------------------------------
      ! Tupel: (SpeziesIndex,-Koeffzien) für alle aktiven Edukte (links)
      ! Tupel: (SpeziesIndex,+Koeffzien) für alle aktiven Produkte (rechts)
      WRITE(Unit,'(*(7X,I5,3X,F6.3))', ADVANCE='NO')                     &
      &                   ( PositionSpeciesAll(ActiveEduct(i)%Species)   &
      &                  ,  -ActiveEduct(i)%Koeff,i=1,nEduct )   &
      &                  ,( PositionSpeciesAll(ActiveProduct(i)%Species) &
      &                  ,   ActiveProduct(i)%Koeff,i=1,nProd )
      WRITE(Unit,*)
      !
      IF (RS(iR)%TypeConstant=='SPECIAL') THEN
        WRITE(Unit,*) TRIM(RS(iR)%Line3)
      ELSE
        WRITE(Unit,*) SIZE(RS(iR)%Constants), RS(iR)%Constants
      END IF
      !
      ! #Reaktionskonstanten, Reaktionskonstanten 1:#

      IF (RS(iR)%Factor /= '') WRITE(Unit,*) 'FACTOR:  ',RS(iR)%Factor

      SELECT CASE (RS(iR)%Factor)
        CASE ('$RO2');   hasRO2   = .TRUE.
        CASE ('$RO2aq'); hasRO2aq = .TRUE.
      END SELECT
      
      IF (CK) WRITE(Unit,*) 'EXTRA1:  ',ADJUSTL(TRIM(RS(iR)%Line2))
      IF (CK) WRITE(Unit,*) 'EXTRA2:  ',ADJUSTL(TRIM(RS(iR)%Line3))
    END DO
  

    CALL PrintFinalReactions( Unit )
    CLOSE(ChemUnit)
  END SUBROUTINE Print_ChemFile


  !
  SUBROUTINE InsertReaction(List,Line,TypeR)
    TYPE(ListReaction_T) :: List
    CHARACTER(*) :: Line(1:4)
    CHARACTER(*) :: TypeR
    !
    INTEGER :: PosColon,PosEqual,PosFactor,PosSpecial
    CHARACTER(LenLine) :: Left,Right
    TYPE(Reaction_T), POINTER :: Reaction
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)

    INTEGER :: i
    !
    IF (ASSOCIATED(List%Start)) THEN
      ALLOCATE(List%End%Next)
      List%End => List%End%Next
    ELSE
      ALLOCATE(List%End)
      List%Start => List%End
    END IF
    List%LenList = List%LenList + 1
    Reaction => List%End
    PosColon = Index(Line(1),':')
    Reaction%Type  = ADJUSTL(Line(1)(PosColon+1:))
    Reaction%Line1 = TRIM(Line(2))
    Reaction%Line3 = TRIM(Line(3))
    PosEqual = Index(Reaction%Line1,' = ')
    IF (PosEqual==0) THEN
      WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A,I0,A)') 'ERROR: Missing seperator " = " in reaction ',fNumber,':  '//TRIM(Line(2)) 
      WRITE(*,*); WRITE(*,*)
      STOP 
    ELSE
      fNumber = fNumber + 1
    END IF
    
    ! extract educts
    Left = Reaction%Line1(:PosEqual-1)
    CALL ExtractSpecies( Left, Reaction%Educt,     &
    &                    Reaction%InActEduct,      &
    &                    Reaction%nInActEd         )  ! geändert
    
    ! extract products
    Right = Reaction%Line1(PosEqual+3:)
    CALL ExtractSpecies( Right, Reaction%Product,  &
    &                    Reaction%InActProduct,    &
    &                    Reaction%nInActPro        )! geändert


    IF ( INDEX(Line(3),'SPECIAL') > 0 ) THEN
      Ducts = [Reaction%Educt(:)%Species , Reaction%Product(:)%Species]
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant,Reaction%Special)
    ELSE
      ! extract constants
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant)
    END IF
    Reaction%Line2  = Line(3)
    PosFactor = INDEX(Line(4),'FACTOR: ')

    IF (PosFactor > 0) THEN
      Reaction%Factor = TRIM(Line(4)(PosFactor+8:)) 
    ELSE
      Reaction%Factor = ''
    END IF
    !
    TypeR = ADJUSTL(Reaction%TypeConstant)
  END SUBROUTINE InsertReaction
  !
  !
  SUBROUTINE ReadUnits
    !
    INTEGER :: Pos
    CHARACTER(LenLine) :: LocLine
    INTEGER :: ios
    CHARACTER(400) :: iom
    LOGICAL :: comments
    !
    REWIND(InputUnit)
    DO 
      READ(InputUnit,'(A400)',IOSTAT=ios,IOMSG=iom) LocLine
      IF ( ios /= 0 ) THEN
        WRITE(*,*) ' Error Reading Units'
        WRITE(*,*) ' Error Message  ::  ',TRIM(iom)
      ELSE
        comments = (ADJUSTL(LocLine(1:1)) == '#'        .OR. &
        &           ADJUSTL(LocLine(1:7)) == 'COMMENT'  .OR. &
        &           LEN_TRIM(LocLine)     == 0 )
        IF ( comments ) CYCLE

        Pos = INDEX(LocLine,'GAS')
        IF ( Pos > 0 ) READ(LocLine(Pos+3:),*) UnitGas 
        Pos = INDEX(LocLine,'AQUA')
        IF ( Pos > 0 ) THEN
          READ(LocLine(Pos+4:),*) UnitAQUA 
          EXIT
        END IF
      END IF
    END DO

  END SUBROUTINE ReadUnits
  !
  !
  SUBROUTINE ExtractConstants(String,Ducts,Constants,Type,SpecialForm)
    CHARACTER(*) :: String
    REAL(dp), ALLOCATABLE :: Constants(:)
    !REAL(dp), POINTER :: Constants(:)
    CHARACTER(*) :: Type
    !
    INTEGER :: PosColon,PosName,PosComment,PosSemiColon
    INTEGER :: i,PosNum1,PosNum2,PosNum3,NumNum,PosElem
    CHARACTER(4)  :: NameNumNum
    CHARACTER(10) :: DummyString
    CHARACTER(LEN(String)) :: LocString
    CHARACTER(LEN(String)) :: LocString1
    CHARACTER(LEN(String)) :: NameConstant
    REAL(dp) :: Dummy
    INTEGER :: is

    ! this is for the new special formula input
    CHARACTER(LenLine)   :: SString
    TYPE(Special_T), OPTIONAL :: SpecialForm
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)
    INTEGER, ALLOCATABLE :: iSortedDucts(:)
    INTEGER :: nvs, cnt, idxDuct, lenDuct
    !
    LocString  = String
    String     = ''
    PosColon   = INDEX(LocString,':')
    Type       = LocString(:PosColon-1)
    LocString  = ADJUSTL(LocString(PosColon+1:))
    PosComment = INDEX(LocString,'#')

    IF ( PosComment > 0 ) LocString = LocString(:PosComment-1)
    
    LocString1 = LocString

    IF ( Type /= 'SPECIAL' ) THEN
      ALLOCATE(Constants(0))
      DO
        PosColon = INDEX(LocString1,':')
        IF ( PosColon > 0 ) THEN
          LocString1 = ADJUSTL(LocString1(PosColon+1:))
          READ( LocString1 , * , IOSTAT=is ) Dummy, DummyString
          PosName   = INDEX(LocString1,TRIM(DummyString))
          Constants = [Constants , Dummy]
          IF ( PosName > 0 ) LocString1 = LocString1(PosName:)
        ELSE
          EXIT
        END IF
      END DO

    ELSE
      ! special rate constatn formula
      IF ( PosColon==0 ) THEN
        WRITE(*,*) '  Missing seperator ":" for TypeConstant SPECIAL'
        STOP
      END IF

      PosSemiColon = Index(LocString,';')
      IF ( PosSemiColon==0 ) THEN
        WRITE(*,*) '  Missing seperator ";" for TypeConstant SPECIAL'
        STOP
      END IF

      ! save formula
      SString = ADJUSTL(LocString1(:PosSemiColon-1))

      ! read number of variables in formula
      READ( LocString1(PosSemiColon+1:) , * , IOSTAT=is ) nvs

      SpecialForm%nVariables = nvs      
      SpecialForm%Formula  = ADJUSTL(SString)
      IF (INDEX(SString,'TEMP')>0) SpecialForm%Temp = .TRUE.

      ALLOCATE(SpecialForm%cVariables(nvs))
      ALLOCATE(iSortedDucts(SIZE(Ducts)) )

      CALL StringSort(Ducts,iSortedDucts)

      cnt = 0
      DO i = SIZE(Ducts), 1 , -1
        idxDuct = INDEX(SString,TRIM(Ducts(i)))

        ! if a species if katalytic ( on both sides of the reaction equation )
        IF ( i > 1 ) THEN
          IF ( Ducts(i)==Ducts(i-1) ) CYCLE
        END IF

        IF ( idxDuct > 0 ) THEN
          cnt = cnt + 1
          SpecialForm%cVariables(cnt) = ADJUSTL(Ducts(i))

          DO
            idxDuct = INDEX(SString,TRIM(Ducts(i)))
            IF ( idxDuct == 0 ) EXIT
            lenDuct = LEN_TRIM(Ducts(i))
            SString(idxDuct:idxDuct+lenDuct) = REPEAT('_',lenDuct)
          END DO

        END IF
      END DO

      IF ( SpecialForm%Temp ) THEN
        SpecialForm%cVariables(cnt+1) = 'TEMP'
      ELSE
        WRITE(*,*) '      Warning: Missing temperature "TEMP" in SPECIAL Formula :: '&
        &           ,TRIM(SpecialForm%Formula)
      END IF
   
    END IF
  END SUBROUTINE ExtractConstants
  !
  !
  SUBROUTINE ExtractSpecies(String,Duct,InActDuct,NumInActDucts)
    ! IN:
    CHARACTER(*) :: String
    ! OUT:
    TYPE(Duct_T), POINTER :: Duct(:)
    TYPE(Duct_T), POINTER :: InActDuct(:)
    INTEGER :: NumInActDucts
    !
    INTEGER :: PosMinus,PosPlus,NumSpec,PosSpecies
    REAL(dp) :: PreFac !NumberSpecies
    CHARACTER(LenLine) :: Species
    CHARACTER(LEN(String)) :: LocString
    INTEGER :: sbL, sbR
    INTEGER :: ios
    LOGICAL :: dummy=.FALSE.

    !
    LocString = String
    NumSpec   = 1

    !WRITE(*,*) ' Current sys String :: ',TRIM(String)

    ! count species
    DO
     PosPlus  = INDEX(LocString,' + ')
     PosMinus = INDEX(LocString,' - ')

     IF ( PosPlus>0 .AND. (PosMinus==0 .OR. PosMinus>PosPlus) ) THEN
       LocString = LocString(PosPlus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosMinus>0 .AND. (PosPlus==0 .OR. PosPlus>PosMinus) ) THEN
       LocString = LocString(PosMinus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT
    END DO
   
    ALLOCATE( Duct(NumSpec) , InActDuct(NumSpec) )
    LocString = String
    NumSpec   = 0
    DO
      PosPlus  = INDEX(LocString,' + ')
      PosMinus = INDEX(LocString,' - ')
      IF ( PosMinus>0 .AND. PosMinus<PosPlus ) THEN
        PreFac = mONE
      ELSE
        PreFac = ONE
      END IF
      IF      ( PosPlus>0  .AND. (PosPlus<PosMinus.OR.PosMinus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosPlus-1))
        LocString = LocString(PosPlus+3:)
      ELSE IF ( PosMinus>0 .AND. (PosMinus<PosPlus.OR.PosPlus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosMinus-1))
        LocString = LocString(PosMinus+3:)
      ELSE
        Species = ADJUSTL(LocString)
      END IF

      PosSpecies = SCAN(Species,SetSpecies)
      NumSpec = NumSpec + 1


      ! check if there is a dummy species e.g. 0.000 with no species
      ! or species: (dummy)
      dummy = SCAN(TRIM(ADJUSTL(Species)) , SetSpecies) == 0 .OR. &
            & INDEX(TRIM(ADJUSTL(Species)) , '(dummy)') /= 0 

      !WRITE(*,*), ' species  ::: ', TRIM(Species)//'   is dummy = ', dummy

      IF (PosSpecies==1) THEN           
        sbL = INDEX(TRIM(Species),'[');  sbR = INDEX(TRIM(Species),']')

        ! check if species if passive
        IF ( sbL==1 .AND. LEN_TRIM(Species)==sbR-sbL+1 ) THEN
          ! works if there's just one InActEduct
          InActDuct(1)%Koeff   = PreFac
          InActDuct(1)%Species = Species
          NumInActDucts        = NumInActDucts + 1
        END IF
        Duct(NumSpec)%Koeff   = PreFac
        Duct(NumSpec)%Species = Species
      ELSE

        IF (.NOT.dummy) THEN
          READ( Species(1:PosSpecies-1) ,*, IOSTAT=ios) Duct(NumSpec)%Koeff
          IF ( ios == 0 ) THEN
            Duct(NumSpec)%Koeff   = PreFac * Duct(NumSpec)%Koeff
            Duct(NumSpec)%Species = Species(PosSpecies:)
          ELSE 
            WRITE(*,*) ' Error reading species and coefficients :: ', Species
            WRITE(*,*) ' IOSTAT = ', ios
          END IF
        ELSE
          Duct(NumSpec)%Koeff   = ZERO
          Duct(NumSpec)%Species = '(dummy)'
          !WRITE(*,*) ' ?? Dummy Species = ',TRIM(Species)
        END IF
      END IF

      ! Syntax check for missing seperators, e.g. "+CO2" insted of "+ CO2" or "CO2+" insted of "CO2 +"
      IF  ( INDEX(TRIM(ADJUSTL(Duct(NumSpec)%Species)),' ',.TRUE.) > 0 ) THEN
        WRITE(*,*); WRITE(*,*)
        WRITE(*,777) 'ERROR: Missing white space in reaction substring: '//TRIM(String)
        WRITE(*,777) '       Species = '//TRIM(Duct(NumSpec)%Species)
        WRITE(*,777) '       Check syntax in '//TRIM(SysFile)//'.sys!'
        WRITE(*,*); WRITE(*,*)
        STOP 
      END IF

      ! if no dummy species was found then add new species to hash table
      CALL InsertSpecies( Duct(NumSpec)%Species, Duct(NumSpec)%Type )

      ! if there are no further species exit the loop
      IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT

    END DO
    777 FORMAT(10X,A)
  END SUBROUTINE ExtractSpecies
  
  !
  SUBROUTINE InsertSpecies(Species,Type)
    CHARACTER(*) :: Species
    CHARACTER(*) :: Type

    INTEGER :: len

    len = LEN_TRIM(Species)

    IF (Species(1:1)=='p') THEN
      CALL InsertHash( ListPartic , TRIM(ADJUSTL(Species)) , ns_PARTI)
      Type = 'Partic'
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      CALL InsertHash( ListAqua   , TRIM(ADJUSTL(Species)) , ns_AQUA)
      Type = 'Aqua'
    ELSE IF (Species(1:1)=='s') THEN
      CALL InsertHash( ListSolid  , TRIM(ADJUSTL(Species)) , ns_SOLID)
      Type = 'Solid'
    ELSE IF ( Species(1:1)=='[' .AND. Species(len:len)==']' &
            &                   .AND. len<maxLENinActDuct   ) THEN
      CALL InsertHash( ListNonReac, TRIM(ADJUSTL(Species)) , ns_KAT)
      Type = 'Inert'
    ELSE IF ( MAXVAL(INDEX(Species(1:1) , ['(','0'])) > 0 ) THEN  ! dummy species
      Type = 'Dummy'
    ELSE
      CALL InsertHash( ListGas ,TRIM(ADJUSTL(Species)) ,     ns_GAS)
      Type = 'Gas'
    END IF
  END SUBROUTINE InsertSpecies
  !
  !
  FUNCTION PositionListSpecies(Species)
    TYPE(Species_T), POINTER :: PositionListSpecies
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    !
    PositionSpecies=0
    PositionListSpecies=>NULL()
    IF (Species(1:1)=='p') THEN
      PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListPartic2(PositionSpecies)
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListAqua2(PositionSpecies)
    ELSE IF (Species(1:1)=='s') THEN
      PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListSolid2(PositionSpecies)
    ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListNonReac2(PositionSpecies)
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListGas2(PositionSpecies)
    END IF
  END FUNCTION PositionListSpecies
  !
  !
  FUNCTION PositionSpeciesGas(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesGas
    ! 
    PositionSpeciesGas=0
    PositionSpeciesGas=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END FUNCTION
  !
  !
  FUNCTION PositionSpeciesCK(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesCK
    ! 
    PositionSpeciesCK=0
    PositionSpeciesCK=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END FUNCTION PositionSpeciesCK
  !
  !
  !
  FUNCTION PositionSpecies(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    ! 
    PositionSpecies=0
    IF (Species(1:1)=='p') THEN
      PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))     
    ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))       &
      &               + ns_GAS
    ELSE IF (Species(1:1)=='s') THEN
      PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))      
    ELSE IF (Species(1:1)=='['.AND.                                  &
    &        Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))  
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
    END IF
  END FUNCTION PositionSpecies
  !
  FUNCTION PositionAtom(Atom) RESULT(Pos)
    CHARACTER(*) :: Atom
    !
    INTEGER :: Pos
    ! 
    Pos = 0
    Pos = GetHash(ListAtoms,TRIM(ADJUSTL(Atom)))
  END FUNCTION PositionAtom
  !
  !
  FUNCTION PositionSpeciesAll(Species) RESULT(Pos)
    CHARACTER(*) :: Species
    !
    INTEGER :: Pos
    INTEGER :: len
   
    Pos = 0
    len = LEN_TRIM(Species)

    !WRITE(*,*) ' number ob species = ',species, len
  
    ! Combustion system
    IF ( Teq ) THEN
      Pos = -1
      Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
    ELSE

    ! tropospheric system
      IF ( Species(1:1) == 'p' ) THEN
        Pos = GetHash(ListPartic,TRIM(ADJUSTL(Species))) 
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA + ns_SOLID         
      
      ! AQUA 
      ELSE IF ( Species(1:1)=='a'.OR.SCAN(Species,'pm')>0 ) THEN
        Pos = GetHash(ListAqua,TRIM(ADJUSTL(Species)))
        IF (Pos>0) Pos = Pos + ns_GAS
 
      ! SOLID
      ELSE IF ( Species(1:1)=='s' ) THEN
        Pos = GetHash(ListSolid,TRIM(ADJUSTL(Species)))      
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA         

      ! NonReac
      ELSE IF ( Species(1:1) == '[' .AND. Species(len:len) == ']' .AND. &
      &        len < maxLENinActDuct) THEN
        Pos = GetHash(ListNonReac,TRIM(ADJUSTL(Species)))    
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI
        
      ! GAS
      ELSE
        Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
      END IF
    END IF
  END FUNCTION PositionSpeciesAll
  !
  !
  SUBROUTINE OpenFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) THEN
      OPEN(UNIT=InputUnit,FILE=TRIM(Filename)//'.'//TRIM(Type),STATUS='UNKNOWN')
    END IF
  END SUBROUTINE OpenFile
  !
  !
  SUBROUTINE CloseFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) CLOSE(UNIT=InputUnit)
  END SUBROUTINE CloseFile
  !
  !
  SUBROUTINE ReadThermoData(FileName)
    CHARACTER(*) :: Filename
    !
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenName) :: Name
    !
    INTEGER :: is,iLine,i
    REAL(8) :: Hf,Gf,Cp
    TYPE(Species_T), POINTER :: Species
    TYPE(Reaction_T), POINTER :: Current
    !
    CALL OpenFile(FileName,'dat')
    iLine=0
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
      IF (ABS(is)>0) THEN
        EXIT
      END IF
      IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
        READ(LocLine,*) Name
        IF (PositionSpecies(Name)>0) THEN
          iLine=iLine+1
        END IF
      END IF
    END DO
    REWIND(InputUnit)
    DO
      READ(InputUnit,'(a400)',IOSTAT=is) LocLine
      IF (ABS(is)>0) THEN
        EXIT
      END IF
      IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
        READ(LocLine,*) Name,Hf,Gf,Cp
        Species=>PositionListSpecies(Name)
        IF (ASSOCIATED(Species)) THEN
          Species%Hf=Hf
          Species%Gf=Gf
          Species%Cp=Cp
        END IF
      END IF
    END DO
    CLOSE(InputUnit)
    !
    Current=>ListRSolid%Start
    DO
      IF (ASSOCIATED(Current)) THEN
        IF (Current%TypeConstant=='DTEMP3') THEN
          Hf=0.0d0
          Gf=0.0d0
          Cp=0.0d0
          DO i=1,SIZE(Current%Educt)
            Species=>PositionListSpecies(Current%Educt(i)%Species)
            Hf=Hf-Current%Educt(i)%Koeff*Species%Hf
            Gf=Gf-Current%Educt(i)%Koeff*Species%Gf
            Cp=Cp-Current%Educt(i)%Koeff*Species%Cp
          END DO
          DO i=1,SIZE(Current%Product)
            Species=>PositionListSpecies(Current%Product(i)%Species)
            Hf=Hf+Current%Product(i)%Koeff*Species%Hf
            Gf=Gf+Current%Product(i)%Koeff*Species%Gf
            Cp=Cp+Current%Product(i)%Koeff*Species%Cp
          END DO
          Hf=Hf*1000.0d0
          Gf=Gf*1000.0d0
          IF (ALLOCATED(Current%Constants)) THEN
            DEALLOCATE(Current%Constants)
          END IF
          ALLOCATE(Current%Constants(3))
          Current%Constants(1)=EXP(-Gf/(RGas*TRef))
          Current%Constants(2)=Cp/RGas
          Current%Constants(3)=-Hf/RGas+Cp*TRef/RGas
        END IF
        Current=>Current%Next
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE ReadThermoData
  !
  !
  SUBROUTINE ReadSystem(FileName)
    CHARACTER(*) :: Filename
    !
    LOGICAL :: Out

    FileName = FileName(:INDEX(FileName,'.')-1)
    !
    CALL InitHashTable(ListAqua,100)
    CALL InitHashTable(ListGas,100)
    CALL InitHashTable(ListSolid,100)
    CALL InitHashTable(ListPartic,100)
    CALL InitHashTable(ListNonReac,100)
    CALL OpenFile(FileName,'spc')
    DO
      CALL ReadSpecies(Out);  IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'spc')
    CALL OpenFile(FileName,'sys')
    CALL ReadUnits
    DO
      CALL ReadReaction(Out); IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'sys')
    ALLOCATE(ListGas2(ns_GAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    ALLOCATE(ListAqua2(ns_AQUA))
    CALL HashTableToList(ListAqua,ListAqua2)
    CALL SortList(ListAqua2)
    CALL ListToHashTable(ListAqua2,ListAqua)
    ALLOCATE(ListSolid2(ns_SOLID))
    CALL HashTableToList(ListSolid,ListSolid2)
    CALL SortList(ListSolid2)
    CALL ListToHashTable(ListSolid2,ListSolid)
    ALLOCATE(ListPartic2(ns_PARTI))
    CALL HashTableToList(ListPartic,ListPartic2)
    CALL SortList(ListPartic2)
    CALL ListToHashTable(ListPartic2,ListPartic)
    ALLOCATE(ListNonReac2(ns_KAT))
    CALL HashTableToList(ListNonReac,ListNonReac2)
    CALL SortList(ListNonReac2)
    CALL ListToHashTable(ListNonReac2,ListNonReac)

  END SUBROUTINE ReadSystem
  !
  !
  SUBROUTINE SortList(List)
    TYPE(Species_T) :: List(:)
    !
    TYPE(Species_T) :: Temp
    INTEGER :: i,j
    !
    DO i=1,SIZE(List)
      DO j=1,SIZE(List)-i
        IF (List(j+1)%Species<List(j)%Species) THEN
          Temp=List(j+1)
          List(j+1)=List(j)
          List(j)=Temp
        END IF
      END DO
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='OHm') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='Hp') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
  END SUBROUTINE SortList
  !
  !
  SUBROUTINE ListToHashTable(List,Table)
    TYPE(Species_T) :: List(:)
    TYPE(hash_tbl_sll) :: Table
    !
    INTEGER :: i
    DO i=1,SIZE(List)
      CALL Table%put(TRIM(ADJUSTL(List(i)%Species)),i)
    END DO
  END SUBROUTINE ListToHashTable
  !
  !
  SUBROUTINE HashTableToList(Table,List)
    TYPE(hash_tbl_sll) :: Table
    TYPE(Species_T) :: List(:)
    !
    INTEGER :: i,j
    TYPE(sllist), POINTER :: child => NULL()
    !
    DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
      IF (ALLOCATED(table%vec(i)%key)) THEN
        DO j=1,SIZE(table%vec(i)%key)
          List(table%vec(i)%Val)%Species(j:j)=table%vec(i)%key(j)
        END DO
      END IF
      Child=>table%vec(i)%Child
      DO
        IF (ASSOCIATED(Child)) THEN
          DO j=1,SIZE(Child%key)
            List(Child%Val)%Species(j:j)=Child%key(j)
          END DO
          Child=>Child%Child
        ELSE
          EXIT
        END IF
      END DO
    END DO
  END SUBROUTINE HashTableToList
 
  

  SUBROUTINE AllListsToArray(ReacStruct,LGas,LHenry,LAqua,LDiss,LSolid,LPartic,LMicro)
    ! OUT:
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStruct(:)
    ! IN:
    TYPE(ListReaction_T) :: LGas, LHenry, LAqua, LDiss  &
    &                     , LSolid, LPartic, LMicro
    
    INTEGER :: i, j, iList, iEq
    INTEGER :: nList

    INTEGER :: icnt(47), icntFAC(10), iHen
   
    ! #Reactions
    neq  = nr_gas + 2*nr_henry + nr_aqua + 2*nr_diss &
    &    + nr_solid + nr_parti + nr_micphys
    nr   = neq
    
    ! #Spezies berechnen
    nspc = ns_GAS + ns_AQUA + ns_SOLID + ns_PARTI
    ns   = nspc
    
    nr_liquid = nr_aqua + 2*nr_diss

    !-----------------------------------------------------------------------
    ! --- Set logicals
    !-----------------------------------------------------------------------
		IF ( ns_GAS   > 0 ) hasGasSpc   = .TRUE.; nPhases = nPhases + 1
		IF ( ns_AQUA  > 0 ) hasAquaSpc  = .TRUE.; nPhases = nPhases + 1
		IF ( ns_SOLID > 0 ) hasSolidSpc = .TRUE.; nPhases = nPhases + 1
		IF ( ns_PARTI > 0 ) hasPartiSpc = .TRUE.; nPhases = nPhases + 1

		IF ( nr_gas    > 0 ) hasGasReac    = .TRUE.
		IF ( nr_aqua   > 0 ) hasAquaReac   = .TRUE.
		IF ( nr_henry  > 0 ) hasHenryReac  = .TRUE.
		IF ( nr_solid  > 0 ) hasSolidReac  = .TRUE.
		IF ( nr_parti  > 0 ) hasPartiReac  = .TRUE.
		IF ( nr_diss   > 0 ) hasDissReac   = .TRUE.
		IF ( nr_liquid > 0 ) hasLiquidReac = .TRUE.

    hasPhotoReac  = (nr_G_photo+nr_A_photo) > 0
		hasFactorReac = nr_FACTOR > 0
    
    nList = 0
    IF (ASSOCIATED(LGas%Start))    nList=nList+1
    IF (ASSOCIATED(LHenry%Start))  nList=nList+1
    IF (ASSOCIATED(LAqua%Start))   nList=nList+1
    IF (ASSOCIATED(LDiss%Start))   nList=nList+1
    IF (ASSOCIATED(LSolid%Start))  nList=nList+1
    IF (ASSOCIATED(LPartic%Start)) nList=nList+1
    IF (ASSOCIATED(LMicro%Start))  nList=nList+1
    ALLOCATE(CompleteReactionList(nList))
    nList = 0
    IF (ASSOCIATED(LGas%Start))   THEN
      nList=nList+1; CompleteReactionList(nList)=LGas
    END IF
    IF (ASSOCIATED(LHenry%Start)) THEN
      !ALLOCATE(isHenry(2*nr_henry,2)); isHenry = 0
      nList=nList+1; CompleteReactionList(nList)=LHenry
    END IF
    IF (ASSOCIATED(LAqua%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LAqua
    END IF
    IF (ASSOCIATED(LDiss%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LDiss
    END IF
    IF (ASSOCIATED(LSolid%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LSolid
    END IF
    IF (ASSOCIATED(LPartic%Start)) THEN
      nList=nList+1; CompleteReactionList(nList)=LPartic
    END IF
    IF (ASSOCIATED(LMicro%Start)) THEN
      nList=nList+1; CompleteReactionList(nList)=LMicro
    END IF

    ALLOCATE( ReacStruct(neq) )
  
 
    i=1
    iHen = 0
    icnt = 0
    icntFAC = 0

    DO iList=1,nList
      Current => CompleteReactionList(iList)%Start
      DO WHILE (ASSOCIATED(Current)) 
        ReacStruct(i)%Type   = Current%Type
        ReacStruct(i)%Line1  = Current%Line1
        ReacStruct(i)%Line2  = Current%Line2
        ReacStruct(i)%Line3  = Current%Line3
        ReacStruct(i)%Factor = Current%Factor
        ReacStruct(i)%TypeConstant = Current%TypeConstant

        CALL Check_ReacParameter( i                    &
        &                       , Current%TypeConstant &
        &                       , Current%Line1        &
        &                       , Current%Constants)
         
        !
        ReacStruct(i)%nActEd  = SIZE(Current%Educt)
        ReacStruct(i)%nActPro = SIZE(Current%Product)
        ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )

        DO j = 1 , ReacStruct(i)%nActEd
          !write(*,*) 'curr educ = ',i,j,TRIM(Current%Educt(j)%Species)
          ReacStruct(i)%Educt(j)%Species  = Current%Educt(j)%Species
          ReacStruct(i)%Educt(j)%Type     = Current%Educt(j)%Type
          ReacStruct(i)%Educt(j)%Koeff    = Current%Educt(j)%Koeff
          ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
        END DO
      
        DO j = 1 , ReacStruct(i)%nActPro
          !write(*,*) 'curr prod = ',i,j,TRIM(Current%Product(j)%Species)
          ReacStruct(i)%Product(j)%Species  = Current%Product(j)%Species
          ReacStruct(i)%Product(j)%Type     = Current%Product(j)%Type
          ReacStruct(i)%Product(j)%Koeff    = Current%Product(j)%Koeff
          ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
        END DO
        
        IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
          j = SIZE(Current%Special%cVariables)
          IF (Current%Special%Temp) THEN
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j-1))
          ELSE
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j))
          END IF

          ReacStruct(i)%Special%nVariables = j
          ReacStruct(i)%Special%Formula  = Current%Special%Formula
          ReacStruct(i)%Special%Temp     = Current%Special%Temp
          DO j = 1,ReacStruct(i)%Special%nVariables
            ReacStruct(i)%Special%cVariables(j) = Current%Special%cVariables(j)
            IF ( Current%Special%cVariables(j) /= 'TEMP' ) THEN 
              ReacStruct(i)%Special%iVariables(j) = PositionSpeciesAll(Current%Special%cVariables(j))
            END IF
          END DO
        ELSE
          ReacStruct(i)%nConst  = SIZE(Current%Constants)
          ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
          ReacStruct(i)%Constants = Current%Constants
        END IF
        !
        ! DAS HIER MUSS ANDERS WERDEN: INDIZES SO ABSPEICHER DASS SIE AUF INDEX DER TEMP3 REAKTIONEN ZEIGEN
        ! iR%iHENRY(icnt(29),1) müsste auf icnt(7) zeigen
        IF ( Current%Type == 'HENRY' ) THEN
          ReacStruct(i)%direction = 'GA'
          ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)
        END IF
        !
        !
        ReacStruct(i)%nInActEd  = Current%nInActEd
        ReacStruct(i)%nInActPro = Current%nInActPro
        ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActEd),      &
                & ReacStruct(i)%InActEductSpc(Current%nInActEd),   & 
                & ReacStruct(i)%InActProduct(Current%nInActPro),   &
                & ReacStruct(i)%InActProductSpc(Current%nInActPro) )

        DO j = 1 , Current%nInActEd
          ReacStruct(i)%InActEduct(j)    = Current%InActEduct(j)%Koeff
          ReacStruct(i)%InActEductSpc(j) = Current%InActEduct(j)%Species
          nFirst_orderKAT = nFirst_orderKAT + 1
        END DO
        DO j = 1 , Current%nInActPro
          ReacStruct(i)%InActProduct(j)    = Current%InActProduct(j)%Koeff    
          ReacStruct(i)%InActProductSpc(j) = Current%InActProduct(j)%Species
        END DO
        
        ReacStruct(i)%SumAqCoef = SUM(Current%Educt%Koeff) - ONE
       
        IF (ReacStruct(i)%nInActEd > 0 ) THEN
          IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
            ReacStruct(i)%SumAqCoef = ReacStruct(i)%SumAqCoef + ONE
          END IF
        END IF
        IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
          IF ( ReacStruct(i)%SumAqCoef > ZERO ) nr_HOaqua = nr_HOaqua + 1
        END IF
        !
        ! for equilibrium reactions save <-- direction
        SELECT CASE (Current%Type)
          CASE ('DISS','HENRY')
            i=i+1
            iEq = INDEX(Current%Line1,' = ')
            ReacStruct(i)%Type   = Current%Type
            ReacStruct(i)%Line1  = TRIM(Current%Line1(iEq+3:))//' = '//TRIM(Current%Line1(:iEq))
            ReacStruct(i)%Line2  = 'reverse reaction'
            ReacStruct(i)%bR     = .TRUE.
            ReacStruct(i)%Line3  = Current%Line3
            ReacStruct(i)%Factor = Current%Factor
            ReacStruct(i)%TypeConstant = Current%TypeConstant
           
            ReacStruct(i)%nActEd  = SIZE(Current%Product)
            ReacStruct(i)%nActPro = SIZE(Current%Educt)
            ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                    & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )

            DO j=1,ReacStruct(i)%nActEd
              ReacStruct(i)%Educt(j)%Species  = Current%Product(j)%Species
              ReacStruct(i)%Educt(j)%Type     = Current%Product(j)%Type
              ReacStruct(i)%Educt(j)%Koeff    = Current%Product(j)%Koeff
              ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
            END DO
            DO j=1,ReacStruct(i)%nActPro
              ReacStruct(i)%Product(j)%Species  = Current%Educt(j)%Species
              ReacStruct(i)%Product(j)%Type     = Current%Educt(j)%Type
              ReacStruct(i)%Product(j)%Koeff    = Current%Educt(j)%Koeff
              ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
            END DO
            
            IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
              ReacStruct(i)%Special%nVariables = ReacStruct(i-1)%Special%nVariables
              ReacStruct(i)%Special%Formula    = ReacStruct(i-1)%Special%Formula
              ReacStruct(i)%Special%Temp       = ReacStruct(i-1)%Special%Temp
              ReacStruct(i)%Special%cVariables = ReacStruct(i-1)%Special%cVariables
              ReacStruct(i)%Special%iVariables = ReacStruct(i-1)%Special%iVariables
            ELSE
              ReacStruct(i)%nConst  = SIZE(Current%Constants)
              ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
              ReacStruct(i)%Constants = Current%Constants
            END IF
            
            IF ( Current%Type == 'HENRY' ) THEN
              ReacStruct(i)%direction = 'AG'
              ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
            END IF
            !
            ReacStruct(i)%nInActEd  = Current%nInActPro
            ReacStruct(i)%nInActPro = Current%nInActEd
            ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActPro),    &
                    & ReacStruct(i)%InActEductSpc(Current%nInActPro), & 
                    & ReacStruct(i)%InActProduct(Current%nInActEd),   &
                    & ReacStruct(i)%InActProductSpc(Current%nInActEd) )

            DO j = 1 , Current%nInActPro
              ReacStruct(i)%InActEduct(j)    = Current%InActProduct(j)%Koeff
              ReacStruct(i)%InActEductSpc(j) = Current%InActProduct(j)%Species
              nFirst_orderKAT = nFirst_orderKAT + 1
            END DO
            DO j = 1 , Current%nInActEd
              ReacStruct(i)%InActProduct(j)    = Current%InActEduct(j)%Koeff    
              ReacStruct(i)%InActProductSpc(j) = Current%InActEduct(j)%Species
            END DO
            
            ReacStruct(i)%SumAqCoef = SUM(Current%Product%Koeff) - ONE
            
            IF (ReacStruct(i)%nInActEd > 0 ) THEN
              IF (TRIM(ADJUSTL(ReacStruct(i)%InActEductSpc(1)))=='[aH2O]') THEN 
                ReacStruct(i)%SumAqCoef = ReacStruct(i)%SumAqCoef + ONE
              END IF
            END IF
            IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
              IF ( ReacStruct(i)%SumAqCoef > ZERO ) nr_HOaqua = nr_HOaqua + 1
            END IF
            !
        END SELECT
        !
        Current=>Current%Next
        i=i+1
      END DO
    END DO
  

  END SUBROUTINE AllListsToArray

 

  SUBROUTINE Check_ReacParameter(iReac,TypeR,Line1,C)
    INTEGER,      INTENT(IN) :: iReac
    CHARACTER(*), INTENT(IN) :: TypeR
    CHARACTER(*), INTENT(IN) :: Line1
    REAL(dp),     INTENT(IN) :: C(:)

    SELECT CASE ( TRIM(TypeR) )
      CASE ('PHOTABC')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('PHOTMCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('PHOTAB')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('CONST')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TEMP1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TEMP2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TEMP3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TEMP4')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROE')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROEF')
        IF ( SIZE(C)<5 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROEQ')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROEQF')
        IF ( SIZE(C)<7 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROEXP')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('TROEMCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC3')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC4')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC1MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC2MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC3MCM')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC4MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC5MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC6MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC7MCM')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC8MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('SPEC9MCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
      CASE ('S4H2O')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('T1H2O')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('ASPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('ASPEC2')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('ASPEC3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DCONST')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DTEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DTEMP2')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DTEMP3')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DTEMP4')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
      CASE ('DTEMP5')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE ('MESKHIDZE')
        IF ( SIZE(C)<7 ) CALL ErrorMSG(iReac,Line1)
      CASE ('PHOTO')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
      CASE ('PHOTO2')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
      CASE ('PHOTO3')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
      CASE ('HOM1')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
      CASE DEFAULT
        WRITE(*,*) ''
        WRITE(*,*) ' Reaction Type unknown:  ',TRIM(TypeR),'  --> check input file'
        WRITE(*,*) ''
    END SELECT

    CONTAINS

      SUBROUTINE ErrorMSG(iR,Line)
        INTEGER      :: iR
        CHARACTER(*) :: Line
        WRITE(*,*);  WRITE(*,*)
        WRITE(*,*) '  ERROR -- > check parameter in reaction: ', iR ,'  ::  '//Line
        WRITE(*,*);  WRITE(*,*)
        STOP
      END SUBROUTINE ErrorMSG

  END SUBROUTINE Check_ReacParameter
  
  
  SUBROUTINE CheckConstants(RS)
    TYPE(ReactionStruct_T) :: RS(:)     ! reaction system
    CHARACTER(15) :: Const_T            ! constant type for reaction i
    !
    INTEGER :: i,j
    !

    DO i = 1,SIZE(RS)
      Const_T = RS(i)%TypeConstant
      DO j=1,nReacTypes  ! 37 comes from mo_reac def_par
        IF ( TRIM(reac_par(j)%name_type) == ADJUSTL(TRIM(Const_T)) ) EXIT
      END DO
      IF ( SIZE(RS(i)%Constants) /= reac_par(j)%n_par ) THEN
        WRITE(*,*) 'ERROR: Wrong number of constants:'
        WRITE(*,*) '----->  reaction:     ',i, '   ', TRIM(RS(i)%Line1)
        WRITE(*,*) '----->  desired #consts: ', reac_par(j)%n_par, j
        WRITE(*,*) '----->  actual  #consts: ', SIZE(RS(i)%Constants)
        WRITE(*,*) '       Check sys-file for syntax errors!'
        STOP
      END IF

    END DO
  END SUBROUTINE CheckConstants


  SUBROUTINE SearchReactions(Species)
    CHARACTER(*) :: Species
    CHARACTER(80) :: tmpSpc
    INTEGER :: iR, jD, uPath
    INTEGER :: cRcnt, pRcnt

    uPath = 13
    tmpSpc = TRIM(ADJUSTL(Species))
    cRcnt = 0
    pRcnt = 0

    OPEN ( UNIT=uPath , FILE='REACTION_PATHS/'//TRIM(tmpSpc)//'_path.txt' , STATUS='REPLACE' )
    WRITE(uPath,*) ' ********************************************************************************************'
    WRITE(uPath,*) '  '
    WRITE(uPath,*) '  Chemical Mechanism ::              ', TRIM(BSP)
    WRITE(uPath,*) '     System contains ::              ', neq , ' reactions'
    WRITE(uPath,*) '                                     ', nspc, ' species'
    WRITE(uPath,*) '  '
    WRITE(uPath,*) '  All reactions including species :: ', TRIM(tmpSpc)
    WRITE(uPath,*) '  '

    DO iR = 1 , neq
      ! Check educts
      DO jD = 1 , ReactionSystem(iR)%nActEd
        IF (TRIM(ReactionSystem(iR)%Educt(jD)%Species) == TRIM(tmpSpc))   cRcnt = cRcnt + 1
      END DO
      ! Check products
      DO jD = 1 , ReactionSystem(iR)%nActPro
        IF (TRIM(ReactionSystem(iR)%Product(jD)%Species) == TRIM(tmpSpc)) pRcnt = pRcnt + 1
      END DO
    END DO

    WRITE(uPath,*) '    + Number of Reactions where ',TRIM(tmpSpc),' is involved: ', cRcnt+pRcnt
    WRITE(uPath,*) '        - Number of consuming Reactions: ', cRcnt
    WRITE(uPath,*) '        - Number of producing Reactions: ', pRcnt

    DO iR = 1 , neq
      ! Check educts
      DO jD = 1 , ReactionSystem(iR)%nActEd
        IF (TRIM(ReactionSystem(iR)%Educt(jD)%Species) == TRIM(tmpSpc)) THEN
          CALL PrintReaction(iR,uPath)
        END IF
      END DO
      ! Check products
      DO jD = 1 , ReactionSystem(iR)%nActPro
        IF (TRIM(ReactionSystem(iR)%Product(jD)%Species) == TRIM(tmpSpc)) THEN
          CALL PrintReaction(iR,uPath)
        END IF
      END DO
    END DO

    CLOSE( UNIT=13 )

    WRITE(*,*) '  All reactions containing ',TRIM(tmpSpc), &
    &          ' saved in REACTION_PATHs/'//TRIM(tmpSpc)//'_path.txt'
  END SUBROUTINE SearchReactions

  SUBROUTINE PrintReaction(iR,Unit)
    INTEGER :: iR
    INTEGER :: Unit

    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) '  Reaction Number   :: ', iR
    WRITE(Unit,*) '  Reaction Class    :: ', TRIM(ReactionSystem(iR)%Type)
    WRITE(Unit,*) '  Constant Type     :: ', TRIM(ReactionSystem(iR)%TypeConstant)
    WRITE(Unit,*) '  Reaction          :: ', TRIM(ReactionSystem(iR)%Line1)
    WRITE(Unit,*) '  Order of Reaction :: ', INT(SUM(ReactionSystem(iR)%Educt%Koeff))
    WRITE(Unit,*) '  Factor            :: ', TRIM(ReactionSystem(iR)%Factor)
    WRITE(Unit,*) '  Constants         :: ', ReactionSystem(iR)%Constants
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintReaction

  SUBROUTINE Permute_Species(perm)
    INTEGER, ALLOCATABLE :: perm(:)
    INTEGER :: j
    INTEGER :: par(4) = [ 99999999,99999998,99999997,99999996 ]
    INTEGER, ALLOCATABLE :: p1(:), p2(:), p3(:), p4(:)
    INTEGER, ALLOCATABLE :: idx(:)

    ! Henry species in the middle
    ALLOCATE(idx(nspc))
    idx(1:ns_GAS)  = par(1)
    idx(iR%iHENRY(:,2)) = par(2)
    idx(ns_GAS+1:) = par(4)
    idx(iR%iHENRY(:,4)) = par(3)

    DO j=1,nspc
      IF     ( idx(j) /= par(1) ) THEN
        p1 = [ p1 , j ]
      ELSEIF ( idx(j) /= par(2) ) THEN
        p2 = [ p2 , j ]
      ELSEIF ( idx(j) /= par(3) ) THEN
        p3 = [ p3 , j ]
      ELSEIF ( idx(j) /= par(4) ) THEN
        p4 = [ p4 , j ]
      END IF
    END DO
    perm = [p1, p2, p3, p4]

  END SUBROUTINE Permute_Species


END MODULE Chemsys_Mod
