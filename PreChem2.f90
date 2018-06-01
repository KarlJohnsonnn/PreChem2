!=======================================
!=======================================
! vorlaeufiges Hauptprogramm fuer Tests
! box model chemical kinetics simulation
!=======================================
!=======================================
!
!
PROGRAM PreChem2
  !
  USE Kind_Mod
  USE Chemsys_Mod
  USE mo_control
  USE mo_reac
  USE fparser
  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  INTEGER  :: i

  ! NetCDF stuff
  REAL(dp) :: StartTimer
  

  !
  !================================================================
  !===                     MAIN Programm
  !================================================================

  WRITE(*,*); WRITE(*,*)
  WRITE(*,*) "                                                            _..._                                                                   "       
  WRITE(*,*) "                                                         .-'_..._''.                                                    .-''-.      "
  WRITE(*,*) "  _________   _...._                   __.....__       .' .'      '.\  .              __.....__      __  __   ___     .' .-.  )     "
  WRITE(*,*) "  \        |.'      '-.            .-''         '.    / .'           .'|          .-''         '.   |  |/  `.'   `.  / .'  / /      "
  WRITE(*,*) "   \        .'```'.    '. .-,.--. /     .-'''"//'"'//"-.  `. . '            <  |         /     .-''"//'"'//"'-.  `. |   .-.  .-.   '(_/   / /       "
  WRITE(*,*) "    \      |       \    \|  .-. |/     /________\   \| |             | |        /     /________\   \|  |  |  |  |  |     / /        "
  WRITE(*,*) "    |     |        |    || |  | ||                  || |             | | .'''-. |                  ||  |  |  |  |  |    / /          "
  WRITE(*,*) "    |      \      /    . | |  | |\    .-------------'. '             | |/.'''. \\    .-------------'|  |  |  |  |  |   . '           "
  WRITE(*,*) "    |     |\`'-.-'   .'  | |  '-  \    '-.____...---. \ '.          .|  /    | | \    '-.____...---.|  |  |  |  |  |  / /    _.-')   "
  WRITE(*,*) "    |     | '-....-'`    | |       `.             .'   '. `._____.-'/| |     | |  `.             .' |__|  |__|  |__|.' '  _.'.-''    "
  WRITE(*,*) "   .'     '.             | |         `''-...... -'       `-.______ / | |     | |    `''-...... -'                  /  /.-'_.'        "
  WRITE(*,*) " '-----------'           |_|                                      `  | '.    | '.                                 /    _.'           " 
  WRITE(*,*) "                                                                     '---'   '---'                               ( _.-'              "

  

  !----------------------------------------------------------------
  ! --- Read run control parameters (which runfile)
  CALL getarg( 1 , FileName0 )             
  IF ( FileName0 == '' ) THEN
    WRITE(*,777,ADVANCE='NO') 'Input Sys-File: '; READ(*,*)   FileName0
  END IF
  SysFile  = TRIM(ADJUSTL(FileName0))
  ChemFile = ADJUSTL(SysFile(:INDEX(SysFile,'.sys')-1)//'.chem')

  !================================================================
  !===                     Initialization
  !================================================================
  !
 

  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  CALL CPU_TIME(StartTimer)

  WRITE(*,777,ADVANCE='NO') 'Reading sys-file .............'

  CALL ReadSystem( SysFile )
  WRITE(*,*) 'done'
  WRITE(*,*)

  !-----------------------------------------------------------------------
  ! --- print reactions and build A, B and (B-A) structure
  !-----------------------------------------------------------------------
  CALL Print_ChemFile( ReactionSystem , ChemFile , ChemUnit , .FALSE. )

  !-----------------------------------------------------------------------
  ! --- initialize fpraser for reactions with special rate formula
  !-----------------------------------------------------------------------
!  IF ( nr_special > 0 ) THEN
!    
!    ! Initialize function parser for n special functions
!    CALL initf( nr_special ) 
!    
!    ! Parse and bytecompile ith function string 
!    DO i = 1,nr_special
!      CALL parsef ( i, ReactionSystem(iR%iSPECIAL(i))%Special%Formula    &
!      &              , ReactionSystem(iR%iSPECIAL(i))%Special%cVariables )
!    END DO
!  END IF
  
  !-----------------------------------------------------------------------
  ! --- Timers
  CALL CPU_TIME(Time_Read)
  Time_Read   = Time_Read - StartTimer ! Stop Input timer
  WRITE(*,'(10X,A,Es8.2, A)') 'CPU Time = ' , Time_Read , ' [sec]'
  WRITE(*,*)
  WRITE(*,777) 'Chem-File location  ::  '//TRIM(ChemFile)
  !-----------------------------------------------------------------------


  WRITE(*,*); WRITE(*,*)

  777  FORMAT(10X,A)

END PROGRAM PreChem2
