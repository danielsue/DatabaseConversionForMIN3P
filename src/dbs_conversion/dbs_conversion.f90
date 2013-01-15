!  dbs_conversion.f90 
!
!  FUNCTIONS:
!  dbs_conversion - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: dbs_conversion
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program dbs_conversion    

    use global
    use logfile,            only :  filePathLog, OpenLogFile, WriteLog, WriteLogSummary
    use dbs_toughreact,     only :  filePathDbsTR, OpenDbsTR, ReadDbsTR, CloseDbsTR
    use convt_tough_min3p,  only :  convert2Min3PDbs
    use dbs_min3p,          only :  OpenAllDbsMin3P, WriteAllDbsMin3P
    use alias,              only :  OpenDbsAlias, ReadDbsAlias
    use name_truncation,    only :  OpenNameTruncation, WriteNameTruncation
    
    implicit none

    ! Variables    
    character(1024)  ::  strBuffer
    logical         ::  file_exists =   .false.
    integer         ::  i
    
    ! Body of dbs_conversion
    
    write(*,2) 
    
2   format (/,'    ------------------------------------------------'   &
     &      //,'                  DATABASE CONVERSION'                  &
     &       /,'                          FOR '                         &
     &       /,'                  TOUGHREACT TO MIN3P'                  &
     &      //,'                       DANYANG SU'                      &
     &      //,'                 EMAIL: DSU@EOS.UBC.CA'                 &
     &       /,'   -------------------------------------------------'//)  
    
    write(*,*) "Type in the file path of toughreact database: "
    
    ! Read input file: toughreact database
100 read(*,"(a)") strBuffer
    
    inquire(file = trim(adjustl(strBuffer)), exist = file_exists)
    
    if (file_exists == .false.) then
        write(*,"(a)") "File does not exist, please retype in the file path of toughreact database: " // trim(adjustl(strBuffer))
        goto 100
    end if
    
    filePathDbsTR = trim(strBuffer)
    
    ! Generate log file path
    i = index(filePathDbsTR, "." , .true.)
    if (i > 0) then
        filePathLog = filePathDbsTR(1:i) // "log"
    else
        filePathLog = trim(filePathDbsTR) // ".log"
    end if    
    
    ! Open log file
    call OpenLogFile

    ! Open alias database if the file exists and open
    call OpenDbsAlias

    call ReadDbsAlias
    
    ! Read database
    call OpenDbsTR

    call ReadDbsTR

    call CloseDbsTR
    
    ! Convert database
    call convert2Min3PDbs
    
    ! Write to min3p database
    call OpenAllDbsMin3P
    
    call WriteAllDbsMin3P
    
    ! Write name truncation
    call OpenNameTruncation
    
    call WriteNameTruncation
    
    ! Write log summary
    call WriteLogSummary
    
    ! Close all files
    call CloseAllFiles
    
    stop

    end program dbs_conversion

