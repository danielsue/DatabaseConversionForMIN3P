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
    use dbs_toughreact,     only :  OpenDbsTR, ReadDbsTR, CloseDbsTR
    use convt_tough_min3p,  only :  convert2Min3PDbs
    use dbs_min3p,          only :  OpenAllDbsMin3P, WriteAllDbsMin3P, temperature_min3p
    use alias,              only :  OpenDbsAlias, ReadDbsAlias
    use name_truncation,    only :  OpenNameTruncation, WriteNameTruncation
    use inputfile,          only :  OpenInp, ReadInp, CloseInp, filePathInp, & 
                                     sourceDatabaseType, targetDatabaseType
    
    implicit none

    ! Variables    
    character(1024)  ::  strBuffer
    logical         ::  file_exists =   .false.
    integer         ::  i
    
    ! Body of dbs_conversion
    
    write(*,2) 
    
2   format (/,'    ------------------------------------------------'   &
     &      //,'                  DATABASE CONVERSION'                  &
     &       /,'                        FOR MIN3P'                      &
     &      //,'                       DANYANG SU'                      &
     &      //,'                 EMAIL: DSU@EOS.UBC.CA'                 &
     &       /,'   -------------------------------------------------'//)  
    
    write(*,*) "Type in the input file path (*.inp): "
    
    ! Read input file: toughreact database
100 read(*,"(a)") strBuffer
    
    inquire(file = trim(adjustl(strBuffer)), exist = file_exists)
    
    if (file_exists == .false.) then
        write(*,"(a)") "File does not exist, please retype in the input file path (*.inp): " // trim(adjustl(strBuffer))
        goto 100
    end if
    
    filePathInp = trim(strBuffer)

    ! Generate log file path
    i = index(filePathInp, "." , .true.)
    if (i > 0) then
        filePathLog = filePathInp(1:i) // "log"
    else
        filePathLog = trim(filePathInp) // ".log"
    end if
    
    ! Open log file
    call OpenLogFile

    ! Open input file
    call OpenInp

    ! Read input file
    call ReadInp

    ! Close input file
    call CloseInp

    ! Open alias database if the file exists and open
    call OpenDbsAlias

    call ReadDbsAlias

    ! Read source database
    select case (trim(sourceDatabaseType))
        
        case ("toughreact")

            ! Read database
            call OpenDbsTR

            call ReadDbsTR

            call CloseDbsTR

        case default

            call WriteLog ("Unknown source database: " // trim(sourceDatabaseType))

            call ErrorHandling

    end select

    ! Conver and write to target database

    select case (trim(targetDatabaseType))

        case ("min3p")
   
            ! Convert database
            call convert2Min3PDbs
    
            ! Write to min3p database
            call OpenAllDbsMin3P
    
            call WriteAllDbsMin3P

        case default

            call WriteLog ("Unknown target database: " // trim(targetDatabaseType))

            call ErrorHandling

    end select
    
    ! Write name truncation
    call OpenNameTruncation
    
    call WriteNameTruncation
    
    ! Write log summary
    call WriteLogSummary
    
    ! Close all files
    call CloseAllFiles
    
    stop

    end program dbs_conversion

