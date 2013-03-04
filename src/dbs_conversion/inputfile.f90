!module of input file

module inputfile

use global, only  : iUnitInp, bOpenInp, ErrorHandling
use logfile, only : WriteLog
use file_utility, only : LowerCase

implicit none

character(1),   parameter   ::  strComment      =   "!"
character(1024)              ::  strBuffer        =   ""
character(1024)              ::  filePathInp      =   ""
integer                      ::  iReadStat       =   1
logical                       ::  bEndOfFile      =   .false.

!Block 1: Database settings
character(128) :: sourceDatabaseType = ""       !Source database type
character(256) :: sourceDatabasePath = ""       !Source database path
character(128) :: targetDatabaseType = ""       !Source database type
character(256) :: targetDatabasePath = ""       !Source database path

!Block 2: Conversion settings
real            :: targetTemperature = -9999     !Target database temperature  
character(60)  :: targetMasterVariable = ""      !Target database master varialbe

contains

    ! Open input file
    subroutine OpenInp
    
        implicit none
        
        integer :: stat 
        
        open(iUnitInp, file = filePathInp,  iostat = stat)
        
        if (stat /= 0) then
            call WriteLog("Error in opening input file file: " // trim(filePathInp))
            call ErrorHandling
        else
            bOpenInp = .true.
            call WriteLog("Open input file success: " // trim(filePathInp))
        end if 
    
    end subroutine OpenInp    
    
    
    ! Read input file
     subroutine ReadInp
     
        implicit none  
        
        call WriteLog ("Begin reading input file")
        
        do while(.not. bEndOfFile)

            call readNextLine

            if (bEndofFile) then
                exit
            end if

            select case (trim(strBuffer))

                case ("source database type")
                    call readNextLine
                    sourceDatabaseType = strBuffer
                    call WriteLog("source database type: " // trim(sourceDatabaseType))
                case ("source database path")
                    call readNextLine
                    sourceDatabasePath = strBuffer
                    call WriteLog("source database path: " // trim(sourceDatabasePath))
                case ("target database type")
                    call readNextLine
                    targetDatabaseType = strBuffer
                    call WriteLog("target database type: " // trim(targetDatabaseType))
                case ("target database path")
                    call readNextLine
                    targetDatabasePath = strBuffer
                    call WriteLog("target database path: " // trim(targetDatabasePath))
                case ("target database temperature")
                    call readNextLine
                    read(strBuffer, *) targetTemperature
                    call WriteLog("target database temperature: " // trim(strBuffer))
                case ("target database master variable")
                    call readNextLine
                    targetMasterVariable = strBuffer
                    call WriteLog("target database master variable: " // trim(targetMasterVariable))
                case default
                    call WriteLog("Unknow instruction: " // trim(strBuffer))
                    call ErrorHandling

            end select

        end do

        call WriteLog ("End reading input file")
     
     end subroutine ReadInp
     
     ! Close input file
     subroutine CloseInp
     
        implicit none
        
        if (bOpenInp) then
            
            close(iUnitInp)
            
        end if
     
     end subroutine CloseInp
     
    
    ! Read next line
    subroutine readNextLine
    
        implicit none
        
        integer :: i

        strBuffer = ""

        do while(.not. bEndOfFile)

            read(iUnitInp, "(a)", iostat = iReadStat) strBuffer
     
            if (iReadStat > 0) then     !Error in reading
            
                call WriteLog("Error in reading input file.")
            
                call ErrorHandling
            
            else if (iReadStat < 0) then    !End of file
            
                bEndOfFile = .true.
                
                exit            
            
            end if

            strBuffer = adjustl(strBuffer)

            if(len_trim(strBuffer) /= 0 ) then
            
                ! Conver all case to lower case
                call LowerCase(strBuffer)

                if (strBuffer(1:1) /= strComment) then
                    exit
                end if

            end if

        end do
    
    endsubroutine readNextLine

end module inputfile