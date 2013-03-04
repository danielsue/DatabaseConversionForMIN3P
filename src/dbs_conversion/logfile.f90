! module of logfile
module logfile

    use global, only : iUnitLog, bOpenLog, ErrorHandling
    
    implicit none
    
    character(260) ::   filePathLog =   ""          !the maximum length that can be used by windows

    integer :: nWarnings = 0
    integer :: nErrors   = 0

    interface WriteLog

        module procedure writeLogStr

        module procedure writeLogInt

        module procedure writeLogReal

        module procedure writeLogIntArray

        module procedure writeLogRealArray
        
        module procedure writeLogRealArray2

    end interface
    
contains
    
    !Open log file
    subroutine OpenLogFile
    
        implicit none
                
        integer :: stat        
       
        open(iUnitLog, file = filePathLog,  iostat = stat, status = "replace")
        
        if (stat /= 0) then            
            call ErrorHandling
        else
            bOpenLog = .true.
            call WriteLog("Open log file success:" // trim(filePathLog))
        end if  
        
    end subroutine OpenLogFile
     
  
    ! Write string to log file
    subroutine writeLogStr(strLog)
    
        implicit none
        
        character(*), intent(in) :: strLog
        
        write(*,"(a)") trim(strLog)
        
        if(bOpenLog) then
            write(iUnitLog, "(a)")  trim(strLog)
        end if
        
    end subroutine writeLogStr

    ! Write integer to log file 
    subroutine writeLogInt(intValue)

        implicit none

        integer, intent(in) :: intValue

        write (*,*) intValue

        if(bOpenLog) then
            write(iUnitLog, *)  intValue
        end if

    end subroutine writeLogInt

    ! Write real  to log file
    subroutine writeLogReal(realValue)

        implicit none

        real, intent(in) :: realValue

        write (*,*) realValue

        if(bOpenLog) then
            write(iUnitLog, *)  realValue
        end if

    end subroutine writeLogReal

    ! Write integer array to log file 
    subroutine writeLogIntArray(intArray)

        implicit none

        integer, allocatable, intent(in) :: intArray(:)

        write (*,*) intArray

        if(bOpenLog) then
            write(iUnitLog, *)  intArray
        end if

    end subroutine writeLogIntArray

    ! Write real array to log file
    subroutine writeLogRealArray(realArray)

        implicit none

        real, allocatable, intent(in) :: realArray(:)

        write (*,*) realArray

        if(bOpenLog) then
            write(iUnitLog, *)  realArray
        end if

    end subroutine writeLogRealArray
    
    ! Write real array to log file
    subroutine writeLogRealArray2(n, realArray)

        implicit none
        integer, intent(in) :: n
        real, intent(in) :: realArray(n)

        write (*,*) realArray

        if(bOpenLog) then
            write(iUnitLog, *)  realArray
        end if

    end subroutine writeLogRealArray2
    
    ! Write log summary
    subroutine WriteLogSummary
    
        implicit none
        
        if(bOpenLog) then
            write(iUnitLog, "(a9,6x,i6)")  "Warnings:", nWarnings
            write(iUnitLog, "(a7,8x,i6)")  "Errors:", nErrors
        end if
        write(*, "(a9,6x,i6)")  "Warnings:", nWarnings
        write(*, "(a7,8x,i6)")  "Errors:", nErrors
        write(*,*)
        write(*,"(a)") "Please see " // trim(filePathLog) // " for detail conversion informations."
    
    end subroutine WriteLogSummary
    
end module logfile