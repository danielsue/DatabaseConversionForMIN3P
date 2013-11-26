!****************************************************************************
!    Geochemitry database conversion program
!    Author: Danyang Su
!    Email: dsu@eos.ubc.ca; danyang.su@gmail.com
!****************************************************************************
! Module of NameTruncation, similar to alias.
! This file is used to collect the name and truncation information during database conversion.
! It help you to decide the which name should be given a better alias.
! You can modify the truncation name and then copy it to alias.dbs. 
! Maximum name length is 60, maximum data amount is 1000. 
! File format: name  truncation
! e.g., abcdefghijkl abcde


module name_truncation

    use global, only : iUnitNameTruncation, bOpenNameTruncation, ErrorHandling
    use logfile, only : WriteLog, nWarnings, nErrors
    use file_utility, only : SetLowerCase

    implicit none
    
    character(19) ::   filePathNameTruncation = "name_truncation.txt"          !the maximum length that can be used by windows
    character(1),parameter :: strComment = "!"
    integer, parameter :: maxNameTruncation = 1000       !Assume the maximum of NameTruncation is 1000, enlarge it if necessary
    character(60)  :: longNames(maxNameTruncation)       
    character(60)  :: truncateNames(maxNameTruncation)
    integer        :: nNameTruncation = 0
    
contains

    !Open NameTruncation database
    subroutine OpenNameTruncation
    
        implicit none
                
        integer :: stat        
        logical ::  file_exists = .false. 

        open(iUnitNameTruncation, file = filePathNameTruncation,  iostat = stat, status = "replace")
        
        if (stat /= 0) then
            call WriteLog("Open/Create name truncation file failed:" // trim(filePathNameTruncation)) 
            call ErrorHandling
        else
            bOpenNameTruncation = .true.
            call WriteLog("Open/Create name truncation file success:" // trim(filePathNameTruncation))
        end if  
        
    end subroutine OpenNameTruncation
    
    !Write name truncation
    subroutine WriteNameTruncation
    
        implicit none        
        integer :: i
        
        if (.not. bOpenNameTruncation) then
            return
        end if
        
        call WriteLog("Write name truncation: " // trim(filePathNameTruncation))
        
        ! write header
        write(iUnitNameTruncation, 100)

100     format( "! This file is used to collect the name and truncation information during database conversion." &
                / "! It help you to decide the which name should be given a better alias."                      &
                / "! You can modify the truncation name and then copy it to alias.dbs."                         &
                / "! Maximum name length is 60, maximum data amount is 1000. "                                  &
                / "! File format: name ; truncation"                                                             &
                / "! e.g., abcdefghijkl abcde" )
        
        do i = 1, nNameTruncation
            write(iUnitNameTruncation, "(a60,a1,1x,a)")  longNames(i), ";", trim(truncateNames(i))
        end do
        
        call WriteLog("End writing name truncation.")
        call WriteLog("Number of truncated names:")
        call WriteLog(nNameTruncation)
    
    end subroutine WriteNameTruncation
   
    !Add NameTruncation
    subroutine AddNameTruncation(strName, strNameTruncation)
    
        implicit none
        
        character(*), intent(in) :: strName, strNameTruncation
        integer :: i
        
        do i = 1, nNameTruncation
            if(trim(longNames(i)) == trim(strName)) then
                return
            end if
        end do
        
        if (checkDuplicateName(strName)) then
            return
        end if
        
        if (checkDuplicateTruncation(strNameTruncation)) then
            nErrors = nErrors + 1
            call WriteLog("Error: Truncate name duplicated, the same truncate name already exist: " // trim(strNameTruncation))            
        end if
        
        if (checkConflictTruncation(strNameTruncation)) then
            nErrors = nErrors + 1
            call WriteLog("Error: Truncate name conflict with the origional name: " // trim(strNameTruncation))            
        end if
           
        if (nNameTruncation == maxNameTruncation) then
            nErrors = nErrors + 1
            call WriteLog("Error: The maximum database of NameTruncation is 1,000 lines." // &
                            "Please modify the source code parameter maxNameTruncation to read more data.")
            call ErrorHandling
        end if 
        
        nNameTruncation = nNameTruncation + 1
        longNames(nNameTruncation) = trim(strName)
        truncateNames(nNameTruncation) = trim(strNameTruncation)

        call WriteLog("New name and truncation added: " // trim(longNames(nNameTruncation)) // "  " // trim(truncateNames(nNameTruncation)))
    
    end subroutine AddNameTruncation
    
    !Get Truncation from name
    function GetTruncationFromName(strName) result(strNameTruncation)
    
        implicit none
        
        character(*), intent(in) :: strName
        character(60) :: strNameTruncation
        integer :: i
        
        strNameTruncation = ""
        
        do i = 1, nNameTruncation
            if(trim(longNames(i)) == trim(strName)) then
                strNameTruncation = trim(truncateNames(i))
                return
            end if
        end do
    
    end function GetTruncationFromName
    
    
    !Get name from name NameTruncation
    function GetNameFromTruncation(strNameTruncation) result(strName)
    
        implicit none
        
        character(*), intent(in) :: strNameTruncation
        character(60) :: strName
        integer :: i
        
        strName = ""
        
        do i = 1, nNameTruncation
            if(trim(truncateNames(i)) == trim(strNameTruncation)) then
                strName = trim(longNames(i))
                return
            end if
        end do
    
    end function GetNameFromTruncation
    
    !Check if there is duplicate truncate name
    function checkDuplicateTruncation(strNameTruncation) result(bFlag)
    
        implicit none
        
        character(*), intent(in) :: strNameTruncation
        logical :: bFlag       
        integer :: i
        
        bFlag = .false.
        do i = 1, nNameTruncation
            if(trim(truncateNames(i)) == trim(strNameTruncation)) then
                bFlag = .true.
                return
            end if
        end do
    
    end function checkDuplicateTruncation
    
    !Check if the truncate name is the same with existed origional names (longNames)
    function checkConflictTruncation(strNameTruncation) result(bFlag)
        
        implicit none
        
        character(*), intent(in) :: strNameTruncation
        logical :: bFlag       
        integer :: i
        
        bFlag = .false.
        do i = 1, nNameTruncation
            if(trim(longNames(i)) == trim(strNameTruncation)) then
                bFlag = .true.
                return
            end if
        end do
    
    end function
    
    !Check if there is duplicate name
    function checkDuplicateName(strName) result(bFlag)
    
        implicit none
        
        character(*), intent(in) :: strName
        logical :: bFlag   
        integer :: i
        
        bFlag = .false.
        do i = 1, nNameTruncation
            if(trim(longNames(i)) == trim(strName)) then
                bFlag = .true.
                return
            end if
        end do
    
    end function checkDuplicateName

end module name_truncation