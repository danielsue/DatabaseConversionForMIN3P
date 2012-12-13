!Module of alias management
    
module alias

    use global, only : iUnitDbsAlias, bOpenDbsAlias, ErrorHandling
    use logfile, only : WriteLog, nWarnings, nErrors
    use file_utility, only : LowerCase

    implicit none
    
    character(9) ::   filePathDbsAlias = "alias.dbs"          !the maximum length that can be used by windows
    character(1),parameter :: strComment = "!"
    integer, parameter :: maxAlias = 1000       !Assume the maximum of alias is 1000, enlarge it if necessary
    character(60)  :: nameOrigional(maxAlias)       
    character(60)  :: nameAlias(maxAlias)
    integer        :: nNameAlias = 0


    
contains

    !Open alias database
    subroutine OpenDbsAlias
    
        implicit none
                
        integer :: stat        
        logical ::  file_exists = .false. 

        inquire(file = trim(filePathDbsAlias), exist = file_exists)

        if (.not.file_exists) then
            call WriteLog("No alias database found: " // trim(filePathDbsAlias))
            return
        end if

        open(iUnitDbsAlias, file = filePathDbsAlias,  iostat = stat)
        
        if (stat /= 0) then
            call WriteLog("Open alias database failed:" // trim(filePathDbsAlias)) 
            call ErrorHandling
        else
            bOpenDbsAlias = .true.
            call WriteLog("Open alias database success:" // trim(filePathDbsAlias))
        end if  
        
    end subroutine OpenDbsAlias
    
    !Read alias database
    subroutine ReadDbsAlias
    
        implicit none

        logical :: bEndOfFile = .false.
        integer :: iReadStat = 0
        integer :: i = 0, j = 0
        character(60) :: strName, strAlias
        character(260) :: strBuffer 
      
        if(.not.bOpenDbsAlias) then 
            return
        end if

        call WriteLog("Start reading alias database: " // trim(filePathDbsAlias))

        do while(.not. bEndOfFile)

            read(iUnitDbsAlias, "(a)", iostat = iReadStat) strBuffer
       
            if (iReadStat > 0) then     !Error in reading
            
                call WriteLog("Error in reading alias databse file:" // trim(filePathDbsAlias))
            
                call ErrorHandling
            
            else if (iReadStat < 0) then    !End of file
            
                bEndOfFile = .true.
                
                exit            
            
            end if

            strBuffer = adjustl(strBuffer)
            
            ! Conver all case to lower case
            call lowerCase(strBuffer)

            if (len_trim(strBuffer) == 0 ) then
                continue
            end if

            if (strBuffer(1:1) /= strComment) then
                i = index(strBuffer," ")
                strName = strBuffer(1:i-1)
                strAlias = trim(adjustl(strBuffer(i:)))
                call AddAlias(strName, strAlias)
            end if

        end do

        call WriteLog("End reading alias database: " // trim(filePathDbsAlias))

    
    end subroutine ReadDbsAlias

   
    !Add alias
    subroutine addAlias(strName, strAlias)
    
        implicit none
        
        character(*), intent(in) :: strName, strAlias
        integer :: i
        
        do i = 1, nNameAlias
            if(trim(nameOrigional(i)) == trim(strName)) then
                nWarnings = nWarnings + 1
                call WriteLog("Warning(ignored): Name duplicated, the same name already exist: " // trim(strName))  
                return
            end if
        end do
        
        do i = 1, nNameAlias
            if(trim(nameOrigional(i)) == trim(strAlias)) then
                nWarnings = nWarnings + 1
                call WriteLog("Warning(ignored): Alias conflict with the origional name: " // trim(strAlias))    
                return
            end if
        end do
        
        do i = 1, nNameAlias
            if(trim(nameAlias(i)) == trim(strAlias)) then
                nWarnings = nWarnings + 1
                call WriteLog("Warning(ignored): Alias duplicated, the same name already exist: " // trim(strAlias))  
                return
            end if
        end do  
           
        if (nNameAlias == maxAlias) then
            nErrors = nErrors + 1
            call WriteLog("The maximum database of alias is 1,000 lines." // &
                            "Please modify the source code parameter maxAlias to read more data.")
            call ErrorHandling
        end if 
        
        nNameAlias = nNameAlias + 1
        nameOrigional(nNameAlias) = trim(strName)
        nameAlias(nNameAlias) = trim(strAlias)

        call WriteLog("Name and alias added: " // trim(adjustl(nameOrigional(nNameAlias))) // "  " // trim(adjustl(nameAlias(nNameAlias))))
    
    end subroutine AddAlias
    
    !Get alias from name
    function GetAliasFromName(strName) result(strAlias)
    
        implicit none
        
        character(*), intent(in) :: strName
        character(60) :: strAlias
        integer :: i
        
        strAlias = ""
        
        do i = 1, nNameAlias
            if(trim(nameOrigional(i)) == trim(strName)) then
                strAlias = trim(nameAlias(i))
                return
            end if
        end do
    
    end function GetAliasFromName
    
    
    !Get name from name alias
    function GetNameFromAlias(strAlias) result(strName)
    
        implicit none
        
        character(*), intent(in) :: strAlias
        character(60) :: strName
        integer :: i
        
        strName = ""
        
        do i = 1, nNameAlias
            if(trim(nameAlias(i)) == trim(strAlias)) then
                strName = trim(nameOrigional(i))
                return
            end if
        end do
    
    end function GetNameFromAlias
    
end module alias