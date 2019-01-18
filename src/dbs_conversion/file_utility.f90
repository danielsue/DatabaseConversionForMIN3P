!****************************************************************************
!    Geochemitry database conversion program
!    Author: Danyang Su
!    Email: dsu@eos.ubc.ca; danyang.su@gmail.com
!****************************************************************************
!module of file and string utility
module file_utility

    use global, only : ErrorHandling
    
    use logfile, only : WriteLog

    implicit none

    contains

    !Check if the number is zero
    function IsZero(dvalue) result(bFlag)
        implicit none
        real, intent(in) :: dvalue
        logical ::bFlag
        bFlag = .false.
        if(abs(dvalue) <= 1.0E-100) then
            bFlag = .true.
        else
            bFlag = .false.
        end if
    end function IsZero

    ! Convert upper case to lower case
    subroutine SetLowerCase(string)

        implicit none

        character(*), intent(INOUT) :: string

        integer :: i,j

        do i = 1,len_trim(string)

            j = ichar(string(i:i))

            if (j >= 65 .AND. j <= 90) then

                j = j + 32
                string(i:i) = char(j)

            end if

        end do

    end subroutine SetLowerCase
    
    ! Get lower case of a string
    function GetLowerCase(n, string) result (lowString)
        
        implicit none
        integer, intent(in) :: n
        character(n), intent(in) :: string
        character(n) :: lowString
        
        lowString = string
        
        call SetLowerCase(lowString)
    
    end function GetLowerCase
    

    ! replace character in string
    subroutine replaceCharacter(str, stra, strb)
    
        implicit none
        
        character(*), intent(inout) :: str
        character(*), intent(in) :: stra, strb
        character(1024) :: tempstr 
        
        integer :: i, j, k, n1, n2, n3
        
        tempstr = ""
        
        n1 = len(str)
        n2 = len(stra)
        n3 = len(strb)        
        j = 1
        k = 0
        do i = 1, n1 - n2 + 1
            if (k>0) then
                k = k - 1
                cycle
            end if
            if(str(i:i+n2-1) == stra(1:n2)) then
                if(n3 > 0) then
                    tempstr(j: j + n3 - 1) = strb(1:n3)
                    j = j + n3
                    k = n2 - 1
                end if
            else                
                tempstr(j:j) = str(i:i)
                j = j + 1
            end if
        end do
        
        str = tempstr
        
    end subroutine replaceCharacter
   
     !Skip n values and return the left value
     !Assume all the values are separated by blank space
     subroutine skipNValues(string, n)
     
        implicit none
        
        character(*), intent(inout) :: string
        integer, intent(in)         :: n
        integer :: i, j
        
        string = adjustl(string)
        
        do i = 1, n
            if(string(1:1) == "'") then
                string = string(2:)
                j = index(string, "'")
                string = adjustl(string(j+1:))
            else
                j = index(string, " ")
                if(j < 1) then
                    string = ""
                    exit
                end if
                string = adjustl(string(j:))
            end if
        end do
        
     end subroutine skipNValues
     
     !count the number of datas in a string
     function countNumberOfDatas(string, separator) result (iNumber)
     
        implicit none
        character(*), intent(in) :: string
        character(*), intent(in) :: separator
        integer :: iNumber        
        integer :: i, j, k  
        character(1024) :: str
        
        str =adjustl(string)
        iNumber = 0
        k = len(separator)
        
        do while(len_trim(str) > 0) 
            j = index(trim(str), separator)
            if( j > 0) then
                iNumber = iNumber + 1
                str = adjustl(str(j+k:))
            else
                iNumber = iNUmber + 1
                exit
            end if
        end do
     
     end function countNumberOfDatas
     
     !Get name from string, string format: 'name' aa bb cc
     subroutine getNameFromString(string, namelen, name)
     
        implicit none
        character(*), intent(in) :: string       
        integer, intent(in) :: namelen 
        character(namelen), intent(out) :: name
        character(1024)             :: tempStr
        integer :: i
        
        i = index(string, "'")
        if (i < 1) then
            call WriteLog("Error detected in converting species name: " // trim(string))
            call ErrorHandling 
        end if
        tempStr = string(i + 1:)

        i = index(tempStr, "'")
        if (i < 2) then
            call WriteLog("Error detected in converting species name: " // trim(string))
            call ErrorHandling 
        end if
        tempStr = adjustl(tempStr(:i-1))
        
        !if (len_trim(tempStr) > nNameLength) then
        !    name = tempStr(1:nNameLength)
        !else
            name = trim(tempStr)
        !end if
        
     end subroutine getNameFromString
     
    !index of first character (a to z) in the string
    function iIndexOfAlphabet(string) result (iIndex)
    
        implicit none    
        character(*), intent(in) :: string  
        integer :: iIndex        
        integer :: i, j, n
        
        iIndex = -1
        n = len(string)
        do i = 1, n
            j = ichar(string(i:i))
            if ((j >= 65 .and. j<= 90) .or. (j>=97 .and. j<=122)) then
                iIndex = i
                exit
            end if
        end do
    
    end function
    
    !check if the string contains only alphabetic
    function bIsAlphabetic(string) result(bFlag)
    
        implicit none
        character(*), intent(in) :: string
        logical :: bFlag
        integer :: i, j, n
        
        bFlag = .true.
        n = len(string)
        do i = 1, n
            j = ichar(string(i:i))
            if (j < 65 .or. (j > 90 .and. j < 97) .or. j > 122) then
                bFlag = .false.
                exit
            end if
        end do
            
    end function bIsAlphabetic
    
    !check if the string contains only alphabetic
    function bIsDigitalPart(string) result(bFlag)
    
        implicit none
        character(*), intent(in) :: string
        logical :: bFlag
        integer :: i, j, n
        
        bFlag = .true.
        n = len(string)
        do i = 1, n
            j = ichar(string(i:i))
            if (j < 43 .or. j > 57) then
                bFlag = .false.
                exit
            end if
        end do
            
    end function bIsDigitalPart

end module