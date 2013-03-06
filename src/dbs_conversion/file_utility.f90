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
        if(abs(dvalue) < 1.0E-100) then
            bFlag = .true.
        else
            bFlag = .false.
        end if
    end function IsZero

    ! Convert upper case to lower case
    subroutine LowerCase(string)

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

    end subroutine lowerCase

    ! replace character in string
    subroutine replaceCharacter(str, stra, strb)
    
        implicit none
        
        character(*), intent(inout) :: str
        character(*), intent(in) :: stra, strb
        character(1024) :: tempstr 
        
        integer :: i, j, k, n1, n2, n3
        
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
                tempstr(j: j + n3 - 1) = strb(1:n3)
                j = j + n3
                k = n2 - 1
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
     
     !Get name from string
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

end module