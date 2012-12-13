!module of file utility
module file_utility

    implicit none

    contains

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

end module