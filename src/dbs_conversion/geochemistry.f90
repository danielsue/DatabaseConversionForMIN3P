!module of geochemistry
module geochemistry

use logfile, only : WriteLog, nErrors, nWarnings
use global, only : ErrorHandling
use file_utility, only : iIndexOfAlphabet, bIsAlphabetic
use elements, only:  getMWT, bMatchElement

implicit none

real ::  aklog_o2aq(20)         !Equilibrium constants (log(K) in base 10) for o2(aq) or h2(aq) at each discrete temperature listed in record
real ::  aklog_h2aq(20)         !Equilibrium constants (log(K) in base 10) for o2(aq) or h2(aq) at each discrete temperature listed in record


!!Sample data                             
!!O2(aq)=2H2O-2H2(aq)								
!!2	H2O	-1	O2(aq)	-2	H2(aq)			
!!logK	100.8872	92.1438	81.88	72.2846	62.594	54.73	48.1708	42.5362
!!								
!!H2(aq)=H2O-0.5O2(aq)								
!!1	H2O	-0.5	O2(aq)	-1	H2(aq)			
!!logK	50.4436	46.0719	40.94	36.1423	31.297	27.365	24.0854	21.2681

real :: r0tol = 0.00000001d0

contains

!switch master variable for secondary species and minerals
subroutine switchMasterVariable(spename, ncp, namelen, nameofstq, stq, ntemp, aklog,targetvar)

    implicit none

    integer :: ncp          !number of the components defining the secondary species
    integer :: namelen      !name length of the components
    character(*) :: spename
    character(namelen) :: nameofstq(20) !name of the components
    real     :: stq(20)     !stoichiometric coefficients of the components
    integer :: ntemp        !number of temperatures
    real     :: aklog(20)   !the equilibrium constants (log(K) in base 10) for the given reaction at each temperature
    character(60) :: targetvar  !Target database master varialbe

    integer :: i, j
    real    :: stq_h2o, stq_o2aq, stq_h2aq
    logical :: b_include_h2o, b_include_o2aq, b_include_h2aq
    

    b_include_h2o  = .false.
    b_include_o2aq = .false.
    b_include_h2aq = .false.

    stq_h2o = 0.0
    stq_o2aq = 0.0
    stq_h2aq = 0.0

    do i = 1, ncp

        if(trim(nameofstq(i)) == "h2o") then
            b_include_h2o = .true.
            stq_h2o = stq(i)
        end if
        
        if(trim(nameofstq(i)) == "o2(aq)") then
            b_include_o2aq = .true.
            stq_o2aq = stq(i)
        end if

        if (trim(nameofstq(i)) == "h2(aq)") then
            b_include_h2aq = .true.
            stq_h2aq = stq(i)
        end if

    end do

    if( (trim(targetvar) == "o2(aq)" .and. b_include_h2aq == .false.) .or.  &
        (trim(targetvar) == "h2(aq)" .and. b_include_o2aq == .false.) ) then
        return
    end if

    if (trim(targetvar) == "h2(aq)") then
        
         call WriteLog("Switch master variable to h2(aq) for "//trim(spename))
         
         call WriteLog ("aklog before converting:")
         
         call WriteLog (ntemp, aklog(:ntemp))
         
!         call WriteLog ("stq_o2aq:")
!         
!         call WriteLog (stq_o2aq)
!         
!         call WriteLog ("aklog_h2aq:")
!         
!         call WriteLog (ntemp, aklog_h2aq (:ntemp))         

         do i = 1, ntemp
            aklog(i) = aklog(i) + aklog_h2aq (i) * 2 * stq_o2aq            
         end do
         
         call WriteLog ("aklog after converting:")
         
         call WriteLog (ntemp, aklog(:ntemp))        
         
         call WriteLog ("Stoichiometry coefficient before converting:")
         call WriteLog ("h2o ", stq_o2aq)
         call WriteLog ("h2(aq) ", stq_h2aq)
         call WriteLog ("o2(aq) ", stq_o2aq)

         stq_h2o= stq_h2o + 2.0 * stq_o2aq         
         stq_h2aq = -2.0 * stq_o2aq
         stq_o2aq = 1.0
         
         call WriteLog ("Stoichiometry coefficient after converting:")
         call WriteLog ("h2o ", stq_o2aq)
         call WriteLog ("h2(aq) ", stq_h2aq)
         call WriteLog ("o2(aq) ", stq_o2aq)

         !replace o2(aq) to h2(aq)
         do i = 1, ncp
             if(trim(nameofstq(i)) == "o2(aq)") then                
                 nameofstq(i) = "h2(aq)"
                 stq(i) = stq_h2aq
                 exit
             end if
         end do

         !check if the reaction include h2o before converting
         if(b_include_h2o) then
            !if the coefficient of h2o is zero, remove it
            if(abs(stq_h2o) < r0tol) then
                 do i = 1, ncp
                     if(trim(nameofstq(i)) == "h2o") then                
                        do j = i, ncp - 1
                            nameofstq(j) = nameofstq(j + 1) 
                            stq(j) = stq(j + 1)
                        end do
                        ncp = ncp - 1
                        exit
                     end if
                 end do
            else
                 do i = 1, ncp
                     if(trim(nameofstq(i)) == "h2o") then                
                         stq(i) = stq_h2o
                         exit
                     end if
                 end do
            end if
         else
            ncp = ncp + 1
            nameofstq(ncp) = "h2o" 
            stq(ncp) = stq_h2o

         end if

    else if ( trim(targetvar) == "o2(aq)") then
        
         call WriteLog("Switch master variable to o2(aq) for "//trim(spename))
         
         call WriteLog ("aklog before converting:")
            
         call WriteLog (ntemp, aklog(:ntemp))
         
!         call WriteLog ("stq_h2aq:")
!         
!         call WriteLog (stq_h2aq)
!         
!         call WriteLog ("aklog_o2aq:")
!         
!         call WriteLog (ntemp, aklog_o2aq (:ntemp))   

         do i = 1, ntemp
            aklog(i) = aklog(i) + aklog_o2aq (i) * 0.5 * stq_h2aq            
         end do         
         
         call WriteLog ("aklog after converting:")
         
         call WriteLog (ntemp, aklog(:ntemp)) 
         
         call WriteLog ("Stoichiometry coefficient before converting:")
         call WriteLog ("h2o ", stq_o2aq)
         call WriteLog ("h2(aq) ", stq_h2aq)
         call WriteLog ("o2(aq) ", stq_o2aq)

         stq_h2o= stq_h2o + 1.0 * stq_h2aq
         stq_o2aq = -0.5 * stq_h2aq
         stq_h2aq = 0.0
         
         call WriteLog ("Stoichiometry coefficient after converting:")
         call WriteLog ("h2o ", stq_o2aq)
         call WriteLog ("h2(aq) ", stq_h2aq)
         call WriteLog ("o2(aq) ", stq_o2aq)

         !replace h2(aq) to o2(aq)
         do i = 1, ncp
             if(trim(nameofstq(i)) == "h2(aq)") then                
                 nameofstq(i) = "o2(aq)"
                 stq(i) = stq_o2aq
                 exit
             end if
         end do

         !check if the reaction include h2o before converting
         if(b_include_h2o) then
            !if the coefficient of h2o is zero, remove it
            if(abs(stq_h2o) < r0tol) then
                 do i = 1, ncp
                     if(trim(nameofstq(i)) == "h2o") then                
                        do j = i, ncp - 1
                            nameofstq(j) = nameofstq(j + 1) 
                            stq(j) = stq(j + 1)
                        end do
                        ncp = ncp - 1
                        exit
                     end if
                 end do
            else
                 do i = 1, ncp
                     if(trim(nameofstq(i)) == "h2o") then                
                         stq(i) = stq_h2o
                         exit
                     end if
                 end do
            end if

         else
            ncp = ncp + 1
            nameofstq(ncp) = "h2o" 
            stq(ncp) = stq_h2o

         end if

    end if

end subroutine switchMasterVariable

     !Estimate charge from the specie name
     !Format: A+, A+++, A-, A---, A+2, A-3 
     subroutine estimateChargeFromName (string, charge)
     
        implicit none
        character(*), intent(in) :: string
        real, intent(out) :: charge
        integer:: i, i2, j
        
        i = index(string, "+")        
        i2 = index(string, "-")
        
        if (i> 0 .and. i2 > 0) then
            goto 9999
        end if
        
        if (i > 0) then
            j = index(string, "+", back = .true.)
            if (i == j) then    ! A+, A+2
                if(j == len_trim(string)) then  !A+
                    charge = 1.0
                else if (j < len_trim(string)) then  !A+2
                    read(string(j:),* , end = 9999, err = 9999) charge
                end if
            else if (i<j) then  !A+++
                charge = j-i+1.0
            end if
        end if
        
        if (i2 > 0) then
            j = index(string, "-", back = .true.)
            if (i2 == j) then    ! A-, A-2
                if(j == len_trim(string)) then  !A-
                    charge = -1.0
                else if (j < len_trim(string)) then  !A-2
                    read(string(j:),* , end = 9999, err = 9999) charge
                end if
            else if (i<j) then  !A+++
                charge = -1.0*(j-i+1.0)
            end if
        end if
        
        if (i < 1 .and. i2 < 1) then
            charge = 0.0
        end if
        
        return
        
9999    call WriteLog("Error detected in estimating charge from specie name")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end subroutine estimateChargeFromName
     
     !Calculate molecular mass from the formula
     !Note: the formula is case sensitive
     function getMolecularMass (string) result (molecularMass)
     
        implicit none
        character(*), intent(in) :: string
        real :: molecularMass
        real :: rTempCoef, rTempInner
        integer:: i, j, k, k2, n
        
        character(1024) :: str
        
        molecularMass = 0
        str = string
        
        if (trim(str) == "e-") then
            return
        end if
        
        !remove the "charge part"        
        i = index(str, "+") 
        if(i == 1) then
            goto 9999
        else if(i > 1) then
            str = str(:i-1)
        end if
        
        i = index(str, "-") 
        if(i ==1 ) then
            goto 9999
        else if(i > 1) then
            str = str(:i-1)
        end if
        
        str = trim(adjustl(str))    
        n = len_trim(str)
        !calculate the molecular mass
        !only alphabetic (a-z), number(2, 3.5) and bracket pair are allowed in the formula
        i = 1
        do while(i <= n)
            !j = iIndexOfAlphabet(str)            
            if(str(i:i) == "(") then
                j = index(str, ")")
                if(j<i+2) then
                    goto 9999
                end if
                rtempInner = calculateMWTSection(str(i+1:j-1))
                str(i:j) = "#"
                
                k = iIndexOfAlphabet(str)
                k2 = index(str,"(")
                if(k > 0 .and. k2 > 0) then
                    k = min(k, k2)
                end if
                
                if(k > j+1) then
                    read(str(j + 1: k - 1),*, end = 9999, err = 9999) rTempCoef
                    molecularMass = molecularMass + rTempInner * rTempCoef
                    str(j + 1: k - 1)="#"
                    i = k
                else if(k == j + 1) then
                    molecularMass = molecularMass + rTempInner
                    i = k
                else
                    if(j == n) then
                        molecularMass = molecularMass + rTempInner
                    else if (j < n) then
                        read(str(j + 1: n),*, end = 9999, err = 9999) rTempCoef
                        molecularMass = molecularMass + rTempInner * rTempCoef
                    end if    
                    exit                                
                end if
             
            else if (bIsAlphabetic(str(i:i))) then
                j = index(str, "(")
                if(j > i) then
                    molecularMass = molecularMass + calculateMWTSection(str(i:j-1))
                    str(i:j-1) = "#"
                    i = j
                else
                    molecularMass = molecularMass + calculateMWTSection(str(i:n))
                    exit
                end if
            else
                goto 9999            
            end if
        end do
        
        
        return
        
9999    call WriteLog("Error detected in the formula form")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end function getMolecularMass 
     
     !calculate the atomic mass of a section, e.g, in the bracket pair
     function calculateMWTSection(string) result(mass)
     
        implicit none
        character(*), intent(in) :: string
        real :: mass        
        real :: rTempCoef, rTempInner
        character(1024) :: str
        integer :: nlen     
        
        integer :: i, k, n        
        mass = 0        
        str = string
        
        n = len_trim(str)
        i = 1
        do while(i <= n)
            if(i < n) then
                !check if the string contains user-specified specie
                if (bMatchElement(str(i:), nlen, rTempInner)) then
                    mass = mass + rTempInner
                    str(i:i+nlen-1) = "#"
                    i = i + nlen
                    cycle
                end if
                
                rTempInner = getMWT(str(i:i+1))
                if(rTempInner < 0) then
                    rTempInner = getMWT(str(i:i))
                    if(rTempInner < 0) then
                        goto 9999
                    else
                        str(i:i)="#"
                        k = iIndexOfAlphabet(str)
                        if(k > i + 1) then
                            read(str(i + 1: k - 1),*, end = 9999, err = 9999) rTempCoef
                            mass = mass + rTempInner * rTempCoef
                            str(i + 1: k - 1)="#"
                            i = k
                         else if(k == i + 1) then
                            mass = mass + rTempInner
                            i = k
                         else
                            if(i == n) then
                                mass = mass + rTempInner
                            else if (i < n) then
                                read(str(i + 1: n),*, end = 9999, err = 9999) rTempCoef
                                mass = mass + rTempInner * rTempCoef
                            end if
                            exit                                
                        end if
                    end if
                else
                    str(i:i+1)="#"
                    k = iIndexOfAlphabet(str)
                    if(k > i + 2) then
                        read(str(i + 2: k - 1),*, end = 9999, err = 9999) rTempCoef
                        mass = mass + rTempInner * rTempCoef
                        str(i + 2: k - 1)="#"
                        i = k
                    else if(k == i + 2) then
                        mass = mass + rTempInner
                        i = k
                    else
                        if(i + 1 == n) then
                            mass = mass + rTempInner
                        else if (i + 1 < n) then
                            read(str(i + 2: n),*, end = 9999, err = 9999) rTempCoef
                            mass = mass + rTempInner * rTempCoef
                        end if
                        exit                                
                    end if                  
                end if 
            else if (i == n) then
                rTempInner = getMWT(str(i:i))
                if(rTempInner < 0) then
                    goto 9999
                else
                    mass = mass + rTempInner
                end if
                exit
            end if
        end do
        
        return
        
9999    call WriteLog("Error detected in the formula form")
        call WriteLog(trim(string))
        call ErrorHandling  
               
     end function calculateMWTSection
     
     !get specie name without charge part, e.g, CO3-2 to CO3, Ca++ to Ca
     function getSpecieNameWithoutCharge(string) result(str)
     
        implicit none
        character(*), intent(in) :: string
        character(1024) :: str
        integer :: i 
        
        str = string
        
        !remove the "charge part"        
        i = index(str, "+") 
        if(i == 1) then
            goto 9999
        else if(i > 1) then
            str = str(:i-1)
        end if
        
        i = index(str, "-") 
        if(i ==1 ) then
            goto 9999
        else if(i > 1) then
            str = str(:i-1)
        end if
        
        return
        
9999    call WriteLog("Error detected in the formula form")
        call WriteLog(trim(string))
        call ErrorHandling 
     
     end function getSpecieNameWithoutCharge
     
     
end module geochemistry