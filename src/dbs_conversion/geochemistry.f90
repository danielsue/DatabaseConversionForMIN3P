!module of geochemistry
module geochemistry

use logfile, only : WriteLog, nErrors, nWarnings
use global, only : ErrorHandling
use file_utility, only : iIndexOfAlphabet, bIsAlphabetic, replaceCharacter
use elements, only:  getMWT, bMatchElement
use sourcedata, only : nNameLength, nSwitchedReactions, switchedReactions

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
    
    
    !switch master variable for secondary species and minerals
    subroutine switchReactionComponent(spename, ncp, namelen, nameofstq, stq, ntemp, aklog, iassociation)

        implicit none

        integer, intent(inout) :: ncp          !number of the components defining the secondary species
        integer, intent(in) :: namelen      !name length of the components
        character(*), intent(in) :: spename
        character(namelen), intent(inout) :: nameofstq(20) !name of the components
        real, intent(inout)     :: stq(20)     !stoichiometric coefficients of the components
        integer, intent(in) :: ntemp        !number of temperatures
        real, intent(inout)     :: aklog(20)   !the equilibrium constants (log(K) in base 10) for the given reaction at each temperature
        integer, intent(in) :: iassociation
    
        integer :: i, j, k, k2
        real :: coeff
        logical :: bFlag
    
        call WriteLog("Switch reaction component for "//trim(spename))
        
        do i = 1, nSwitchedReactions
            do j = 1, ncp
                if(trim(switchedReactions(i)%Name) == trim(nameofstq(j))) then
                    call WriteLog("switch "//trim(switchedReactions(i)%Name))
                    coeff = stq(j) * iassociation * switchedReactions(i)%iAssociation
                    !recalculate logK
                    do k = 1, ntemp
                        aklog(k) = aklog(k) + coeff * switchedReactions(i)%AKLOG(k)
                    end do
                    !recalculate stoichiometric coefficient                    
                    do k = 1, switchedReactions(i)%NCP
                        bFlag = .true.
                        do k2 = 1, ncp
                            if(trim(nameofstq(k2)) == trim(switchedReactions(i)%NameOfSTQ(k))) then
                                stq(k2) = stq(k2) + stq(j)*switchedReactions(i)%STQ(k)                                
                                bFlag = .false.
                                exit
                            end if
                        end do
                        if(bFlag) then
                            ncp = ncp + 1
                            nameofstq(ncp) = trim(switchedReactions(i)%NameOfSTQ(k))
                            stq(ncp) = stq(j)*switchedReactions(i)%STQ(k)
                        end if                        
                    end do
                    !remove the switched component
                    nameofstq(j) = nameofstq(ncp)
                    stq(j) = stq(ncp)
                    ncp = ncp - 1
                    !remove other component if the stoichiometric coefficient is zero
                    k = 1
                    do while(k <= ncp)
                        if (abs(stq(k)) < r0tol) then
                            nameofstq(k) = nameofstq(ncp)
                            stq(k) = stq(ncp)
                            ncp = ncp - 1 
                        else
                            k = k + 1
                        end if
                    end do
                    !output to log file
                    call WriteLog("logK")
                    call WriteLog(ntemp, aklog)
                    do k = 1, ncp
                        call WriteLog("  "//trim(nameofstq(k)), stq(k))
                    end do
                    exit
                end if
            end do 
            
        end do
    
    end subroutine switchReactionComponent

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
     
     !calculate molecular weight from the components
     function calculateMWTFromComponents(ncp, stq, nameofstq) result(mwt)
     
        implicit none
        
        integer, intent(in) :: ncp
        real, intent(in) :: stq(20)
        character(nNameLength), intent(in) :: nameofstq(20)
        real :: mwt, tempmwt
        
        integer :: i
        
        mwt = 0.0d0
        
        do i = 1, ncp
            
            !tempmwt = getMWTFromMasterSpecies(trim(nameofstq(i)))            
            !if(tempmwt < 0) then
            tempmwt = getMolecularMass(trim(nameofstq(i)))
            !end if
            
            mwt = mwt + stq(i) * tempmwt
            
        end do
     
     end function calculateMWTFromComponents
     
     !analyse reaction equation
     !Format 1 association reaction: +1.000Al+3                 +1.000SO4-2                = AlSO4+
     !Format 2 association reaction: +1.000K+                   +1.000H2O                  -1.000H+                   = KOH
     !Format 3 neither association nor dissociation: Feii+2 + H2O = FeiiOH+ + H+
     !Format 4 disassociation reaction: SiO2 + 1OH- + 1H2O =  SiO(OH)3-
     !Format 5 disassociation reaction: SiO2 + 1O H- + 1 H2O =  SiO(OH)3-
     !Note: the stoichiometric coefficient of the master species should be 1.0
     subroutine analyseReactionEquation(string, namelen, bRightFirst, masterSpecieName, nSpeciesOut, NameOfSTQ, STQ, bCollectMasterName)
     
        implicit none
        character(*), intent(in) :: string      !reaction equation
        integer, intent(in) :: namelen          !maximum length of the specie name
        logical, intent(in) :: bRightFirst      !If true, the right first specie after "=" is the master speice name and the reaction is association reaction;
                                                !otherwise, the left first specie is the master specie name and the reaction is disassociation reaction.
        character(namelen), intent(out) :: masterSpecieName !name of master specie                                        
        integer, intent(out) :: nSpeciesOut     !Number of species in the reaction
        character(namelen), intent(out) :: NameOfSTQ(20)    !name of all the species in the reaction, the first one is the master specie 
        real, intent(out) :: STQ(20)            !stoichiometric coefficients of all the species in the reaction. 
                                                !Note: the reaction will be converted into association or disassociation reaction.
        logical, optional, intent(in) :: bCollectMasterName !check if the reaction string start with the master specie name.
                                                            !SiO2 ; SiO2 + 1OH- + 1H2O =  SiO(OH)3-	; log_K   	1.475988; for this reaction, bCollectMasterName = .true.
                                                            !SiO2 + 1OH- + 1H2O =  SiO(OH)3-	; log_K   	1.475988; for this reaction, bCollectMasterName = .false.
        
        character(1024) :: tempStr
        character(1024) :: tempStr2
        integer :: indexEq                      !index of "="
        integer :: i, j, k, ii   
        character(2) :: labels(2)
        
        tempStr = ""
        tempStr2 = ""
        
        labels = (/" +", " -"/)
                                                
        masterSpecieName = ""                                       
        nSpeciesOut = 0
        NameOfSTQ = ""
        STQ = 0
        
        indexEq = index(string, "=")         
       
        !Association reaction, this is for secondary species
        if(bRightFirst) then
            if((.not. present(bCollectMasterName)) .or. (present(bCollectMasterName) .and. bCollectMasterName)) then
                !master specie
                call getSpecieName(string, namelen, masterSpecieName, bRightFirst)            
                !add the first one
                tempStr = trim(adjustl(string(:indexEq - 1)))
                call WriteLog("Analyse: " // trim(tempStr))
                nSpeciesOut = nSpeciesOut + 1                
                j = iIndexOfAlphabet(tempStr)  
           
                if (j>1) then
                    tempStr2 = tempStr(:j-1)
                    call replaceCharacter(tempStr2, " ", "")
                    call WriteLog("Read coefficient from: " // trim(tempStr2))
                    read(tempStr2, *, end = 9999, err = 9999) STQ(nSpeciesOut)
                    call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))                
                else if (j == 1) then
                    STQ(nSpeciesOut) = 1.0
                    call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                else
                    nErrors = nErrors + 1
                    call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                    call WriteLog(trim(string))
                    call ErrorHandling
                end if
            end if
            
            !!left of "="
            do k = 1, 2                
                tempStr = string(:indexEq - 1)            
                i = index(tempStr, labels(k))                
                do while(i>0) 
                    nSpeciesOut = nSpeciesOut + 1  
                    tempStr = tempStr(i+1: )
                    call WriteLog("Analyse: " // trim(tempStr))
                    j = iIndexOfAlphabet(tempStr)
                    if (j>1) then
                        if(trim(adjustl(tempStr(:j-1))) == "+") then
                            STQ(nSpeciesOut) = 1.0
                        else if(trim(adjustl(tempStr(:j-1))) == "-") then
                            STQ(nSpeciesOut) = -1.0
                        else
                            tempStr2 = tempStr(:j-1)
                            call replaceCharacter(tempStr2, " ", "") 
                            read(tempStr2,* , end = 9999, err = 9999) STQ(nSpeciesOut)
                        end if
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else if (j == 1) then
                        STQ(nSpeciesOut) = 1.0
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else
                        nErrors = nErrors + 1
                        call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                        call WriteLog(trim(string))
                        call ErrorHandling
                    end if
                    
                    i = index(tempStr, labels(k))
                end do
            end do
            
            !!right of "=", the sign of stoichiometric coefficient should be reversed for Association reaction
            ii = index(string, ";")
            if (ii < 1) then
                ii = len_trim(string)
            end if
            
            do k = 1, 2                
                tempStr = string(indexEq+1: ii - 1)            
                i = index(tempStr, labels(k))                
                do while(i>0) 
                    nSpeciesOut = nSpeciesOut + 1  
                    tempStr = tempStr(i+1: )
                    call WriteLog("Analyse: " // trim(tempStr))
                    j = iIndexOfAlphabet(tempStr)
                    if (j>1) then
                        if(trim(adjustl(tempStr(:j-1))) == "+") then
                            STQ(nSpeciesOut) = -1.0
                        else if(trim(adjustl(tempStr(:j-1))) == "-") then
                            STQ(nSpeciesOut) = 1.0
                        else
                            tempStr2 = tempStr(:j-1)
                            call replaceCharacter(tempStr2, " ", "") 
                            read(tempStr2,* , end = 9999, err = 9999) STQ(nSpeciesOut)
                            STQ(nSpeciesOut) = -STQ(nSpeciesOut)
                        end if
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else if (j == 1) then
                        STQ(nSpeciesOut) = -1.0
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else
                        nErrors = nErrors + 1
                        call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                        call WriteLog(trim(string))
                        call ErrorHandling
                    end if
                    
                    i = index(tempStr, labels(k))
                end do
            end do
            
        !Disassociation reaction, for gases and minerals
        else
            !get name
            ii = 0
            if((.not. present(bCollectMasterName)) .or. (present(bCollectMasterName) .and. bCollectMasterName)) then
                ii = index(string, ";") 
                if (ii < indexEq) then
                    masterSpecieName = trim(adjustl((string(:ii-1))))
                else
                    nErrors = nErrors + 1
                    call WriteLog("Error: Unrecognized specie name. The species name should be in the head of the data followed by the reaction equation.")
                    call WriteLog(trim(string))
                    call ErrorHandling
                end if    
            end if
            !!left of "="
            do k = 1, 2                
                tempStr = trim(adjustl(string(ii + 1:indexEq - 1)))   
                
                i = index(tempStr, labels(k))                
                do while(i>0) 
                    nSpeciesOut = nSpeciesOut + 1  
                    tempStr = tempStr(i+1: )
                    call WriteLog("Analyse: " // trim(tempStr))
                    j = iIndexOfAlphabet(tempStr)
                    if (j>1) then
                        if(trim(adjustl(tempStr(:j-1))) == "+") then
                            STQ(nSpeciesOut) = -1.0
                        else if(trim(adjustl(tempStr(:j-1))) == "-") then
                            STQ(nSpeciesOut) = 1.0
                        else
                            tempStr2 = tempStr(:j-1)
                            call replaceCharacter(tempStr2, " ", "") 
                            read(tempStr2,* , end = 9999, err = 9999) STQ(nSpeciesOut)
                            STQ(nSpeciesOut) = -STQ(nSpeciesOut)
                        end if
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else if (j == 1) then
                        STQ(nSpeciesOut) = -1.0
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else
                        nErrors = nErrors + 1
                        call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                        call WriteLog(trim(string))
                        call ErrorHandling
                    end if
                    
                    i = index(tempStr, labels(k))
                end do
            end do
            
            !!right of "=", the sign of stoichiometric coefficient should be reversed for Association reaction
            ii = ii + index(string(ii+1:), ";")
            if (ii < 1) then
                ii = len_trim(string)
            end if
            
            !add the first one
            tempStr = trim(adjustl(string(indexEq+1: ii - 1)))
            call WriteLog("Analyse: " // trim(tempStr))
            nSpeciesOut = nSpeciesOut + 1                
            j = iIndexOfAlphabet(tempStr)  
           
            if (j>1) then
                tempStr2 = tempStr(:j-1)
                call replaceCharacter(tempStr2, " ", "")
                call WriteLog("Read coefficient from: " // trim(tempStr2))
                read(tempStr2, *, end = 9999, err = 9999) STQ(nSpeciesOut)
                call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))                
            else if (j == 1) then
                STQ(nSpeciesOut) = 1.0
                call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
            else
                nErrors = nErrors + 1
                call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                call WriteLog(trim(string))
                call ErrorHandling
            end if
            
            
            do k = 1, 2                
                tempStr = trim(adjustl(string(indexEq+1: ii - 1)))
                
                i = index(tempStr, labels(k))                
                do while(i>0) 
                    nSpeciesOut = nSpeciesOut + 1  
                    tempStr = tempStr(i+1: )
                    call WriteLog("Analyse: " // trim(tempStr))
                    j = iIndexOfAlphabet(tempStr)
                    if (j>1) then
                        if(trim(adjustl(tempStr(:j-1))) == "+") then
                            STQ(nSpeciesOut) = 1.0
                        else if(trim(adjustl(tempStr(:j-1))) == "-") then
                            STQ(nSpeciesOut) = -1.0
                        else
                            tempStr2 = tempStr(:j-1)
                            call replaceCharacter(tempStr2, " ", "") 
                            read(tempStr2,* , end = 9999, err = 9999) STQ(nSpeciesOut)
                        end if
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else if (j == 1) then
                        STQ(nSpeciesOut) = 1.0
                        call getSpecieName(tempStr(j:), namelen, NameOfSTQ(nSpeciesOut))
                    else
                        nErrors = nErrors + 1
                        call WriteLog("Error: Unrecognized specie name, the specie name should start with alphabet")
                        call WriteLog(trim(string))
                        call ErrorHandling
                    end if
                    
                    i = index(tempStr, labels(k))
                end do                
            end do
            
        end if
        
        return
        
9999    call WriteLog("Error detected in analysing reaction equation")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end subroutine analyseReactionEquation
     
     !Get name from string, seperated by blank space, string format: name aa bb cc
     !NO3- = NO3-                 : bRightFirst, NO3-; default, NO3-
     !Ca+2 + H2O = CaOH+  + H+    : bRightFirst, CaOH+; default, Ca+2
     !Format 1: 
     !+1.000Al+3   +1.000SO4-2  = AlSO4+ ;  -gamma   4.00 0.064	; log_K	3.900111 ; -analytical_expression -59.74759294081 0	2048.61497066133 22.9451885920	0
     !Format 2: Format 2 will be converted into format 1 in the program
     !+1.000Al+3   +1.000SO4-2  = AlSO4+ 
     !-gamma   4.00 0.064	
     !log_K	3.900111 
     !-analytical_expression -59.74759294081 0	2048.61497066133 22.9451885920	0

     subroutine getSpecieName(string, namelen, name, bRightFirst)
     
        implicit none
        character(*), intent(in) :: string       
        integer, intent(in) :: namelen 
        character(namelen), intent(out) :: name
        character(1024)             :: tempStr
        logical, optional :: bRightFirst
        integer :: i, j, k
        
        tempStr = trim(adjustl(string))
        
        if (present(bRightFirst)) then 
            if (bRightFirst) then
                i = index(tempStr, "=") 
                if (i > 1) then
                    tempStr = adjustl(tempStr(i+1:))
                    i = index(trim(tempStr), " ")
                    j = index(trim(tempStr), ";")
                    if(i < 1) then
                        i = len_trim(tempStr)
                    else
                        i = i -1
                    end if
                    if(j < 1) then
                        j = len_trim(tempStr)
                    else
                        j = j - 1
                    end if
                    k = min(i, j)
                    name = trim(tempStr(:k))
                else
                    call WriteLog("Error detected in getting specie name from string")
                    call WriteLog(trim(string))
                    call ErrorHandling                     
                end if
            else
                i = index(tempStr, " ")
                if (i < 1) then
                    name = trim(tempStr)            
                else
                    name = trim(tempStr(1:i))
                end if  
            end if
        else
            i = index(tempStr, " ")
            if (i < 1) then
                name = trim(tempStr)            
            else
                name = trim(tempStr(1:i))
            end if     
        end if
     
     end subroutine getSpecieName
     
     
end module geochemistry