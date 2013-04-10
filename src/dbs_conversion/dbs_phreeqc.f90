!Module to read and store Phreeqc database
    
module dbs_phreeqc 

    use global, only : iUnitDbsPhreeqc, bOpenDbsPhreeqc, ErrorHandling
    
    use logfile, only : WriteLog, nWarnings, nErrors

    use file_utility, only : SetLowerCase, GetLowerCase, IsZero, getNameFromString, skipNValues, replaceCharacter, countNumberOfDatas
    
    use alias, only : GetNameFromAlias

    use inputfile, only : sourceDatabasePath

    use geochemistry, only : aklog_o2aq, aklog_h2aq, estimateChargeFromName, getMolecularMass, getSpecieNameWithoutCharge
    
    use elements, only : InitElementList, userDefinedElementList, bIncludeElement, reorderUserDefinedElements

    use sourcedata
    
    implicit none
    
    character(1024)      ::      filePathDbsPhreeqc   !file path for Phreeqc database
    
    character(1),   parameter   ::  strComment      =  "#"
    character(9),   parameter   ::  strSection      = "#section:"
    character(13),  parameter   ::  strSectionEnd   = "#section: end"
    character(33),  parameter   ::  strSectionMaster = "#section: solution master species"
    character(31),  parameter   ::  strSectionPMatchMaster = "#section: pmatch master species"
    character(41),  parameter   ::  strSectionPmatchSecondMaster = "#section: pmatch secondary master species"
    character(15),  parameter   ::  strSectionGases = "#section: gases"
    character(18),  parameter   ::  strSectionMinerals = "#section: minerals"

    character(1),   parameter   ::  strData         =   "'"
    character(4),   parameter   ::  strSection1     =   "null"
    character(3),   parameter   ::  strSection2     =   "end"
    character(6),   parameter   ::  strRedoxReaction = "o2(aq)"

    integer, parameter          ::  maxStrBuffers = 100000
    character(1024)              ::  strBuffers(maxStrBuffers)       !Assume the maximum number of each data section 
    character(1024)              ::  strBuffer        =   ""
 
    integer                      ::  strBuffersLineIndex(maxStrBuffers)
    integer                      ::  iReadStat       =   1
    logical                      ::  bEndOfFile      =   .false.
    
    integer :: currentLine = 0        
    
    type(TypePrimarySpecies), allocatable :: masterSpecies(:)

    
contains

    ! Open Phreeqc database
    subroutine OpenDbsPhreeqc
    
        implicit none
        
        integer :: stat 

        filePathDbsPhreeqc = trim(sourceDatabasePath)
        
        open(iUnitDbsPhreeqc, file = filePathDbsPhreeqc,  iostat = stat)
        
        if (stat /= 0) then
            call WriteLog("Error in opening Phreeqc databse file: " // trim(filePathDbsPhreeqc))
            call ErrorHandling
        else
            bOpenDbsPhreeqc = .true.
            call WriteLog("Open Phreeqc database success: " // trim(filePathDbsPhreeqc))
        end if 
    
    end subroutine    
    
    
    ! Read Phreeqc database
     subroutine ReadDbsPhreeqc
     
        implicit none  
        
        call WriteLog ("Start initializing element and atomic mass list")
        
        call InitElementList
        
        call WriteLog ("End initializing element and atomic mass list")
        
        call WriteLog ("Start reading Phreeqc database")
        
        call readTemperature
        
        call readSpecies
        
        call readAqueousSpecies
        
        call readGases
        
        call readMinerals

        call WriteLog ("End reading Phreeqc database")
        
        call WriteLog ("Start converting specie name into lower case")
        
        call SetLowerCaseAll
        
        call WriteLog("End converting specie name into lower case")
     
     end subroutine ReadDbsPhreeqc
     
     ! Close Phreeqc database
     subroutine CloseDbsPhreeqc
     
        implicit none
        
        if (bOpenDbsPhreeqc) then
            
            close(iUnitDbsPhreeqc)
            
        end if
     
     end subroutine CloseDbsPhreeqc

  
    
    ! Read next line
    function bReadNextLine(bSetLowerCase) result(bFlag) 
    
        implicit none
        logical, optional :: bSetLowerCase
        logical :: bFlag
        integer :: i
        
        
        bFlag = .false.
        
        do while(.not. bEndOfFile)

            read(iUnitDbsPhreeqc, "(a)", iostat = iReadStat) strBuffer
            
            currentLine = currentLine + 1
       
            if (iReadStat > 0) then     !Error in reading
            
                call WriteLog("Error in reading Phreeqc databse file.")
            
                call ErrorHandling
            
            else if (iReadStat < 0) then    !End of file
            
                bEndOfFile = .true.
               
                exit            
            
            end if

            strBuffer = adjustl(strBuffer)

            if(len_trim(strBuffer) /= 0 ) then
            
                ! Conver all case to lower case  
                if(present(bSetLowerCase)) then
                    if(bSetLowerCase) then
                        call SetLowerCase(strBuffer)
                    end if
                end if

                if(GetLowerCase(9, strBuffer(1:9)) == strSection) then
                    call replaceCharacter(strBuffer, achar(9), " ")
                    bFlag = .true.
                    exit
                end if

                if (GetLowerCase(1,strBuffer(1:1)) /= strComment) then
                    call replaceCharacter(strBuffer, achar(9), " ")
                    i = index(strBuffer, strComment)
                    if (i > 0) then
                        strBuffer = strBuffer(:i-1)
                    end if
                    
                    bFlag = .true.                    
                    exit                    
                end if

            end if

        end do
    
    end function bReadNextLine
    
    !
    !Check if the line is the indication of master species or reactions
    !Line 0:  SOLUTION_SPECIES 
    !Line 1a: SO4-2 = SO4-2                         !iType = 1, This is master species than can be converted into min3p database
    !Line 2a:      log_k     0.0                    !iType = 0
    !Line 5a:      -gamma    5.0     -0.04
    !Line 1b: SO4-2 + 9H+ + 8e- = HS- + 4H2O        !iType = 2, This can not be converted into min3p database at present, skipped
    !Line 2b:      log_k     33.652
    !Line 3b:      delta_h   -40.14
    !Line 5b:      -gamma    3.5     0.0
    !
    function iTypeOfData(str, strEndLabel) result (iType)
    
        implicit none
        
        character(*), intent(in) :: str
        character(*), intent(in) :: strEndLabel
        integer :: iType        !0 - unknown, 1 - master specie, 2 - reaction
        integer :: i, j
        
        if (trim(str) == trim(strEndLabel)) then
            iType = -1
            return
        end if  

        iType = 0 
        
        i = index(str, "=")
        
        if (i > 0) then
            j = index(str, strComment)
            if(j > i) then
                if (trim(adjustl(str(1:i-1))) == trim(adjustl(str(i+1:j-1)))) then
                    iType = 1       !master species, left == right
                else
                    iType = 2       !reaction equation
                end if
            else
                if (trim(adjustl(str(1:i-1))) == trim(adjustl(str(i+1:)))) then
                    iType = 1       !master species, left == right
                else
                    iType = 2       !
                end if
            end if

        end if
 
    end function iTypeOfData

    !check if the string contains keyword
    function iTypeKeyWord(string, strEndLabel) result (iType)
    
        implicit none
        
        
        character(*), intent(in) :: string
        character(*), intent(in) :: strEndLabel
        integer :: iType        !0 - unknown, 1 - master specie, 2 - reaction
        integer :: i, j
        
        character(1024) :: str
        str = string
        
        call SetLowerCase(str)
        
        if (trim(str) == trim(strEndLabel)) then
            iType = -1
            return
        end if  

        iType = 0 
        
        if (index(str,"log_k")>0) then
            iType = 1           !log_k data
        end if

        if (index(str,"-gamma")>0) then
            iType = 2           !delta_h data
        end if

        if (index(str,"delta_h")>0) then
            iType = 3          !delta_h data
        end if

        if (index(str,"-analytic")>0) then
            iType = 4           !delta_h data
        end if

        if (index(str,"-no_check")>0) then
            iType = 5           !delta_h data
        end if      

 
    end function iTypeKeyWord  
    
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
     
     !Get Debye-Hukel a, Debye-Hukel b from string, label -gamma
     subroutine getDH_AB(string, dha, dhb)
     
        implicit none
        character(*), intent(in) :: string
        real, intent(out) :: dha, dhb
        integer :: i
        
        character(1024) :: str
        str = string
        call SetLowerCase(str)
        
        i = index(str, "-gamma")
        
        write(*,*) "index of -gamma", i
        
        if (i > 0) then
            read(str(i + 6: ),*, end = 9999, err = 9999) dha, dhb
        else
            dha = 0
            dhb = 0
        end if
        
        return
        
9999    call WriteLog("Error detected in getting Deebye-huckel parameters from the follow string")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end subroutine getDH_AB
     
     !Get equilibrium constant from the string
     subroutine getLogK(string, logK)
     
        implicit none
        character(*), intent(in) :: string
        real, intent(out) :: logK
        integer :: i
        
        character(1024) :: str
        str = string
        call SetLowerCase(str)
        
        i = index(str, "log_k")
        
        if (i > 0) then
            read(str(i + 5: ),*, end = 9999, err = 9999) logK
        else
            logK = 0            
        end if
        
        return
        
9999    call WriteLog("Error detected in getting logK from the follow string")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end subroutine getLogK
     
     
    
     !Get equilibrium constant from the string
     subroutine getEnthalpyChange(string, enthalpyChange)
     
        implicit none
        character(*), intent(in) :: string
        real, intent(out) :: enthalpyChange
        integer :: i
        
        character(1024) :: str
        str = string
        call SetLowerCase(str)
        
        i = index(str, "delta_h")
        
        if (i > 0) then
            read(str(i + 7: ),*, end = 9999, err = 9999) enthalpyChange
        else
            enthalpyChange = 0            
        end if
        
        return
        
9999    call WriteLog("Error detected in getting enthalpy change from the follow string")
        call WriteLog(trim(string))
        call ErrorHandling  
     
     end subroutine getEnthalpyChange
     
     !analyse reaction equation
     !Format 1 association reaction: +1.000Al+3                 +1.000SO4-2                = AlSO4+
     !Format 2 association reaction: +1.000K+                   +1.000H2O                  -1.000H+                   = KOH
     !Format 3 neither association nor dissociation: Feii+2 + H2O = FeiiOH+ + H+
     !Format 4 disassociation reaction: SiO2 + 1OH- + 1H2O =  SiO(OH)3-
     !Format 5 disassociation reaction: SiO2 + 1O H- + 1 H2O =  SiO(OH)3-
     !Note: the stoichiometric coefficient of the master species should be 1.0
     subroutine analyseReactionEquation(string, namelen, bRightFirst, masterSpecieName, nSpeciesOut, NameOfSTQ, STQ)
     
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
            ii = index(string, ";") 
            if (ii < indexEq) then
                masterSpecieName = trim(adjustl((string(:ii-1))))
            else
                nErrors = nErrors + 1
                call WriteLog("Error: Unrecognized phase (gas or mineral) name. The phase name should be in the head of the data followed by the reaction equation.")
                call WriteLog(trim(string))
                call ErrorHandling
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

    
     ! Read temperature points
     subroutine readTemperature

        implicit none
               
        nTemperature = 1                

     end subroutine readTemperature
    
    
     ! Read species
     subroutine readSpecies

        implicit none

        integer :: i, j, k
        character (nNameLength) :: strName, strName2
        logical :: bFindSpecies
        integer :: nTempSpecies = 0

        nSpecies = 0        
        bEndOfFile = .false. 
        
        ! read SOLUTION_SPECIES data
        rewind iUnitDbsPhreeqc        
        bFindSpecies = .false.
        do while (.not. bEndOfFile)
            if(.not. bReadNextLine()) then
                exit
            end if
            if (trim(adjustl(strbuffer)) == strSectionPMatchMaster) then
                bFindSpecies = .true.
                exit
            end if
        end do
        
        if (bFindSpecies) then
            
            do while (.not. bEndOfFile)
                
                if (nSpecies == maxStrBuffers) then
                    nErrors = nErrors + 1
                    call WriteLog("Error: The maximum database of Phreeqc is 100,000 lines for each section." // & 
                                    "Please modify the source code parameter maxStrBuffers to read more data.")
                    call ErrorHandling
                end if
                
                if(.not. bReadNextLine()) then
                    exit
                end if                
                !Check if it is a master species
                if (iTypeOfData(strbuffer, strSectionEnd) == 1 ) then                   
                    nSpecies = nSpecies + 1
                    strBuffers(nSpecies) = strBuffer
                    strBuffersLineIndex(nSpecies) = currentLine                
                else if (iTypeOfData(strbuffer, strSectionEnd) == 0 ) then
                    strBuffers(nSpecies) = trim(strBuffers(nSpecies)) // " ; " // trim(strBuffer)
                else if (iTypeOfData(strbuffer, strSectionEnd) == -1) then 
                    exit
                end if
            end do
            
            if (nSpecies > 0) then
                if(allocated(species)) then
                    deallocate(species)
                end if
                allocate(species(nSpecies))
                do i = 1, nSpecies
                    call WriteLog("Species data: " // trim (strBuffers(i)))
                    call getSpecieName(strBuffers(i), nNameLength, species(i)%Name)
                    
                    !Check if this name is used as an alias
                    strName2 = GetLowerCase(len(species(i)%Name), species(i)%Name)
                    strName = GetNameFromAlias(strName2)
                    if (trim(strName) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(strName2))
                        call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(strName2) // " are the same")
                    end if
                    
                    !Read Debye-Hukel a, Debye-Hukel b 
                    call getDH_AB(strBuffers(i), species(i)%DHA, species(i)%DHB)
                    
                    !Set charge of the species
                    call estimateChargeFromName (species(i)%Name, species(i)%Z)
                    
                end do
            end if
        
        end if
        
!read SOLUTION_MASTER_SPECIES data
!read alkalinity factor and molecular weight
        
!!#  elemen     species        alk   gfw_formula element_gfw atomic 
!!#                                                          number
!!#
!!H               H+          -1.0        H         1.008   #   1       !!h+, OK
!!H(0)            H2           0.0        H                 #    
!!H(1)            H+          -1.0        H                 #    
!!E               e-           0.0        0.0        0.0    #           !!e-, OK
!!O               H2O          0.0        O        15.999   #   8       !!h2o, OK
!!O(0)            O2           0.0        O                 #        
!!O(-2)           H2O          0.0        O                 #    
!!C              HCO3-         1.0        C        12.011   #   6       !!hco3-, ERROR, the value is not correct, should take 61.016 instead
!!C(+4)          HCO3-         1.0      HCO3-               #    
!!C(-4)           CH4          0.0       CH4                #    
!!Alkalinity     HCO3-         1.0      HCO3-      61.016   #           !!Note: molecular for hco3- should be manually set, the program only take the first one, 12.011

        rewind iUnitDbsPhreeqc        
        bFindSpecies = .false.
        do while (.not. bEndOfFile)
            if(.not. bReadNextLine()) then
                exit
            end if
            if (trim(adjustl(strbuffer)) == strSectionMaster) then
                bFindSpecies = .true.
                exit
            end if
        end do
        
        if (bFindSpecies) then
            
            do while (.not. bEndOfFile)
                
                if (nTempSpecies == maxStrBuffers) then
                    nErrors = nErrors + 1
                    call WriteLog("Error: The maximum database of Phreeqc is 100,000 lines for each section." // & 
                                    "Please modify the source code parameter maxStrBuffers to read more data.")
                    call ErrorHandling
                end if
                
                if(.not. bReadNextLine()) then
                    exit
                end if                
                !Check if it is a master species
                if (iTypeOfData(strbuffer, strSectionEnd) == 0) then                   
                    nTempSpecies = nTempSpecies + 1
                    strBuffers(nTempSpecies) = strBuffer
                    strBuffersLineIndex(nTempSpecies) = currentLine
                else if (iTypeOfData(strbuffer, strSectionEnd) == -1) then 
                    exit
                end if
                
            end do
            
            if (nTempSpecies > 0) then
                
                if (allocated(masterSpecies)) then
                    deallocate(masterSpecies)
                end if
                allocate(masterSpecies(nTempSpecies))
                
                j = 0                
                do i = 1, nTempSpecies
                    call getSpecieName(strBuffers(i), nNameLength, masterSpecies(i)%Name)
                    call skipNValues(strBuffers(i),2)
                    read(strBuffers(i),*, end = 9999, err = 9999) masterSpecies(i)%AlkFac
                    call skipNValues(strBuffers(i),2)
                    if(len_trim(strBuffers(i)) > 0) then
                        read(strBuffers(i),*, end = 9999, err = 9999) masterSpecies(i)%MWT
                        if(.not. bIncludeElement(getSpecieNameWithoutCharge(masterSpecies(i)%Name))) then
                            j = j + 1
                        end if
                    else
                        masterSpecies(i)%MWT = -1.0d0
                    end if                    
                end do

                if (j > 0) then
                    if (allocated(userDefinedElementList)) then
                        deallocate(userDefinedElementList)
                    end if
                    allocate(userDefinedElementList(j))   
                    k = 0
                    do i = 1, nTempSpecies  
                        if(masterSpecies(i)%MWT > 0) then
                            if(.not. bIncludeElement(getSpecieNameWithoutCharge(masterSpecies(i)%Name))) then
                                k = k + 1
                                userDefinedElementList(k)%Name = trim(getSpecieNameWithoutCharge(masterSpecies(i)%Name))
                                userDefinedElementList(k)%MWT = masterSpecies(i)%MWT
                                call WriteLog("User-defined specie added: "//trim(userDefinedElementList(k)%Name), userDefinedElementList(k)%MWT)
                            end if
                        end if
                    end do
                    call reorderUserDefinedElements
                    call WriteLog("Reordered user-defined species")
                    do i =1, size(userDefinedElementList,1)
                        call WriteLog(trim(userDefinedElementList(i)%Name), userDefinedElementList(i)%MWT)
                    end do
                end if
                
                do i = 1, nSpecies                    
                    do j = 1, nTempSpecies                        
                        if (trim(species(i)%Name) == trim(masterSpecies(j)%Name)) then
                            species(i)%AlkFac = masterSpecies(j)%AlkFac
                            exit
                        end if                    
                    end do                    
                end do
                
            end if
        
        end if
        
        if (nSpecies > 0) then
            do i = 1, nSpecies 
                
                !species(i)%MWT = getMWTFromMasterSpecies(species(i)%Name)
                !if(species(i)%MWT < 0) then
                    !calculate molecular weight
                    species(i)%MWT = getMolecularMass(trim(species(i)%Name))
                    call WriteLog("Calculate molecular weight for " // trim(species(i)%Name),  species(i)%MWT)
                !end if
            end do
        end if
        

        call WriteLog("Number of Species: ")
        call WriteLog(nSpecies)
        call WriteLog("Read species success")
        
	    return

9999    call WriteLog("Error detected in reading/converting data in the follow string")
        call WriteLog(trim(strBuffer))
        call ErrorHandling    

     end subroutine readSpecies
     
     ! Extract molecular weight from the solution master species section
     function getMWTFromMasterSpecies(string) result(mwt)
     
        implicit none
        character(*), intent(in) :: string
        real :: mwt
        
        integer :: i, n
        
        n = size(masterSpecies, 1)
        
        mwt = -1.0d0
        
        do i = 1, n
            if(trim(masterSpecies(i)%Name) == trim(string)) then
                mwt = masterSpecies(i)%MWT
                exit
            end if
        end do
        
     
     end function getMWTFromMasterSpecies

     ! Read aqueous species
     subroutine readAqueousSpecies
        
        implicit none
        
        integer :: i, j, k
        character(nNameLength) :: strName
        logical :: bFindSpecies
        real :: rTempMWT

        i = 0
        nAqueousSpecies = 0
        nRedoxReaction = 0
        
        ! read SOLUTION_SPECIES data
        bEndOfFile = .false.
        rewind iUnitDbsPhreeqc        
        bFindSpecies = .false.
        do while (.not. bEndOfFile)
            if(.not. bReadNextLine()) then
                exit
            end if
            if (trim(adjustl(strbuffer)) == strSectionPmatchSecondMaster) then
                bFindSpecies = .true.
                exit
            end if
        end do
        
        if (bFindSpecies) then
            
            do while (.not. bEndOfFile)
                
                if (nAqueousSpecies == maxStrBuffers) then
                    nErrors = nErrors + 1
                    call WriteLog("Error: The maximum database of Phreeqc is 100,000 lines for each section." // & 
                                    "Please modify the source code parameter maxStrBuffers to read more data.")
                    call ErrorHandling
                end if
                
                if(.not. bReadNextLine()) then
                    exit
                end if               

                if (iTypeOfData(strbuffer, strSectionEnd) == 1  .or. iTypeOfData(strbuffer, strSectionEnd) == 2 ) then                   
                    nAqueousSpecies = nAqueousSpecies + 1
                    strBuffers(nAqueousSpecies) = strBuffer
                    strBuffersLineIndex(nAqueousSpecies) = currentLine 
                else if (iTypeOfData(strbuffer, strSectionEnd) == 0) then     
                    strBuffers(nAqueousSpecies) = trim(strBuffers(nAqueousSpecies)) // " ; " // trim(strBuffer)
                else if (iTypeOfData(strbuffer, strSectionEnd) == -1) then 
                    exit
                end if
                
            end do
            
            if (nAqueousSpecies > 0) then
                if(allocated(aqueousSpecies)) then
                    deallocate(aqueousSpecies)
                end if
                allocate(aqueousSpecies(nAqueousSpecies))
                do i = 1, nAqueousSpecies 
                    
                    call WriteLog("Aqueous number:", i)
                    call WriteLog(strBuffers(i))
                    
                    call analyseReactionEquation(strBuffers(i), nNameLength, .true., aqueousSpecies(i)%Name, aqueousSpecies(i)%NCP, aqueousSpecies(i)%NameOfSTQ, aqueousSpecies(i)%STQ)
                    aqueousSpecies(i)%iAssociation = 1
                    
                    call getDH_AB(strBuffers(i), aqueousSpecies(i)%DHA, aqueousSpecies(i)%DHB)                    
                    call getLogK(strBuffers(i), aqueousSpecies(i)%AKLOG(1))
                    call getEnthalpyChange(strBuffers(i), aqueousSpecies(i)%EnthalpyChange)
                    call estimateChargeFromName (aqueousSpecies(i)%Name, aqueousSpecies(i)%Z)
                    
                    !calculate molecular weight                  
                    aqueousSpecies(i)%MWT = calculateMWTFromComponents(aqueousSpecies(i)%NCP, aqueousSpecies(i)%STQ, aqueousSpecies(i)%NameOfSTQ)
                    
                    !rTempMWT = getMWTFromMasterSpecies(trim(aqueousSpecies(i)%Name))
                    !if(rTempMWT < 0) then
                    rTempMWT = getMolecularMass(trim(aqueousSpecies(i)%Name))
                    !end if
                    
                    if(abs(aqueousSpecies(i)%MWT - rTempMWT) > 1.0d-3) then
                        nWarnings = nWarnings + 1                        
                        call WriteLog("Warning: Check the calculated molecular weight for : "// trim(aqueousSpecies(i)%Name))
                    end if                                        
                    
                    call WriteLog("master specie name: " // trim(adjustl(aqueousSpecies(i)%Name)))
                    call WriteLog("Debye-Hukel A: ", aqueousSpecies(i)%DHA)
                    call WriteLog("Debye-Hukel B: ", aqueousSpecies(i)%DHB)                    
                    call WriteLog("LogK: ", aqueousSpecies(i)%AKLOG(1))
                    call WriteLog("Enthalpy change: ", aqueousSpecies(i)%EnthalpyChange)
                    call WriteLog("Molecular weight: ", aqueousSpecies(i)%MWT)
                    
                    
                    do j = 1, aqueousSpecies(i)%NCP
                        call WriteLog(aqueousSpecies(i)%NameOfSTQ(j), aqueousSpecies(i)%STQ(j))
                    end do
                    
                end do
            end if
        
        end if

        call WriteLog("Number of aqueous species: ")
        call WriteLog(nAqueousSpecies)
        call WriteLog("Number of redox reactions:")
        call WriteLog(nRedoxReaction)
        call WriteLog("Read aqueous species success")
        
	    return

9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling

     end subroutine readAqueousSpecies
     

     ! Read gases
     subroutine readGases
        
        implicit none
        
        integer :: i, j, k, n
        character(nNameLength) :: strName
        logical :: bFindGases
        character(1024) :: strBufferTemp

        i = 0
        n = 0
        nGases = 0
        strBufferTemp = ""
        
        ! read SOLUTION_SPECIES data
        bEndOfFile = .false.
        rewind iUnitDbsPhreeqc        
        bFindGases = .false.
        do while (.not. bEndOfFile)
            if(.not. bReadNextLine()) then
                exit
            end if
            if (trim(adjustl(strbuffer)) == strSectionGases) then
                bFindGases = .true.
                exit
            end if
        end do
        
        if (bFindGases) then
            
            do while (.not. bEndOfFile)
                
                if (nGases == maxStrBuffers) then
                    nErrors = nErrors + 1
                    call WriteLog("Error: The maximum database of Phreeqc is 100,000 lines for each section." // & 
                                    "Please modify the source code parameter maxStrBuffers to read more data.")
                    call ErrorHandling
                end if
                
                if(.not. bReadNextLine()) then
                    exit
                end if               
            
                if (iTypeOfData(strbuffer, strSectionEnd) == 1 .or. iTypeOfData(strbuffer, strSectionEnd) == 2) then
                    nGases = nGases + 1
                    if(len_trim(strBufferTemp) > 0) then
                        strBuffers(nGases) = trim(strBufferTemp) // " ; " //  trim(strBuffer)
                    else
                        strBuffers(nGases) = trim(strBuffer)
                    end if
                    strBuffersLineIndex(nGases) = currentLine
              
                else if (iTypeOfData(strbuffer, strSectionEnd) == -1 ) then 
                    exit
                else
                    if(iTypeKeyWord(strbuffer, strSectionEnd) > 0) then
                        if(nGases> 0) then
                            strBuffers(nGases) = trim(strBuffers(nGases)) // " ; " // trim(strBuffer)
                        end if
                        strBufferTemp = ""
                    else
                        strBufferTemp = trim(strBuffer)
                    end if
                end if
                
            end do
            
            if (nGases > 0) then
                if(allocated(gases)) then
                    deallocate(gases)
                end if
                allocate(gases(nGases))
                do i = 1, nGases 
                    
                    call WriteLog("Gases number:", i)
                    call WriteLog(strBuffers(i))
                    
                    call analyseReactionEquation(strBuffers(i), nNameLength, .false., gases(i)%Name, gases(i)%NCP, gases(i)%NameOfSTQ, gases(i)%STQ)
                    gases(i)%iAssociation = -1
                    
                    call getLogK(strBuffers(i), gases(i)%AKLOG(1))
                    
                    call getEnthalpyChange(strBuffers(i), gases(i)%EnthalpyChange)
                    
                    !calculate molecular weight                  
                    gases(i)%MWT = calculateMWTFromComponents(gases(i)%NCP, gases(i)%STQ, gases(i)%NameOfSTQ)                      
                    
                    call WriteLog("gases name: " // trim(adjustl(gases(i)%Name)))
                    call WriteLog("LogK: ", gases(i)%AKLOG(1))
                    call WriteLog("Enthalpy change: ", gases(i)%EnthalpyChange)
                    call WriteLog("Molecular weight: ", gases(i)%MWT)
                    
                    do j = 1, gases(i)%NCP
                        call WriteLog(gases(i)%NameOfSTQ(j), gases(i)%STQ(j))
                    end do
                    
                end do
            end if
            
        end if

        call WriteLog("Number of gases: ")
        call WriteLog(nGases)
        call WriteLog("Read gases success")
        
        return
        
9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling

     end subroutine readGases
     

     ! Read minerals
     subroutine readMinerals
        
        implicit none
        
        integer :: i, j, k, n
        character(nNameLength) :: strName
        logical :: bFindMinerals
        character(1024) :: strBufferTemp

        i = 0
        n = 0
        nMinerals = 0
        strBufferTemp = ""
        
        ! read SOLUTION_SPECIES data
        bEndOfFile = .false.
        rewind iUnitDbsPhreeqc        
        bFindMinerals = .false.
        do while (.not. bEndOfFile)
            if(.not. bReadNextLine()) then
                exit
            end if
            if (trim(adjustl(strbuffer)) == strSectionMinerals) then
                bFindMinerals = .true.
                exit
            end if
        end do
        
        if (bFindMinerals) then
            
            do while (.not. bEndOfFile)
                
                if (nMinerals == maxStrBuffers) then
                    nErrors = nErrors + 1
                    call WriteLog("Error: The maximum database of Phreeqc is 100,000 lines for each section." // & 
                                    "Please modify the source code parameter maxStrBuffers to read more data.")
                    call ErrorHandling
                end if
                
                if(.not. bReadNextLine()) then
                    exit
                end if               
            
                if (iTypeOfData(strbuffer, strSectionEnd) == 1 .or. iTypeOfData(strbuffer, strSectionEnd) == 2) then
                    nMinerals = nMinerals + 1
                    if(len_trim(strBufferTemp) > 0) then
                        strBuffers(nMinerals) = trim(strBufferTemp) // " ; " //  trim(strBuffer)
                    else
                        strBuffers(nMinerals) = trim(strBuffer)
                    end if
                    strBuffersLineIndex(nMinerals) = currentLine
              
                else if (iTypeOfData(strbuffer, strSectionEnd) == -1 ) then 
                    exit
                else
                    if(iTypeKeyWord(strbuffer, strSectionEnd) > 0) then
                        if(nMinerals> 0) then
                            strBuffers(nMinerals) = trim(strBuffers(nMinerals)) // " ; " // trim(strBuffer)
                        end if
                        strBufferTemp = ""
                    else
                        strBufferTemp = trim(strBuffer)
                    end if
                end if
                
            end do
            
            if (nMinerals > 0) then
                if(allocated(minerals)) then
                    deallocate(minerals)
                end if
                allocate(minerals(nMinerals))
                do i = 1, nMinerals 
                    
                    call WriteLog("Minerals number:", i)
                    call WriteLog(strBuffers(i))
                    
                    call analyseReactionEquation(strBuffers(i), nNameLength, .false., minerals(i)%Name, minerals(i)%NCP, minerals(i)%NameOfSTQ, minerals(i)%STQ)
                    minerals(i)%iAssociation = -1
                    
                    call getLogK(strBuffers(i), minerals(i)%AKLOG(1))
                    
                    call getEnthalpyChange(strBuffers(i), minerals(i)%EnthalpyChange)
                    
                    !calculate molecular weight                  
                    minerals(i)%MWT = calculateMWTFromComponents(minerals(i)%NCP, minerals(i)%STQ, minerals(i)%NameOfSTQ)  
                    
                    call WriteLog("minerals name: " // trim(adjustl(minerals(i)%Name)))
                    call WriteLog("LogK: ", minerals(i)%AKLOG(1))
                    call WriteLog("Enthalpy change: ", minerals(i)%EnthalpyChange)
                    call WriteLog("Molecular weight: ", minerals(i)%MWT)
                    
                    do j = 1, minerals(i)%NCP
                        call WriteLog(minerals(i)%NameOfSTQ(j), minerals(i)%STQ(j))
                    end do
                    
                end do
            end if
            
        end if

        call WriteLog("Number of minerals: ")
        call WriteLog(nMinerals)
        call WriteLog("Read minerals success")
        
        return
        
9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling

     end subroutine readMinerals  
     
     !SetLowerCase all the species and phases
     subroutine SetLowerCaseAll
     
        implicit none
        
        integer :: i , j 
        
        do i = 1, nSpecies
            call SetLowerCase(species(i)%Name)
        end do
        
        do i = 1, nAqueousSpecies
            call SetLowerCase(aqueousSpecies(i)%Name)
            do j = 1, aqueousSpecies(i)%NCP
                call SetLowerCase(aqueousSpecies(i)%NameOfSTQ(j))
            end do
        end do
        
        do i = 1, nGases
            call SetLowerCase(gases(i)%Name)
            do j = 1, gases(i)%NCP
                call SetLowerCase(gases(i)%NameOfSTQ(j))
            end do
        end do  
        
        do i = 1, nMinerals
            call SetLowerCase(minerals(i)%Name)
            do j = 1, minerals(i)%NCP
                call SetLowerCase(minerals(i)%NameOfSTQ(j))
            end do
        end do 
     
     end subroutine
    
end module