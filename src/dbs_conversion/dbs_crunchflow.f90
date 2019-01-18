!****************************************************************************
!    Geochemitry database conversion program
!    Author: Danyang Su
!    Email: dsu@eos.ubc.ca; danyang.su@gmail.com
!****************************************************************************
!Module to read and store CRUNCHFLOW database
    
module dbs_crunchflow 

    use global, only : iUnitDbsCF, bOpenDbsCF, ErrorHandling
    
    use logfile, only : WriteLog, nWarnings, nErrors

    use file_utility, only : SetLowerCase, IsZero, getNameFromString, skipNValues, replaceCharacter
    
    use alias, only : GetNameFromAlias

    use inputfile, only : sourceDatabasePath

    use geochemistry,    only :     aklog_o2aq, aklog_h2aq

    use sourcedata
    
    implicit none    
        
    character(1024)      ::      filePathDbsCF   !file path for CRUNCHFLOW database
    
    character(1),   parameter   ::  strComment      =   "#"
    character(1),   parameter   ::  strData         =   "'"

    character(14),  parameter   ::  strEndPrimary   =   "end of primary"
    character(16),  parameter   ::  strEndSecondary =   "end of secondary"
    character(12),  parameter   ::  strEndGases     =   "end of gases"
    character(15),  parameter   ::  strEndMinerals  =   "end of minerals"
    character(26),  parameter   ::  strBeginSurfComplex = "begin surface complexation"
    character(27),  parameter   ::  strEndSurfComplex   = "end of surface complexation"
    character(18),  parameter   ::  strTemperature  =   "temperature points"
    character(6),   parameter   ::  strRedoxReaction =  "o2(aq)"
    character(18),  parameter   ::  strDebyeHuckelAdh = "'debye-huckel adh'"
    character(18),  parameter   ::  strDebyeHuckelBdh = "'debye-huckel bdh'"
    character(18),  parameter   ::  strDebyeHuckelBdt = "'debye-huckel bdt'"
    

    integer, parameter          ::  maxStrBuffers = 100000
    character(1024)              ::  strBuffers(maxStrBuffers)       !Assume the maximum number of each data section 
    character(1024)              ::  strBuffer        =   ""
    integer                     ::  strBuffersLineIndex(maxStrBuffers)
    integer                      ::  iReadStat       =   1
    logical                      ::  bEndOfFile      =   .false.   
    
    integer :: currentLine = 0          

    
contains

    ! Open CRUNCHFLOW database
    subroutine OpenDbsCF
    
        implicit none
        
        integer :: stat 

        filePathDbsCF = trim(sourceDatabasePath)
        
        open(iUnitDbsCF, file = filePathDbsCF,  iostat = stat)
        
        if (stat /= 0) then
            call WriteLog("Error in opening CRUNCHFLOW databse file: " // trim(filePathDbsCF))
            call ErrorHandling
        else
            bOpenDbsCF = .true.
            call WriteLog("Open CRUNCHFLOW database success: " // trim(filePathDbsCF))
        end if 
    
    end subroutine OpenDbsCF
    
    
    ! Read CRUNCHFLOW database
     subroutine ReadDbsCF
     
        implicit none  
        
        call WriteLog ("Begin reading CRUNCHFLOW database")        
      
        call readTemperature 
        
        call readDebyeHuckelParams

        call readSpecies
        
        call readAqueousSpecies
        
        call readGases
        
        call readMinerals
        
        call readSurfaceComplexes

        call WriteLog ("End reading CRUNCHFLOW database")
     
     end subroutine ReadDbsCF
     
     ! Close CRUNCHFLOW database
     subroutine CloseDbsCF
     
        implicit none
        
        if (bOpenDbsCF) then
            
            close(iUnitDbsCF)
            
        end if
     
     end subroutine CloseDbsCF
     
    
    ! Read next line
    subroutine readNextLine
    
        implicit none
        
        integer :: i

        do while(.not. bEndOfFile)

            read(iUnitDbsCF, "(a)", iostat = iReadStat) strBuffer
            
            currentLine = currentLine + 1
       
            if (iReadStat > 0) then     !Error in reading
            
                call WriteLog("Error in reading CRUNCHFLOW databse file.")
            
                call ErrorHandling
            
            else if (iReadStat < 0) then    !End of file
            
                bEndOfFile = .true.
                
                exit            
            
            end if

            strBuffer = adjustl(strBuffer)

            if(len_trim(strBuffer) /= 0 ) then
            
                ! Conver all case to lower case
                call SetLowerCase(strBuffer)

                if (strBuffer(1:1) /= strComment) then
                    call replaceCharacter(strBuffer, achar(9), " ")
                    exit
                end if

            end if

        end do
    
    endsubroutine readNextLine
    
   
    ! Read Debye-Huckel Parameters
    subroutine readDebyeHuckelParams

        implicit none

        logical :: bFlagAdh = .false.
        logical :: bFlagBdh = .false.
        logical :: bFlagBdt = .false.

        do while (.not. bEndOfFile)
            call readNextLine
            if (index(strBuffer,strDebyeHuckelAdh) == 1) then
                bFlagAdh = .true.
            else if (index(strBuffer,strDebyeHuckelBdh) == 1) then
                bFlagBdh = .true.
            else if (index(strBuffer,strDebyeHuckelBdt) == 1) then
                bFlagBdt = .true.
            end if              
            
            if(bFlagAdh .and. bFlagBdh .and. bFlagBdt) then
                exit
            end if
        end do

        if (bFlagAdh == .false. .or. bFlagBdh == .false. .or. bFlagBdt == .false.) then
            call WriteLog("Error in reading Debye-Huckel Parameters adh, bdh, or bdt.")
            call ErrorHandling
        else
            call WriteLog("Read Debye-Huckel Parameters success.")
        end if

    end subroutine readDebyeHuckelParams
   
     ! Read temperature points
     subroutine readTemperature

        implicit none
        
        logical :: bFlag = .true.
       
        character(1024) :: strTemp

        do while (.not. bEndOfFile)
            call readNextLine
            if (index(strBuffer,strTemperature) > 0) then 
              
                call skipNValues(strBuffer,1)
                
                read(strBuffer, *, end = 9999, err = 9999)  nTemperature
                
                write(*,*) nTemperature

                if (allocated(temperature)) then
                    deallocate(temperature)
                end if
            
                allocate(temperature(nTemperature))
            
                read(strBuffer, *, end = 9999, err = 9999)  nTemperature, temperature(1:nTemperature)
            
                bFlag = .false.
                exit
                
            end if
        end do

        if (bFlag) then
            call WriteLog("Error in reading temperature points, label not found or data error:" // strTemperature)
            call ErrorHandling
        else
            call WriteLog("Read temperature points success.")
            call WriteLog("Number of termperature points:")
            call WriteLog(nTemperature)
            call WriteLog("Temperatures:")
            call WriteLog(temperature)
        end if
        
	    return
        
9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(currentLine)
        call ErrorHandling

     end subroutine readTemperature

     ! Read species
     subroutine readSpecies

        implicit none

        integer :: i

        character (nNameLength) :: strName

        nSpecies = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, nNameLength, strName)
            if (trim(strName) == strEndPrimary) then                
                exit
            end if
            if (nSpecies == maxStrBuffers) then
                nErrors = nErrors + 1
                call WriteLog("Error: The maximum database of CRUNCHFLOW is 100,000 lines for each section." // & 
                                "Please modify the source code parameter maxStrBuffers to read more data.")
                call ErrorHandling
            end if
            nSpecies = nSpecies + 1
            strBuffers(nSpecies) = strBuffer
            strBuffersLineIndex(nSpecies) = currentLine
        end do

        if (nSpecies > 0) then
            if(allocated(species)) then
                deallocate(species)
            end if
            allocate(species(nSpecies))
            do i = 1, nSpecies
                call getNameFromString(strBuffers(i), nNameLength, species(i)%Name)
                !Check if this name is used as an alias
                strName = GetNameFromAlias(species(i)%Name)
                if (trim(strName) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(species(i)%Name))
                    call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(species(i)%Name) // " are the same")
                end if
                call skipNValues(strBuffers(i), 1)
                read(strBuffers(i),* , end = 9999, err = 9999) species(i)%A0, species(i)%Z, species(i)%MWT
            end do
        end if

        call WriteLog("Number of Species: ")
        call WriteLog(nSpecies)
        call WriteLog("Read species success")
        
	    return

9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(i))
        call ErrorHandling    

     end subroutine readSpecies    

     
!     !Allocate STQ and NAM for aqueous species
!     subroutine allocateOneAqueousSpecies(i,n, ntmp)
!        
!        implicit none
!        integer, intent(in) :: i, n, ntmp
!        
!        if(allocated(aqueousSpecies(i)%STQ)) then            
!            deallocate(aqueousSpecies(i)%STQ)
!        end if
!        allocate(aqueousSpecies(i)%STQ(n))
!        
!        if(allocated(aqueousSpecies(i)%NameOfSTQ)) then            
!            deallocate(aqueousSpecies(i)%NameOfSTQ)
!        end if
!        allocate(aqueousSpecies(i)%NameOfSTQ(n)) 
!        
!        if(allocated(aqueousSpecies(i)%AKLOG)) then            
!            deallocate(aqueousSpecies(i)%AKLOG)
!        end if
!        allocate(aqueousSpecies(i)%AKLOG(ntmp)) 
!            
!     end subroutine allocateOneAqueousSpecies

     ! Read aqueous species
     subroutine readAqueousSpecies
        
        implicit none
        
        integer :: i, j, k
        character(nNameLength) :: strName
        character(nNameLength) :: tempNames(3)

        i = 0
        nAqueousSpecies = 0
        nRedoxReaction = 0
        
        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, nNameLength, strName)
            if (trim(strName) == strEndSecondary) then                
                exit
            end if
            nAqueousSpecies = nAqueousSpecies + 1
            strBuffers(nAqueousSpecies) = strBuffer
            strBuffersLineIndex(nAqueousSpecies) = currentLine
        end do
       
        if (nAqueousSpecies > 0) then
            if(allocated(aqueousSpecies)) then
                deallocate(aqueousSpecies)
            end if
            allocate(aqueousSpecies(nAqueousSpecies))
            
            aqueousSpecies(1:nAqueousSpecies)%iAssociation = -1
            
            do i = 1, nAqueousSpecies                
                !read components
                j = i
                call getNameFromString(strBuffers(j), nNameLength, aqueousSpecies(i)%Name)
                !Check if this name is used as an alias
                strName = GetNameFromAlias(aqueousSpecies(i)%Name)
                if (trim(strName) /= "") then
                    nErrors = nErrors + 1                    
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(aqueousSpecies(i)%Name))
                    call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(aqueousSpecies(i)%Name) // " are the same")
                end if
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*, end = 9999, err = 9999) aqueousSpecies(i)%NCP
                call skipNValues(strBuffers(j),1) 
                do k = 1, aqueousSpecies(i)%NCP, 1
                    read(strBuffers(j),*, end = 9999, err = 9999) aqueousSpecies(i)%STQ(k)
                    call getNameFromString(strBuffers(j), nNameLength, aqueousSpecies(i)%NameOfSTQ(k))
                    !Check if this name is used as an alias                  
                    strName = GetNameFromAlias(aqueousSpecies(i)%NameOfSTQ(k))
                    if (trim(strName) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(aqueousSpecies(i)%NameOfSTQ(k)))
                        call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(aqueousSpecies(i)%NameOfSTQ(k)) // " are the same")
                    end if
                    
                    call skipNValues(strBuffers(j),2)
                    !write(*,*) trim(aqueousSpecies(i)%NameOfSTQ(m)), aqueousSpecies(i)%STQ(m)
                end do               
                
                !read equilibrium constants, molecular weight
                read(strBuffers(j),*, end = 9999, err = 9999) aqueousSpecies(i)%AKLOG(1:nTemperature), &
                    aqueousSpecies(i)%A0, aqueousSpecies(i)%Z, aqueousSpecies(i)%MWT
   

                if(trim(aqueousSpecies(i)%Name) == "h2(aq)") then
                    aklog_h2aq = 0
                    do k = 1, nTemperature
                        aklog_h2aq(k) = aqueousSpecies(i)%AKLOG(k)
                        aklog_o2aq(k) = 2.0 * aklog_h2aq(k)
                    end do
                    !aklog_h2aq(1:nTemperature) = aqueousSpecies(i)%AKLOG(1:nTemperature)
                    call WriteLog ("h2(aq) AKLOG")
                    call WriteLog (nTemperature, aklog_h2aq(:nTemperature))
                else if(trim(aqueousSpecies(i)%Name) == "o2(aq)") then
                    do k = 1, nTemperature
                        aklog_o2aq(k) = aqueousSpecies(i)%AKLOG(k)
                        aklog_h2aq(k) = 0.5 * aklog_o2aq(k)
                    end do
                    !aklog_o2aq(1:nTemperature) = aqueousSpecies(i)%AKLOG(1:nTemperature)
                    call WriteLog ("o2(aq) AKLOG")
                    call WriteLog (nTemperature, aklog_o2aq(:nTemperature))
                end if

            end do
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
     
!     !Allocate STQ and NAM for mineral
!     subroutine allocateOneMineral(i,n, ntmp)
!        
!        implicit none
!        integer, intent(in) :: i, n, ntmp
!        
!        if(allocated(minerals(i)%STQ)) then            
!            deallocate(minerals(i)%STQ)
!        end if
!        allocate(minerals(i)%STQ(n))
!        
!        if(allocated(minerals(i)%NameOfSTQ)) then            
!            deallocate(minerals(i)%NameOfSTQ)
!        end if
!        allocate(minerals(i)%NameOfSTQ(n)) 
!        
!        if(allocated(minerals(i)%AKLOG)) then            
!            deallocate(minerals(i)%AKLOG)
!        end if
!        allocate(minerals(i)%AKLOG(ntmp)) 
!            
!     end subroutine allocateOneMineral
     
     ! Read minerals
     subroutine readMinerals
        
        implicit none
        
        integer :: i, j, k
        character(nNameLength) :: strName
        character(nNameLength) :: tempNames(3)

        i = 0
        nMinerals = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, nNameLength, strName)
            if (trim(strName) == strEndMinerals) then                
                exit
            end if
            nMinerals = nMinerals + 1
            strBuffers(nMinerals) = strBuffer
            strBuffersLineIndex(nMinerals) = currentLine
        end do
        
        if (nMinerals > 0) then
            if(allocated(minerals)) then
                deallocate(minerals)
            end if
            allocate(minerals(nMinerals))
            
            minerals(1:nMinerals)%iAssociation = -1
            
            do i = 1, nMinerals                
                !read first line
                j = i
                call getNameFromString(strBuffers(j), nNameLength,minerals(i)%Name)
                
                !Check if this name is used as an alias
                strName = GetNameFromAlias(minerals(i)%Name)
                if (trim(strName) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(minerals(i)%Name))
                    call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(minerals(i)%Name) // " are the same")
                end if                
                
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*, end = 9999, err = 9999) minerals(i)%VMIN, minerals(i)%NCP
                
                call skipNValues(strBuffers(j),2)  
                
                do k = 1, minerals(i)%NCP, 1
                    read(strBuffers(j),*, end = 9999, err = 9999) minerals(i)%STQ(k)
                    call getNameFromString(strBuffers(j), nNameLength, minerals(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    strName = GetNameFromAlias(minerals(i)%NameOfSTQ(k))
                    if (trim(strName) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(minerals(i)%NameOfSTQ(k)))
                        call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(minerals(i)%NameOfSTQ(k)) // " are the same")
                    end if
                    
                    call skipNValues(strBuffers(j),2)
                    !write(*,*) trim(minerals(i)%NameOfSTQ(k)), minerals(i)%STQ(k)
                end do

                read(strBuffers(j),*, end = 9999, err = 9999) minerals(i)%AKLOG(1:nTemperature), minerals(i)%MWT
                
            end do
        end if

        call WriteLog("Number of minerals: ")
        call WriteLog(nMinerals)
        call WriteLog("Read minerals success")
        
        return
        
9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling

     end subroutine readMinerals
     
!      !Allocate STQ and NAM for gas
!     subroutine allocateOneGas(i,n, ntmp)
!        
!        implicit none
!        integer, intent(in) :: i, n, ntmp
!        
!        if(allocated(gases(i)%STQ)) then            
!            deallocate(gases(i)%STQ)
!        end if
!        allocate(gases(i)%STQ(n))
!        
!        if(allocated(gases(i)%NameOfSTQ)) then            
!            deallocate(gases(i)%NameOfSTQ)
!        end if
!        allocate(gases(i)%NameOfSTQ(n)) 
!        
!        if(allocated(gases(i)%AKLOG)) then            
!            deallocate(gases(i)%AKLOG)
!        end if
!        allocate(gases(i)%AKLOG(ntmp)) 
!            
!     end subroutine allocateOneGas
     
     ! Read gases
     subroutine readGases
     
        implicit none
        
        integer :: i, j, k
        character(nNameLength) :: strName
        character(nNameLength) :: tempNames(3)

        i = 0
        nGases = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, nNameLength, strName)
            if (trim(strName) == strEndGases) then                
                exit
            end if
            nGases = nGases + 1
            strBuffers(nGases) = strBuffer
            strBuffersLineIndex(nGases) = currentLine
        end do
        
        if (nGases > 0) then
            if(allocated(gases)) then
                deallocate(gases)
            end if
            allocate(gases(nGases))
            
            gases(1:nGases)%iAssociation = -1
            
            do i = 1, nGases                
                !read first line
                j = i
                call getNameFromString(strBuffers(j), nNameLength, gases(i)%Name)
                
                !Check if this name is used as an alias
                strName = GetNameFromAlias(gases(i)%Name)
                if (trim(strName) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(gases(i)%Name))
                    call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(gases(i)%Name) // " are the same")
                end if
                
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*, end = 9999, err = 9999) gases(i)%VMIN, gases(i)%NCP                
        
                call skipNValues(strBuffers(j),2)  
                
                do k = 1, gases(i)%NCP, 1
                    read(strBuffers(j),*, end = 9999, err = 9999) gases(i)%STQ(k)
                    call getNameFromString(strBuffers(j), nNameLength, gases(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    strName = GetNameFromAlias(gases(i)%NameOfSTQ(k))
                    if (trim(strName) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(gases(i)%NameOfSTQ(k)))
                        call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(gases(i)%NameOfSTQ(k)) // " are the same")
                    end if
                    
                    call skipNValues(strBuffers(j),2)
                    !write(*,*) trim(gases(i)%NameOfSTQ(k)), gases(i)%STQ(k)
                end do
                read(strBuffers(j),*, end = 9999, err = 9999) gases(i)%AKLOG(1:nTemperature), gases(i)%MWT                
            end do
        end if

        call WriteLog("Number of gases: ")
        call WriteLog(nGases)
        call WriteLog("Read gases success")
 
		return

9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling
        
     end subroutine readGases
     
!      !Allocate STQ and NAM for surface complexes
!     subroutine allocateOneSurfaceComplexes(i,n, ntmp)
!        
!        implicit none
!        integer, intent(in) :: i, n, ntmp
!        
!        if(allocated(surfaceComplexes(i)%STQ)) then            
!            deallocate(surfaceComplexes(i)%STQ)
!        end if
!        allocate(surfaceComplexes(i)%STQ(n))
!        
!        if(allocated(surfaceComplexes(i)%NameOfSTQ)) then            
!            deallocate(surfaceComplexes(i)%NameOfSTQ)
!        end if
!        allocate(surfaceComplexes(i)%NameOfSTQ(n)) 
!        
!        if(allocated(surfaceComplexes(i)%AKLOG)) then            
!            deallocate(surfaceComplexes(i)%AKLOG)
!        end if
!        allocate(surfaceComplexes(i)%AKLOG(ntmp)) 
!            
!     end subroutine allocateOneSurfaceComplexes
     
     ! Read surface complexes
     subroutine readSurfaceComplexes
     
        implicit none
        
        integer :: i, j, k
        character(nNameLength) :: strName
        character(nNameLength) :: tempNames(3)
        logical :: bBeginSurfComplex = .false.

        i = 0
        nSurfaceComplexes = 0

        do while (.not. bEndOfFile)
            call readNextLine
            
            if(.not. bBeginSurfComplex) then
                if (index(strBuffer, strBeginSurfComplex) == 1) then                
                    bBeginSurfComplex = .true.
                    cycle
                end if
            end if
           
            if (index(strBuffer, strEndSurfComplex) == 1) then                
                exit
            end if
            
            if (bBeginSurfComplex) then
                nSurfaceComplexes = nSurfaceComplexes + 1
                strBuffers(nSurfaceComplexes) = strBuffer
                strBuffersLineIndex(nSurfaceComplexes) = currentLine   
            end if
        end do
               
        if (nSurfaceComplexes > 0) then
            if(allocated(surfaceComplexes)) then
                deallocate(surfaceComplexes)
            end if
            allocate(surfaceComplexes(nSurfaceComplexes))
            
            surfaceComplexes(1:nSurfaceComplexes)%iAssociation = -1
            
            do i = 1, nSurfaceComplexes                
                !read first line
                j = i
                call getNameFromString(strBuffers(j), nNameLength,surfaceComplexes(i)%Name)
                
                !Check if this name is used as an alias
                strName = GetNameFromAlias(surfaceComplexes(i)%Name)
                if (trim(strName) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(surfaceComplexes(i)%Name))
                    call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(surfaceComplexes(i)%Name) // " are the same")
                end if
                
                
                call skipNValues(strBuffers(j),1)
                
                read(strBuffers(j),*, end = 9999, err = 9999) surfaceComplexes(i)%NCP

                call skipNValues(strBuffers(j),1) 
                
                do k = 1, surfaceComplexes(i)%NCP, 1
                    read(strBuffers(j),*, end = 9999, err = 9999) surfaceComplexes(i)%STQ(k)
                    call getNameFromString(strBuffers(j), nNameLength, surfaceComplexes(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    strName = GetNameFromAlias(surfaceComplexes(i)%NameOfSTQ(k))
                    if (trim(strName) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "// trim(strName) // "  ->  " //trim(surfaceComplexes(i)%NameOfSTQ(k)))
                        call WriteLog("Ignore this error if "// trim(strName) // " and " //trim(surfaceComplexes(i)%NameOfSTQ(k)) // " are the same")
                    end if                     
                    
                    call skipNValues(strBuffers(j),2)
                    !write(*,*) trim(surfaceComplexes(i)%NameOfSTQ(k)), surfaceComplexes(i)%STQ(k)
                end do

                read(strBuffers(j),*, end = 9999, err = 9999) surfaceComplexes(i)%AKLOG(1:nTemperature)
            end do
        end if

        call WriteLog("Number of surface complexes: ")
        call WriteLog(nSurfaceComplexes)
        call WriteLog("Read surface complexes success")
        
        return

9999    call WriteLog("Error detected in reading/converting data in the follow line number")
        call WriteLog(strBuffersLineIndex(j))
        call ErrorHandling
        
     end subroutine readSurfaceComplexes
    
end module