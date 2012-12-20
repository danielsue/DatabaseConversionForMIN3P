!Module to read and store toughreact database
    
module dbs_toughreact 

    use global, only : iUnitDbsTR, bOpenDbsTR, ErrorHandling
    
    use logfile, only : WriteLog, nWarnings, nErrors

    use file_utility, only : LowerCase
    
    use alias, only : GetNameFromAlias
    
    implicit none
    
    integer, parameter:: nMNLTR = 60
    
    !
    type TypePrimarySpecies
        
        character(nMNLTR)   ::      Name    !name or chemical formula of aqueous basis species, in quotes (truncated after 20 characters)
        
        real            ::      A0      !Ion effective or hydrated radius used to compute the Debye-Hückel a0 parameter
        
        real            ::      Z       !the ion electric charge
        
        real            ::      MWT     !Molecular weight of the aqueous species (g/mol)
        
    end type TypePrimarySpecies   
     
    !
    type TypeDerivedSpecies
        
        character(nMNLTR)   ::      Name    !chemical formula of secondary species, in quotes (truncated after 20 characters).
        
        real            ::      A0      !Ion effective or hydrated radius used to compute the Debye-Hückel a0 parameter
        
        real            ::      Z       !the ion electric charge
        
        real            ::      MWT     !Molecular weight of the aqueous species (g/mol)
        
        integer         ::      NCP     !number of basis species defining the secondary species
        
        real,allocatable          ::  STQ(:)          !stoichiometric coefficients of component (basis) species NAM included in the dissociation
                                                        !reaction of the derived species (negative and positive values for reactants and products,
                                                        !respectively). The derived species is always assumed to have a stoichiometric coefficient of -1.0,
                                                        !which is not included in STQ.
                                                        
        
        character(nMNLTR),allocatable ::  NameOfSTQ(:)    !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).
                                                         
        real,allocatable          ::  AKLOG(:)         !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)         !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)         !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ΔV
                                                         !(in cm3/mol) for the reaction as a function of temperature (average ΔV over the pressure interval
                                                         !P0 to P), with ΔV = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists

        logical                     ::  bRedoxReaction = .false.   !To seperate if it is complexation reaction or redox reaction
                                                                     !If NameOfSTQ contain "o2(aq)", it is redox reaction
                                                                     !...........IS IT CORRECT?.....................
            
    end type TypeDerivedSpecies   
                                             
    
    !
    type TypeMineral
        
        character(nMNLTR)   ::      Name    !name or chemical formula of a mineral, in quotes (truncated after 20 characters).
        
        real            ::      MWT     !molecular weight (g/mol)
        
        real            ::      VMIN    !molar volume (cm3/mole)
        
        integer         ::      NCP    !the number of component species defining the mineral.
        
        real,allocatable            ::  STQ(:)      !contains the stoichiometric coefficient of basis species NAM in the dissociation (hydrolysis)
                                                    !reaction of the mineral (negative and positive values for reactants and products, respectively). The
                                                    !mineral species is always assumed to have a stoichiometric coefficient of -1.0, which is not
                                                    !included in STQ.
                                                    
        character(nMNLTR),allocatable   ::  NameOfSTQ(:)  !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).                                              
                                                    
        real,allocatable            ::  AKLOG(:)         !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)        !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)         !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ΔV
                                                         !(in cm3/mol) for the reaction as a function of temperature (average ΔV over the pressure interval
                                                         !P0 to P), with ΔV = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists  
    
    end type TypeMineral
    
    !
    type TypeGas
        
        character(nMNLTR)   ::      Name    !name or chemical formula of a gas species, in quotes (truncated after 20 characters).
        
        real            ::      MWT    !molecular weight (g/mol)
        
        real            ::      DMDIAM  !molecular diameter (m) used to calculate gas diffusion coefficient
        
        integer         ::      NCP    !the number of basis species defining the gas.
        
        real,allocatable            ::  STQ(:)     !contains the stoichiometric coefficient of component species NAM in the dissociation reaction
                                                    !of the gas (negative and positive values for reactants and products, respectively). The gas is
                                                    !always assumed to have a stoichiometric coefficient of -1.0, which is not included in STQ.
                                                    
        character(nMNLTR),allocatable   ::  NameOfSTQ(:)    !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).  
                                                         
        real,allocatable            ::  AKLOG(:)         !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)         !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)         !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ΔV
                                                         !(in cm3/mol) for the reaction as a function of temperature (average ΔV over the pressure interval
                                                         !P0 to P), with ΔV = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists                                                
        
    end type TypeGas
    
    !
    type TypeSurfaceComplexes
        
        character(nMNLTR)   ::      Name    !name or chemical formula of a gas species, in quotes (truncated after 20 characters).
        
        real            ::      Z       !electric charge of surface complex.
        
        integer         ::      NCP     !the number of basis species defining the surface complex.
        
        real,allocatable            ::  STQ(:)      !contains the stoichiometric coefficient of component species NAM in the dissociation reaction
                                                    !of the gas (negative and positive values for reactants and products, respectively). The gas is
                                                    !always assumed to have a stoichiometric coefficient of -1.0, which is not included in STQ.
                                                    
        character(nMNLTR),allocatable   ::  NameOfSTQ(:)  !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).  
                                                         
        real,allocatable            ::  AKLOG(:)         !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)         !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.                                             
        
    end type TypeSurfaceComplexes
    
    
    character(260)      ::      filePathDbsTR   !file path for toughreact database
    
    character(1),   parameter   ::  strComment      =   "#"
    character(1),   parameter   ::  strData         =   "'"
    character(14),  parameter   ::  strEndOfHeader  =   "!end-of-header"
    character(4),   parameter   ::  strSection1     =   "null"
    character(3),   parameter   ::  strSection2     =   "end"
    character(18),  parameter   ::  strTemperature  =   "temperature points"
    character(6),   parameter   ::  strRedoxReaction = "o2(aq)"

    integer, parameter          ::  maxStrBuffers = 10000
    character(256)              ::  strBuffers(maxStrBuffers)       !Assume the maximum number of each data section 
    character(256)              ::  strBuffer        =   ""
    integer                      ::  iReadStat       =   1
    logical                      ::  bEndOfFile      =   .false.

    real,allocatable           ::  temperature(:)       ! temperature points
    integer                    ::  nTemperature = 0    ! number of temperature points

    type(TypePrimarySpecies), allocatable :: species(:)             ! species
    integer                               :: nSpecies = 0           ! number of species
    
    type(TypeDerivedSpecies), allocatable :: aqueousSpecies(:)      ! aqueous species
    integer                               :: nAqueousSpecies = 0    ! number of aqueous species
    integer                               :: nRedoxReaction  = 0
    
    type(TypeMineral), allocatable        :: minerals(:)            ! minerals
    integer                               :: nMinerals = 0          ! number of minerals
    
    type(TypeGas), allocatable            :: gases(:)               ! gases
    integer                               :: nGases = 0             ! number of gases    
    
    type(TypeSurfaceComplexes),allocatable  ::  surfaceComplexes(:)    ! surfaceComplexes
    integer                                 ::  nSurfaceComplexes = 0  ! number of surface complexes    

    
contains
    
    ! Open toughreact database
    subroutine OpenDbsTR
    
        implicit none
        
        integer :: stat 
        
        open(iUnitDbsTR, file = filePathDbsTR,  iostat = stat)
        
        if (stat /= 0) then
            call WriteLog("Error in opening TOUGHREACT databse file: " // trim(filePathDbsTR))
            call ErrorHandling
        else
            bOpenDbsTR = .true.
            call WriteLog("Open TOUGHREACT database success: " // trim(filePathDbsTR))
        end if 
    
    end subroutine    
    
    
    ! Read toughreact database
     subroutine ReadDbsTR
     
        implicit none  
        
        call WriteLog ("Begin reading TOUGHREACT database")
        
        call readHead
        
        call readTemperature     

        call readSpecies
        
        call readAqueousSpecies
        
        call readMinerals
        
        call readGases
        
        call readSurfaceComplexes

        call WriteLog ("End reading TOUGHREACT database")
     
     end subroutine ReadDbsTR
     
     ! Close toughreact database
     subroutine CloseDbsTR
     
        implicit none
        
        if (bOpenDbsTR) then
            
            close(iUnitDbsTR)
            
        end if
     
     end subroutine CloseDbsTR
     
    
    ! Read next line
    subroutine readNextLine
    
        implicit none
        
        integer :: i

        do while(.not. bEndOfFile)

            read(iUnitDbsTR, "(a)", iostat = iReadStat) strBuffer
       
            if (iReadStat > 0) then     !Error in reading
            
                call WriteLog("Error in reading toughreact databse file.")
            
                call ErrorHandling
            
            else if (iReadStat < 0) then    !End of file
            
                bEndOfFile = .true.
                
                exit            
            
            end if

            strBuffer = adjustl(strBuffer)

            if(len_trim(strBuffer) /= 0 ) then
            
                ! Conver all case to lower case
                call lowerCase(strBuffer)

                if (strBuffer(1:1) /= strComment) then
                    exit
                end if

            end if

        end do
    
    endsubroutine readNextLine
   
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
     subroutine getNameFromString(string, name)
     
        implicit none
        character(*), intent(in) :: string        
        character(nMNLTR), intent(out) :: name
        character(256)             :: tempStr
        integer :: i
        
        i = index(string, "'")
        tempStr = string(i + 1:)
        i = index(tempStr, "'")
        tempStr = adjustl(tempStr(:i-1))
        
        !if (len_trim(tempStr) > nMNLTR) then
        !    name = tempStr(1:nMNLTR)
        !else
            name = trim(tempStr)
        !end if
        
     end subroutine getNameFromString
    
   
    ! Read head
    subroutine readHead

        implicit none

        logical :: bFlag = .true.

        do while (.not. bEndOfFile)
            call readNextLine
            if (index(strBuffer,strEndOfHeader) == 1) then
                bFlag = .false.
                exit
            end if        
        end do

        if (bFlag) then
            call WriteLog("Error in finding lable: " // strEndOfheader)
            call ErrorHandling
        else
            call WriteLog("Read head success.")
        end if

    end subroutine readHead

    !Check if the number is zero
    function isZero(dvalue) result(bFlag)
        implicit none
        real, intent(in) :: dvalue
        logical ::bFlag
        bFlag = .false.
        if(abs(dvalue) < 1.0E-100) then
            bFlag = .true.
        else
            bFlag = .false.
        end if
    end function isZero
    
     ! Read temperature points
     subroutine readTemperature

        implicit none
        
        logical :: bFlag = .true.
       
        character(260) :: strTemp

        do while (.not. bEndOfFile)
            call readNextLine
            if (index(strBuffer,strTemperature) > 0) then                   
              
                call skipNValues(strBuffer,1)
                
                if (index(strTemperature, " ") > 0) then
                    call skipNValues(strBuffer,1)
                end if
            
                write(*,*) trim(strBuffer)
                
                read(strBuffer, *)  nTemperature

                if (allocated(temperature)) then
                    deallocate(temperature)
                end if
            
                allocate(temperature(nTemperature))
            
                read(strBuffer, *)  nTemperature, temperature
            
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
            !call WriteLog("Temperatures:")
            !call WriteLog(temperature)
        end if

     end subroutine readTemperature

     ! Read species
     subroutine readSpecies

        implicit none

        integer :: i

        character (nMNLTR) :: strName

        nSpecies = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, strName)
            if (trim(strName) == strSection1 .or. trim(strName) == strSection2) then                
                exit
            end if
            if (nSpecies == maxStrBuffers) then
                nErrors = nErrors + 1
                call WriteLog("The maximum database of TOUGHREACT is 10,000 lines." // & 
                                "Please modify the source code parameter maxStrBuffers to read more data.")
                call ErrorHandling
            end if
            nSpecies = nSpecies + 1
            strBuffers(nSpecies) = strBuffer
        end do

        if (nSpecies > 0) then
            if(allocated(species)) then
                deallocate(species)
            end if
            allocate(species(nSpecies))
            do i = 1, nSpecies
                call getNameFromString(strBuffers(i),species(i)%Name)
                !Check if this name is used as an alias
                if (GetNameFromAlias(species(i)%Name) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(species(i)%Name))
                end if
                call skipNValues(strBuffers(i), 1)
                read(strBuffers(i),*) species(i)%A0, species(i)%Z, species(i)%MWT
            end do
        end if

        call WriteLog("Number of Species: ")
        call WriteLog(nSpecies)
        call WriteLog("Read species success")

     end subroutine readSpecies    

     
     !Allocate STQ and NAM for aqueous species
     subroutine allocateOneAqueousSpecies(i,n, ntmp)
        
        implicit none
        integer, intent(in) :: i, n, ntmp
        
        if(allocated(aqueousSpecies(i)%STQ)) then            
            deallocate(aqueousSpecies(i)%STQ)
        end if
        allocate(aqueousSpecies(i)%STQ(n))
        
        if(allocated(aqueousSpecies(i)%NameOfSTQ)) then            
            deallocate(aqueousSpecies(i)%NameOfSTQ)
        end if
        allocate(aqueousSpecies(i)%NameOfSTQ(n)) 
        
        if(allocated(aqueousSpecies(i)%AKLOG)) then            
            deallocate(aqueousSpecies(i)%AKLOG)
        end if
        allocate(aqueousSpecies(i)%AKLOG(ntmp)) 
            
     end subroutine allocateOneAqueousSpecies

     ! Read aqueous species
     subroutine readAqueousSpecies
        
        implicit none
        
        integer :: i, j, k
        character(nMNLTR) :: strName
        character(nMNLTR) :: tempNames(3)

        i = 0
        nAqueousSpecies = 0
        nRedoxReaction = 0
        
        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, strName)
            if (trim(strName) == strSection1 .or. trim(strName) == strSection2) then                
                exit
            end if
            nAqueousSpecies = nAqueousSpecies + 1
            strBuffers(nAqueousSpecies) = strBuffer

            !check if the three names are the same
            i = i + 1
            tempNames(i) = trim(strName)
            if(i == 3) then
                if(tempNames(1) /= tempNames(2) .or. tempNames(1) /= tempNames(3)) then
                    call WriteLog("Error detected in reading aqueous species. Three names do not match: " &
                    // trim(tempNames(1)) // ", " // trim(tempNames(2)) // ", " // trim(tempNames(3)))
                    call ErrorHandling
                else
                    call WriteLog("Read in data: "//trim(tempNames(1)))
                end if
                i = 0
            end if
        end do

        if(mod(nAqueousSpecies,3) /= 0) then
            call WriteLog("Error detected in reading aqueous species. Number of data lines is not 3x.")
            call ErrorHandling
        end if

        nAqueousSpecies = nAqueousSpecies / 3
        
        if (nAqueousSpecies > 0) then
            if(allocated(aqueousSpecies)) then
                deallocate(aqueousSpecies)
            end if
            allocate(aqueousSpecies(nAqueousSpecies))
            do i = 1, nAqueousSpecies                
                !read first line
                j = 3 * i - 2
                call getNameFromString(strBuffers(j),aqueousSpecies(i)%Name)
                !Check if this name is used as an alias
                if (GetNameFromAlias(aqueousSpecies(i)%Name) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(species(i)%Name))
                end if
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) aqueousSpecies(i)%MWT, aqueousSpecies(i)%A0, aqueousSpecies(i)%Z, aqueousSpecies(i)%NCP
                call allocateOneAqueousSpecies(i,aqueousSpecies(i)%NCP, nTemperature)
                call skipNValues(strBuffers(j),4)            
                do k = 1, aqueousSpecies(i)%NCP, 1
                    read(strBuffers(j),*) aqueousSpecies(i)%STQ(k)
                    call getNameFromString(strBuffers(j),aqueousSpecies(i)%NameOfSTQ(k))
                    !Check if this name is used as an alias
                    if (GetNameFromAlias(aqueousSpecies(i)%NameOfSTQ(k)) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(aqueousSpecies(i)%NameOfSTQ(k)))
                    end if
                    if(k < aqueousSpecies(i)%NCP) then
                        call skipNValues(strBuffers(j),2)
                    end if
                    !write(*,*) trim(aqueousSpecies(i)%NameOfSTQ(m)), aqueousSpecies(i)%STQ(m)
                end do

                !Check if it is complexation reaction or redox reaction
                !................IS IT CORRECT?........................
                !do k = 1, aqueousSpecies(i)%NCP
                !    if(trim(aqueousSpecies(i)%NameOfSTQ(k)) == strRedoxReaction) then
                !        aqueousSpecies(i)%bRedoxReaction = .true.
                !        nRedoxReaction = nRedoxReaction + 1
                !        exit
                !    end if
                !end do

                !read the second line
                j = 3 * i - 1
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) aqueousSpecies(i)%AKLOG
                !write(*,*) aqueousSpecies(i)%AKLOG
                !read the third line
                j = 3 * i
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) aqueousSpecies(i)%AKCOE
                call skipNValues(strBuffers(j),5)
                
                if (len_trim(strBuffers(j)) > 0) then
                    read(strBuffers(j),*) aqueousSpecies(i)%AKCOP(1)
                    if(.not. isZero(aqueousSpecies(i)%AKCOP(1))) then
                        read(strBuffers(j),*) aqueousSpecies(i)%AKCOP
                        aqueousSpecies(i)%bAKCOP = .true.
                    end if
                end if   
                !write(*,*) aqueousSpecies(i)%AKCOE
            end do
        end if

        call WriteLog("Number of aqueous species: ")
        call WriteLog(nAqueousSpecies)
        call WriteLog("Number of redox reactions:")
        call WriteLog(nRedoxReaction)
        call WriteLog("Read aqueous species success")

     end subroutine readAqueousSpecies
     
     !Allocate STQ and NAM for mineral
     subroutine allocateOneMineral(i,n, ntmp)
        
        implicit none
        integer, intent(in) :: i, n, ntmp
        
        if(allocated(minerals(i)%STQ)) then            
            deallocate(minerals(i)%STQ)
        end if
        allocate(minerals(i)%STQ(n))
        
        if(allocated(minerals(i)%NameOfSTQ)) then            
            deallocate(minerals(i)%NameOfSTQ)
        end if
        allocate(minerals(i)%NameOfSTQ(n)) 
        
        if(allocated(minerals(i)%AKLOG)) then            
            deallocate(minerals(i)%AKLOG)
        end if
        allocate(minerals(i)%AKLOG(ntmp)) 
            
     end subroutine allocateOneMineral
     
     ! Read minerals
     subroutine readMinerals
        
        implicit none
        
        integer :: i, j, k
        character(nMNLTR) :: strName
        character(nMNLTR) :: tempNames(3)

        i = 0
        nMinerals = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, strName)
            if (trim(strName) == strSection1 .or. trim(strName) == strSection2) then                
                exit
            end if
            nMinerals = nMinerals + 1
            strBuffers(nMinerals) = strBuffer

          !check if the three names are the same
            i = i + 1
            tempNames(i) = trim(strName)
            if(i == 3) then
                if(tempNames(1) /= tempNames(2) .or. tempNames(1) /= tempNames(3)) then
                    call WriteLog("Error detected in reading minerals. Three names do not match: " &
                    // trim(tempNames(1)) // ", " // trim(tempNames(2)) // ", " // trim(tempNames(3)))
                    call ErrorHandling
                else
                    call WriteLog("Read in data: "//trim(tempNames(1)))
                end if
                i = 0
            end if
        end do

        if(mod(nMinerals,3) /= 0) then
            call WriteLog("Error detected in reading minerals. Number of data lines is not 3x.")
            call ErrorHandling
        end if


        nMinerals = nMinerals / 3
        
        if (nMinerals > 0) then
            if(allocated(minerals)) then
                deallocate(minerals)
            end if
            allocate(minerals(nMinerals))
            do i = 1, nMinerals                
                !read first line
                j = 3 * i - 2
                call getNameFromString(strBuffers(j),minerals(i)%Name)
                
                !Check if this name is used as an alias
                if (GetNameFromAlias(minerals(i)%Name) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(minerals(i)%Name))
                end if
                
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) minerals(i)%MWT, minerals(i)%VMIN, minerals(i)%NCP
                call allocateOneMineral(i,minerals(i)%NCP, nTemperature)
                call skipNValues(strBuffers(j),3)            
                do k = 1, minerals(i)%NCP, 1
                    read(strBuffers(j),*) minerals(i)%STQ(k)
                    call getNameFromString(strBuffers(j),minerals(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    if (GetNameFromAlias(minerals(i)%NameOfSTQ(k)) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(minerals(i)%NameOfSTQ(k)))
                    end if
                    
                    if(k < minerals(i)%NCP) then
                        call skipNValues(strBuffers(j),2)
                    end if
                    !write(*,*) trim(minerals(i)%NameOfSTQ(k)), minerals(i)%STQ(k)
                end do
                !read the second line
                j = 3 * i - 1
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) minerals(i)%AKLOG
                !write(*,*) minerals(i)%AKLOG
                !read the third line
                j = 3 * i
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) minerals(i)%AKCOE
                call skipNValues(strBuffers(j),5)
   
                if (len_trim(strBuffers(j)) > 0) then
                    read(strBuffers(j),*) minerals(i)%AKCOP(1)
                    if(.not. isZero(minerals(i)%AKCOP(1))) then
                        read(strBuffers(j),*) minerals(i)%AKCOP
                        minerals(i)%bAKCOP = .true.
                    end if
                end if 

                !write(*,*) minerals(i)%AKCOE
            end do
        end if

        call WriteLog("Number of minerals: ")
        call WriteLog(nMinerals)
        call WriteLog("Read minerals success")

     end subroutine readMinerals
     
      !Allocate STQ and NAM for gas
     subroutine allocateOneGas(i,n, ntmp)
        
        implicit none
        integer, intent(in) :: i, n, ntmp
        
        if(allocated(gases(i)%STQ)) then            
            deallocate(gases(i)%STQ)
        end if
        allocate(gases(i)%STQ(n))
        
        if(allocated(gases(i)%NameOfSTQ)) then            
            deallocate(gases(i)%NameOfSTQ)
        end if
        allocate(gases(i)%NameOfSTQ(n)) 
        
        if(allocated(gases(i)%AKLOG)) then            
            deallocate(gases(i)%AKLOG)
        end if
        allocate(gases(i)%AKLOG(ntmp)) 
            
     end subroutine allocateOneGas
     
     ! Read gases
     subroutine readGases
     
        implicit none
        
        integer :: i, j, k
        character(nMNLTR) :: strName
        character(nMNLTR) :: tempNames(3)

        i = 0
        nGases = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, strName)
            if (trim(strName) == strSection1 .or. trim(strName) == strSection2) then                
                exit
            end if
            nGases = nGases + 1
            strBuffers(nGases) = strBuffer
           !check if the three names are the same
            i = i + 1
            tempNames(i) = trim(strName)
            if(i == 3) then
                if(tempNames(1) /= tempNames(2) .or. tempNames(1) /= tempNames(3)) then
                    call WriteLog("Error detected in reading gases. Three names do not match: " &
                    // trim(tempNames(1)) // ", " // trim(tempNames(2)) // ", " // trim(tempNames(3)))
                    call ErrorHandling
                else
                    call WriteLog("Read in data: "//trim(tempNames(1)))
                end if
                i = 0
            end if
        end do

        if(mod(nGases,3) /= 0) then
            call WriteLog("Error detected in reading gases. Number of data lines is not 3x.")
            call ErrorHandling
        end if

        nGases = nGases / 3
        
        if (nGases > 0) then
            if(allocated(gases)) then
                deallocate(gases)
            end if
            allocate(gases(nGases))
            do i = 1, nGases                
                !read first line
                j = 3 * i - 2
                call getNameFromString(strBuffers(j),gases(i)%Name)
                
                !Check if this name is used as an alias
                if (GetNameFromAlias(gases(i)%Name) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(gases(i)%Name))
                end if
                
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) gases(i)%MWT, gases(i)%DMDIAM, gases(i)%NCP
                call allocateOneGas(i,gases(i)%NCP, nTemperature)
                call skipNValues(strBuffers(j),3)            
                do k = 1, gases(i)%NCP, 1
                    read(strBuffers(j),*) gases(i)%STQ(k)
                    call getNameFromString(strBuffers(j),gases(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    if (GetNameFromAlias(gases(i)%NameOfSTQ(k)) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(gases(i)%NameOfSTQ(k)))
                    end if     
                    
                    if(k < gases(i)%NCP) then
                        call skipNValues(strBuffers(j),2)
                    end if
                    !write(*,*) trim(gases(i)%NameOfSTQ(k)), gases(i)%STQ(k)
                end do
                !read the second line
                j = 3 * i - 1
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) gases(i)%AKLOG
                !write(*,*) gases(i)%AKLOG
                !read the third line
                j = 3 * i
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) gases(i)%AKCOE
                call skipNValues(strBuffers(j),5)
                if (len_trim(strBuffers(j)) > 0) then
                    read(strBuffers(j),*) gases(i)%AKCOP(1)
                    if(.not. isZero(gases(i)%AKCOP(1))) then
                        read(strBuffers(j),*) gases(i)%AKCOP
                        gases(i)%bAKCOP = .true.
                    end if
                end if   
                !write(*,*) gases(i)%AKCOE
            end do
        end if

        call WriteLog("Number of gases: ")
        call WriteLog(nGases)
        call WriteLog("Read gases success")
     
     end subroutine readGases
     
      !Allocate STQ and NAM for surface complexes
     subroutine allocateOneSurfaceComplexes(i,n, ntmp)
        
        implicit none
        integer, intent(in) :: i, n, ntmp
        
        if(allocated(surfaceComplexes(i)%STQ)) then            
            deallocate(surfaceComplexes(i)%STQ)
        end if
        allocate(surfaceComplexes(i)%STQ(n))
        
        if(allocated(surfaceComplexes(i)%NameOfSTQ)) then            
            deallocate(surfaceComplexes(i)%NameOfSTQ)
        end if
        allocate(surfaceComplexes(i)%NameOfSTQ(n)) 
        
        if(allocated(surfaceComplexes(i)%AKLOG)) then            
            deallocate(surfaceComplexes(i)%AKLOG)
        end if
        allocate(surfaceComplexes(i)%AKLOG(ntmp)) 
            
     end subroutine allocateOneSurfaceComplexes
     
     ! Read surface complexes
     subroutine readSurfaceComplexes
     
        implicit none
        
        integer :: i, j, k
        character(nMNLTR) :: strName
        character(nMNLTR) :: tempNames(3)

        i = 0
        nSurfaceComplexes = 0

        do while (.not. bEndOfFile)
            call readNextLine
            call getNameFromString(strBuffer, strName)
            if (trim(strName) == strSection1 .or. trim(strName) == strSection2) then                
                exit
            end if
            nSurfaceComplexes = nSurfaceComplexes + 1
            strBuffers(nSurfaceComplexes) = strBuffer
           !check if the three names are the same
            i = i + 1
            tempNames(i) = trim(strName)
            if(i == 3) then
                if(tempNames(1) /= tempNames(2) .or. tempNames(1) /= tempNames(3)) then
                    call WriteLog("Error detected in reading surface complexes. Three names do not match: " &
                    // trim(tempNames(1)) // ", " // trim(tempNames(2)) // ", " // trim(tempNames(3)))
                    call ErrorHandling
                else
                    call WriteLog("Read in data: "//trim(tempNames(1)))
                end if
                i = 0
            end if
        end do

        if(mod(nSurfaceComplexes,3) /= 0) then
            call WriteLog("Error detected in reading surface complexes. Number of data lines is not 3x.")
            call ErrorHandling
        end if

        nSurfaceComplexes = nSurfaceComplexes / 3
        
        if (nSurfaceComplexes > 0) then
            if(allocated(surfaceComplexes)) then
                deallocate(surfaceComplexes)
            end if
            allocate(surfaceComplexes(nSurfaceComplexes))
            do i = 1, nSurfaceComplexes                
                !read first line
                j = 3 * i - 2
                call getNameFromString(strBuffers(j),surfaceComplexes(i)%Name)
                
                !Check if this name is used as an alias
                if (GetNameFromAlias(surfaceComplexes(i)%Name) /= "") then
                    nErrors = nErrors + 1
                    call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(surfaceComplexes(i)%Name))
                end if
                
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) surfaceComplexes(i)%Z, surfaceComplexes(i)%NCP
                call allocateOneSurfaceComplexes(i,surfaceComplexes(i)%NCP, nTemperature)
                call skipNValues(strBuffers(j),2) 
                do k = 1, surfaceComplexes(i)%NCP, 1
                    read(strBuffers(j),*) surfaceComplexes(i)%STQ(k)
                    call getNameFromString(strBuffers(j),surfaceComplexes(i)%NameOfSTQ(k))
                    
                    !Check if this name is used as an alias
                    if (GetNameFromAlias(surfaceComplexes(i)%NameOfSTQ(k)) /= "") then
                        nErrors = nErrors + 1
                        call WriteLog("Error: This name has already been used as an alias in alias.dbs: "//trim(surfaceComplexes(i)%NameOfSTQ(k)))
                    end if
                    
                    if(k < surfaceComplexes(i)%NCP) then
                        call skipNValues(strBuffers(j),2)
                    end if
                    !write(*,*) trim(surfaceComplexes(i)%NameOfSTQ(k)), surfaceComplexes(i)%STQ(k)
                end do
                !read the second line
                j = 3 * i - 1
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) surfaceComplexes(i)%AKLOG
                !write(*,*) surfaceComplexes(i)%AKLOG
                !read the third line
                j = 3 * i
                call skipNValues(strBuffers(j),1)
                read(strBuffers(j),*) surfaceComplexes(i)%AKCOE                
                !write(*,*) surfaceComplexes(i)%AKCOE
            end do
        end if

        call WriteLog("Number of surface complexes: ")
        call WriteLog(nSurfaceComplexes)
        call WriteLog("Read surface complexes success")
     
     end subroutine readSurfaceComplexes
    
end module