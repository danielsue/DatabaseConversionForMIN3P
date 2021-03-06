!****************************************************************************
!    Geochemitry database conversion program
!    Author: Danyang Su
!    Email: dsu@eos.ubc.ca; danyang.su@gmail.com
!****************************************************************************
! Module of min3p database
    
module dbs_min3p
        
    use global
    
    use logfile, only : WriteLog
    
    use inputfile, only : targetDatabasePath, IsInSpecifiedExport, bSpecifiedExport, sourceDatabaseType
    
    implicit none
    
    integer, parameter:: nMNLmin3P = 12         !Please note: if you enlarge the maximum length of min3p, you should modify the output format.
    
    real                ::     temperature_min3p = 25       !Default temperature 
    character(60)      ::     masterVariable_min3p = ""     !Database master varialbe
    
    ! define species data structure
    type typeMin3PSpecies        
        character(nMNLmin3P)   ::      Name = ""       !name or chemical formula of aqueous basis species       
        real            ::      Z = 0.0d0       !the ion electric charge        
        real            ::      DHA = 0.0d0     !Debye-Huckel constants a        
        real            ::      DHB = 0.0d0     !Debye-Huckel constants b
        real            ::      MWT = 0.0d0     !Molecular weight of the aqueous species (g/mol)
        real            ::      AlkFac = 0.0d0  !Alkalinity factor    
    end type typeMin3PSpecies

    ! define aqueous data structure
    type typeMin3PAqueousSpecies

        character(nMNLmin3P)   ::      Name = ""                  !Name of the aqueous species
        real            ::      EnthalpyChange = 0.0d0     !Enthalpy change
        real            ::      AKLOG_exp = 0.0d0          !Equilibrium constant
        real            ::      AKLOG(20)                  !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                            !temperature listed in record
        real            ::      Z = 0.0d0                  !charge of the species
        real            ::      DHA = 0.0d0                !Debye-Huckel constants a        
        real            ::      DHB = 0.0d0                !Debye-Huckel constants b
        real            ::      AlkFac = 0.0d0             !Alkalinity factor
        real            ::      MWT = 0.0d0                !Molecular weight of the aqueous species (g/mol)
        integer         ::      NCP = 0                    !number of basis species defining the secondary species
        real            ::      STQ(20) = 0.0d0             !stoichiometric coefficients of the components comprising the aqueous species                                                        
        character(nMNLmin3P)   ::      NameOfSTQ(20) = ""          !name of the reactant or product
        
        integer         ::      iAssociation = 1           !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation

    end type typeMin3PAqueousSpecies
    
    ! define gas data structure
    type typeMin3PGas        
        
        character(nMNLmin3P)   ::      Name = ""                 !Name of the aqueous species
        real            ::      EnthalpyChange = 0.0d0     !Enthalpy change
        real            ::      AKLOG_exp = 0.0d0          !Equilibrium constant
        real            ::      AKLOG(20)                  !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                            !temperature listed in record
        real            ::      MWT = 0.0d0                !Molecular weight of the aqueous species (g/mol)
        integer         ::      NCP = 0                    !number of basis species defining the secondary species
        real            ::      STQ(20) = 0.0d0             !stoichiometric coefficients of the components                                                       
        character(nMNLmin3P)   ::      NameOfSTQ(20) = ""          !name of the reactant or product  
        
        integer         ::      iAssociation = 1           !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation
        
    end type typeMin3PGas
    
    ! define miniral data structure
    type typeMin3PMineral
        
        character(nMNLmin3P)   ::      Name = ""                 !Name of the mineral phase
        real            ::      AKLOG_exp = 0.0d0          !Equilibrium constant
        real            ::      AKLOG(20)                  !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                            !temperature listed in record
        real            ::      EnthalpyChange = 0.0d0     !Enthalpy change
        real            ::      MWT = 0.0d0               !Molecular weight of the mineral (g/mol)
        real            ::      Density = 0.0d0           !Density, g/cm^3
        integer         ::      NCP = 0                    !number of basis species defining the secondary species
        real            ::      STQ(20) = 0.0d0             !stoichiometric coefficients of the components                                                       
        character(nMNLmin3P)   ::      NameOfSTQ(20) = ""          !name of the reactant or product 
        
        integer         ::      iAssociation = 1           !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation
        
    end type typeMin3PMineral   

    
    type(typeMin3PSpecies), allocatable :: min3PSpecies(:)
    integer :: nMin3PSpecies = 0

    type(typeMin3PAqueousSpecies), allocatable :: min3PComplexReactions(:)
    integer :: nMin3PComplexReactions = 0

    type(typeMin3PAqueousSpecies), allocatable :: min3PRedoxReactions(:)
    integer :: nMin3PRedoxReactions = 0
    
    type(typeMin3PGas), allocatable :: min3PGases(:)
    integer :: nMin3PGases = 0
    
    type(typeMin3PMineral), allocatable :: min3PMinerals(:)
    integer :: nMin3PMinerals
    
    character(8), parameter  ::  filePathDbsMin3PComp      =   "comp.dbs"       !comp.dbs for min3p
    character(11), parameter ::  filePathDbsMin3PComplex   =   "complex.dbs"    !complex.dbs for min3p
    character(9), parameter  ::  filePathDbsMin3PGases     =   "gases.dbs"      !gases.dbs for min3p
    character(9), parameter  ::  filePathDbsMin3PRedox     =   "redox.dbs"      !redox.dbs for min3p
    character(11), parameter ::  filePathDbsMin3PMinerals  =   "mineral.dbs"    !mineral.dbs for min3p
    
    !!**************************************************************************************
    !! IMPORTANT NOTE ON THE SIGH OF EQUILIBRIUM CONSTANT
    !! 1. Database of MIN3P uses association reactions
    !! 2. Database of Toughreact or Crunchflow uses dissociation reaction (sign turn around)
    !! 3. Database of Phreeqc has both association and dissociation reactions
    !!    -Secondary species: association (the same sign)
    !!    -Minerals (Phases): dissociation (sign turn around)
    !!    -Gases (Phases): dissociation (sign turn around)
    !!**************************************************************************************
    !!MIN3P uses association reactions (formation of dependent/secondary species). 
    !!(the sign turn around, if reactions are expressed as dissociation reactions (not done in MIN3P for complexes, gases and minerals).
    !!
    !!For complexes:
    !!------------------
    !!Example:
    !!
    !!mgoh+            15.9520  -11.4400                 1.00 6.50  .00  41.3190    .00
    !!      3   mg+2           1.000 h2o            1.000 h+1           -1.000
    !!
    !!This stands for the association reaction:  mg+2   + h2o - h+1  -> mgoh+
    !!
    !!K = products/reactants
    !!K = [MgOH+] [H+] [Mg2+]^-1
    !![MgOH+] = K x  [H+]^-1 [Mg2+]
    !![MgOH+] = 10^- 11.4400  x  [H+]^-1 [Mg2+]
    !!
    !!
    !!For gases:
    !!-------------
    !!Example:
    !!
    !!ch4(g)             3.373   2.8894                                  16.0432  !modified 9/99
    !!      1   ch4(aq)        1.000
    !!
    !!This stands for the association reaction: ch4(aq) -> ch4(g)
    !!K = pCH4/{CH4(aq)}
    !!pCH4 = K {CH4(aq)}
    !!pCH4 = 10^2.8894 {CH4(aq)}
    !!
    !!
    !!for minerals:
    !! --------------
    !!Example:
    !!
    !!'calcite'
    !!'surface'
    !!100.0894    2.7100
    !!2  'ca+2'    1.000  'co3-2'    1.000
    !!'reversible'    8.4750  2.5850
    !!
    !!Stands for the association reaction: Ca2+ + CO32- ?calcite    K_a
    !!K_a = 1/{Ca2+}{CO32-}
    !!10^8.4750  = 1/{Ca2+}{CO32-}
    !! 
    !!
    !!This is consistent with the above for complexes and gases, but is a bit backward in terms of the standard definition:
    !!SR = IAP/Ksp = {Ca2+}{CO32-}/Ksp
    !!Which is valid for a dissociation reaction
    !!Calcite -> Ca2+ + CO32-     Ksp
    
    
contains

    ! Open all min3p database
    subroutine OpenAllDbsMin3P
    
        use ifport
        USE IFPOSIX
    
        implicit none
        integer :: opendirid, ierror
        logical ::  path_exists = .false. 
        
        ierror = -1
        
        !check and create folder if necessary
        
        if(len_trim(targetDatabasePath) > 0) then
            !inquire(file = trim(targetDatabasePath), exist = path_exists)      !do not use this, not correct
            CALL PXFOPENDIR (trim(targetDatabasePath),len_trim(targetDatabasePath),opendirid,ierror) 
            
            if (ierror /= 0) then
                path_exists = makedirqq(trim(targetDatabasePath))
                if (path_exists) then
                    call WriteLog("Create output folder for min3p database: success")
                else
                    call WriteLog("Create output folder for min3p database: failed")
                    call ErrorHandling
                end if
            end if
        end if
        
        if (len_trim(targetDatabasePath) == 0) then
            call openDbsMin3P(iUnitDbsMin3PComp, filePathDbsMin3PComp)
            call openDbsMin3P(iUnitDbsMin3PComplex, filePathDbsMin3PComplex)
            call openDbsMin3P(iUnitDbsMin3PRedox, filePathDbsMin3PRedox)
            call openDbsMin3P(iUnitDbsMin3PGases, filePathDbsMin3PGases)
            call openDbsMin3P(iUnitDbsMin3PMinerals, filePathDbsMin3PMinerals) 
        else
            call openDbsMin3P(iUnitDbsMin3PComp, trim(targetDatabasePath) // "\" // filePathDbsMin3PComp)
            call openDbsMin3P(iUnitDbsMin3PComplex, trim(targetDatabasePath) // "\" // filePathDbsMin3PComplex)
            call openDbsMin3P(iUnitDbsMin3PRedox, trim(targetDatabasePath) // "\" // filePathDbsMin3PRedox)
            call openDbsMin3P(iUnitDbsMin3PGases, trim(targetDatabasePath) // "\" // filePathDbsMin3PGases)
            call openDbsMin3P(iUnitDbsMin3PMinerals, trim(targetDatabasePath) // "\" // filePathDbsMin3PMinerals)
        end if
    
    end subroutine OpenAllDbsMin3P
    
    ! Write all min3p database
    subroutine WriteAllDbsMin3P
    
        implicit none
        
        call writeDbsMin3PComp
        call writeDbsMin3PComplex
        call writeDbsMin3PRedox
        call writeDbsMin3PGases
        call writeDbsMin3PMinerals
    
    end subroutine WriteAllDbsMin3P
   

    ! Open min3p database comp.dbs
    subroutine openDbsMin3P(iUnit,filePath)
    
        implicit none
        
        integer, intent(in) :: iUnit
        character(*), intent(in) :: filePath
        
        integer :: stat 
        
        open(iUnit, file = trim(filePath),  iostat = stat, status = "replace")
        
        if (stat /= 0) then
            call WriteLog("Error in opening/creating min3p databse file: " // trim(filePath))
            call ErrorHandling
        else
            bOpenDbsMin3PComp = .true.
            call WriteLog("Open min3p database success: " // trim(filePath))
        end if
    
    end subroutine openDbsMin3P    
    
    ! Write min3p database comp.dbs
    subroutine writeDbsMin3PComp
    
        implicit none
        
        integer ::  i
        
        !write(iUnitDbsMin3PComp, 99) "!This database is converted from " // trim(sourceDatabaseType)
        !if(trim(sourceDatabaseType) == "phreeqc") then
        !    write(iUnitDbsMin3PComp, 99) "!IMPORTANT NOTE: PLEASE CHECK IF THE MOLECULAR WEIGHT OF EACH COMPONENT IS CORRECT"
        !end if
        do i = 1, nMin3PSpecies
            if (bSpecifiedExport) then
                if (.not. IsInSpecifiedExport(min3PSpecies(i)%Name)) then
                    cycle
                end if
            end if
            write(iUnitDbsMin3PComp, 100) min3PSpecies(i)%Name, min3PSpecies(i)%Z, min3PSpecies(i)%DHA,  &
                                            min3PSpecies(i)%DHB, min3PSpecies(i)%MWT, min3PSpecies(i)%AlkFac
        end do
        
        call WriteLog("Write min3p database success: " // trim(filePathDbsMin3PComp)) 

99      format(a)        
100     format(a12,f4.1,4x,f5.2,f5.2,8x,f11.5,f7.2)
    
    end subroutine writeDbsMin3PComp

    
    ! Write min3p database complex.dbs
    subroutine writeDbsMin3PComplex
    
        implicit none
        
        integer ::  i, j
        
        !write(iUnitDbsMin3PComplex, 199) "!This database is converted from " // trim(sourceDatabaseType)
        do i = 1, nMin3PComplexReactions
            if (bSpecifiedExport) then
                if (.not. IsInSpecifiedExport(min3PComplexReactions(i)%Name)) then
                    cycle
                end if
            end if
            write(iUnitDbsMin3PComplex, 200) min3PComplexReactions(i)%Name, min3PComplexReactions(i)%EnthalpyChange, &
                                                min3PComplexReactions(i)%AKLOG_exp, min3PComplexReactions(i)%Z, &
                                                min3PComplexReactions(i)%DHA, min3PComplexReactions(i)%DHB, &
                                                min3PComplexReactions(i)%MWT, min3PComplexReactions(i)%AlkFac
            write(iUnitDbsMin3PComplex, 201) min3PComplexReactions(i)%NCP, ((min3PComplexReactions(i)%NameOfSTQ(j),min3PComplexReactions(i)%STQ(j)), &
                                            j = 1, min3PComplexReactions(i)%NCP)
        end do
        
        call WriteLog("Write min3p database success: " // trim(filePathDbsMin3PComplex)) 

199     format(a)        
200     format(a12,2x,2f10.4,16x,3f5.2,f9.4,f7.2)
201     format(6x,i1,3x,20(a12,1x,f7.3,1x))     
    
    end subroutine writeDbsMin3PComplex   

    
    ! Write min3p database redox.dbs
    subroutine writeDbsMin3PRedox
    
        implicit none
        
        integer ::  i, j
        
        !write(iUnitDbsMin3PRedox, 300) "!This database is converted from " // trim(sourceDatabaseType)
        do i = 1, nMin3PRedoxReactions
            if (bSpecifiedExport) then
                if (.not. IsInSpecifiedExport(min3PRedoxReactions(i)%Name)) then
                    cycle
                end if
            end if
            
            write(iUnitDbsMin3PRedox, 300) "! "
            write(iUnitDbsMin3PRedox, 300) "! "// trim(min3PRedoxReactions(i)%Name)
            write(iUnitDbsMin3PRedox, 300) "! equilibrium constant is calculated from other database"
            write(iUnitDbsMin3PRedox, 300) "! enthalpy change is set to be zero if not available"            
            write(iUnitDbsMin3PRedox, 300) "'"//trim(min3PRedoxReactions(i)%Name)//"'"
            write(iUnitDbsMin3PRedox, 301) min3PRedoxReactions(i)%NCP, (("'"//trim(min3PRedoxReactions(i)%NameOfSTQ(j))//"'",min3PRedoxReactions(i)%STQ(j)), &
                                            j = 1, min3PRedoxReactions(i)%NCP)
            write(iUnitDbsMin3PRedox, 302) min3PRedoxReactions(i)%AKLOG_exp, min3PRedoxReactions(i)%EnthalpyChange
        end do
        
        call WriteLog("Write min3p database success: " // trim(filePathDbsMin3PRedox)) 
        
        !Unlike for the database files comp.dbs, complex.dbs, gases.dbs and
        !sorption.dbs, the input in redox.dbs is format-free. Each line starting with ! is
        !considered a comment line.
        
300     format(a)
301     format(i1,3x,20(a,1x,f7.3,1x))       
302     format("'equilibrium'", 4x, 2f8.4)         
    
    end subroutine writeDbsMin3PRedox    
 
    ! Write min3p database gases.dbs
    subroutine writeDbsMin3PGases
    
        implicit none
        
        integer ::  i, j
        !write(iUnitDbsMin3PGases, 399) "!This database is converted from " // trim(sourceDatabaseType)
        do i = 1, nmin3PGases
            if (bSpecifiedExport) then 
                if (.not. IsInSpecifiedExport(min3PGases(i)%Name)) then
                    cycle
                end if
            end if
            write(iUnitDbsMin3PGases, 400) min3PGases(i)%Name,min3PGases(i)%EnthalpyChange, min3PGases(i)%AKLOG_exp, min3PGases(i)%MWT
            write(iUnitDbsMin3PGases, 401) min3PGases(i)%NCP, ((min3PGases(i)%NameOfSTQ(j),min3PGases(i)%STQ(j)),j = 1, min3PGases(i)%NCP)
        end do
        
        call WriteLog("Write min3p database success: " // trim(filePathDbsMin3PGases))         
399     format(a)        
400     format(a12,2x,2f10.4,31x,f9.4)
401     format(6x,i1,3x,20(a12,1x,f7.3,1x))        
    
    end subroutine writeDbsMin3PGases
    
    ! Write min3p database minerals.dbs
    !How to know from TOUGHREACT if the mineral is reversible or reversible 
    !and the parallel reaction pathways, which are needed in MIN3P? 
    !At present, exported to as �surface?and 'reversible'
    !?'uo3(c)'
    !?'surface'
    !? 286.0272    8.0435
    !? 3  'h+1'   -2.000  'uo2+2'    1.000  'h2o'    1.000
    !?'reversible'   -7.7190  19.3150
    subroutine writeDbsMin3PMinerals
    
        implicit none
        
        integer ::  i, j
        !write(iUnitDbsMin3PMinerals, 500) "!This database is converted from " // trim(sourceDatabaseType)
        do i = 1, nMin3PMinerals
            if (bSpecifiedExport) then 
                if (.not. IsInSpecifiedExport(min3PMinerals(i)%Name)) then
                    cycle
                end if
            end if

            write(iUnitDbsMin3PMinerals, 500) "! "
            write(iUnitDbsMin3PMinerals, 500) "! " // trim(min3PMinerals(i)%Name)            
            write(iUnitDbsMin3PMinerals, 500) "! equilibrium constant and formula weight are calculated from other database"
            write(iUnitDbsMin3PMinerals, 500) "! enthalpy change is set to 0 if not available"
            if(trim(sourceDatabaseType) == "phreeqc") then
                write(iUnitDbsMin3PMinerals, 500) "! Lack density, set this value to  1.0 instead"
            end if
            write(iUnitDbsMin3PMinerals, 500) "'"//trim(min3PMinerals(i)%Name)//"'"
            write(iUnitDbsMin3PMinerals, 500) "'surface'"
            write(iUnitDbsMin3PMinerals, 501) min3PMinerals(i)%MWT, min3PMinerals(i)%Density            
            write(iUnitDbsMin3PMinerals, 502) min3PMinerals(i)%NCP, (("'"//trim(min3PMinerals(i)%NameOfSTQ(j))//"'",min3PMinerals(i)%STQ(j)),j = 1, min3PMinerals(i)%NCP)
            write(iUnitDbsMin3PMinerals, 503) "'reversible'",min3PMinerals(i)%AKLOG_exp, min3PMinerals(i)%EnthalpyChange
        end do
        
        call WriteLog("Write min3p database success: " // trim(filePathDbsMin3PMinerals))         
        
500     format(a)
501     format(2(f10.4,2x))
502     format(i2,2x,20(a,1x,f7.3,1x))      !change from (i1,3x,20(a,1x,f7.3,1x))
503     format(a, 4x, 2(f10.4,1x))        
    
    end subroutine writeDbsMin3PMinerals
    
end module dbs_min3p