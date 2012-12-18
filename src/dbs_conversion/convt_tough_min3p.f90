! module of toughreact database to min3p converter

module convt_tough_min3p

    use dbs_toughreact, only :     temperature, nTemperature, species, nSpecies, aqueousSpecies,            &
                                    nAqueousSpecies, minerals, nMinerals, gases, nGases, surfaceComplexes,  &
                                    nSurfaceComplexes, nRedoxReaction, nMNLTR
    
    use dbs_min3p,      only :     min3PSpecies, nMin3PSpecies, min3PComplexReactions,                      &
                                    nMin3PComplexReactions, min3PRedoxReactions, nMin3PRedoxReactions,      &
                                    min3PGases,nMin3PGases, min3PMinerals, nMin3PMinerals,nMNLmin3P
    
    use logfile,        only :     WriteLog, nWarnings, nErrors
    
    use name_truncation, only :     AddNameTruncation
    use alias, only :               GetAliasFromName
    
contains    

    ! Convert toughreact database to min3p
    subroutine convert2Min3PDbs
    
        implicit none
        
        call convert2Min3PSpecies

        call convert2Min3PAqueousSpecies
        
        call convert2Min3PGases
        
        call convert2Min3PMinerals
    
    end subroutine convert2Min3PDbs

    ! Convert species
    subroutine convert2Min3PSpecies
    
        implicit none
        
        integer :: i        
        
        nMin3PSpecies = nSpecies
        if(allocated(min3PSpecies)) then
            deallocate(min3PSpecies)
        end if
        
        if (nMin3PSpecies > 0) then
            allocate(min3PSpecies(nSpecies))
        end if
        
        call WriteLog("Start converting species")
        
        do i = 1, nSpecies
            !Convert name
            !Length of min3p species name is 12 maximum while toughreact species name is 20,
            !Deal this first if there exists the same name after truncation
            !   .........................................            
            
            
            if(len_trim(species(i)%Name)>nMNLmin3P) then
                min3PSpecies(i)%Name = GetAliasFromName(species(i)%Name)       !get alias for long name
                if(trim(min3PSpecies(i)%Name) == "") then
                    min3PSpecies(i)%Name = species(i)%Name(1:nMNLmin3P)
                    call WriteLog("Name truncation: " // trim(species(i)%Name) // " to " // species(i)%Name(1:nMNLmin3P)) 
                    call AddNameTruncation(trim(species(i)%Name), species(i)%Name(1:nMNLmin3P))
                else
                    call WriteLog("Alias used for : " // trim(species(i)%Name) // " to " // trim(min3PSpecies(i)%Name))
                end if
            else
                min3PSpecies(i)%Name = GetAliasFromName(species(i)%Name)
                if(trim(min3PSpecies(i)%Name) == "") then
                    min3PSpecies(i)%Name = species(i)%Name
                else
                    call WriteLog("Alias used for : " // trim(species(i)%Name) // " to " // trim(min3PSpecies(i)%Name))
                end if    
            end if
            
            !Convert charge
            min3PSpecies(i)%Z = species(i)%Z
            
            !Convert Debye-Huckel constants
            !Ion effective or hydrated radius used to compute the Debye-Huckel a0 parameter (see Appendix H
            !of Toughreact_V2_User_Guide.pdf for details).
            !DHA = 2(r + 1.91|Z|)/(|Z| + 1) for anions      r is the effective ionic radius, referred to %A0
            !DHA = 2(r + 1.81|Z|)/(|Z| + 1) for cations      r is the effective ionic radius, referred to %A0
            !   .........................................
            if (species(i)%Z < 0) then
                min3PSpecies(i)%DHA = 2.0*(species(i)%A0 + 1.91 * abs(species(i)%Z)) / (abs(species(i)%Z) + 1.0)
            else if (species(i)%Z > 0) then
                min3PSpecies(i)%DHA = 2.0*(species(i)%A0 + 1.81 * abs(species(i)%Z)) / (abs(species(i)%Z) + 1.0)
            else                                                            !!When Z = 0, what is the formulation
                min3PSpecies(i)%DHA = species(i)%A0
            end if
                
            min3PSpecies(i)%DHB = 0
            
            !Convert the gram formula weight, molecular weight
            min3PSpecies(i)%MWT = species(i)%MWT
            
            !Calculate alkalinity factor
            if(trim(min3PSpecies(i)%Name) == "co3-2" .or. trim(min3PSpecies(i)%Name) == "co3--") then
                min3PSpecies(i)%AlkFac = 2.0d0
            else if(trim(min3PSpecies(i)%Name) == "hco3-" .or. trim(min3PSpecies(i)%Name) == "hco3-1") then
                min3PSpecies(i)%AlkFac = 1.0d0
            else if (trim(min3PSpecies(i)%Name) == "h+" .or. trim(min3PSpecies(i)%Name) == "h+1") then
                min3PSpecies(i)%AlkFac = -1.0d0
            end if
            
            call WriteLog("Convert species:" // trim(min3PSpecies(i)%Name))
            
        end do
        
        call WriteLog("End converting species")
    
    end subroutine convert2Min3PSpecies

    !Convert aqueous species
    !Aqueous species in TOUGHREACT include Complexation Reactions and Redox Reactions in min3p, judged by o2(aq)?
    !For some complex.dbs, it also includes o2(aq) ?
    !At present, all the reactant or product that contain o2(aq) are treated as redox reactions and the others are treated as the complexation
    subroutine convert2Min3PAqueousSpecies

        implicit none

        integer :: i, j, k , m, n

        nMin3PRedoxReactions = nRedoxReaction
        nMin3PComplexReactions = nAqueousSpecies - nRedoxReaction

        if (allocated(min3PRedoxReactions)) then
            deallocate(min3PRedoxReactions)
        end if
        
        if(nMin3PRedoxReactions>0) then
            allocate(min3PRedoxReactions(nMin3PRedoxReactions))
        end if
        
        if (allocated(min3PComplexReactions)) then
            deallocate(min3PComplexReactions)
        end if  
        
        if (nMin3PComplexReactions > 0) then
            allocate(min3PComplexReactions(nMin3PComplexReactions))
        end if

        j = 0
        k = 0
        
        call WriteLog("Start converting aqueous species: complexation reactions and redox reactions")
        
        do i = 1, nAqueousSpecies

            if(aqueousSpecies(i)%bRedoxReaction) then       !Redox reaction
                !Convert name
                !Length of min3p species name is 12 maximum while toughreact species name is 20,
                !Deal this first if there exists the same name after truncation
                !   .........................................
                j = j + 1

                if(len_trim(aqueousSpecies(i)%Name)>nMNLmin3P) then
                    min3PRedoxReactions(j)%Name = GetAliasFromName(aqueousSpecies(i)%Name)
                    if (trim(min3PRedoxReactions(j)%Name) == "") then                        
                        min3PRedoxReactions(j)%Name = aqueousSpecies(i)%Name(1:nMNLmin3P)
                        call WriteLog("Name truncation: " // trim(aqueousSpecies(i)%Name) // " to " // aqueousSpecies(i)%Name(1:nMNLmin3P)) 
                        call AddNameTruncation(trim(aqueousSpecies(i)%Name),aqueousSpecies(i)%Name(1:nMNLmin3P))
                    else
                        call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%Name) // " to " // trim(min3PRedoxReactions(j)%Name))
                    end if
                else
                    min3PRedoxReactions(j)%Name = GetAliasFromName(aqueousSpecies(i)%Name)
                    if(trim(min3PRedoxReactions(j)%Name) == "") then
                        min3PRedoxReactions(j)%Name = aqueousSpecies(i)%Name
                    else
                        call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%Name) // " to " // trim(min3PRedoxReactions(j)%Name))
                    end if
                end if
                
                min3PRedoxReactions(j)%EnthalpyChange = 0.0d0
                !The second temperature in TOUGHREACT 25C is used. 
                !Should modify if 25C does not exist or in the different position.
                min3PRedoxReactions(j)%AKLOG = -aqueousSpecies(i)%AKLOG(2)      !Min3P use the reverse equilibrium constant
                min3PRedoxReactions(j)%Z = aqueousSpecies(i)%Z
                !Convert Debye-Huckel constants
                !Ion effective or hydrated radius used to compute the Debye-Huckel a0 parameter (see Appendix H
                !of Toughreact_V2_User_Guide.pdf for details).
                !DHA = 2(r + 1.91|Z|)/(|Z| + 1) for anions      r is the effective ionic radius, referred to %A0
                !DHA = 2(r + 1.81|Z|)/(|Z| + 1) for cations      r is the effective ionic radius, referred to %A0
                !   .........................................
                if (aqueousSpecies(i)%Z < 0) then
                    min3PRedoxReactions(j)%DHA = 2.0*(aqueousSpecies(i)%A0 + 1.91 * abs(aqueousSpecies(i)%Z)) / (abs(aqueousSpecies(i)%Z) + 1.0)
                else if (aqueousSpecies(i)%Z > 0) then
                    min3PRedoxReactions(j)%DHA = 2.0*(aqueousSpecies(i)%A0 + 1.81 * abs(aqueousSpecies(i)%Z)) / (abs(aqueousSpecies(i)%Z) + 1.0)
                else                                                            !!When Z = 0, what is the formulation
                    min3PRedoxReactions(j)%DHA = aqueousSpecies(i)%A0
                end if
                min3PRedoxReactions(j)%DHB = 0
                min3PRedoxReactions(j)%MWT = aqueousSpecies(i)%MWT

                !Calculate alkalinity factor
                !..............................?
                if(trim(min3PRedoxReactions(j)%Name) == "co3-2" .or. &
                   trim(min3PRedoxReactions(j)%Name) == "co3--") then
                    min3PRedoxReactions(j)%AlkFac = 2.0d0
                else if(trim(min3PRedoxReactions(j)%Name) == "hco3-" .or. &
                         trim(min3PRedoxReactions(j)%Name) == "hco3-1") then
                    min3PRedoxReactions(j)%AlkFac = 1.0d0
                else if (trim(min3PRedoxReactions(j)%Name) == "h+" .or. &
                          trim(min3PRedoxReactions(j)%Name) == "h+1") then
                    min3PRedoxReactions(j)%AlkFac = -1.0d0
                end if

                !Some of the aqueous/minerals in TOUGHREACT have more than 5 reactant or product, 
                !but min3p can only support 5, at present, the program only export 5
                !change to support 20
                !.....................................
                min3PRedoxReactions(j)%NCP = min(aqueousSpecies(i)%NCP,20)
                
                m = min3PRedoxReactions(j)%NCP                
                              
                if (aqueousSpecies(i)%NCP > 20) then
                    call WriteLog("Warning:"// "Number of components exceed 20, only 20 components are converted.") 
                    nWarnings = nWarnings + 1
                end if
                min3PRedoxReactions(j)%STQ(1:m) = aqueousSpecies(i)%STQ(1:m)
                
                do n = 1, m
                    if(len_trim(aqueousSpecies(i)%NameOfSTQ(n)) > nMNLmin3P) then                        
                        min3PRedoxReactions(j)%NameOfSTQ(n) = GetAliasFromName(aqueousSpecies(i)%NameOfSTQ(n))                        
                        if(trim(min3PRedoxReactions(j)%NameOfSTQ(n)) == "") then
                            min3PRedoxReactions(j)%NameOfSTQ(n) = aqueousSpecies(i)%NameOfSTQ(n)(1:nMNLmin3P)     !no quotes                        
                            call WriteLog("Name truncation: " // trim(aqueousSpecies(i)%NameOfSTQ(n)) &
                                        // " to " // trim(min3PRedoxReactions(j)%NameOfSTQ(n)))                         
                            call AddNameTruncation( trim(aqueousSpecies(i)%NameOfSTQ(n)), aqueousSpecies(i)%NameOfSTQ(n)(1:nMNLmin3P))
                        else
                            call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%NameOfSTQ(n)) // " to " // trim(min3PRedoxReactions(j)%NameOfSTQ(n)))
                        end if
                    else
                        min3PRedoxReactions(j)%NameOfSTQ(n) = GetAliasFromName(aqueousSpecies(i)%NameOfSTQ(n))         
                        if(trim(min3PRedoxReactions(j)%NameOfSTQ(n)) == "" ) then
                            min3PRedoxReactions(j)%NameOfSTQ(n) = aqueousSpecies(i)%NameOfSTQ(n)                  !no quotes
                        else
                            call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%NameOfSTQ(n)) // " to " // trim(min3PRedoxReactions(j)%NameOfSTQ(n)))
                        end if
                    end if
                end do
                
                call WriteLog("Convert redox reaction: "// trim(min3PRedoxReactions(j)%Name))

            else 
                k = k + 1
                !Complex reaction
                !Convert name
                !Length of min3p species name is 12 maximum while toughreact species name is 20,
                !Deal this first if there exists the same name after truncation
                !   .........................................
                
                if(len_trim(aqueousSpecies(i)%Name)>nMNLmin3P) then
                    min3PComplexReactions(k)%Name = GetAliasFromName(aqueousSpecies(i)%Name)
                    if(trim(min3PComplexReactions(k)%Name) == "") then
                        min3PComplexReactions(k)%Name = aqueousSpecies(i)%Name(1:nMNLmin3P)
                        call WriteLog("Name truncation: " // trim(aqueousSpecies(i)%Name) // " to " // aqueousSpecies(i)%Name(1:nMNLmin3P)) 
                        call AddNameTruncation(trim(aqueousSpecies(i)%Name), aqueousSpecies(i)%Name(1:nMNLmin3P))
                    else
                        call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%Name) // " to " // trim(min3PComplexReactions(k)%Name))
                    end if
                else
                    min3PComplexReactions(k)%Name = GetAliasFromName(aqueousSpecies(i)%Name)
                    if(trim(min3PComplexReactions(k)%Name) == "") then
                        min3PComplexReactions(k)%Name = aqueousSpecies(i)%Name
                    else
                        call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%Name) // " to " // trim(min3PComplexReactions(k)%Name))
                    end if
                end if
                
                min3PComplexReactions(k)%EnthalpyChange = 0.0d0
                !The second temperature in TOUGHREACT 25C is used. 
                !Should modify if 25C does not exist or in the different position.
                min3PComplexReactions(k)%AKLOG = -aqueousSpecies(i)%AKLOG(2)   
                min3PComplexReactions(k)%Z = aqueousSpecies(i)%Z
                !Convert Debye-Huckel constants
                !Ion effective or hydrated radius used to compute the Debye-Huckel a0 parameter (see Appendix H
                !of Toughreact_V2_User_Guide.pdf for details).
                !DHA = 2(r + 1.91|Z|)/(|Z| + 1) for anions      r is the effective ionic radius, referred to %A0
                !DHA = 2(r + 1.81|Z|)/(|Z| + 1) for cations      r is the effective ionic radius, referred to %A0
                !   .........................................
                if (aqueousSpecies(i)%Z < 0) then
                    min3PComplexReactions(k)%DHA = 2.0*(aqueousSpecies(i)%A0 + 1.91 * abs(aqueousSpecies(i)%Z)) / (abs(aqueousSpecies(i)%Z) + 1.0)
                else if(aqueousSpecies(i)%Z > 0) then
                    min3PComplexReactions(k)%DHA = 2.0*(aqueousSpecies(i)%A0 + 1.81 * abs(aqueousSpecies(i)%Z)) / (abs(aqueousSpecies(i)%Z) + 1.0)
                else                                                            !!When Z = 0, what is the formulation
                    min3PComplexReactions(k)%DHA = aqueousSpecies(i)%A0
                end if
                min3PComplexReactions(k)%DHB = 0
                min3PComplexReactions(k)%MWT = aqueousSpecies(i)%MWT

                !Calculate alkalinity factor
                !..............................?
                if(trim(min3PComplexReactions(k)%Name) == "co3-2" .or. &
                   trim(min3PComplexReactions(k)%Name) == "co3--") then
                    min3PComplexReactions(k)%AlkFac = 2.0d0
                else if(trim(min3PComplexReactions(k)%Name) == "hco3-" .or. &
                         trim(min3PComplexReactions(k)%Name) == "hco3-1") then
                    min3PComplexReactions(k)%AlkFac = 1.0d0
                else if (trim(min3PComplexReactions(k)%Name) == "h+" .or. &
                          trim(min3PComplexReactions(k)%Name) == "h+1") then
                    min3PComplexReactions(k)%AlkFac = -1.0d0
                end if

                !
                min3PComplexReactions(k)%NCP = min(aqueousSpecies(i)%NCP,20)
               
                m = min3PComplexReactions(k)%NCP
                
                if (aqueousSpecies(i)%NCP > 20) then
                    call WriteLog("Warning:"// "Number of components exceed 20, only 20 components are converted.") 
                    nWarnings = nWarnings + 1
                end if

                min3PComplexReactions(k)%STQ(1:m) = aqueousSpecies(i)%STQ(1:m)                
                
                do n = 1, m
                    if(len_trim(aqueousSpecies(i)%NameOfSTQ(n)) > nMNLmin3P) then
                        min3PComplexReactions(k)%NameOfSTQ(n) = GetAliasFromName(aqueousSpecies(i)%NameOfSTQ(n))
                        if(trim(min3PComplexReactions(k)%NameOfSTQ(n)) == "") then
                            min3PComplexReactions(k)%NameOfSTQ(n) = aqueousSpecies(i)%NameOfSTQ(n)(1:nMNLmin3P)     !no quotes                        
                            call WriteLog("Name truncation: " // trim(aqueousSpecies(i)%NameOfSTQ(n)) &
                                        // " to " // aqueousSpecies(i)%NameOfSTQ(n)(1:nMNLmin3P))                         
                            call AddNameTruncation(trim(aqueousSpecies(i)%NameOfSTQ(n)), aqueousSpecies(i)%NameOfSTQ(n)(1:nMNLmin3P))
                        else
                            call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%NameOfSTQ(n)) // " to " // trim(min3PComplexReactions(k)%NameOfSTQ(n)))
                        end if
                    else
                        min3PComplexReactions(k)%NameOfSTQ(n) = GetAliasFromName(aqueousSpecies(i)%NameOfSTQ(n))
                        if(trim(min3PComplexReactions(k)%NameOfSTQ(n)) == "") then 
                            min3PComplexReactions(k)%NameOfSTQ(n) = aqueousSpecies(i)%NameOfSTQ(n)                  !no quotes
                        else
                            call WriteLog("Alias used for : " // trim(aqueousSpecies(i)%NameOfSTQ(n)) // " to " // trim(min3PComplexReactions(k)%NameOfSTQ(n)))
                        end if
                    end if                    
                end do
                
                call WriteLog("Convert complexation reaction: "// trim(min3PComplexReactions(k)%Name))
                
            end if
        
        end do
        
        call WriteLog("End converting aqueous species: complexation reactions and redox reactions")

    end subroutine convert2Min3PAqueousSpecies
    
    !Convert to min3p gases
    subroutine convert2Min3PGases
    
        integer :: i, j

        nMin3PGases = nGases

        if (allocated(min3PGases)) then
            deallocate(min3PGases)
        end if
        
        if (nMin3PGases > 0) then
            allocate(min3PGases(nMin3PGases))        
        end if
        
        call WriteLog("Start converting gases")
        
        do i = 1, nGases

            !Convert name
            !Length of min3p species name is 12 maximum while toughreact species name is 20,
            !Deal this first if there exists the same name after truncation
            !   .........................................
            if(len_trim(gases(i)%Name)>nMNLmin3P) then
                min3PGases(i)%Name = GetAliasFromName(gases(i)%Name)
                if(trim(min3PGases(i)%Name)=="") then                    
                    min3PGases(i)%Name = gases(i)%Name(1:nMNLmin3P)
                    call WriteLog("Name truncation: " // trim(gases(i)%Name) // " to " // gases(i)%Name(1:nMNLmin3P)) 
                    call AddNameTruncation(trim(gases(i)%Name), gases(i)%Name(1:nMNLmin3P))
                else
                    call WriteLog("Alias used for : " // trim(gases(i)%Name) // " to " // trim(min3PGases(i)%Name))
                end if
            else
                min3PGases(i)%Name = GetAliasFromName(gases(i)%Name)
                if(trim(min3PGases(i)%Name) == "") then
                    min3PGases(i)%Name = gases(i)%Name
                else
                    call WriteLog("Alias used for : " // trim(gases(i)%Name) // " to " // trim(min3PGases(i)%Name))
                end if
            end if
            
            min3PGases(i)%EnthalpyChange = 0.0d0
            !The second temperature in TOUGHREACT 25C is used. 
            !Should modify if 25C does not exist or in the different position.
            min3PGases(i)%AKLOG = -gases(i)%AKLOG(2)   

            !Molecular weight
            min3PGases(i)%MWT = gases(i)%MWT

            !
            min3PGases(i)%NCP = min(gases(i)%NCP,20)
                
            j = min3PGases(i)%NCP
               
            if (gases(i)%NCP > 20) then
                call WriteLog("Warning:"// "Number of components exceed 20, only 20 components are converted.") 
                nWarnings = nWarnings + 1
            end if

            min3PGases(i)%STQ(1:j) = gases(i)%STQ(1:j)
               
            do k = 1, j
                if(len_trim(gases(i)%NameOfSTQ(k)) > nMNLmin3P) then
                    min3PGases(i)%NameOfSTQ(k) = GetAliasFromName(gases(i)%NameOfSTQ(k))
                    if(trim(min3PGases(i)%NameOfSTQ(k)) == "") then
                        min3PGases(i)%NameOfSTQ(k) = gases(i)%NameOfSTQ(k)(1:nMNLmin3P)     !no quotes         
                        call WriteLog("Name truncation: " // trim(gases(i)%NameOfSTQ(k)) &
                                    // " to " // trim(min3PGases(i)%NameOfSTQ(k))) 
                        call AddNameTruncation(trim(gases(i)%NameOfSTQ(k)), gases(i)%NameOfSTQ(k)(1:nMNLmin3P))
                    else
                        call WriteLog("Alias used for : " // trim(gases(i)%NameOfSTQ(k)) // " to " // trim(min3PGases(i)%NameOfSTQ(k)))
                    end if
                else
                    min3PGases(i)%NameOfSTQ(k) = GetAliasFromName(gases(i)%NameOfSTQ(k))
                    if(trim(min3PGases(i)%NameOfSTQ(k)) == "") then
                        min3PGases(i)%NameOfSTQ(k) = gases(i)%NameOfSTQ(k)
                    else
                        call WriteLog("Alias used for : " // trim(gases(i)%NameOfSTQ(k)) // " to " // trim(min3PGases(i)%NameOfSTQ(k)))
                    end if
                end if
            end do
                
            call WriteLog("Convert gas: "// trim(min3PGases(i)%Name))
        
        end do
        
    end subroutine convert2Min3PGases
    
    !Convert to min3p minerals
    subroutine convert2Min3PMinerals
    
        integer :: i, j

        nMin3PMinerals = nMinerals

        if (allocated(min3PMinerals)) then
            deallocate(min3PMinerals)
        end if
        
        if (nMin3PMinerals > 0) then
            allocate(min3PMinerals(nMin3PMinerals))        
        end if
        
        call WriteLog("Start converting minerals")
        
        do i = 1, nMinerals

            !Convert name
            !Length of min3p species name is 12 maximum while toughreact species name is 20,
            !Deal this first if there exists the same name after truncation
            !   .........................................
            
            if(len_trim(minerals(i)%Name)>nMNLmin3P) then
                min3PMinerals(i)%Name = GetAliasFromName(minerals(i)%Name)
                if(trim(min3PMinerals(i)%Name) == "") then
                    min3PMinerals(i)%Name = minerals(i)%Name(1:nMNLmin3P)
                    call WriteLog("Name truncation: " // trim(minerals(i)%Name) // " to " // minerals(i)%Name(1:nMNLmin3P)) 
                    call AddNameTruncation(trim(minerals(i)%Name), minerals(i)%Name(1:nMNLmin3P))
                else
                    call WriteLog("Alias used for : " // trim(minerals(i)%Name) // " to " // trim(min3PMinerals(i)%Name))
                end if
            else
                min3PMinerals(i)%Name = GetAliasFromName(minerals(i)%Name)
                if(trim(min3PMinerals(i)%Name) == "") then
                    min3PMinerals(i)%Name = minerals(i)%Name
                else
                    call WriteLog("Alias used for : " // trim(minerals(i)%Name) // " to " // trim(min3PMinerals(i)%Name))
                end if
            end if
            
            min3PMinerals(i)%EnthalpyChange = 0.0d0
            !The second temperature in TOUGHREACT 25C is used. 
            !Should modify if 25C does not exist or in the different position.
            min3PMinerals(i)%AKLOG = -minerals(i)%AKLOG(2)   

            !Molecular weight
            min3PMinerals(i)%MWT = minerals(i)%MWT
            
            !Density
            min3PMinerals(i)%Density = minerals(i)%MWT / minerals(i)%VMIN

            !Number of components
            min3PMinerals(i)%NCP = min(minerals(i)%NCP,20)
                
            j = min3PMinerals(i)%NCP
               
            if (minerals(i)%NCP > 20) then
                call WriteLog("Warning:"// "Number of components exceed 20, only 20 components are converted.") 
                nWarnings = nWarnings + 1
            end if

            min3PMinerals(i)%STQ(1:j) = minerals(i)%STQ(1:j)
            
            do k = 1, j
                if(len_trim(minerals(i)%NameOfSTQ(k)) > nMNLmin3P) then
                    min3PMinerals(i)%NameOfSTQ(k) = GetAliasFromName(minerals(i)%NameOfSTQ(k))
                    if(trim(min3PMinerals(i)%NameOfSTQ(k)) == "") then
                        min3PMinerals(i)%NameOfSTQ(k) = minerals(i)%NameOfSTQ(k)(1:nMNLmin3P)     !no quotes
                        call WriteLog("Name truncation: " // trim(minerals(i)%NameOfSTQ(k)) &
                                    // " to " // trim(min3PMinerals(i)%NameOfSTQ(k))) 
                        call AddNameTruncation(trim(minerals(i)%NameOfSTQ(k)), minerals(i)%NameOfSTQ(k)(1:nMNLmin3P))
                    else
                        call WriteLog("Alias used for : " // trim(minerals(i)%NameOfSTQ(k)) // " to " // trim(min3PMinerals(i)%NameOfSTQ(k)))
                    end if
                else
                    min3PMinerals(i)%NameOfSTQ(k) = GetAliasFromName(minerals(i)%NameOfSTQ(k))
                    if(trim(min3PMinerals(i)%NameOfSTQ(k)) == "") then
                        min3PMinerals(i)%NameOfSTQ(k) = minerals(i)%NameOfSTQ(k)                  !no quotes
                    else
                        call WriteLog("Alias used for : " // trim(minerals(i)%NameOfSTQ(k)) // " to " // trim(min3PMinerals(i)%NameOfSTQ(k)))
                    end if
                end if
            end do
                
            call WriteLog("Convert minerals: "// trim(min3PMinerals(i)%Name))
        
        end do
        
    end subroutine convert2Min3PMinerals

end module convt_tough_min3p