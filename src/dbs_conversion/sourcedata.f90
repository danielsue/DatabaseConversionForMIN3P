module sourcedata

    integer, parameter                      :: nNameLength = 60

    !
    type TypePrimarySpecies
        
        character(nNameLength)   ::      Name    !name or chemical formula of aqueous basis species, in quotes (truncated after 20 characters), case-insensitive
        
        character(nNameLength)   ::      Formula !Formula, case-sensitive
        
        real            ::      A0      !Ion effective or hydrated radius used to compute the Debye-Hukel a0 parameter        
        real            ::      DHA     !Debye-Hukel a
        real            ::      DHB     !Debye-Hukel b 
        
        real            ::      Z       !the ion electric charge
        
        real            ::      MWT     !Molecular weight of the aqueous species (g/mol)
        
        real            ::      AlkFac = 0.0d0  !Alkalinity factor  
        
    end type TypePrimarySpecies   
     
    !
    type TypeDerivedSpecies
        
        character(nNameLength)   ::      Name    !chemical formula of secondary species, in quotes (truncated after 20 characters).
        character(nNameLength)   ::      Formula !Formula, case-sensitive
        
        real            ::      A0      !Ion effective or hydrated radius used to compute the Debye-Hukel a0 parameter
        real            ::      DHA     !Debye-Hukel a
        real            ::      DHB     !Debye-Hukel b 
        
        real            ::      Z       !the ion electric charge
        
        real            ::      MWT     !Molecular weight of the aqueous species (g/mol)
        
        integer         ::      NCP     !number of basis species defining the secondary species
        
        real            ::  STQ(20)     !stoichiometric coefficients of component (basis) species NAM included in the dissociation
                                        !reaction of the derived species (negative and positive values for reactants and products,
                                        !respectively). The derived species is always assumed to have a stoichiometric coefficient of -1.0,
                                        !which is not included in STQ.
                                                        
        
        character(nNameLength)      ::  NameOfSTQ(20)    !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).
                                                         
        real                        ::  AKLOG(20)        !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)         !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)         !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ?V
                                                         !(in cm3/mol) for the reaction as a function of temperature (average ?V over the pressure interval
                                                         !P0 to P), with ?V = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists

        logical                     ::  bRedoxReaction = .false.   !To seperate if it is complexation reaction or redox reaction
                                                                   !If NameOfSTQ contain "o2(aq)", it is redox reaction
                                                                   !...........IS IT CORRECT?.....................
                                                                     
        integer         ::      iAssociation = 1         !Indication of reaction direction, if == 1, association, else if == -1, dissociation
        real            ::      EnthalpyChange = 0.0d0     !Enthalpy change
            
    end type TypeDerivedSpecies   
                                             
    
    !
    type TypeMineral
        
        character(nNameLength)   ::      Name    !name or chemical formula of a mineral, in quotes (truncated after 20 characters).
        character(nNameLength)   ::      Formula !Formula, case-sensitive
        
        real            ::      MWT     !molecular weight (g/mol)
        
        real            ::      VMIN    !molar volume (cm3/mole)
        
        integer         ::      NCP     !the number of component species defining the mineral.
        
        real            ::     STQ(20)  !contains the stoichiometric coefficient of basis species NAM in the dissociation (hydrolysis)
                                        !reaction of the mineral (negative and positive values for reactants and products, respectively). The
                                        !mineral species is always assumed to have a stoichiometric coefficient of -1.0, which is not
                                        !included in STQ.
                                                    
        character(nNameLength)      ::  NameOfSTQ(20)    !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                         !basis species).                                              
                                                    
        real                        ::  AKLOG(20)        !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                         !temperature listed in record
                                                         
        real                        ::  AKCOE(5)         !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                         !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                         !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)         !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ?V
                                                         !(in cm3/mol) for the reaction as a function of temperature (average ?V over the pressure interval
                                                         !P0 to P), with ?V = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists 
        
        integer         ::      iAssociation = 1         !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation
        
        real            ::      EnthalpyChange = 0.0d0   !Enthalpy change
    
    end type TypeMineral
    
    !
    type TypeGas
        
        character(nNameLength)   ::      Name    !name or chemical formula of a gas species, in quotes (truncated after 20 characters).
        character(nNameLength)   ::      Formula !Formula, case-sensitive
        
        real            ::      MWT      !molecular weight (g/mol)
        
        real            ::      VMIN     !molecular volume
        
        real            ::      DMDIAM   !molecular diameter (m) used to calculate gas diffusion coefficient
        
        integer         ::      NCP      !the number of basis species defining the gas.
        
        real            ::      STQ(20)  !contains the stoichiometric coefficient of component species NAM in the dissociation reaction
                                         !of the gas (negative and positive values for reactants and products, respectively). The gas is
                                         !always assumed to have a stoichiometric coefficient of -1.0, which is not included in STQ.
                                                    
        character(nNameLength)      ::  NameOfSTQ(20)  !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                       !basis species).  
                                                         
        real                        ::  AKLOG(20)      !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                       !temperature listed in record
                                                         
        real                        ::  AKCOE(5)       !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                       !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                       !where Tk is absolute temperature (K), and ln stands for natural logarithm.
                                                         
        real                        ::  AKCOP(5)       !(optional) contains regression coefficients a, b, c, d, and e to calculate the volume change ?V
                                                       !(in cm3/mol) for the reaction as a function of temperature (average ?V over the pressure interval
                                                       !P0 to P), with ?V = a + b*Tk + c*Tk^2 + d/Tk + e/Tk^2.
                                                         
        logical                     ::  bAKCOP = .false. !Indicate if AKCOP exists     
        
        integer         ::      iAssociation = 1         !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation
        
        real            ::      EnthalpyChange = 0.0d0   !Enthalpy change
        
    end type TypeGas
    
    !
    type TypeSurfaceComplexes
        
        character(nNameLength)   ::      Name    !name or chemical formula of a gas species, in quotes (truncated after 20 characters).
        character(nNameLength)   ::      Formula !Formula, case-sensitive
        
        real            ::      Z          !electric charge of surface complex.
        
        integer         ::      NCP        !the number of basis species defining the surface complex.
        
        real            ::      STQ(20)    !contains the stoichiometric coefficient of component species NAM in the dissociation reaction
                                           !of the gas (negative and positive values for reactants and products, respectively). The gas is
                                           !always assumed to have a stoichiometric coefficient of -1.0, which is not included in STQ.
                                                    
        character(nNameLength)      ::  NameOfSTQ(20)  !name of the reactant or product, in quotes (truncated after 20 characters; must match one of the
                                                       !basis species).  
                                                         
        real                        ::  AKLOG(20)      !contains the equilibrium constants (log(K) in base 10) for the given reaction at each discrete
                                                       !temperature listed in record
                                                         
        real                        ::  AKCOE(5)       !contains regression coefficients a, b, c, d, and e to calculate log10(K) as a function of
                                                       !temperature (at a reference pressure P0) with log10(K)T,P0 = a*ln(Tk) + b + c*Tk + d/Tk + e/Tk^2,
                                                       !where Tk is absolute temperature (K), and ln stands for natural logarithm.  
                                                       
        integer         ::      iAssociation = 1       !Indication of reaction direction, if iAssociation = 1, association, else if iAssociation = -1, dissociation
        
    end type TypeSurfaceComplexes

    real,allocatable                         ::  temperature(:)         ! temperature points
    integer                                  ::  nTemperature = 0       ! number of temperature points
                                                                      
    type(TypePrimarySpecies), allocatable    :: species(:)              ! species
    integer                                  :: nSpecies = 0            ! number of species
                                                                      
    type(TypeDerivedSpecies), allocatable    :: aqueousSpecies(:)       ! aqueous species
    integer                                  :: nAqueousSpecies = 0     ! number of aqueous species
    integer                                  :: nRedoxReaction  = 0   
                                                                      
    type(TypeMineral), allocatable           :: minerals(:)             ! minerals
    integer                                  :: nMinerals = 0           ! number of minerals
                                                                      
    type(TypeGas), allocatable               :: gases(:)                ! gases
    integer                                  :: nGases = 0              ! number of gases    
    
    type(TypeSurfaceComplexes),allocatable   ::  surfaceComplexes(:)    ! surfaceComplexes
    integer                                  ::  nSurfaceComplexes = 0  ! number of surface complexes 

end module sourcedata