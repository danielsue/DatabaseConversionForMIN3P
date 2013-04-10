!Module of global parameters
!The parameters are by default public and can be accessed from outside of the module 
    
module global
    
    implicit none
    
    integer, parameter  ::  iUnit                   =   10 
    integer, parameter  ::  iUnitLog                =   iUnit + 1                   !logfile
    integer, parameter  ::  iUnitDbsTR              =   iUnitLog + 1                !database of toughreact
    integer, parameter  ::  iUnitDbsMin3PComp       =   iUnitDbsTR + 1              !comp.dbs for min3p
    integer, parameter  ::  iUnitDbsMin3PComplex    =   iUnitDbsMin3PComp + 1       !complex.dbs for min3p
    integer, parameter  ::  iUnitDbsMin3PGases      =   iUnitDbsMin3PComplex + 1    !gases.dbs for min3p
    integer, parameter  ::  iUnitDbsMin3PRedox      =   iUnitDbsMin3PGases + 1      !redox.dbs for min3p
    integer, parameter  ::  iUnitDbsMin3PMinerals   =   iUnitDbsMin3PRedox + 1      !mineral.dbs for min3p
    integer, parameter  ::  iUnitDbsAlias           =   iUnitDbsMin3PMinerals + 1   !alias.dbs
    integer, parameter  ::  iUnitNameTruncation     =   iUnitDbsAlias + 1           !name_truncation.txt
    integer, parameter  ::  iUnitInp                =   iUnitNameTruncation + 1     !input file for dbs conversion
    integer, parameter  ::  iUnitDbsCF              =   iUnitInp + 1                !database of crunch flow
    integer, parameter  ::  iUnitDbsPhreeqc         =   iUnitDbsCF + 1              !database of Phreeqc
    
    logical             ::  bOpenLog                =   .false.
    logical             ::  bOpenDbsTR              =   .false. 
    logical             ::  bOpenDbsMin3PComp       =   .false.
    logical             ::  bOpenDbsMin3PComplex    =   .false.
    logical             ::  bOpenDbsMin3PGases      =   .false.
    logical             ::  bOpenDbsMin3PRedox      =   .false.
    logical             ::  bOpenDbsMin3PMinerals   =   .false.    
    logical             ::  bOpenDbsAlias           =   .false.
    logical             ::  bOpenNameTruncation     =   .false.
    logical             ::  bOpenInp                =   .false.
    logical             ::  bOpenDbsCF              =   .false.
    logical             ::  bOpenDbsPhreeqc         =   .false.
    
 
contains

    !Error handling
    subroutine ErrorHandling
    
        implicit none        
        
        call CloseAllFiles
        
        stop
    
    end subroutine ErrorHandling
    
    !Close all files
    subroutine CloseAllFiles
    
        implicit none
        
        if (bOpenLog) then
            close(iUnitLog)
        end if
        
        if (bOpenDbsTR) then
            close(iUnitDbsTR)
        end if
        
        if (bOpenDbsMin3PComp) then
            close(iUnitDbsMin3PComp)
        end if
        
        if (bOpenDbsMin3PComplex) then
            close(iUnitDbsMin3PComplex)
        end if
        
        if (bOpenDbsMin3PGases) then
            close(iUnitDbsMin3PGases)
        end if
        
        if (bOpenDbsMin3PRedox) then
            close(iUnitDbsMin3PRedox)
        end if
        
        if (bOpenDbsMin3PMinerals) then
            close(iUnitDbsMin3PMinerals)
        end if
        
        if (bOpenDbsAlias) then
            close(iUnitDbsAlias)
        end if
        
        if (bOpenNameTruncation) then
            close(iUnitNameTruncation)
        end if

        if (bOpenInp) then
            close(iUnitInp)
        end if
        
        if (bOpenDbsCF) then
            close(iUnitDbsCF)
        end if
        
        if (bOpenDbsPhreeqc) then
            close(iUnitDbsPhreeqc)
        end if
    
    end subroutine CloseAllFiles
    
end module global