!module of geochemistry
module geochemistry

use logfile, only : WriteLog

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
        
         call WriteLog("Switch master variable h2(aq) for "//trim(spename))
         
         call WriteLog ("aklog before converting:")
         
         call WriteLog (ntemp, aklog(:ntemp))
         
         call WriteLog ("stq_o2aq:")
         
         call WriteLog (stq_o2aq)
         
         call WriteLog ("aklog_h2aq:")
         
         call WriteLog (ntemp, aklog_h2aq (:ntemp))         

         do i = 1, ntemp
            aklog(i) = aklog(i) + aklog_h2aq (i) * 2 * stq_o2aq            
         end do
         
         call WriteLog ("aklog after converting:")
         
         call WriteLog (ntemp, aklog(:ntemp))        
         

         stq_h2o= stq_h2o + 2.0 * stq_o2aq
         stq_o2aq = stq_o2aq - 1.0 * stq_o2aq
         stq_h2aq = stq_h2aq - 2.0 * stq_o2aq

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
         
         call WriteLog ("stq_h2aq:")
         
         call WriteLog (stq_h2aq)
         
         call WriteLog ("aklog_o2aq:")
         
         call WriteLog (ntemp, aklog_o2aq (:ntemp))   

         do i = 1, ntemp
            aklog(i) = aklog(i) + aklog_o2aq (i) * 0.5 * stq_h2aq            
         end do         
         
         call WriteLog ("aklog after converting:")
         
         call WriteLog (ntemp, aklog(:ntemp)) 

         stq_h2o= stq_h2o + 1.0 * stq_h2aq
         stq_o2aq = stq_o2aq - 0.5 * stq_h2aq
         stq_h2aq = stq_h2aq - 1.0 * stq_h2aq

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



end module geochemistry