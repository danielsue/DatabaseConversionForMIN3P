!****************************************************************************
!    Geochemitry database conversion program
!    Author: Danyang Su
!    Email: dsu@eos.ubc.ca; danyang.su@gmail.com
!****************************************************************************
!module of elements
module elements

    use logfile, only : WriteLog
    
    use sourcedata, only : nNameLength

    implicit none
    
    type ElementType
        
        character(nNameLength) :: Name        !Element name is case sensitive
        
        real :: MWT
        
    end type
    
    type(ElementType) :: elementList(109)
    type(ElementType), allocatable :: userDefinedElementList(:)
    
    contains
    
    !initialize element-mass list
    subroutine InitElementList
    
        implicit none
        
        integer :: i
        
        elementList(1)%Name="H"
        elementList(2)%Name="He"
        elementList(3)%Name="Li"
        elementList(4)%Name="Be"
        elementList(5)%Name="B"
        elementList(6)%Name="C"
        elementList(7)%Name="N"
        elementList(8)%Name="O"
        elementList(9)%Name="F"
        elementList(10)%Name="Ne"
        elementList(11)%Name="Na"
        elementList(12)%Name="Mg"
        elementList(13)%Name="Al"
        elementList(14)%Name="Si"
        elementList(15)%Name="P"
        elementList(16)%Name="S"
        elementList(17)%Name="Cl"
        elementList(18)%Name="Ar"
        elementList(19)%Name="K"
        elementList(20)%Name="Ca"
        elementList(21)%Name="Sc"
        elementList(22)%Name="Ti"
        elementList(23)%Name="V"
        elementList(24)%Name="Cr"
        elementList(25)%Name="Mn"
        elementList(26)%Name="Fe"
        elementList(27)%Name="Co"
        elementList(28)%Name="Ni"
        elementList(29)%Name="Cu"
        elementList(30)%Name="Zn"
        elementList(31)%Name="Ga"
        elementList(32)%Name="Ge"
        elementList(33)%Name="As"
        elementList(34)%Name="Se"
        elementList(35)%Name="Br"
        elementList(36)%Name="Kr"
        elementList(37)%Name="Rb"
        elementList(38)%Name="Sr"
        elementList(39)%Name="Y"
        elementList(40)%Name="Zr"
        elementList(41)%Name="Nb"
        elementList(42)%Name="Mo"
        elementList(43)%Name="Tc"
        elementList(44)%Name="Ru"
        elementList(45)%Name="Rh"
        elementList(46)%Name="Pd"
        elementList(47)%Name="Ag"
        elementList(48)%Name="Cd"
        elementList(49)%Name="In"
        elementList(50)%Name="Sn"
        elementList(51)%Name="Sb"
        elementList(52)%Name="Te"
        elementList(53)%Name="I"
        elementList(54)%Name="Xe"
        elementList(55)%Name="Cs"
        elementList(56)%Name="Ba"
        elementList(57)%Name="La"
        elementList(58)%Name="Ce"
        elementList(59)%Name="Pr"
        elementList(60)%Name="Nd"
        elementList(61)%Name="Pm"
        elementList(62)%Name="Sm"
        elementList(63)%Name="Eu"
        elementList(64)%Name="Gd"
        elementList(65)%Name="Tb"
        elementList(66)%Name="Dy"
        elementList(67)%Name="Ho"
        elementList(68)%Name="Er"
        elementList(69)%Name="Tm"
        elementList(70)%Name="Yb"
        elementList(71)%Name="Lu"
        elementList(72)%Name="Hf"
        elementList(73)%Name="Ta"
        elementList(74)%Name="W"
        elementList(75)%Name="Re"
        elementList(76)%Name="Os"
        elementList(77)%Name="Ir"
        elementList(78)%Name="Pt"
        elementList(79)%Name="Au"
        elementList(80)%Name="Hg"
        elementList(81)%Name="Tl"
        elementList(82)%Name="Pb"
        elementList(83)%Name="Bi"
        elementList(84)%Name="Po"
        elementList(85)%Name="At"
        elementList(86)%Name="Rn"
        elementList(87)%Name="Fr"
        elementList(88)%Name="Ra"
        elementList(89)%Name="Ac"
        elementList(90)%Name="Th"
        elementList(91)%Name="Pa"
        elementList(92)%Name="U"
        elementList(93)%Name="Np"
        elementList(94)%Name="Pu"
        elementList(95)%Name="Am"
        elementList(96)%Name="Cm"
        elementList(97)%Name="Bk"
        elementList(98)%Name="Cf"
        elementList(99)%Name="Es"
        elementList(100)%Name="Fm"
        elementList(101)%Name="Md"
        elementList(102)%Name="No"
        elementList(103)%Name="Lr"
        elementList(104)%Name="Rf"
        elementList(105)%Name="Db"
        elementList(106)%Name="Sg"
        elementList(107)%Name="Bh"
        elementList(108)%Name="Hs"
        elementList(109)%Name="Mt"
        
        elementList(1)%MWT=1.0079
        elementList(2)%MWT=4.0026
        elementList(3)%MWT=6.941
        elementList(4)%MWT=9.0122
        elementList(5)%MWT=10.811
        elementList(6)%MWT=12.0107
        elementList(7)%MWT=14.0067
        elementList(8)%MWT=15.9994
        elementList(9)%MWT=18.9984
        elementList(10)%MWT=20.1797
        elementList(11)%MWT=22.9897
        elementList(12)%MWT=24.305
        elementList(13)%MWT=26.9815
        elementList(14)%MWT=28.0855
        elementList(15)%MWT=30.9738
        elementList(16)%MWT=32.065
        elementList(17)%MWT=35.453
        elementList(18)%MWT=39.948
        elementList(19)%MWT=39.0983
        elementList(20)%MWT=40.078
        elementList(21)%MWT=44.9559
        elementList(22)%MWT=47.867
        elementList(23)%MWT=50.9415
        elementList(24)%MWT=51.9961
        elementList(25)%MWT=54.938
        elementList(26)%MWT=55.845
        elementList(27)%MWT=58.9332
        elementList(28)%MWT=58.6934
        elementList(29)%MWT=63.546
        elementList(30)%MWT=65.39
        elementList(31)%MWT=69.723
        elementList(32)%MWT=72.64
        elementList(33)%MWT=74.9216
        elementList(34)%MWT=78.96
        elementList(35)%MWT=79.904
        elementList(36)%MWT=83.8
        elementList(37)%MWT=85.4678
        elementList(38)%MWT=87.62
        elementList(39)%MWT=88.9059
        elementList(40)%MWT=91.224
        elementList(41)%MWT=92.9064
        elementList(42)%MWT=95.94
        elementList(43)%MWT=98
        elementList(44)%MWT=101.07
        elementList(45)%MWT=102.9055
        elementList(46)%MWT=106.42
        elementList(47)%MWT=107.8682
        elementList(48)%MWT=112.411
        elementList(49)%MWT=114.818
        elementList(50)%MWT=118.71
        elementList(51)%MWT=121.76
        elementList(52)%MWT=127.6
        elementList(53)%MWT=126.9045
        elementList(54)%MWT=131.293
        elementList(55)%MWT=132.9055
        elementList(56)%MWT=137.327
        elementList(57)%MWT=138.9055
        elementList(58)%MWT=140.116
        elementList(59)%MWT=140.9077
        elementList(60)%MWT=144.24
        elementList(61)%MWT=145
        elementList(62)%MWT=150.36
        elementList(63)%MWT=151.964
        elementList(64)%MWT=157.25
        elementList(65)%MWT=158.9253
        elementList(66)%MWT=162.5
        elementList(67)%MWT=164.9303
        elementList(68)%MWT=167.259
        elementList(69)%MWT=168.9342
        elementList(70)%MWT=173.04
        elementList(71)%MWT=174.967
        elementList(72)%MWT=178.49
        elementList(73)%MWT=180.9479
        elementList(74)%MWT=183.84
        elementList(75)%MWT=186.207
        elementList(76)%MWT=190.23
        elementList(77)%MWT=192.217
        elementList(78)%MWT=195.078
        elementList(79)%MWT=196.9665
        elementList(80)%MWT=200.59
        elementList(81)%MWT=204.3833
        elementList(82)%MWT=207.2
        elementList(83)%MWT=208.9804
        elementList(84)%MWT=209
        elementList(85)%MWT=210
        elementList(86)%MWT=222
        elementList(87)%MWT=223
        elementList(88)%MWT=226
        elementList(89)%MWT=227
        elementList(90)%MWT=232.0381
        elementList(91)%MWT=231.0359
        elementList(92)%MWT=238.0289
        elementList(93)%MWT=237
        elementList(94)%MWT=244
        elementList(95)%MWT=243
        elementList(96)%MWT=247
        elementList(97)%MWT=247
        elementList(98)%MWT=251
        elementList(99)%MWT=252
        elementList(100)%MWT=257
        elementList(101)%MWT=258
        elementList(102)%MWT=259
        elementList(103)%MWT=262
        elementList(104)%MWT=261
        elementList(105)%MWT=262
        elementList(106)%MWT=266
        elementList(107)%MWT=264
        elementList(108)%MWT=277
        elementList(109)%MWT=268
        
        do i = 1, size(elementList, 1)
            call WriteLog(elementList(i)%Name, elementList(i)%MWT)
        end do
    
     end subroutine
    
     !get the atomic mass
     function getMWT(string) result(MWT)
     
        implicit none
        
        character(*), intent(in) :: string
        real :: MWT
        
        integer :: i, n
        
        !seperate string
        n = size(elementList, 1)
        MWT = -1.0d0
        do i = 1, n
            if(trim(elementList(i)%Name) == trim(string)) then
                MWT = elementList(i)%MWT
                exit
            end if
        end do        
     
     end function getMWT
     
     
     !get the atomic mass
     function getUserDefinedMWT(string) result(MWT)
     
        implicit none
        
        character(*), intent(in) :: string
        real :: MWT
        
        integer :: i, n
        
        !seperate string
        n = size(userDefinedElementList, 1)
        MWT = -1.0d0
        do i = 1, n
            if(trim(userDefinedElementList(i)%Name) == trim(string)) then
                MWT = userDefinedElementList(i)%MWT
                exit
            end if
        end do        
     
     end function getUserDefinedMWT
     
     !check if the specified element is included in the element list
     function bIncludeElement(string) result(bFlag)
     
        implicit none
        character(*), intent(in) :: string
        logical :: bFlag
        integer :: i, n
        
        !seperate string
        n = size(elementList, 1) 
        bFlag = .false.
        
        do i = 1, n
            if(trim(elementList(i)%Name) == trim(string)) then
                bFlag = .true.
                exit
            end if
        end do 
     
     end function bIncludeElement
     
     !check if the user-defined element name is included in the specified string
     function bMatchElement(string, nlen, mwt) result(bFlag)
     
        implicit none
        character(*), intent(in) :: string
        integer, intent(out) :: nlen
        real, intent(out) :: mwt
        logical :: bFlag
        integer :: i, n
        
        !seperate string
        n = size(userDefinedElementList, 1) 
        bFlag = .false.
        nlen = 0
        mwt = -1.0d0
        
        do i = 1, n
            if(index(string, trim(userDefinedElementList(i)%Name)) == 1) then
                bFlag = .true.
                nlen = len_trim(userDefinedElementList(i)%Name)
                mwt = userDefinedElementList(i)%MWT
                exit
            end if
        end do 
     
     end function bMatchElement
     
     !Reorder user defined element in decreasing order of name length .
     !Make name best match
     subroutine reorderUserDefinedElements()
     
        implicit none
        integer :: i , j, n
        type(ElementType) :: switchElement
        n = size(userDefinedElementList, 1)
        
        do i =1, n - 1
            do j = i + 1, n
                if(len_trim(userDefinedElementList(j)%Name) > len_trim(userDefinedElementList(i)%Name)) then
                    switchElement = userDefinedElementList(i)
                    userDefinedElementList(i) = userDefinedElementList(j)
                    userDefinedElementList(j) = switchElement
                end if
            end do
        end do
        
     
     end subroutine reorderUserDefinedElements
     

end module elements