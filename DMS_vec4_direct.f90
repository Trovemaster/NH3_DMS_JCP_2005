      program DIPOLE_XY3
!    

!    This program computes the dipole moment of XY3 molecules using the FUNCTION DMS, see below.
!    The dipole moment \vec{mu} is given by its projections \bar{mu}_i (i=1,2,3)
!    on the three molecular bonds X-Y1, X-Y2, and X-Y3:
!    \bar{mu}_i = (\vec{mu}\cdot \vec{r_i}/|r_i| ),
!    where (a\cdot b) denotes a scalar product of the vectors a and b and
!    r_i is a vector, coinciding with the bond X-Y_i and pointing to the atom Y_i.
!    For the details on the definition of \bar{mu}_i see the paper.
!    Although the function DMS is based the geometrically defined coordinates
!    r1,r2,r3,alpha23,alpha13,alpha12, one can also use the Cartesian coordinates in order
!    to define the molecular geometry. Thus, in order to run this program one needs
!    an input file with a list of molecular geometries, see the input file examples.
!    The molecular geometries can be given either as nine Cartesian coordinates
!    (coords='CARTESIAN') or 6 geometrically defined coordinates. The Cartesian coordinates
!    must be given in the following order (e.g. for NH3):
!    x_N, y_N, z_N, x_H1, y_H1, z_H1, x_H2, y_H2, z_H2, x_H2, y_H3, z_H3. 
!    For the geometrically defined coordinates please use the following order
!    r1,r2,r3,alpha23,alpha13,alpha12. The Cartesian coordinates are assumed in Angstroms,
!    as well as the bond lengths. The inter-bond angles must be given in degrees.
!    The line positions in the input file are important.
!    The input files contain weight factors for each dipole moment parameter, which 
!    are coming from the fitting procedure are redundant here. 
!    The relative positions of the weights are important in the input
!    file, as well as positions of the parameter label and value. They must stay  
!    within the right sections, see either the input file examples or
!    the corresponding READ command. We hope the program is simple and self-explanatory.
!    The Lapack routine dgelss is required for coords='CARTESIAN'.
!
!    15.09.2008: The instability of the calculations of the dipole moment in the Cartesian coordinates 
!    representation appearing for geometries close to the planarity have been fixed by introducing 
!    an extra projection  of the dipole moment to the trisector-vector 'b' (see the code).
!    Now the three Cartesinan componets are computed from the four projections (to r1,r2,r3,b) of the dipole 
!    moment by solving a 4x3 system of linear equations with the Lapack rouitne  dgelss. 
!    The input file thus also contains parameters defining the projection of DMS to 'b' 
!    which did not appear in J.Chem.Phys. 122, 104317 (2005).
!    

! Here we go:
!
! Parameters:  
 ! 
 ! number of parameters  
 !
 integer,parameter          ::  parmax1=125,parmax2 = 112, parmax = parmax1+parmax2
 !
 ! where the output goes to 
 ! 
 integer,parameter          :: f_inp=5, f_out=6
 !
 ! correspondence between the dipole terms and expansion orders
 !
 integer,parameter          ::  DMS_orders(0:4) = (/1,11,24,56,126/)
 !
 ! constants 
 ! 
 double precision,parameter ::pi=3.141592653589793
 !
 !
 integer           :: npts,i
 ! 
 character(len=300):: title(4),longlabel ! input title of the job from the first four input lines 
 character(len=10) :: parnam(parmax) ! Names of the parameters, as they appear in the input file  
 character(len=10) :: label
 !
 double precision  :: local(6)
 integer           :: alloc
 logical           :: yes
 integer           :: ivartmp,ivar(parmax),order,rank0,ierror,ix,iy,iz,n,x_singular
 double precision  :: paramtmp,DMS,charge,EC,Xshift,VCPRM,mu(3),DMS_A
 double precision  :: r1,r2,r3,a1,a2,a3,dip(4),edip(4,1),wspace(50),tmat(4,3),dip_cart(3),dip_A
 double precision  :: x(4),y(4),z(4),x0,y0,z0,tsing(3),bmat(3),f_t,amat(3,3),mass(4),total_mass
 double precision  :: v12(3),v23(3),v31(3),n3(3)
 double precision  :: levich(3,3,3)    ! epsil - antisymmetric tensor
 character(len=80) :: coords
 !
 ! some matrices, to be allocated  
 !
 double precision, allocatable :: param(:)
 !
 ! Here we go!
 !
 ! Array's allocation:
 !
 allocate (param(parmax),stat=alloc)
 if (alloc/=0) then 
  write(f_out,"('parmax - out of memory')") 
 stop
 endif
! 
!  ELEMENTARY CHARGE IN COULOMBS
!
      EC=1.60217733D-19
      VCPRM=8.854187817D-12
!
! ***** ELEMENTARY CHARGE IN ELECTROSTATIC UNITS (E.S.U.)
!
      EC=EC/SQRT(4.0D-9*PI*VCPRM)
!
! ***** ELEMENTARY CHARGE IN DEBYE/ANGSTROM =
!
!       1.0E8 * DEBYE/CM = 1.0E-10 * E.S.U. * CM/CM
!
      EC=1.0D10*EC
 !
 levich = 0 
 levich(1,2,3) = 1.0d0
 levich(1,3,2) =-1.0d0
 levich(2,1,3) =-1.0d0
 levich(2,3,1) = 1.0d0
 levich(3,1,2) = 1.0d0
 levich(3,2,1) =-1.0d0

 !
 ! Input data from file: 
 !
 call skiplines(f_inp,1)
 !
 !  input the job title 
 !
 do i=1,4
   read  (f_inp,"(a80)") title(i)
 enddo
 !
 !  output the job title 
 !
 write (f_out,"(3(132('*')/),16('*'),100x,16('*'))") 
 write (f_out,"(4(16('*'),10x,a80,10x,16('*')/),16('*'),100x,16('*')/3(132('*')/))") (title(i), i=1,4)
 !
 call skiplines(f_inp,3)
 !
 ! There are two options availible for the input data:
 ! 1) coords = 'CARTESIAN', when the cartesian coordinates are used. They must 
 !    be given in the following order:
 !    x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3 
 !    where x0,y0,z0 are the coordinates of the central atom (e.g. N) 
 !          x1,y1,z1 are the coordinates of H1
 !          x2,y2,z2 are the coordinates of H2
 !          x3,y3,z3 are the coordinates of H3
 ! 2) coords = 'INTERNAL', for the geometrically defined coordinates r1,r2,r3,alpha23,alpha13,alpah12 
 !    The bond length must be given in Angstrom, bond angles in degrees.
 !
 read (f_inp,"(a80)") coords
 !
 coords = adjustl(coords)
 !
 call skiplines(f_inp,3)
 !
 ! Molecular charge
 !
 read (f_inp,*) charge
 write (f_out,"(f10.0,' <= Molecular charge is ')") charge

 !
 call skiplines(f_inp,4)


 !
 ! input masses and re equilibrium constant
 !
 read (f_inp,*) mass(:)
 total_mass = sum(mass(:)) 

 !
 ! Input for the dipole moment parameters 
 !
 call skiplines(f_inp,3)
 do i = 1,parmax
   read (f_inp,"(a300)") longlabel
   read (longlabel,"(a10,i4,d18.8)") label,ivartmp,paramtmp
   if (ivartmp/=-1) then 
     parnam(i) = label 
     ivar(i)   = ivartmp
     param(i)  = paramtmp
   else 
     write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i,parmax
     stop 'Too few parameters'
   endif 
 enddo
 ! 
 read (f_inp,"(a300)") longlabel
 read (longlabel,"(a10,i4,d18.8)") label,ivartmp,paramtmp
 if (ivartmp/=-1) then 
   write(f_out,"('Wrong number of parameters (',I6,'), has to be (',I6,')')") i,parmax
   stop 'Too many parameters'
 endif 
 !
 ! Output of the potential parameters 
 ! 
 write (f_out,"(//' Dipole moment expansion coefficients:')")
 do i = 1,int(parmax/3)*3,3 
 write (f_out,"(3('  ',a9,'=',f16.3))") parnam(i),param(i),     &
                                        parnam(i+1),param(i+1), &
                                        parnam(i+2),param(i+2)
 enddo 
 write (f_out,"(3('  ',a9,'=',f16.3))") (parnam(int(parmax/3)*3+i), &
                                         param(int(parmax/3)*3+i),i=1,mod(parmax,3))
 !
 ! Check ivar and count the maximal order used 
 !
 order = 0
 do i=1,4
 ivartmp = sum(ivar(DMS_orders(i-1):DMS_orders(i)-1))
 if (ivartmp/=0) then
    order = order +float(1)
 endif 
 enddo 
 !
 select case (trim(coords))
 !
 case default
   !
   write (f_out,"('Bad coordinates type',a12)") trim(coords)
   stop 'bad coodrs'
   !
 case ('CARTESIAN','cartesian') 
   !
   write (f_out,"(/'    r1        r2        r3       alpha1    alpha2    alpha3   ',&
                   '    Mu_x          Mu_y           Mu_z            Mu_x(ab in)   Mu2_y(ab in)  Mu3_z(ab in)',&
                   '   dMu_x         dMu_y          dMu_z            Mu1           Mu2           Mu3')") 
   !
 case ('INTERNAL','internal')
   !
   write (f_out,"(/'    r1        r2        r3       alpha1    alpha2    alpha3       Mu1           Mu2           Mu3')") 
   ! 
 end select 
 !
 npts = 0
 yes = .true. 
 do while( yes )
    !
    npts=npts+1
    !
    read (f_inp,"(a300)") longlabel
    label  = adjustl(longlabel)
    if ( trim(longlabel)/='end'.and.trim(longlabel)/='END') then
       !
       select case (trim(coords))
       !
       case default
         !
         write (f_out,"('Bad coordinates type',a12)") trim(coords)
         stop 'bad coodrs'
         !
       case ('CARTESIAN','cartesian') 
         !
         !
         read (longlabel,*) x(4),y(4),z(4),x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),mu(1),mu(2),mu(3)
         !
         x0  =(sum(x(:)*mass(:)))/total_mass
         y0  =(sum(y(:)*mass(:)))/total_mass
         z0  =(sum(z(:)*mass(:)))/total_mass
         !
         x = x - x(4)
         y = y - y(4)
         z = z - z(4)
         !
         r1 = sqrt(x(1)**2+y(1)**2+z(1)**2)
         r2 = sqrt(x(2)**2+y(2)**2+z(2)**2)
         r3 = sqrt(x(3)**2+y(3)**2+z(3)**2)
         !
         a1 = acos((x(2)*x(3)+y(2)*y(3)+z(2)*z(3))/(r2*r3))
         a2 = acos((x(1)*x(3)+y(1)*y(3)+z(1)*z(3))/(r1*r3))
         a3 = acos((x(1)*x(2)+y(1)*y(2)+z(1)*z(2))/(r1*r2))
         !
         local(1)=r1
         local(2)=r2
         local(3)=r3
         !
         ! Running interbond angles
         !
         local(4)=a1
         local(5)=a2
         local(6)=a3
         ! 
         ! Value of the dipole moment fucntion at the current geometry 
         !
         dip(1) = DMS(1,order,parmax1,param(1:parmax1),local)
         dip(2) = DMS(2,order,parmax1,param(1:parmax1),local)
         dip(3) = DMS(3,order,parmax1,param(1:parmax1),local)
         dip(4) = DMS_A(parmax2,param(parmax1+1:parmax),local)
         !
         edip(1:4,1)= dip(1:4)
         !
         tmat(1,1) = x(1)/r1
         tmat(1,2) = y(1)/r1
         tmat(1,3) = z(1)/r1
         !
         tmat(2,1) = x(2)/r2
         tmat(2,2) = y(2)/r2
         tmat(2,3) = z(2)/r2
         !
         tmat(3,1) = x(3)/r3
         tmat(3,2) = y(3)/r3
         tmat(3,3) = z(3)/r3
         !
         call vector_product(tmat(1,:)/r1,tmat(2,:)/r2,v12)
         call vector_product(tmat(2,:)/r2,tmat(3,:)/r3,v23)
         call vector_product(tmat(3,:)/r3,tmat(1,:)/r1,v31)
         !
         n3 = v12 + v23 + v31
         !
         n3 = n3/sqrt(sum(n3(:)**2))
         !
         tmat(4, :) = n3(:)
         !
         lspace = size(wspace)
         !
         call dgelss(4,3,1,tmat,4,edip,4,tsing,1.D-16,rank0,wspace,size(wspace),ierror)
         !
         dip_cart(1:3) = edip(1:3,1) + EC*x0*charge
         !
         write (f_out,"(3f10.6,3f10.4,3(3f14.8,5x,3f14.8))") & 
              r1,r2,r3,a1*180.d0/pi,a2*180.d0/pi,a3*180.d0/pi,&
              dip_cart(1:3),mu(1:3),mu(1:3)-dip_cart(1:3),dip(1),dip(2),dip(3)
         !
         continue
         !
         !
       case ('INTERNAL','internal')
         !
         read (longlabel,*) r1,r2,r3,a1,a2,a3
         !
         a1=a1*pi/180.d0
         a2=a2*pi/180.d0
         a3=a3*pi/180.d0
         !
         local(1)=r1
         local(2)=r2
         local(3)=r3
         !
         ! Running interbond angles
         !
         local(4)=a1
         local(5)=a2
         local(6)=a3
         ! 
         ! Value of the dipole moment fucntion at the current geometry 
         !
         dip(1) = DMS(1,order,parmax1,param(1:parmax1),local)
         dip(2) = DMS(2,order,parmax1,param(1:parmax1),local)
         dip(3) = DMS(3,order,parmax1,param(1:parmax1),local)
         !
         !dip = dip - EC*Xshift*charge
         !
         write (f_out,"(3f10.6,3f10.4,3f14.8)") & 
                r1,r2,r3,a1*180.d0/pi,a2*180.d0/pi,a3*180.d0/pi,&
                dip(1),dip(2),dip(3)
         !
       end select 
       !
       !
      else
      !
      yes = .false.
      !
    endif
    !
 enddo 
 !
 ! Number of ab initio points  
 npts=npts-1
 ! 
end program DIPOLE_XY3


double precision function DMS(ix,order,parmax,param,local)

 integer,intent(in)          ::  ix,parmax,order
 double precision,intent(in) ::  param(parmax)
 double precision,intent(in) ::  local(6)

 double precision         ::  r14,r24,r34,alpha1,alpha2,alpha3
 double precision         ::  xi1,xi2,xi3,xi4,xi5,xi6
 double precision         ::  alphae

 double precision         ::  alphaedg

 double precision         ::  re14 ,beta,t1,t2,t3,t4

 double precision          &
    F5,F4,F3,F1,F56,F55,F46,F44,                                     &
    F36,F35,F34,F33,F23,F16,F14,F13,                                 &
    F11,F556,F555,F466,F456,F445,F444,F344,                          &
    F334,F333,F266,F256,F255,F246,F245,F235,                         &
    F234,F226,F225,F223,F156,F155,F146,F144,                         &
    F136,F135,F133,F124,F123,F115,F114,F112,                         &
    F111,F6666,F5666,F5566,F4556,F4555,F4466,F4456,                  &
    F4445,F4444,F3666,F3566,F3556,F3555,F3445,F3335,                 &
    F2466,F2456,F2455,F2445,F2444,F2366,F2356,F2345,                 &
    F2344,F2336,F2335,F2334,F2333,F2266,F2256,F2255,                 &
    F2246,F2245,F2244,F2233,F2225,F2224,F2222,F1666,                 &
    F1566,F1466,F1456,F1445,F1444,F1366,F1355,F1335,                 &
    F1334,F1256,F1246,F1245,F1244,F1235,F1234,F1233,                 &
    F1225,F1222,F1156,F1155,F1146,F1144,F1136,F1126,                 &
    F1124,F1123,F1122,F1115,F1114,F1112,F1111

  double precision s1,s2,s3,s4,s5,pi,gamma,delta,mu0
  double precision ta1,ta2,ta3


  !-------------------------------


    pi=3.141592653589793

    alphaedg     = param( 1)
    alphae=alphaedg*pi/180.d0

    re14      =  param(  2)
    beta      =  param(  3)**2
    gamma     =  param(  4)**2
    delta     =  param(  5)
    mu0       =  param(  6)
    F1      =  param(    7)
    F3      =  param(    8)
    F4      =  param(    9)
    F5      =  param(   10)

   if (order>=2) then
      F11     =  param(   11)
      F13     =  param(   12)
      F14     =  param(   13)
      F16     =  param(   14)
      F23     =  param(   15)
      F33     =  param(   16)
      F34     =  param(   17)
      F35     =  param(   18)
      F36     =  param(   19)
      F44     =  param(   20)
      F46     =  param(   21)
      F55     =  param(   22)
      F56     =  param(   23)
   endif
   if (order>=3) then
      F111    =  param(   24)
      F112    =  param(   25)
      F114    =  param(   26)
      F115    =  param(   27)
      F123    =  param(   28)
      F124    =  param(   29)
      F133    =  param(   30)
      F135    =  param(   31)
      F136    =  param(   32)
      F144    =  param(   33)
      F146    =  param(   34)
      F155    =  param(   35)
      F156    =  param(   36)
      F223    =  param(   37)
      F225    =  param(   38)
      F226    =  param(   39)
      F234    =  param(   40)
      F235    =  param(   41)
      F245    =  param(   42)
      F246    =  param(   43)
      F255    =  param(   44)
      F256    =  param(   45)
      F266    =  param(   46)
      F333    =  param(   47)
      F334    =  param(   48)
      F344    =  param(   49)
      F444    =  param(   50)
      F445    =  param(   51)
      F456    =  param(   52)
      F466    =  param(   53)
      F555    =  param(   54)
      F556    =  param(   55)
   endif
   if (order>=4) then
      F1111   =  param(   56)
      F1112   =  param(   57)
      F1114   =  param(   58)
      F1115   =  param(   59)
      F1122   =  param(   60)
      F1123   =  param(   61)
      F1124   =  param(   62)
      F1126   =  param(   63)
      F1136   =  param(   64)
      F1144   =  param(   65)
      F1146   =  param(   66)
      F1155   =  param(   67)
      F1156   =  param(   68)
      F1222   =  param(   69)
      F1225   =  param(   70)
      F1233   =  param(   71)
      F1234   =  param(   72)
      F1235   =  param(   73)
      F1244   =  param(   74)
      F1245   =  param(   75)
      F1246   =  param(   76)
      F1256   =  param(   77)
      F1334   =  param(   78)
      F1335   =  param(   79)
      F1355   =  param(   80)
      F1366   =  param(   81)
      F1444   =  param(   82)
      F1445   =  param(   83)
      F1456   =  param(   84)
      F1466   =  param(   85)
      F1566   =  param(   86)
      F1666   =  param(   87)
      F2222   =  param(   88)
      F2224   =  param(   89)
      F2225   =  param(   90)
      F2233   =  param(   91)
      F2244   =  param(   92)
      F2245   =  param(   93)
      F2246   =  param(   94)
      F2255   =  param(   95)
      F2256   =  param(   96)
      F2266   =  param(   97)
      F2333   =  param(   98)
      F2334   =  param(   99)
      F2335   =  param(  100)
      F2336   =  param(  101)
      F2344   =  param(  102)
      F2345   =  param(  103)
      F2356   =  param(  104)
      F2366   =  param(  105)
      F2444   =  param(  106)
      F2445   =  param(  107)
      F2455   =  param(  108)
      F2456   =  param(  109)
      F2466   =  param(  110)
      F3335   =  param(  111)
      F3445   =  param(  112)
      F3555   =  param(  113)
      F3556   =  param(  114)
      F3566   =  param(  115)
      F3666   =  param(  116)
      F4444   =  param(  117)
      F4445   =  param(  118)
      F4456   =  param(  119)
      F4466   =  param(  120)
      F4555   =  param(  121)
      F4556   =  param(  122)
      F5566   =  param(  123)
      F5666   =  param(  124)
      F6666   =  param(  125)
   endif


   r14    = local(1) ; r24    = local(2) ; r34    = local(3) 
   alpha1 = local(4) ; alpha2 = local(5) ; alpha3 = local(6)

     !
     ta1 = cos(alpha1)-cos(alphae)
     ta2 = cos(alpha2)-cos(alphae)
     ta3 = cos(alpha3)-cos(alphae)
     !

     select case ( ix )
     case default
       write (6,"(' dip. order component',i)") ix
       stop 'dip. order component'
     case (1)
       xi1=(r14-re14) *exp(-beta*(r14-re14)**2)
       xi2=(r24-re14) *exp(-beta*(r24-re14)**2)
       xi3=(r34-re14) *exp(-beta*(r34-re14)**2)

       xi4=ta1
       xi5=ta2
       xi6=ta3
       !
     case (2)
       xi1=(r24-re14) *exp(-beta*(r24-re14)**2)
       xi2=(r34-re14) *exp(-beta*(r34-re14)**2)
       xi3=(r14-re14) *exp(-beta*(r14-re14)**2)

       xi4=ta2
       xi5=ta3
       xi6=ta1
       ! 
     case (3)
       xi1=(r34-re14) *exp(-beta*(r34-re14)**2)
	   xi2=(r14-re14) *exp(-beta*(r14-re14)**2)
       xi3=(r24-re14) *exp(-beta*(r24-re14)**2)

       xi4=ta3
       xi5=ta1
       xi6=ta2
       ! 
     end select


     t1=0 ; t2=0 ; t3=0 ; t4=0 
     if (order>=1) then 
      t1 = F1*xi1+(xi2+xi3)*F3+(xi6+xi5)*F5+F4*xi4
     endif 

     if (order>=2) then 
    t2 = F11*xi1**2+F14*xi1*xi4+F23*xi2*xi3+(xi1*xi3+xi1*xi2)*F13+(xi1*xi5+xi1*xi6)*F16+&
        (xi2*xi4+xi3*xi4)*F34+(xi3**2+xi2**2)*F33+(xi3*xi6+xi2*xi5)*F36+(xi2*xi6+xi3*xi5)*F35+&
        F44*xi4**2+(xi4*xi6+xi4*xi5)*F46+(xi6**2+xi5**2)*F55+F56*xi5*xi6
     endif                                                              
                                                                        
     if (order>=3) then                                                 
        s1 = (xi1**2*xi3+xi1**2*xi2)*F112+(xi3*xi4*xi6+xi2*xi4*xi5)*F245+&
        (xi2*xi5*xi6+xi3*xi6*xi5)*F256+(xi1**2*xi5+xi1**2*xi6)*F115+&
        (xi3*xi5**2+xi2*xi6**2)*F266+(xi2**2*xi6+xi3**2*xi5)*F226+&
        (xi1*xi4*xi6+xi1*xi4*xi5)*F146+(xi1*xi6**2+xi1*xi5**2)*F155+&
        (xi2*xi5**2+xi3*xi6**2)*F255+(xi1*xi2*xi5+xi1*xi3*xi6)*F136+&
        (xi1*xi3**2+xi1*xi2**2)*F133+F144*xi1*xi4**2+F114*xi1**2*xi4+&
        F234*xi2*xi3*xi4+F123*xi1*xi2*xi3+F456*xi4*xi5*xi6
        t3 = s1+(xi1*xi3*xi4+xi1*xi2*xi4)*F124+(xi6**3+xi5**3)*F555+&
        (xi1*xi3*xi5+xi1*xi2*xi6)*F135+(xi3*xi4*xi5+xi2*xi4*xi6)*F246+&
        (xi2**2*xi4+xi3**2*xi4)*F334+(xi3**2*xi6+xi2**2*xi5)*F225+&
        (xi3*xi4**2+xi2*xi4**2)*F344+F156*xi1*xi5*xi6+(xi3*xi2*xi6+xi2*xi3*xi5)*F235+&
        (xi3**3+xi2**3)*F333+F444*xi4**3+(xi4*xi5**2+xi4*xi6**2)*F466+&
        (xi3**2*xi2+xi2**2*xi3)*F223+F111*xi1**3+(xi5**2*xi6+xi6**2*xi5)*F556+&
        (xi4**2*xi6+xi4**2*xi5)*F445
     endif                                                              
                                                                        
     if (order>=4) then                                                       
        s2 = (xi1*xi2**3+xi1*xi3**3)*F1222+ (xi3*xi2**3+xi2*xi3**3)*F2333+ &
        (xi3*xi5**3+xi2*xi6**3)*F3555+                                     &
        (xi1*xi2*xi6**2+xi1*xi3*xi5**2)*F1355+                             &
        (xi3*xi6**3+xi2*xi5**3)*F3666+                                     &
        (xi2*xi4*xi5*xi6+xi3*xi4*xi6*xi5)*F2456+                           &
        (xi2**3*xi6+xi3**3*xi5)*F3335+                                     &
        (xi3*xi4**2*xi6+xi2*xi4**2*xi5)*F2445+                             &
        (xi1*xi6*xi5**2+xi1*xi5*xi6**2)*F1566+                             &
        (xi4**2*xi5**2+xi4**2*xi6**2)*F4466+                               &
        (xi3*xi4*xi6**2+xi2*xi4*xi5**2)*F2455+ (xi5**4+xi6**4)*F6666+      &
        (xi1*xi3*xi6*xi5+xi1*xi2*xi5*xi6)*F1256+                           &
        (xi1*xi3*xi2*xi6+xi1*xi2*xi3*xi5)*F1235+                           &
        (xi1**2*xi3*xi6+xi1**2*xi2*xi5)*F1136+                             &
        (xi3**3*xi4+xi2**3*xi4)*F2224+ (xi4*xi5**3+xi4*xi6**3)*F4555       
        
        s1 = s2+ (xi2*xi3*xi4*xi5+xi3*xi2*xi4*xi6)*F2345+                  &
        (xi2*xi3**2*xi4+xi3*xi2**2*xi4)*F2334+                             &
        (xi1*xi3*xi2**2+xi1*xi2*xi3**2)*F1233+                             &
        (xi1*xi3*xi4*xi5+xi1*xi2*xi4*xi6)*F1246+                           &
        (xi1*xi2*xi4*xi5+xi1*xi3*xi4*xi6)*F1245+                           &
        (xi2**2*xi6**2+xi3**2*xi5**2)*F2266+                               &
        (xi1*xi3**2*xi6+xi1*xi2**2*xi5)*F1225+                             &
        (xi2**2*xi4*xi6+xi3**2*xi4*xi5)*F2246 +F1444*xi1*xi4**3            &
        +F5566*xi5**2*xi6**2 +F2233*xi2**2*xi3**2 +F1114*xi1**3*xi4+       &
        (xi3*xi4**2*xi5+xi2*xi4**2*xi6)*F3445+                             &
        (xi1**2*xi3*xi4+xi1**2*xi2*xi4)*F1124+                             &
        (xi1**3*xi6+xi1**3*xi5)*F1115+ (xi3**4+xi2**4)*F2222+              &
        (xi2**2*xi4**2+xi3**2*xi4**2)*F2244 +F4444*xi4**4 

        s2 = s1          &
        +F1144*xi1**2*xi4**2+ (xi2*xi6*xi5**2+xi3*xi5*xi6**2)*F3566+       &
        (xi1*xi3**2*xi5+xi1*xi2**2*xi6)*F1335+                             &
        (xi6*xi5**3+xi5*xi6**3)*F5666+                                     &
        (xi1*xi4*xi6**2+xi1*xi4*xi5**2)*F1466 +F2356*xi2*xi3*xi5*xi6       &
        +F1234*xi1*xi2*xi3*xi4+ (xi1*xi2*xi5**2+xi1*xi3*xi6**2)*F1366+     &
        (xi1**2*xi4*xi6+xi1**2*xi4*xi5)*F1146+                             &
        (xi1*xi4**2*xi5+xi1*xi4**2*xi6)*F1445+                             &
        (xi3*xi2*xi5**2+xi2*xi3*xi6**2)*F2366+                             &
        (xi3*xi5**2*xi6+xi2*xi6**2*xi5)*F3556 +F1456*xi1*xi4*xi5*xi6+      &
        (xi2*xi3**2*xi6+xi3*xi2**2*xi5)*F2336 +F1111*xi1**4+               &
        (xi4*xi5**2*xi6+xi4*xi6**2*xi5)*F4556+                             &
        (xi3*xi2**2*xi6+xi2*xi3**2*xi5)*F2335 

        t4 = s2+                     &
        (xi2*xi4**3+xi3*xi4**3)*F2444+ (xi1*xi5**3+xi1*xi6**3)*F1666+      &
        (xi1*xi2**2*xi4+xi1*xi3**2*xi4)*F1334+                             &
        (xi2**2*xi4*xi5+xi3**2*xi4*xi6)*F2245+                             &
        (xi1**2*xi3**2+xi1**2*xi2**2)*F1122 +F4456*xi4**2*xi5*xi6+         &
        (xi1**2*xi5**2+xi1**2*xi6**2)*F1155+                               &
        (xi2*xi4*xi6**2+xi3*xi4*xi5**2)*F2466+                             &
        (xi1*xi2*xi4**2+xi1*xi3*xi4**2)*F1244+                             &
        (xi2**3*xi5+xi3**3*xi6)*F2225+ (xi4**3*xi5+xi4**3*xi6)*F4445       &
        +F1156*xi1**2*xi5*xi6+ (xi1**2*xi2*xi6+xi1**2*xi3*xi5)*F1126+      &
        (xi1**3*xi3+xi1**3*xi2)*F1112+                                     &
        (xi2**2*xi5*xi6+xi3**2*xi6*xi5)*F2256 +F2344*xi2*xi3*xi4**2+       &
        (xi2**2*xi5**2+xi3**2*xi6**2)*F2255 +F1123*xi1**2*xi2*xi3           
     endif                                                                  

     DMS = (mu0+t1+t2+t3+t4)

end function DMS



 double precision function DMS_A(parmax,param,local)

 integer,intent(in)          ::  parmax
 double precision,intent(in) ::  param(parmax)
 double precision,intent(in) ::  local(6)

 double precision            ::  r14,r24,r34,alpha1,alpha2,alpha3
 double precision            ::  y1,y2,y3,y4,y5,alpha,rho

 double precision            ::  v,v0,rhoe,pi,v1,v2,v3,v4,v5,v6

 double precision            ::  sinrho,coro,cosrho,beta,de,b0

                        
 double precision                                        &
      fea    ,fea1  ,                                    &
      fea11  ,fea12  ,fea14  ,fea44  ,                   &
      fea111 ,fea112 ,fea114 ,fea123 ,                   &
      fea124 ,fea144 ,fea155 ,fea455 ,                   &
      fea1111,fea1112,fea1114,fea1122,                   &
      fea1123,fea1124,fea1125,fea1144,                   &
      fea1155,fea1244,fea1255,fea1444,                   &
      fea1455,fea4444
                    
 double precision :: &      
      Rhoedg,re14,aa1,ve  ,f0a,  &
      f1a,f2a,f3a,f4a,f5a,f6a,f7a,f8a, &
      f1a1,f2a1,f3a1,f4a1,f5a1,f6a1,  &
      f0a11,f1a11,f2a11,f3a11,f4a11, &
      f0a12,f1a12,f2a12,f3a12,f4a12, &
      f0a14,f1a14,f2a14,f3a14,f4a14, &
      f0a44,f1a44,f2a44,f3a44,f4a44, &
      f0a111,f1a111,f2a111,f3a111  , &
      f0a112,f1a112,f2a112,f3a112  , &
      f0a114,f1a114,f2a114,f3a114  , &
      f0a123,f1a123,f2a123,f3a123  , &
      f0a124,f1a124,f2a124,f3a124  , &
      f0a144,f1a144,f2a144,f3a144  , &
      f0a155,f1a155,f2a155,f3a155  , &
      f0a455,f1a455,f2a455,f3a455  , &
      f0a1111,f1a1111,f2a1111      , &
      f0a1112,f1a1112,f2a1112      , &
      f0a1114,f1a1114,f2a1114      , &
      f0a1122,f1a1122,f2a1122      , &
      f0a1123,f1a1123,f2a1123      , &
      f0a1124,f1a1124,f2a1124      , &
      f0a1125,f1a1125,f2a1125      , &
      f0a1144,f1a1144,f2a1144      , &
      f0a1155,f1a1155,f2a1155      , &
      f0a1244,f1a1244,f2a1244      , &
      f0a1255,f1a1255,f2a1255      , &
      f0a1444,f1a1444,f2a1444      , &
      f0a1455,f1a1455,f2a1455      , &
      f0a4444,f1a4444,f2a4444

  !-------------------------------


   pi=3.141592653589793

   rhoedg     = param(  1)
   rhoe=pi*rhoedg/1.8d+02 

   re14       = param( 2)
   b0         = param( 3)**2
   de         = param( 4)

   f1a        = param(  5)
   f2a        = param(  6)
   f3a        = param(  7)
   f4a        = param(  8)
   f5a        = param(  9)
   f6a        = param( 10)
   f7a        = param( 11)
   !
   f0a        = param( 12)
   f1a1       = param( 13)
   f2a1       = param( 14)
   f3a1       = param( 15)
   f4a1       = param( 16)
   f5a1       = param( 17)
   f6a1       = param( 18)
   f0a11      = param( 19)
   f1a11      = param( 20)
   f2a11      = param( 21)
   f3a11      = param( 22)
   f4a11      = param( 23)
   f0a12      = param( 24)
   f1a12      = param( 25)
   f2a12      = param( 26)
   f3a12      = param( 27)
   f4a12      = param( 28)
   f0a14      = param( 29)
   f1a14      = param( 30)
   f2a14      = param( 31)
   f3a14      = param( 32)
   f4a14      = param( 33)
   f0a44      = param( 34)
   f1a44      = param( 35)
   f2a44      = param( 36)
   f3a44      = param( 37)
   f4a44      = param( 38)
   f0a111     = param( 39)
   f1a111     = param( 40)
   f2a111     = param( 41)
   f3a111     = param( 42)
   f0a112     = param( 43)
   f1a112     = param( 44)
   f2a112     = param( 45)
   f3a112     = param( 46)
   f0a114     = param( 47)
   f1a114     = param( 48)
   f2a114     = param( 49)
   f3a114     = param( 50)
   f0a123     = param( 51)
   f1a123     = param( 52)
   f2a123     = param( 53)
   f3a123     = param( 54)
   f0a124     = param( 55)
   f1a124     = param( 56)
   f2a124     = param( 57)
   f3a124     = param( 58)
   f0a144     = param( 59)
   f1a144     = param( 60)
   f2a144     = param( 61)
   f3a144     = param( 62)
   f0a155     = param( 63)
   f1a155     = param( 64)
   f2a155     = param( 65)
   f3a155     = param( 66)
   f0a455     = param( 67)
   f1a455     = param( 68)
   f2a455     = param( 69)
   f3a455     = param( 70)
   f0a1111    = param( 71)
   f1a1111    = param( 72)
   f2a1111    = param( 73)
   f0a1112    = param( 74)
   f1a1112    = param( 75)
   f2a1112    = param( 76)
   f0a1114    = param( 77)
   f1a1114    = param( 78)
   f2a1114    = param( 79)
   f0a1122    = param( 80)
   f1a1122    = param( 81)
   f2a1122    = param( 82)
   f0a1123    = param( 83)
   f1a1123    = param( 84)
   f2a1123    = param( 85)
   f0a1124    = param( 86)
   f1a1124    = param( 87)
   f2a1124    = param( 88)
   f0a1125    = param( 89)
   f1a1125    = param( 90)
   f2a1125    = param( 91)
   f0a1144    = param( 92)
   f1a1144    = param( 93)
   f2a1144    = param( 94)
   f0a1155    = param( 95)
   f1a1155    = param( 96)
   f2a1155    = param( 97)
   f0a1244    = param( 98)
   f1a1244    = param( 99)
   f2a1244    = param(100)
   f0a1255    = param(101)
   f1a1255    = param(102)
   f2a1255    = param(103)
   f0a1444    = param(104)
   f1a1444    = param(105)
   f2a1444    = param(106)
   f0a1455    = param(107)
   f1a1455    = param(108)
   f2a1455    = param(109)
   f0a4444    = param(110)
   f1a4444    = param(111)
   f2a4444    = param(112)

   r14    = local(1) ;  r24    = local(2) ;  r34    = local(3)
   alpha1 = local(4) ;  alpha2 = local(5) ;  alpha3 = local(6)

   pi=3.141592653589793
   rhoe=pi*rhoedg/1.8d+02

   y4=(2.d0*alpha1-alpha2-alpha3)/sqrt(6.d0)
   y5=(alpha2-alpha3)/sqrt(2.d0)

   alpha=(alpha1+alpha2+alpha3)/3.d0
   rho=pi-dasin(2.d0*sin(alpha*0.5d0)/dsqrt(3.d0))

   if ( 2.d0*sin(alpha*0.5d0)/dsqrt(3.d0).ge.1.0d0 ) then 
     sinrho=1.d0 
   else 
     sinrho = 2.d0*sin(alpha*0.5d0)/dsqrt(3.d0)
   endif 

   cosrho =-sqrt(1.d0-sinrho**2)
   !
   coro=(sin(rhoe)-sinrho)
   !
   y1=1.0d0*(r14-re14) *exp(-b0*(r14-re14)**2)
   y2=1.0d0*(r24-re14) *exp(-b0*(r24-re14)**2)
   y3=1.0d0*(r34-re14) *exp(-b0*(r34-re14)**2)
   !
   v0=de+f1a*coro+f2a*coro**2+f3a*coro**3+f4a*coro**4+f5a*coro**5 &
         +f6a*coro**6+f7a*coro**7

   fea1= f0a + f1a1*coro+f2a1*coro**2+f3a1*coro**3+f4a1*coro**4+f5a1*coro**5+f6a1*coro**6

   fea11=   f0a11+f1a11*coro+f2a11*coro**2+f3a11*coro**3+f4a11*coro**4
   fea12=   f0a12+f1a12*coro+f2a12*coro**2+f3a12*coro**3+f4a12*coro**4
   fea14=   f0a14+f1a14*coro+f2a14*coro**2+f3a14*coro**3+f4a14*coro**4
   fea44=   f0a44+f1a44*coro+f2a44*coro**2+f3a44*coro**3+f4a44*coro**4
 
   fea111= f0a111+f1a111*coro+f2a111*coro**2+f3a111*coro**3
   fea112= f0a112+f1a112*coro+f2a112*coro**2+f3a112*coro**3
   fea114= f0a114+f1a114*coro+f2a114*coro**2+f3a114*coro**3
   fea123= f0a123+f1a123*coro+f2a123*coro**2+f3a123*coro**3
   fea124= f0a124+f1a124*coro+f2a124*coro**2+f3a124*coro**3
   fea144= f0a144+f1a144*coro+f2a144*coro**2+f3a144*coro**3
   fea155= f0a155+f1a155*coro+f2a155*coro**2+f3a155*coro**3
   fea455= f0a455+f1a455*coro+f2a455*coro**2+f3a455*coro**3

   fea1111= f0a1111+f1a1111*coro+f2a1111*coro**2
   fea1112= f0a1112+f1a1112*coro+f2a1112*coro**2
   fea1114= f0a1114+f1a1114*coro+f2a1114*coro**2
   fea1122= f0a1122+f1a1122*coro+f2a1122*coro**2
   fea1123= f0a1123+f1a1123*coro+f2a1123*coro**2
   fea1124= f0a1124+f1a1124*coro+f2a1124*coro**2
   fea1125= f0a1125+f1a1125*coro+f2a1125*coro**2
   fea1144= f0a1144+f1a1144*coro+f2a1144*coro**2
   fea1155= f0a1155+f1a1155*coro+f2a1155*coro**2
   fea1244= f0a1244+f1a1244*coro+f2a1244*coro**2
   fea1255= f0a1255+f1a1255*coro+f2a1255*coro**2
   fea1444= f0a1444+f1a1444*coro+f2a1444*coro**2
   fea1455= f0a1455+f1a1455*coro+f2a1455*coro**2
   fea4444= f0a4444+f1a4444*coro+f2a4444*coro**2

   fea44444 = f0a44444  + f1a44444 *coro+ f2a44444 *coro**2
   fea33455 = f0a33455  + f1a33455 *coro+ f2a33455 *coro**2
   fea33445 = f0a33445  + f1a33445 *coro+ f2a33445 *coro**2
   fea33345 = f0a33345  + f1a33345 *coro+ f2a33345 *coro**2
   fea33344 = f0a33344  + f1a33344 *coro+ f2a33344 *coro**2
   fea33334 = f0a33334  + f1a33334 *coro+ f2a33334 *coro**2
   fea33333 = f0a33333  + f1a33333 *coro+ f2a33333 *coro**2
   fea25555 = f0a25555  + f1a25555 *coro+ f2a25555 *coro**2
   fea24455 = f0a24455  + f1a24455 *coro+ f2a24455 *coro**2
   fea24445 = f0a24445  + f1a24445 *coro+ f2a24445 *coro**2
   fea23333 = f0a23333  + f1a23333 *coro+ f2a23333 *coro**2
   fea13455 = f0a13455  + f1a13455 *coro+ f2a13455 *coro**2
   fea13445 = f0a13445  + f1a13445 *coro+ f2a13445 *coro**2
   fea13345 = f0a13345  + f1a13345 *coro+ f2a13345 *coro**2
   fea12355 = f0a12355  + f1a12355 *coro+ f2a12355 *coro**2
   fea11334 = f0a11334  + f1a11334 *coro+ f2a11334 *coro**2
   fea11333 = f0a11333  + f1a11333 *coro+ f2a11333 *coro**2
   fea11255 = f0a11255  + f1a11255 *coro+ f2a11255 *coro**2
   fea11245 = f0a11245  + f1a11245 *coro+ f2a11245 *coro**2
   fea11234 = f0a11234  + f1a11234 *coro+ f2a11234 *coro**2
   fea11233 = f0a11233  + f1a11233 *coro+ f2a11233 *coro**2
   fea11135 = f0a11135  + f1a11135 *coro+ f2a11135 *coro**2
   fea11134 = f0a11134  + f1a11134 *coro+ f2a11134 *coro**2
   fea11123 = f0a11123  + f1a11123 *coro+ f2a11123 *coro**2
   fea555555= f0a555555 + f1a555555*coro+ f2a555555*coro**2
   fea444444= f0a444444 + f1a444444*coro+ f2a444444*coro**2
   fea335555= f0a335555 + f1a335555*coro+ f2a335555*coro**2
   fea334455= f0a334455 + f1a334455*coro+ f2a334455*coro**2
   fea334445= f0a334445 + f1a334445*coro+ f2a334445*coro**2
   fea333555= f0a333555 + f1a333555*coro+ f2a333555*coro**2
   fea333333= f0a333333 + f1a333333*coro+ f2a333333*coro**2
   fea244555= f0a244555 + f1a244555*coro+ f2a244555*coro**2
   fea244455= f0a244455 + f1a244455*coro+ f2a244455*coro**2
   fea233445= f0a233445 + f1a233445*coro+ f2a233445*coro**2
   fea233444= f0a233444 + f1a233444*coro+ f2a233444*coro**2
   fea233345= f0a233345 + f1a233345*coro+ f2a233345*coro**2
   fea233344= f0a233344 + f1a233344*coro+ f2a233344*coro**2
   fea233335= f0a233335 + f1a233335*coro+ f2a233335*coro**2
   fea223355= f0a223355 + f1a223355*coro+ f2a223355*coro**2
   fea222335= f0a222335 + f1a222335*coro+ f2a222335*coro**2
   fea222334= f0a222334 + f1a222334*coro+ f2a222334*coro**2
   fea222333= f0a222333 + f1a222333*coro+ f2a222333*coro**2
   fea222255= f0a222255 + f1a222255*coro+ f2a222255*coro**2
   fea222245= f0a222245 + f1a222245*coro+ f2a222245*coro**2
   fea222233= f0a222233 + f1a222233*coro+ f2a222233*coro**2
   fea222224= f0a222224 + f1a222224*coro+ f2a222224*coro**2
   fea145555= f0a145555 + f1a145555*coro+ f2a145555*coro**2
   fea134444= f0a134444 + f1a134444*coro+ f2a134444*coro**2
   fea133444= f0a133444 + f1a133444*coro+ f2a133444*coro**2
   fea133345= f0a133345 + f1a133345*coro+ f2a133345*coro**2
   fea133334= f0a133334 + f1a133334*coro+ f2a133334*coro**2
   fea133333= f0a133333 + f1a133333*coro+ f2a133333*coro**2
   fea124555= f0a124555 + f1a124555*coro+ f2a124555*coro**2
   fea124455= f0a124455 + f1a124455*coro+ f2a124455*coro**2
   fea123455= f0a123455 + f1a123455*coro+ f2a123455*coro**2
   fea123345= f0a123345 + f1a123345*coro+ f2a123345*coro**2
   fea113555= f0a113555 + f1a113555*coro+ f2a113555*coro**2
   fea113345= f0a113345 + f1a113345*coro+ f2a113345*coro**2
   fea112355= f0a112355 + f1a112355*coro+ f2a112355*coro**2
   fea112335= f0a112335 + f1a112335*coro+ f2a112335*coro**2
   fea112233= f0a112233 + f1a112233*coro+ f2a112233*coro**2
   fea111444= f0a111444 + f1a111444*coro+ f2a111444*coro**2
   fea111234= f0a111234 + f1a111234*coro+ f2a111234*coro**2
   fea111233= f0a111233 + f1a111233*coro+ f2a111233*coro**2
   fea111123= f0a111123 + f1a111123*coro+ f2a111123*coro**2



   v1 = (y3+y2+y1)*fea1

   v2 = (y2*y3+y1*y3+y1*y2)*fea12                                                                 &
    +(y2**2+y3**2+y1**2)*fea11                                                                    &
    +(-sqrt(3.0d0)*y3*y5/2.0d0-y3*y4/2.0d0+y1*y4+sqrt(3.0d0)*y2*y5/2.0d0-y2*y4/2.0d0)*fea14 &
    +(y5**2+y4**2)*fea44

   v3 = (y1*y3*y4+y1*y2*y4-2.0d0*y2*y3*y4+sqrt(3.0d0)*y1*y2*y5-sqrt(3.0d0)*y1*y3*y5)*fea124   &
    +(3.0d0/4.0d0*y3*y4**2-sqrt(3.0d0)*y3*y4*y5/2.0d0+y1*y5**2+y2*y5**2/4.0d0               & 
    +3.0d0/4.0d0*y2*y4**2+sqrt(3.0d0)*y2*y4*y5/2.0d0+y3*y5**2/4.0d0)*fea155                 &
    +(y2*y3**2+y1*y3**2+y1**2*y3+y1*y2**2+y2**2*y3+y1**2*y2)*fea112+                             &
    (-y4**3/3.0d0+y4*y5**2)*fea455+fea123*y1*y2*y3                                              &
    +(y1*y4**2+3.0d0/4.0d0*y3*y5**2+3.0d0/4.0d0*y2*y5**2+y2*y4**2/4.0d0                     &
    -sqrt(3.0d0)*y2*y4*y5/2.0d0+sqrt(3.0d0)*y3*y4*y5/2.0d0+y3*y4**2/4.0d0)*fea144           &
    +(y3**3+y2**3+y1**3)*fea111+(-y2**2*y4/2.0d0-y3**2*y4/2.0d0+sqrt(3.0d0)*y2**2*y5/2.0d0   & 
    +y1**2*y4-sqrt(3.0d0)*y3**2*y5/2.0d0)*fea114
    !
   s2 = (y4**4+y5**4+2.0d0*y4**2*y5**2)*fea4444+(3.0d0/8.0d0*sqrt(3.0d0)*&
    y2*y5**3-3.0d0/8.0d0*sqrt(3.0d0)*y3*y4**2*y5-3.0d0/8.0d0*sqrt(3.0d0)*y3*&
    y5**3-9.0d0/8.0d0*y2*y4*y5**2-y3*y4**3/8.0d0-y2*y4**3/8.0d0-9.0d0/8.0d0*&
    y3*y4*y5**2+y1*y4**3+3.0d0/8.0d0*sqrt(3.0d0)*y2*y4**2*y5)*fea1444 &
    +(3.0d0/4.0d0*y2**2*y4**2+3.0d0/4.0d0*y3**2*y4**2+y1**2*y5**2+y3**2*y5**2/4.0d0 &
    -sqrt(3.0d0)*y3**2*y4*y5/2.0d0+sqrt(3.0d0)*y2**2*y4*y5/2.0d0+y2**2&
    *y5**2/4.0d0)*fea1155 
    s1 = s2+(y3**2*y4**2/4.0d0+3.0d0/4.0d0*y3**2*y5**2+y1**2*y4**2+y2**2*&
    y4**2/4.0d0+sqrt(3.0d0)*y3**2*y4*y5/2.0d0-sqrt(3.0d0)*y2**2*y4*y5/2.0d0&
    +3.0d0/4.0d0*y2**2*y5**2)*fea1144+(y1**3*y4+sqrt(3.0d0)*y2**3*y5/2.0d0&
    -sqrt(3.0d0)*y3**3*y5/2.0d0-y2**3*y4/2.0d0-y3**3*y4/2.0d0)*fea1114&
    +(y2**4+y1**4+y3**4)*fea1111+(sqrt(3.0d0)*y1*y3*y4*y5+3.0d0/2.0d0*y2*y3*y5**2&
    -y2*y3*y4**2/2.0d0+y1*y2*y4**2-sqrt(3.0d0)*y1*y2*y4*y5+y1*y3*y4**2)*fea1244 
    !
   s2 = s1+(y1*y3*y5**2+y1*y2*y5**2-sqrt(3.0d0)*y1*y3*y4*y5-y2*y3*y5**&
    2/2.0d0+3.0d0/2.0d0*y2*y3*y4**2+sqrt(3.0d0)*y1*y2*y4*y5)*fea1255+&
    (-y1*y3**2*y4/2.0d0+y1**2*y3*y4-sqrt(3.0d0)*y1*y3**2*y5/2.0d0-sqrt(3.0d0)*y2&
    *y3**2*y5/2.0d0+y1**2*y2*y4+sqrt(3.0d0)*y2**2*y3*y5/2.0d0-y2**2*y3*y4&
    /2.0d0+sqrt(3.0d0)*y1*y2**2*y5/2.0d0-y2*y3**2*y4/2.0d0-y1*y2**2*y4/2.0d0&
    )*fea1124+(y1**2*y2*y5+sqrt(3.0d0)*y1*y3**2*y4/2.0d0+sqrt(3.0d0)*y1*&
    y2**2*y4/2.0d0-sqrt(3.0d0)*y2*y3**2*y4/2.0d0-sqrt(3.0d0)*y2**2*y3*y4/2.0d0&
    -y2**2*y3*y5/2.0d0+y2*y3**2*y5/2.0d0-y1*y3**2*y5/2.0d0+y1*y2**2*y5&
    /2.0d0-y1**2*y3*y5)*fea1125 
    !
   v4 = s2+(y2*y3**3+y1**3*y3+y1**3*y2+y1*y2**3+y1*y3**3+y2**3*y3)*fea1112+&
    (y2**2*y3**2+y1**2*y3**2+y1**2*y2**2)*fea1122+(y1*y2**2*y3&
    +y1**2*y2*y3+y1*y2*y3**2)*fea1123+(5.0d0/8.0d0*y2*y4*y5**2+sqrt(3.0d0)*&
    y2*y5**3/8.0d0-sqrt(3.0d0)*y3*y4**2*y5/8.0d0+sqrt(3.0d0)*y2*y4**2*y5/8.0d0&
    -3.0d0/8.0d0*y2*y4**3+y1*y4*y5**2-sqrt(3.0d0)*y3*y5**3/8.0d0&
    +5.0d0/8.0d0*y3*y4*y5**2-3.0d0/8.0d0*y3*y4**3)*fea1455

   
   
   DMS_A = (v0+v1+v2+v3+v4+v5+v6 )*cosrho ! (cos(rhoe)-cosrho)

end function DMS_A



  subroutine MLlinur(dimen,npar,coeff,constant,solution,error)

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  real(8),intent(in)  :: coeff(npar,npar),constant(npar)
  real(8),intent(out) :: solution(npar)
  real(8)          :: a0(npar,npar)
  real(8)          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine MLlinur


  subroutine vector_product(v1,v2,v)
    !
    double precision,intent(in) :: v1(3),v2(3)
    double precision :: v(3)
    !
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    !
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    !
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
    !
  end subroutine vector_product

!
!   Skip n lines  in the input file 
!
  subroutine skiplines( inpunit,n )
  integer,intent(in) :: n,inpunit
  character(len=80)  :: label
  integer            :: i0

    do i0=1,n
       read  (inpunit,"(a80)") label 
    enddo

  end subroutine skiplines

