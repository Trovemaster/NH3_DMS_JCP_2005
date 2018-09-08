# NH3_DMS_JCP_2005: DMS of NH3 from JCP 122, 104317 (2005)

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


