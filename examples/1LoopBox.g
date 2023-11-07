#! @Chunk 1LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1", "k1..2,k4" );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1, k2, k4 ]>
s12 := 2*k1*k2;
#! 2*k1*k2
SetAbbreviation( s12, "s12" );
s14 := 2*k1*k4;
#! 2*k1*k4
SetAbbreviation( s14, "s14" );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );
#! [ k1^2, k2^2, k4^2 ]
rel2 := [ (k1+k2+k4)^2 ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD,
        [ l1^2, l1*k1, l1*k2, l1*k4, s12, s14 ] );
SetPropagators( LD, -[ l1^2, (l1-k1)^2, (l1-k1-k2)^2, (l1+k4)^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s12, s14 ] );
R := RingOfLoopDiagram( LD );
#! Q[d,s12,s14][D1,D2,D3,D4]
ibps := MatrixOfIBPRelations( LD );
#! <A 4 x 1 matrix over a residue class ring>
ibp1 := ibps[1,1];
#! |[ -a2*D1*D2_-s12*a3*D3_-a3*D1*D3_-a4*D1*D4_+d-2*a1-a2-a3-a4 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ -a3 ]|, |[ D1*D3_ ]| ],
#!   [ |[ -a4 ]|, |[ D1*D4_ ]| ],
#!   [ |[ -s12*a3 ]|, |[ D3_ ]| ],
#!   [ |[ d-2*a1-a2-a3-a4 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp1 );
#! Q[d,s12,s14][a1,a2,a3,a4]<D1,D1_,D2,D2_,D3,D3_,D4,D4_>/( D4*D4_-1, D3*D3_-1,\
#!   D2*D2_-1, D1*D1_-1 )
E12 := PairOfMatricesOfLoopDiagramInPropagators( LD );
#! [ <A 4 x 4 matrix over an external ring>,
#!   <A 4 x 4 matrix over an external ring> ]
Display( E12[1] );
#! 2*D1,     D1-D2,     -s12+D2-D3,-D1+D4,
#! D1+D2,    D1-D2,     D2-D3,     s14-D1+D4,
#! s12+D1+D3,s12+D1-D2, D2-D3,     -s12-D1+D4,
#! D1+D4,    -s14+D1-D2,s14+D2-D3, -D1+D4
S := SyzygiesOfColumns( E12 );
#! <A non-zero 4 x 8 matrix over an external ring>
Sred := ReducedBasisOfColumnModule( BasisOfColumnModule( S ) );
#! <A non-zero 4 x 6 matrix over an external ring>
Display( Sred );
#! D2-D4,D1-D3,s12*D4+2*D3*D4-2*D4^2,     s14*D3-2*D3^2+2*D3*D4,-D3*D4^2,D3^2*D4,
#! D4,   -D1,  -s12*D4+D2*D4-D3*D4+2*D4^2,-D1*D3-D3*D4,         D3*D4^2, 0,      
#! 0,    -D1,  D1*D4+D2*D4,               -2*D1*D3,             0,       D1*D3*D4,
#! D2,   0,    2*D2*D4,                   -D1*D3-D2*D3,         D2*D3*D4,0       
Display( EntriesOfHomalgMatrixAsListList( CertainColumns( Sred, [1 .. 3] ) ) );
#! [ [ D2-D4, D1-D3, s12*D4+2*D3*D4-2*D4^2 ],
#!   [ D4, -D1, -s12*D4+D2*D4-D3*D4+2*D4^2 ],
#!   [ 0, -D1, D1*D4+D2*D4 ],
#!   [ D2, 0, 2*D2*D4 ] ]
Sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD );
#! |[ -s14*a2+s14*a4+d*D2-a1*D2-a2*D2-a3*D2-a4*D2-d*D4+a1*D4+a2*D4+a3*D4+a4*D4 ]|
ViewList( DecomposeInMonomials( Sibp1 ) );
#! [ [ |[ d-a1-a2-a3-a4 ]|, |[ D2 ]| ],
#!   [ |[ -d+a1+a2+a3+a4 ]|, |[ D4 ]| ],
#!   [ |[ -s14*a2+s14*a4 ]|, |[ 1 ]| ] ]
sibp1 := IBPRelation( CertainColumns( S, [ 1 ] ), LD, [ 1, 1, 1, 1 ] );
#! |[ d*D2-d*D4-4*D2+4*D4 ]|
ViewList( DecomposeInMonomials( sibp1 ) );
#! [ [ |[ d-4 ]|, |[ D2 ]| ],
#!   [ |[ -d+4 ]|, |[ D4 ]| ] ]
Sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD );
#! |[ -s12*a1+s12*a3+d*D1-a1*D1-a2*D1-a3*D1-a4*D1-d*D3+a1*D3+a2*D3+a3*D3+a4*D3 ]|
ViewList( DecomposeInMonomials( Sibp2 ) );
#! [ [ |[ d-a1-a2-a3-a4 ]|, |[ D1 ]| ],
#!   [ |[ -d+a1+a2+a3+a4 ]|, |[ D3 ]| ],
#!   [ |[ -s12*a1+s12*a3 ]|, |[ 1 ]| ] ]
sibp2 := IBPRelation( CertainColumns( S, [ 2 ] ), LD, [ 1, 1, 1, 1 ] );
#! |[ d*D1-d*D3-4*D1+4*D3 ]|
ViewList( DecomposeInMonomials( sibp2 ) );
#! [ [ |[ d-4 ]|, |[ D1 ]| ],
#!   [ |[ -d+4 ]|, |[ D3 ]| ] ]
Sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD );;
ViewList( DecomposeInMonomials( Sibp3 ) );
#! [ [ |[ 2*d-2*a1-2*a2-2*a3-2*a4+2 ]|, |[ D3*D4 ]| ],
#!   [ |[ -2*d+2*a1+2*a2+2*a3+2*a4-2 ]|, |[ D4^2 ]| ],
#!   [ |[ -s14*a4+s14 ]|, |[ D1 ]| ],
#!   [ |[ -s12*a4+s12 ]|, |[ D2 ]| ],
#!   [ |[ -s14*a4+s14 ]|, |[ D3 ]| ],
#!   [ |[ d*s12-2*s12*a2-2*s14*a2-2*s12*a3-s12*a4+2*s14*a4+s12-2*s14 ]|, |[ D4 ]| ],
#!   [ |[ -s12*s14*a4+s12*s14 ]|, |[ 1 ]| ] ]
sibp3 := IBPRelation( CertainColumns( S, [ 3 ] ), LD, [ 1, 1, 1, 1 ] );
#! |[ d*s12*D4+2*d*D3*D4-2*d*D4^2-4*s12*D4-2*s14*D4-6*D3*D4+6*D4^2 ]|
ViewList( DecomposeInMonomials( sibp3 ) );
#! [ [ |[ 2*d-6 ]|, |[ D3*D4 ]| ],
#!   [ |[ -2*d+6 ]|, |[ D4^2 ]| ],
#!   [ |[ d*s12-4*s12-2*s14 ]|, |[ D4 ]| ] ]
bas := BasisOfIBPRelations( LD );
#! <A non-zero 28 x 1 matrix over a residue class ring>
Sbas := BasisOfSpecialIBPRelations( LD );
#! <A non-zero 28 x 1 matrix over a residue class ring>
bas = Sbas;
#! true
SymanzikPolynomials( LD );
#! [ z1+z2+z3+z4, -s12*z1*z3-s14*z2*z4 ]
SymanzikPolynomials( LD, [ 1, 2, 3, 4 ] );
#! [ z1+z2+z3+z4, -s12*z1*z3-s14*z2*z4 ]
SymanzikPolynomials( LD, [ 1, 2, 3 ] );
#! [ z1+z2+z3, -s12*z1*z3 ]
SymanzikPolynomials( LD, [ 1, 2 ] );
#! [ z1+z2, 0 ]
SymanzikPolynomials( LD, [ 1 ] );
#! [ z1, 0 ]
SymanzikPolynomials( LD, [ ] );
#! [ 0, 0 ]
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 4 matrix over an external ring>
Display( gen );
#! D3*D4,D1*D4,D2*D3,D1*D2
gen2 := GeneratorsOfScalelessSectors( LD, [ 2, 2, 2, 2 ] );
#! <An unevaluated 1 x 4 matrix over an external ring>
Display( gen2 );
#! D1*D2*D3^2*D4^2,D1^2*D2*D3*D4^2,D1*D2^2*D3^2*D4,D1^2*D2^2*D3*D4
prel2 := ColumnReversedMatrixOfCoefficientsOfParametricIBPs( LD, 2 );
#! [ <A non-zero 7 x 11 matrix over an external ring>,
#!   [ D1_*D3_, D1_*D4_, D2_*D4_, D3_*D4_, D1_, D2_, D3_, D4_, 1, D1, D2 ] ]
#! @EndExample

Y := RationalDoubleShiftAlgebra( R );
mibps := Y * ibps;
mbas := BasisOfRows( mibps );

lhs1 := HomalgMatrix( "[a1*D1_]", 1, 1, Y );
R1 := lhs1 - DecideZeroRows( lhs1, mbas );
lhs2 := HomalgMatrix( "[a2*D2_]", 1, 1, Y );
R2 := lhs2 - DecideZeroRows( lhs2, mbas );
lhs3 := HomalgMatrix( "[a3*D3_]", 1, 1, Y );
R3 := lhs3 - DecideZeroRows( lhs3, mbas );
lhs4 := HomalgMatrix( "[a4*D4_]", 1, 1, Y );
R4 := lhs4 - DecideZeroRows( lhs4, mbas );

Ris := UnionOfRows( [ R1, R2, R3, R4 ] );

b1 := HomalgMatrix( "[(d-2*(a1+a2+1))*(d-2*(a1+a4+1))*s12*s14*a1*D1_]", 1, 1, Ypol );
RHS1 := HomalgMatrix( "[ -2 * (d-2*(a1+a2+a4))*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1)) * D3 + 4 * (a3-1)*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1))*D4+(d-2*(a1+a3+a4))*(d-(a1+a2+a3+a4+1))*(d-2*(a1+a2+1))*s14 ]", 1, 1, Ypol );
NF1 := b1 - RHS1;
Assert( 0, IsZero( DecideZeroRows( NF1, ibps ) ) );
b2 := HomalgMatrix( "[(d-2*(a1+a2+1))*(d-2*(a2+a3+1))*s12*s14*a2*D2_]", 1, 1, Ypol );
RHS2 := HomalgMatrix( "[ 4 * (d-(a1+a2+a3+a4))*(a4-1)*(d-(a1+a2+a3+a4+1)) * D3 -2 * (d-2*(a1+a2+a3))*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1))*D4+(d-2*(a2+a3+a4))*(d-(a1+a2+a3+a4+1))*(d-2*(a1+a2+1))*s12 ]", 1, 1, Ypol );
NF2 := b2 - RHS2;
Assert( 0, IsZero( DecideZeroRows( NF2, ibps ) ) );
b3 := HomalgMatrix( "[(d-2*(a2+a3+1))*(d-2*(a3+a4+1))*s12*s14*a3*D3_]", 1, 1, Ypol );
RHS3 := HomalgMatrix( "[ -2 * (d-2*(a2+a3+a4))*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1)) * D3 + 4 * (a1-1)*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1))*D4+(d-2*(a1+a3+a4))*(d-(a1+a2+a3+a4+1))*(d-2*(a2+a3+1))*s14 - 2*(a1-a3)*(d-2*(a2+a3+a4))*(d-(a1+a2+a3+a4+1)) * s12 ]", 1, 1, Ypol );
NF3 := b3 - RHS3;
Assert( 0, IsZero( DecideZeroRows( NF3, ibps ) ) );
b4 := HomalgMatrix( "[(d-2*(a1+a4+1))*(d-2*(a3+a4+1))*s12*s14*a4*D4_]", 1, 1, Ypol );
RHS4 := HomalgMatrix( "[ 4 * (d-(a1+a2+a3+a4))*(a2-1)*(d-(a1+a2+a3+a4+1)) * D3 -2 * (d-2*(a1+a3+a4))*(d-(a1+a2+a3+a4))*(d-(a1+a2+a3+a4+1))*D4+(d-2*(a2+a3+a4))*(d-(a1+a2+a3+a4+1))*(d-2*(a1+a4+1))*s12 - 2 * (a2-a4) * (d-2*(a1+a3+a4)) * (d-(a1+a2+a3+a4+1)) * s14 ]", 1, 1, Ypol );
NF4 := b4 - RHS4;
Assert( 0, IsZero( DecideZeroRows( NF4, ibps ) ) );

nf := UnionOfRows( NF1, NF2, NF3, NF4 );
Assert( 0, IsZero( DecideZeroRows( nf, ibps ) ) );

tau := RightDivide( nf, ibps );

# DecideZeroRows( ibps, nf );
# BasisOfRows( nf );
# #c := RightDivide( mibps, nf );
# Error, the external CAS Maple (which should be running with PID 3441074) seems to have died!
# The last error was:
# maple: fatal error, lost connection to kernel

#From: Jan Piclum <piclum@physik.uni-siegen.de>
#Subject: 1LoopBox
#Date: 30. November 2020 at 16:03:03 CET
#To: Mohamed Barakat <mohamed.barakat@uni-siegen.de>
#Cc: Tobias Huber <huber@physik.uni-siegen.de>, Robin Brüser <Robin.Brueser@uni-siegen.de>

#Hallo Mohamed,

#hier ist die FIRE-Reduktion für die 1LoopBox:

#D1_ -> (-5 + D)/s12 - (4*(-5 + D)*(-3 + D)*D1*D3)/((-6 + D)*s12*s14^2)
#D2_ -> (-5 + D)/s14 - (4*(-5 + D)*(-3 + D)*D2*D4)/((-6 + D)*s12^2*s14)
#D3_ -> (-5 + D)/s12 - (4*(-5 + D)*(-3 + D)*D1*D3)/((-6 + D)*s12*s14^2)
#D4_ -> (-5 + D)/s14 - (4*(-5 + D)*(-3 + D)*D2*D4)/((-6 + D)*s12^2*s14)

#Viele Grüße
#Jan.

# b := "D1_" / Y;
# b := "D2_" / Y;
# b := "D1*D2" / Y; ## scaleless
# NormalFormWrtInitialIntegral( b / Y, mbas );
# Die GB ist nur fuer die Reduktion wichtig, sie würde aber sehr allgemeine Reduktionen erlauben
# Wir sind aber nur an folgende Reduktionen vor dem Einsetzen der a_i's interessiert: Reduziere nur die Di's, Di_'s und die Ni's

#y := AmbientRing( ReversedDoubleShiftAlgebra( RingOfLoopDiagram( LD ) ) );
#SortedList( List( MatrixOfCoefficientsOfIBPs( BasisOfRows( ibps ) )[3], d -> d / y ), function( a, b ) return a = LeadingMonomial( a + b ); end );

Qa := FieldOfCoefficientsOfLoopDiagramInHecke( LD );
prel2 := ColumnReversedMatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa : homalgIOMode := "d" );

Q := CoefficientsRing( AmbientRing( Y ) );
m := Q * prel2[1];
b := BasisOfRows( m );
homalgDisplay( [ "map(factor,", b, "):" ] );

# : subset := [ 1, 3, 5, 13, 15, 16, 17, 18, 20, 21, 22, 25, 26, 27, 30, 31, 32, 34, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 54, 55, 57, 59, 73, 75, 76, 77, 78, 79, 80, 84, 86, 91, 95, 97, 99, 101, 102, 103, 121, 124, 125, 126, 127, 134, 137, 138, 140, 145, 146, 147, 148, 150, 151, 152, 153, 154, 155, 157, 159, 161, 162, 163, 165, 167, 168 ]; # -> 1:54
# : subset := [ 1, 3, 5, 13, 15, 16, 17, 18, 20, 21, 22, 25, 26, 27, 30, 31, 32, 34, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 54, 55, 57, 59, 73, 75, 76, 77, 78, 79, 80, 84, 86, 91, 95, 97, 99, 101, 102, 103, 121, 124, 125, 126, 127, 134, 137, 138, 140, 145, 146, 147, 148, 150, 151, 152, 153, 154, 155, 157, 159, 161, 162, 163, 165, 167, 168, 169, 170, 171, 172, 174, 175, 176, 177, 178, 179 ]; # -> 4:26
# : subset := "all" -> 3:58
