#! @Chunk 1LoopBox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k1", "p1..2,p4" );
#! <A loop diagram with loop momenta [ k1 ] & external momenta [ p1, p2, p4 ]>
s12 := 2*p1*p2;
#! 2*p1*p2
SetAbbreviation( s12, "s12" );
s14 := 2*p1*p4;
#! 2*p1*p4
SetAbbreviation( s14, "s14" );
rel1 := List( ExternalMomenta( LD ), p -> p^2 );
#! [ p1^2, p2^2, p4^2 ]
rel2 := [ (p1+p2+p4)^2 ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD,
        [ k1^2, k1*p1, k1*p2, k1*p4, s12, s14 ] );
SetPropagators( LD, -[ k1^2, (k1-p1)^2, (k1-p1-p2)^2, (k1+p4)^2 ] );
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
#! <A non-zero 4 x 12 matrix over an external ring>
Sred := ReducedBasisOfColumnModule( S );
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
#! @EndExample

Y := RationalDoubleShiftAlgebra( R );
mibps := Y * ibps;

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


#homalgDisplay( [ "subs(a1=1,a2=1,a3=1,a4=1,", EvalRingElement( DecideZeroRows( b, mibps )[1,1] ), ")" ], Q );
# Die GB ist nur wichtig fuer die Reduktion, aber sie würde sehr allgemeine Reduktionen erlauben, wir sind aber nur an den Reduktionen interessiert.
# Reduziere nur die Di's, Di_'s und die Ni's

d := DimensionOfCoefficientsVector( LD ); R := RingOfLoopDiagram( LD );
jan := HomalgMatrix( [ 1,0,0,0,  1,-1,0,0, 1,-1,-1,0,  1,0,0,1 ], d, d, R );
jbps := MatrixOfIBPRelations( jan, LD );
mjbps := P * jbps;
