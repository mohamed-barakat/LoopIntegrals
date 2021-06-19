#! @Chunk 1LoopDoubleMassiveBubble

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1", "k1" : masses := "m" );
#! <A loop diagram with loop momenta [ l1 ] & external momenta [ k1 ] &
#!  masses [ m ]>
s := k1*k1;
#! k1^2
SetAbbreviation( s, "s" );
rel := [ ];;
SetRelationsOfExternalMomenta( LD, rel );
SetIndependentLorentzInvariants( LD, [ ] );
SetPropagators( LD, -[ l1^2 + m^2, ( l1 + k1 )^2 + m^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
R := RingOfLoopDiagram( LD );
#! Q[m,d,s][D1,D2]
ibps := MatrixOfIBPRelations( LD );
#! <A 2 x 1 matrix over a residue class ring>
ibp1 := ibps[1,1];
#! |[ -2*m^2*a1*D1_-2*m^2*a2*D2_-s*a2*D2_-a2*D1*D2_+d-2*a1-a2 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ -2*m^2*a1 ]|, |[ D1_ ]| ],
#!   [ |[ -2*m^2*a2-s*a2 ]|, |[ D2_ ]| ],
#!   [ |[ d-2*a1-a2 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp1 );
#! Q[m,d,s][a1,a2]<D1,D1_,D2,D2_>/( D2*D2_-1, D1*D1_-1 )
bas := BasisOfRows( ibps );
#! <A non-zero 6 x 1 matrix over a residue class ring>
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 1 matrix over an external ring>
Display( gen );
#! D1*D2
gen2 := GeneratorsOfScalelessSectors( LD, [ 2, 2 ] );
#! <An unevaluated 1 x 1 matrix over an external ring>
Display( gen2 );
#! D1^2*D2^2
#! @EndExample

Q := HomalgFieldOfRationalsInMaple();
P := Q * List( Indeterminates( BaseRing( BaseRing( Y ) ) ), String );
P := P * List( RelativeIndeterminateCoordinatesOfDoubleShiftAlgebra( Y ), String );
P := DoubleShiftAlgebra( P, List( IndeterminateShiftsOfDoubleShiftAlgebra( Y ), String ) : pairs := true, steps := -1 );
mibps := P * ibps;
