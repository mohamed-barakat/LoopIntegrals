#! @Chunk OffShellKite

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p" : masses := "m" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p ]>
s := p^2;;
SetAbbreviation( s, "s" );
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2, l1*p, l2*p, s ] );
SetPropagators( LD, -[ l1^2, (l1+p)^2-m^2, l2^2-m^2, (l2-p)^2, (l1+l2)^2-m^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
#! @EndExample

ibps := MatrixOfIBPRelations( LD );
Y := HomalgRing( ibps );
Q := HomalgFieldOfRationalsInMaple();
P := Q * List( Indeterminates( BaseRing( BaseRing( Y ) ) ), String );
P := P * List( RelativeIndeterminateCoordinatesOfDoubleShiftAlgebra( Y ), String );
P := DoubleShiftAlgebra( P, List( IndeterminateShiftsOfDoubleShiftAlgebra( Y ), String ) : pairs := true, steps := -1 );
mibps := P * ibps;

# prel2 := ParametricIBPs( LD, 2 );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );

## prel2 := ParametricIBPs( LD, 2 );
## Manually killed after 488G 2170h
