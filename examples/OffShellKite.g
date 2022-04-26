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
Ypol := HomalgRing( ibps );

R := RingOfLoopDiagram( LD );
#! Q[m,d,s][D1,D2,D3,D4,D5]
Y := RationalDoubleShiftAlgebra( R );
mibps := Y * ibps;

Q := CoefficientsRing( AmbientRing( Y ) );

# prel2 := ParametricIBPs( LD, 2 );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );

## prel2 := ParametricIBPs( LD, 2 );
## Manually killed after 488G 2170h
