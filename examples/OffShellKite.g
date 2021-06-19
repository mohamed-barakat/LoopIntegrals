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

Q := HomalgFieldOfRationalsInMaple();

# prel2 := ParametricIBPs( LD, 2 );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );
