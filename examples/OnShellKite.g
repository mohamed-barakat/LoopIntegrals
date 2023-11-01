#! @Chunk OnShellKite

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "p" );
#! <A loop diagram with loop momenta [ l1, l2 ] & external momenta [ p ]>
s := p^2;;
SetAbbreviation( s, "s" );
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2, l1*p, l2*p, s ] );
SetPropagators( LD, -[ l1^2, (l1+p)^2-p^2, (l2+p)^2-p^2, l2^2, (l1+l2+p)^2-p^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ s ] );
#! @EndExample

R := RingOfLoopDiagram( LD );
Y := RationalDoubleShiftAlgebra( R );
A := BaseRing( Y );
Qa := HomalgRingOfIntegersInOscar( JoinStringsWithSeparator( List( Indeterminates( A ), String ) ) );
#Qa := HomalgFieldOfRationalsInOscar( JoinStringsWithSeparator( List( Indeterminates( A ), String ) ) );
prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa );

Q := CoefficientsRing( AmbientRing( Y ) );

prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2 );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );

#ibps := MatrixOfIBPRelations( LD );
#<A 6 x 1 matrix over a residue class ring>
#BasisOfRows( ibps ); # killed after 870 hours with 62MB memory consumption
