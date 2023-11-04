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

Qa := FieldOfCoefficientsOfLoopDiagramInHecke( LD );
prel2 := ColumnReversedMatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa : homalgIOMode := "d" );
#! [ <A non-zero 5 x 10 matrix over an external ring>,
#!   [ D1_, D2_, D3_, D4_, D5_, 1, D3, D2, D1, D5 ] ]

## computing prel2 takes 30 minutes and 23GB

Q := CoefficientsRing( AmbientRing( Y ) );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );

#ibps := MatrixOfIBPRelations( LD );
#<A 6 x 1 matrix over a residue class ring>
#BasisOfRows( ibps ); # killed after 870 hours with 62MB memory consumption
