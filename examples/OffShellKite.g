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

## prel2 := ParametricIBPs( LD, 2 );
## finished after 328G 327h with the last RowEchelonForm step
## m;
## <A non-zero 21 x 46 matrix over an external ring>

## mat[5]{range} =
## [ D1_*D2_, D2_^2, D1_*D3_, D2_*D3_, D3_^2, D1_*D4_, D2_*D4_, D3_*D4_, D4_^2, D1_*D5_, D2_*D5_, D3_*D5_, D4_*D5_, D1_, D2_, D3_, D4_, D5_, 1, D2, D5, D1, D3, D4, D3*D4, D5^2, D4*D5, D3*D5, D2*D5, D1*D5, D4^2, D2*D4, D1*D4, D3^2, D2*D3, D1*D3, D2^2, D1*D2, D1^2, D4*D5^2, D3*D5^2, D2*D5^2, D1*D5^2, D4^2*D5, D3*D4*D5, D2*D4*D5 ]
