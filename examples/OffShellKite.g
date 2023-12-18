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

Qa := FieldOfCoefficientsOfLoopDiagramInHecke( LD );
#prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa : homalgIOMode := "d", reduced_basis_of_relative_syzygies := true );

## RowEchelonForm 588 x 546 : Z(m,d,s,a1,a2,a3,a4,a5)

Q := CoefficientsRing( AmbientRing( Y ) );
# m := Q * prel2[1];
# b := BasisOfRows( m );
# homalgDisplay( [ "map(factor,", b, ")" ] );

#prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa : homalgIOMode := "d", basis_of_relative_syzygies := false );

## RowEchelonForm 672 x 546 : Z(m,d,s,a1,a2,a3,a4,a5)
## runs out of memory within a day

#prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2, Qa : homalgIOMode := "d", basis_of_relative_syzygies := true );

## RowEchelonForm 756 x 796 : Z(m,d,s,a1,a2,a3,a4,a5)

## prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2 );
## Old Hecke: Manually killed after 488G 2170h

## prel2 := MatrixOfCoefficientsOfParametricIBPs( LD, 2 );
## New Hecke: finished after 342G/156G 327h with the last RowEchelonForm step
## m;
## <A non-zero 21 x 46 matrix over an external ring>

## prel2[2] :=
##   [ D2*D4*D5, D3*D4*D5, D4^2*D5, D1*D5^2, D2*D5^2, D3*D5^2, D4*D5^2,
##      D1^2, D1*D2, D2^2, D1*D3, D2*D3, D3^2, D1*D4, D2*D4, D4^2, D1*D5, D2*D5, D3*D5, D4*D5, D5^2, D3*D4,
##      D4, D3, D1, D5, D2, 1,
##      D5_, D4_, D3_, D2_, D1_, D4_*D5_, D3_*D5_, D2_*D5_, D1_*D5_, D4_^2, D3_*D4_, D2_*D4_, D1_*D4_, D3_^2, D2_*D3_, D1_*D3_, D2_^2, D1_*D2_ ]
