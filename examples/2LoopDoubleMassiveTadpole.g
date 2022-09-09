#! @Chunk 2LoopDoubleMassiveTadpole

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "" : masses := "m" );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [  ] & masses [ m ]>
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2 ] );
SetPropagators( LD, -[ l1^2, l2^2 - m^2, (l1+l2)^2 - m^2 ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ ] );
R := RingOfLoopDiagram( LD );
#! Q[m,D][D1,D2,D3]
ibps := MatrixOfIBPRelations( LD );
#! <A 4 x 1 matrix over a residue class ring>
Ypol := HomalgRing( ibps );
#! Q[m,D][a1,a2,a3]<D1,D1_,D2,D2_,D3,D3_>/( D3*D3_-1, D2*D2_-1, D1*D1_-1 )
BasisOfRows( ibps );
#! <A non-zero 18 x 1 matrix over a residue class ring>
#! @EndExample

Y := RationalDoubleShiftAlgebra( R );
#! Q(m,D)(a1,a2,a3)<D1,D1_,D2,D2_,D3,D3_>/( D1*D1_-1, D2*D2_-1, D3*D3_-1 )
mibps := Y * ibps;
#! <A 4 x 1 matrix over a residue class ring>
mbas := BasisOfRows( mibps );
#! <A non-zero 6 x 1 matrix over a residue class ring>

lhs1 := HomalgMatrix( "[a1*D1_]", 1, 1, Y );
R1 := lhs1 - DecideZeroRows( lhs1, mbas );
lhs2 := HomalgMatrix( "[a2*D2_]", 1, 1, Y );
R2 := lhs2 - DecideZeroRows( lhs2, mbas );
lhs3 := HomalgMatrix( "[a3*D3_]", 1, 1, Y );
R3 := lhs3 - DecideZeroRows( lhs3, mbas );

Ris := UnionOfRows( [ R1, R2, R3 ] );

b1 := HomalgMatrix( "[2*m^2*(d-2*a1-2)*(-a3-a2-2*a1+d)*(d-2*a1-a2-a3-1)*a1*D1_]", 1, 1, Ypol );
RHS1 := HomalgMatrix( "[ -a1*(d-2*a1-2*a3)*(d-2*a1-2*a2)*(d-a1-a2-a3) ]", 1, 1, Ypol );
NF1 := b1 - RHS1;
Assert( 0, IsZero( DecideZeroRows( NF1, ibps ) ) );
b2 := HomalgMatrix( "[2*m^2*(-a3-a2-2*a1+d)*a2*D2_]", 1, 1, Ypol );
RHS2 := HomalgMatrix( "[ -(d-2*a1-2*a2)*(d-a1-a2-a3) ]", 1, 1, Ypol );
NF2 := b2 - RHS2;
Assert( 0, IsZero( DecideZeroRows( NF2, ibps ) ) );
b3 := HomalgMatrix( "[2*m^2*(-a3-a2-2*a1+d)*a3*D3_]", 1, 1, Ypol );
RHS3 := HomalgMatrix( "[ -(d-2*a1-2*a3)*(d-a1-a2-a3) ]", 1, 1, Ypol );
NF3 := b3 - RHS3;
Assert( 0, IsZero( DecideZeroRows( NF3, ibps ) ) );

nf := UnionOfRows( NF1, NF2, NF3 );
Assert( 0, IsZero( DecideZeroRows( nf, ibps ) ) );
