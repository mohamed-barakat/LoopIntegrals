#! @Chunk 3LoopBanana

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..3", "" : masses := "m1,m3" );
#! <A loop diagram with loop momenta [ l1, l2, l3 ] &
#!  external momenta [  ] & masses [ m1, m3 ]>
SetRelationsOfExternalMomenta( LD, [ ] );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l3^2, l1*l2, l1*l3, l2*l3 ] );
SetPropagators( LD, -[ l1^2 - m1^2, l2^2 - m1^2, l3^2 - m3^2, (l1+l2+l3)^2 - m3^2 ] );
SetNumerators( LD, -[ (l1+l2)^2, (l1+l3)^2 ] );
SetExtraLorentzInvariants( LD, [ ] );
ibps := MatrixOfIBPRelations( LD );
#! <A 9 x 1 matrix over a residue class ring>
Y := HomalgRing( ibps );
#! Q[m1,m3,D][a1,a2,a3,a4,a5,a6]<D1,D1_,D2,D2_,D3,D3_,D4,D4_,N5,N5_,N6,N6_>/
#!  ( N6*N6_-1, N5*N5_-1, D4*D4_-1, D3*D3_-1, D2*D2_-1, D1*D1_-1 )
#! @EndExample

Q := HomalgFieldOfRationalsInMaple();
P := Q * List( Indeterminates( BaseRing( BaseRing( Y ) ) ), String );
P := P * List( RelativeIndeterminateCoordinatesOfDoubleShiftAlgebra( Y ), String );
P := DoubleShiftAlgebra( P, List( IndeterminateShiftsOfDoubleShiftAlgebra( Y ), String ) : pairs := true, steps := -1 );
mibps := P * ibps;
