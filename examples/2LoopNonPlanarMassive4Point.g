#! @Chunk 2LoopNonPlanarMassive4Point

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "k1..3" : masses := "mt2" );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [ k1, k2, k3 ] & masses [ mt2 ] >
mW2 := k3^2;
#! k3^2
SetAbbreviation( mW2, "mW2" );
t := mW2 - 2 * k1 * k3;
#! -2*k1*k3+k3^2
SetAbbreviation( t, "t" );
u := mW2 - 2 * k2 * k3;
#! -2*k2*k3+k3^2
SetAbbreviation( u, "u" );
rel1 := List( ExternalMomenta( LD ){[1 .. 2]}, k -> k^2 );
#! [ k1^2, k2^2 ]
rel2 := [ ( k1 + k2 - k3 )^2 - mt2, 2 * k1 * k2 - mt2 - mW2 + t + u ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2 ] );
SetPropagators( LD, [ l1^2 - mt2, l2^2 - mt2, ( l1 + l2 )^2, ( l1 + k1 )^2 - mt2, ( l1 + l2 + k1 + k2 - k3 )^2 - mt2, ( l2 + k2 - k3 )^2, ( l2 - k3 )^2 ] );
SetNumerators( LD, [ ( l1 - k1 + k3 )^2, ( l2 - k1 - k2 + k3 )^2 ] );
SetExtraLorentzInvariants( LD, [ t, u, mW2 ] );
ibps := MatrixOfIBPRelations( LD );
#! <A 10 x 1 matrix over a residue class ring>
#! @EndExample


