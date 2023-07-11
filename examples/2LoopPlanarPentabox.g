#! @Chunk 2LoopPlanarPentabox

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..2", "k1..4" );
#! <A loop diagram with loop momenta [ l1, l2 ] &
#!  external momenta [ k1, k2, k3, k4 ]>
s12 := 2 * k1 * k2;;
SetAbbreviation( s12, "s12" );
s23 := 2 * k2 * k3;;
SetAbbreviation( s23, "s23" );
s34 := 2 * k3 * k4;;
SetAbbreviation( s34, "s34" );
s15 := 2 * k2 * k4 + s23 + s34;
SetAbbreviation( s15, "s15" );
s45 := 2 * k1 * k3 + s12 + s23;
SetAbbreviation( s45, "s45" );
rel1 := List( ExternalMomenta( LD ), k -> k^2 );
#! [ k1^2, k2^2, k3^2, k4^2 ]
rel2 := [ ( k1 + k2 + k3 + k4 )^2, 2 * k1 * k4 + s15 - s23 + s45 ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l1*l2 ] );
SetPropagators( LD, [ l1^2, ( l1 - k1 )^2, ( l1 - k1 - k2 )^2, ( l1 - k1 - k2 - k3 )^2, l2^2, ( l2 + k1 + k2 + k3 + k4 )^2, ( l2 + k1 + k2 + k3 )^2, ( l1 + l2 )^2 ] );
SetNumerators( LD, [ ( l1 - k1 - k2 - k3 - k4 )^2, ( l2 + k1 )^2, ( l2 + k2 )^2 ] );
SetExtraLorentzInvariants( LD, [ s12, s23, s34, s15, s45 ] );
ibps := MatrixOfIBPRelations( LD );
#! <A 12 x 1 matrix over a residue class ring>
#! @EndExample


