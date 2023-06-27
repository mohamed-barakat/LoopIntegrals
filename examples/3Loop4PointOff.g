#! @Chunk 3Loop4PointOff

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "l1..3", "k1..2,k4" );
#! <A loop diagram with loop momenta [ l1, l2, l3 ] &
#!  external momenta [ k1, k2, k4 ]>
s := 2*k1*k2;
#! 2*k1*k2
SetAbbreviation( s, "s" );
t := -s-2*k2*k4;
#! -2*k1*k4-2*k2*k4
SetAbbreviation( t, "t" );
mm := k4^2;
#! k4^2
SetAbbreviation( mm, "mm" );
rel1 := List( ExternalMomenta( LD ){[1 .. 2]}, k -> k^2 );
#! [ k1^2, k2^2 ]
rel2 := [ (k1+k2+k4)^2, 2*k1*k4 - t + mm ];;
SetRelationsOfExternalMomenta( LD, Concatenation( rel1, rel2 ) );
SetIndependentLorentzInvariants( LD, [ l1^2, l2^2, l3^2, l1*l2, l1*l3, l2*l3 ] );
SetPropagators( LD, [ l1^2, (l1-k1)^2, (l1-k1-k2)^2, l3^2, (l3+k1+k2)^2, (l1+l3)^2, (l2-l3)^2, l2^2, (l2-k4)^2, (l2+k1+k2)^2 ] );
SetNumerators( LD, [ (l1+k4)^2, (l2+k1)^2, (l3+k1)^2, (l3+k4)^2, (l1+l2)^2 ] );
SetExtraLorentzInvariants( LD, [ s, t, mm ] );
ibps := MatrixOfIBPRelations( LD );
#! <A 18 x 1 matrix over a residue class ring>
# sibps := MatrixOfSpecialIBPRelations( LD );
# <A 825 x 1 matrix over a residue class ring>
#! @EndExample


## sibps took 293h and 21GB

## [ [ 2*D1, -D1-D8+N15, -D1-D4+D6, D1-D2, s+D2-D3, -mm-D1+N11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
##   [ D1+D2, -D1-N12+N15, -D1+D6-N13, D1-D2, D2-D3, -t-D1+N11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
##   [ -s+D1+D3, s-D1-D10+N15, s-D1-D5+D6, -s+D1-D2, D2-D3, s-D1+N11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
##   [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -D1-D4+D6, D4-D7+D8, 2*D4, -D4+N13, -s+D5-N13, -mm-D4+N14 ],
##   [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s-D3-D4+D6, -s+D4-D7+D10, -s+D4+D5, s-D4+N13, D5-N13, -s-2*mm-D4+N14 ],
##   [ D1-D4+D6, -D1+D4-D7+N15, -D1+D4+D6, D1-D2-D4+N13, D2-D3+D5-N13, -2*mm-D1-D4+N11+N14, 0, 0, 0, 0, 0, 0, D1-D4+D6, -D1+D4-D7+N15, -D1+D4+D6, D1-D2-D4+N13,
##     D2-D3+D5-N13, -2*mm-D1-D4+N11+N14 ],
##   [ 0, 0, 0, 0, 0, 0, D4-D6-D8+N15, -D4+D7+D8, -D4-D7+D8, D4-D8+N12-N13, -D5+D10-N12+N13, 2*mm+D4+D8-D9-N14, -D4+D6+D8-N15,
##     D4-D7-D8, D4+D7-D8, -D4+D8-N12+N13, D5-D10+N12-N13, -2*mm-D4-D8+D9+N14 ],
##   [ 0, 0, 0, 0, 0, 0, -D1-D8+N15, 2*D8, D4-D7+D8, -D8+N12, -s+D10-N12, mm+D8-D9, 0, 0, 0, 0, 0, 0 ],
##   [ 0, 0, 0, 0, 0, 0, mm-D8-N11+N15, -mm+D8+D9, mm+2*D4-D7+D8-N14, -t+mm-D8+N12, t+D10-N12, -mm+D8-D9, 0, 0, 0, 0, 0, 0 ],
##   [ 0, 0, 0, 0, 0, 0, s-D3-D8+N15, -s+D8+D10, -s+D5-D7+D8, s-D8+N12, D10-N12, -s+D8-D9, 0, 0, 0, 0, 0, 0 ] ]
