#! @Chunk 1LoopMassiveOnShell

LoadPackage( "LoopIntegrals" );

#! @Example
LD := LoopDiagram( "k", "q" : masses := "m" );
#! <A loop diagram with loop momenta [ k ] & external momenta [ q ] &
#!  masses [ m ]>
rel := [ q^2 - m^2 ];
#! [ q^2+(-m^2) ]
SetRelationsOfExternalMomenta( LD, rel );
SetIndependentLorentzInvariants( LD, [  ] );
SetPropagators( LD, -[ k^2, k^2 + 2*q*k ] );
SetNumerators( LD, -[ ] );
SetExtraLorentzInvariants( LD, [ ] );
R := RingOfLoopDiagram( LD );
#! Q[m,d][D1,D2]
ibps := MatrixOfIBPRelations( LD );
#! <A 2 x 1 matrix over a residue class ring>
Display( ibps );
#! -a2*D1*D2_+d-2*a1-a2,
#! 2*m^2*a2*D2_-a1*D1_*D2+a2*D1*D2_+a1-a2
#! 
#! modulo [ D2*D2_-1, D1*D1_-1 ]
ibp1 := ibps[1,1];
#! |[ -a2*D1*D2_+d-2*a1-a2 ]|
ViewList( DecomposeInMonomials( ibp1 ) );
#! [ [ |[ -a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ d-2*a1-a2 ]|, |[ 1 ]| ] ]
ibp2 := ibps[2,1];
#! |[ 2*m^2*a2*D2_-a1*D1_*D2+a2*D1*D2_+a1-a2 ]|
ViewList( DecomposeInMonomials( ibp2 ) );
#! [ [ |[ -a1 ]|, |[ D1_*D2 ]| ],
#!   [ |[ a2 ]|, |[ D1*D2_ ]| ],
#!   [ |[ 2*m^2*a2 ]|, |[ D2_ ]| ],
#!   [ |[ a1-a2 ]|, |[ 1 ]| ] ]
ibp2_red := DecideZero( ibp2, ibps{[ 1 ]} );
#! |[ 2*m^2*a2*D2_-a1*D1_*D2+d-a1-2*a2 ]|
ViewList( DecomposeInMonomials( ibp2_red ) );
#! [ [ |[ -a1 ]|, |[ D1_*D2 ]| ],
#!   [ |[ 2*m^2*a2 ]|, |[ D2_ ]| ],
#!   [ |[ d-a1-2*a2 ]|, |[ 1 ]| ] ]
Ypol := HomalgRing( ibp1 );
#! Q[m,d][a1,a2]<D1,D1_,D2,D2_>/( D2*D2_-1, D1*D1_-1 )
gen := GeneratorsOfScalelessSectors( LD );
#! <A 1 x 1 matrix over an external ring>
Display( gen );
#! D2
#! @EndExample

# Below we use the Minkowski coordinates

F := BaseRing( R );
S := F["x1..3"];
ExportVariables( S );
subs := RingMap( [ m, d, x1, x2+m, x3+m ], S, S );
G := Pullback( subs, HomalgMatrix( "[ 2*x1, x2-x1-x3, x2-x1-x3, 2*x3 ]", 2, 2, S ) );
B := HomalgMatrix( [ Determinant( G ) ], 1, 1, S );
Bi := Diff( HomalgMatrix( RelativeIndeterminatesOfPolynomialRing( S ), 3, 1, S ), B );
M := HomalgMatrix( [ B[1,1], x1 * Bi[1,1], x2 * Bi[2,1], x3 * Bi[3,1] ], 4, 1, S );
