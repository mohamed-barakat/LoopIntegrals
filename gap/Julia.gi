##
InstallOtherMethod( LoopDiagram,
        [ IsJuliaObject, IsJuliaObject, IsInt ], 10001,
        
  function( L, K, dim )
    local LD;
    
    LD := LoopDiagram( JuliaToGAP( IsString, L ), JuliaToGAP( IsString, K ), dim );
    
    Perform( LoopMomenta( LD ), function( a ) JuliaEvalString( Concatenation( Name( a ), " = GAP.Globals.", Name( a ) ) ); end );
    Perform( ExternalMomenta( LD ), function( a ) JuliaEvalString( Concatenation( Name( a ), " = GAP.Globals.", Name( a ) ) ); end );
    
    return LD;
    
end );

##
InstallOtherMethod( LoopDiagram,
        [ IsJuliaObject, IsJuliaObject ], 10001,
        
  function( L, K )
    
    return LoopDiagram( L, K, LOOP_INTEGRALS.Dimension );
    
end );

##
InstallOtherMethod( SetAbbreviation,
        [ IsHomalgRingElement, IsJuliaObject ],
        
  function( xy, str )
    
    xy!.Abbreviation := JuliaToGAP( IsString, str );
    
end );

##
InstallMethod( SetRelationsOfMomenta,
        [ IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( LD, rel )
    
    LD!.RelationsOfMomenta := JuliaToGAP( IsList, rel );
    SetFilterObj( LD, HasRelationsOfMomenta );
    
end );

##
InstallMethod( SetIndependentLorentzInvariants,
        [ IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( LD, inv )
    
    LD!.IndependentLorentzInvariants := JuliaToGAP( IsList, inv );
    SetFilterObj( LD, HasIndependentLorentzInvariants );
    
end );

##
InstallMethod( SetPropagators,
        [ IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( LD, Ds )
    
    LD!.Propagators := JuliaToGAP( IsList, Ds );
    SetFilterObj( LD, HasPropagators );
    
end );

##
InstallMethod( SetNumerators,
        [ IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( LD, Ns )
    
    LD!.Numerators := JuliaToGAP( IsList, Ns );
    SetFilterObj( LD, HasNumerators );
    
end );

##
InstallMethod( SetExtraLorentzInvariants,
        [ IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( LD, inv )
    
    LD!.ExtraLorentzInvariants := JuliaToGAP( IsList, inv );
    SetFilterObj( LD, HasExtraLorentzInvariants );
    
end );

##
InstallOtherMethod( IBPRelation,
        [ IsHomalgMatrix, IsLoopDiagram, IsJuliaObject ], 10000000001,
        
  function( row, LD, vec )
    
    return IBPRelation( row, LD, JuliaToGAP( IsList, vec ) );
    
end );
