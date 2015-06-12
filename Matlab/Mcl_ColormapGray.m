function out = Mcl_ColormapGray( nrows )

%	function out = Mcl_ColormapGray( nrows )
% Returns a vector of length nrows and 3 colomns suitable as input for
%	matlab's COLORMAP funciton.
% The map transitions linearly from black to white.

if( nrows < 8 )
	error('nrows must be >= 8');
end

out = repmat( (0.5:nrows)'/nrows, [1 3] );
