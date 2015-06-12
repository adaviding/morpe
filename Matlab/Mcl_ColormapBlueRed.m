function out = Mcl_ColormapBlueRed( nrows )

%	function out = Mcl_ColormapBlueRed( nrows )
% Returns a vector of length nrows and 3 colomns suitable as input for
%	matlab's COLORMAP funciton.
% The map transitions linearly from blue to red with gray in the middle.

if( nrows < 8 )
	error('nrows must be >= 8');
end

out = repmat( (0.5:nrows)'/nrows, [1 3] );
out(:,3) = 1-out(:,3);
out(:,2) = min(out(:,1), out(:,3));