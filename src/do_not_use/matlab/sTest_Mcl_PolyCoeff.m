Rank = int32(3);
Ndims = int32(5);

disp('------------------------------Testing Mcl_PolyCoeff');
disp(' ');
disp(['The properties of a polynomial defined by [Rank, Ndims]=[' num2str(Rank) ',' num2str(Ndims) '] are as follows:']);
disp(' ');
if false
	Ncoeff = Mcl_PolyCoeff(Rank,Ndims);
	disp(['Ncoeff=' num2str(Ncoeff) ' is the number of coefficients in the polynomial (not including the additive constant that is required in some formulations).']);
elseif false
	[Ncoeff, N] = Mcl_PolyCoeff(Rank,Ndims);
	disp(['Ncoeff=' num2str(Ncoeff) ' is the number of coefficients in the polynomial (not including the additive constant that is required in some formulations).']);
	disp(' ');
	disp(['N=[' num2str(reshape(N, [1, numel(N)])) '] is the number of coefficients at each rank (1 to Rank)']);
	disp(['sum(N)=' num2str(sum(N)) ' should be the same as Ncoeff.']);
else
	[Ncoeff, N, CoeffDims] = Mcl_PolyCoeff(Rank,Ndims);
	disp(['Ncoeff=' num2str(Ncoeff) ' is the number of coefficients in the polynomial (not including the additive constant that is required in some formulations).']);
	disp(' ');
	disp(['N=[' num2str(reshape(N, [1, numel(N)])) '] is the number of coefficients at each rank (1 to Rank)']);
	disp(['sum(N)=' num2str(sum(N)) ' should be the same as Ncoeff.']);
	disp(' ');
	disp('CoeffDims (int32 2D array: Ncoeff*Rank), gives the dimensional expansion of each polynomial coefficient.  Dimensions are 0-based, -1 implies empty.');
	disp(CoeffDims);
	ct = 0;
	for iCol=1:(Rank-1)
		ct = ct + sum(CoeffDims(:,iCol+1)>CoeffDims(:,iCol));
	end
	if ct>0
		error('The CoeffDims parameter has invalid values.  Mcl_PolyCoeff failed the test.');
	end
end
disp(' ');
disp('------------------------------Finished Testing Mcl_PolyCoeff');