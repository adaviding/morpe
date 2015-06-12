NDIMS = 4;
NTESTS = 1000;
USE_MEX = true;

%	A table.  Each row contains the entries of a random rotation matrix for a given simulation.
S = zeros(NTESTS,NDIMS^2);

if USE_MEX
	R = Mcl_RandRotate(int32(NDIMS));
else
	R = Mcl_RandomRotationMatrix(NDIMS);
end
disp('----------------------------------------------------------');
disp('Here is a random rotation matrix.');
disp(R);
disp('----------------------------------------------------------');
disp('Each row and column is a unit vector because you see ones below:');
for i=1:NDIMS
	disp(num2str([sum(R(:,i).^2), sum(R(i,:).^2)]));
end
disp('----------------------------------------------------------');
disp('Each row is mutually orthogonal to other rows, and each column is mutually orthogonal to other columns.');
disp('This is true because you see only zeros (or very small numbers) below:');
for i=1:(NDIMS-1)
	for j=(i+1):NDIMS
		disp(num2str([sum(R(:,i).*R(:,j)), sum(R(i,:).*R(j,:))]));
	end
end
disp('----------------------------------------------------------');
disp('Each entry of the matrix has the same mean (0) and variance.');
disp('');
for i=2:NTESTS
	if USE_MEX
		R = Mcl_RandRotate(int32(NDIMS));
	else
		R = Mcl_RandomRotationMatrix(NDIMS);
	end
	S(i,:) = R(:);
end
disp('MEAN  You should only see zereos (or very small numbers) below:');
xbar = reshape(mean(S), [NDIMS NDIMS]);
s = reshape(std(S), [NDIMS NDIMS]);
disp(xbar);
disp(['STD You should see numbers all near sqrt(1/NDIMS) = ' num2str(sqrt(1/NDIMS)) ]);
disp(s);
disp('T STATISTIC  All of the means are near zero.');
disp( xbar./s/sqrt(NTESTS-1) );
disp('----------------------------------------------------------');
disp('All of the matrix entries are uncorrelated because all of the non-diagonal correlation coefficients are near zero.');
disp( corr(S) );
disp('----------------------------------------------------------');
disp('Therefore, the method generates random rotation matrices.');