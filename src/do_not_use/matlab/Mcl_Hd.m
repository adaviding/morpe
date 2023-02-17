function out = Mcl_Hd(dPrime)

% function out = Mcl_Hd(dPrime)
%	Returns the conditional entropy as a function of d'.

h = dPrime/2;
dz = 0.01;
z = -(4+h):dz:(4+h);

a = normpdf(z-h);
b = normpdf(z+h);
ab = a+b;

ind = find(ab>0);
a=a(ind);
b=b(ind);
ab=ab(ind);

out = -sum(0.5*ab.*log2( max(a,b)./ab )) * dz;