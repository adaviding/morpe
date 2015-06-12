function out = Mcl_RandVec(setSize)

%	This function takes the set [1, 2, ..., setSize] and randomizes its order.
%	This function is called with the syntax Mcl_RandVec(setSize).
%	This function returns a vector of length setSize with entries in the range [1, setSize].

pickFrom = [1:setSize];
out = zeros(1,setSize);

index = 0;

pickFromSize = setSize;

while pickFromSize > 0
   
   index = index + 1;
   
   randomIndex = ceil(rand*pickFromSize);
      
   out(index) = pickFrom(randomIndex);
   
   if randomIndex ~= 1
      tempPickFrom = pickFrom(1:(randomIndex-1));
      if randomIndex ~= pickFromSize
         tempPickFrom = [tempPickFrom, pickFrom((randomIndex+1):pickFromSize)];
      end
   elseif pickFromSize > 1
      tempPickFrom = pickFrom((randomIndex+1):pickFromSize);
   else
      tempPickFrom = [];
   end
   
   pickFrom = tempPickFrom;
      
   pickFromSize = pickFromSize - 1;
   
end