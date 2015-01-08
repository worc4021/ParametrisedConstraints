function retVal = altIsContained(A1,b1,A2,b2)
% Checks whether A1*x<=b1 is contained in A2*x<=b2.
% 
if and(~isempty(A1),~isempty(A2))
    V = vertexCompute(A1,b1);
    retVal = 1;
    for i = 1:size(V,2)
        retVal = and(retVal,CheckLinCons(V(i,:)',A2,b2));
    end
else
    retVal = 0;
end