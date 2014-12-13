function out = isContained(A1,b1,A2,b2)
% Checks whether A1*x<=b1 is contained in A2*x<=b2.
% 


% [A2,b2] = linReduce(A2,b2,0);
if length(b2)>=1
out = 1;

    for i = 1:length(b2)
        [x,fval] = cplexlp(-A2(i,:),A1,b1);
        if isempty(fval)
            out = 0;
        elseif -fval - 1e-9 > b2(i)
            out = 0;
        end
    end
else
    out = 0;
end