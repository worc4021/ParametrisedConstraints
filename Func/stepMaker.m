function [stepsize,idx] = stepMaker(wMin,uMin,lambdaMin,kappaMin)


if isempty(wMin)
    wMin = inf;
end
if isempty(uMin)
    uMin = inf;
end
if isempty(lambdaMin)
    lambdaMin = inf;
end
if isempty(kappaMin)
    kappaMin = inf;
end

[stepsize,idx] = min([wMin,uMin,lambdaMin,kappaMin]);