function [uIdx,wIdx,stepsize] = linesearch(theta0,thetaE,actU,actW,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s)

uIdx = reshape(actU',numel(actU),1);
wIdx = reshape(actW',numel(actW),1);


[wVar,wCon,xpVar,xpCon,lambdaVar,lambdaCon,lambdaD,muVar,muCon,...
    muD,zetaVar,zetaCon,zetaD,uVar,uCon,kappaVar,kappaCon,kappaD,...
    rhoVar,rhoCon,rhoD] = forwardsim(wM,wm,xpM,xpm,lambdaM,lambdam,...
    lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,...
    kappaZ,rhoK,rhok,rhoZ,s);

wCandNum = kron(eye(s.N),s.G)*(wVar*theta0 + wCon) - ones(s.N*size(s.G,1),1) - repmat([zeros(size(s.G,1),s.nX+s.nAlpha),ones(size(s.G,1),s.nBeta)]*theta0,s.N,1);
wCandDen = kron(eye(s.N),s.G)*wVar*(thetaE-theta0) - repmat([zeros(size(s.G,1),s.nX+s.nAlpha),ones(size(s.G,1),s.nBeta)]*(thetaE-theta0),s.N,1);

uCandNum = kron(eye(s.N),s.F)*(uVar*theta0 + uCon) - ones(s.N*size(s.F,1),1) - repmat([zeros(size(s.F,1),s.nX),ones(size(s.F,1),s.nAlpha),zeros(size(s.F,1),s.nBeta)]*theta0,s.N,1);
uCandDen = kron(eye(s.N),s.F)*uVar*(thetaE-theta0) - repmat([zeros(size(s.F,1),s.nX),ones(size(s.F,1),s.nAlpha),zeros(size(s.F,1),s.nBeta)]*(thetaE-theta0),s.N,1);

lambdaCandNum = lambdaVar*theta0+lambdaCon;
lambdaCandDen = lambdaVar*(thetaE-theta0);

kappaCandNum = kappaVar*theta0+kappaCon;
kappaCandDen = kappaVar*(thetaE-theta0);

wNEQ = find(wCandDen>s.tol);
uNEQ = find(uCandDen>s.tol);
lambdaNEQ = find(and(lambdaCandDen<-s.tol,lambdaCandNum>s.tol));
kappaNEQ = find(and(kappaCandNum>s.tol,kappaCandDen<-s.tol));

[wCand,wCandIDX] = min(-wCandNum(wNEQ)./wCandDen(wNEQ));
[uCand,uCandIDX] = min(-uCandNum(uNEQ)./uCandDen(uNEQ));
[lambdaCand,lambdaCandIDX] = min(-lambdaCandNum(lambdaNEQ)./lambdaCandDen(lambdaNEQ));
[kappaCand,kappaCandIDX] = min(-kappaCandNum(kappaNEQ)./kappaCandDen(kappaNEQ));

[stepsize,sel] = stepMaker(wCand,uCand,lambdaCand,kappaCand);

if sel==1
    wIdx(wNEQ(wCandIDX)) = 1;
elseif sel==2
    uIdx(uNEQ(uCandIDX)) = 1;
elseif sel==3
    wIdx(lambdaNEQ(lambdaCandIDX)) = 0;
elseif sel==4
    uIdx(kappaNEQ(kappaCandIDX)) = 0;
else
    stepsize = 1e10;
end

uIdx = reshape(uIdx,numel(uIdx)/s.N,s.N)';
wIdx = reshape(wIdx,numel(wIdx)/s.N,s.N)';
