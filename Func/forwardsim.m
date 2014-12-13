function [wVar,wCon,xpVar,xpCon,lambdaVar,lambdaCon,lambdaD,muVar,...
    muCon,muD,zetaVar,zetaCon,zetaD,uVar,uCon,kappaVar,kappaCon,kappaD,...
    rhoVar,rhoCon,rhoD] = forwardsim(wM,wm,xpM,xpm,lambdaM,lambdam,...
    lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,...
    rhoK,rhok,rhoZ,s)

N = size(wM,3);


nVar = s.nX + s.nAlpha + s.nBeta;

wVar = zeros(s.nW*N,nVar);
wCon = zeros(s.nW*N,1);

xpVar = zeros(s.nX*(N+1),nVar);
xpCon = zeros(s.nX*(N+1),1);

lambdaVar = zeros(s.nLambda,nVar);
lambdaCon = zeros(s.nLambda,1);
lambdaD = zeros(s.nLambda,s.ndelta);

muVar = zeros(s.N*s.nMu,nVar);
muCon = zeros(s.N*s.nMu,1);
muD = zeros(s.N*s.nMu,s.ndelta);


zetaVar = zeros(s.nZeta*(N+1),nVar);
zetaCon = zeros(s.nZeta*(N+1),1);
zetaD = zeros(s.nZeta*(N+1),s.ndelta);

uVar = zeros(s.nU*N,nVar);
uCon = zeros(s.nU*N,1);

kappaVar = zeros(s.nKappa*N,nVar);
kappaCon = zeros(s.nKappa*N,1);
kappaD = zeros(s.nKappa*N,s.ndelta);

rhoVar = zeros(s.nRho*N,nVar);
rhoCon = zeros(s.nRho*N,1);
rhoD = zeros(s.nRho*N,s.ndelta);

xpVar(1:s.nX,:) = [eye(s.nX),zeros(s.nAlpha+s.nBeta)];

zetaVar(1:s.nZeta,:) = zetaM(:,[1:s.nX,s.nX+s.nU+1:end],1);
zetaCon(1:s.nZeta) = zetam(:,:,1);
zetaD(1:s.nZeta,:) = zetaZ(:,:,1);


for i = 1:N
    T = [xpVar((i-1)*s.nX+1:i*s.nX,:);
        zeros(s.nAlpha+s.nBeta,s.nX),eye(s.nAlpha+s.nBeta)];
    
    uVar((i-1)*s.nU+1:i*s.nU,:) = uK(:,:,i)*T;% + kappaZ*zetaVar((i-1)*s.nZeta+1:i*s.nZeta);
    uCon((i-1)*s.nU+1:i*s.nU) = uK(:,1:s.nX,i)*xpCon((i-1)*s.nX+1:i*s.nX) + uk(:,:,i);
    
    kappaVar((i-1)*s.nKappa+1:i*s.nKappa,:) = kappaK(:,:,i)*T ...
        + kappaZ(:,:,i)*zetaVar((i-1)*s.nZeta+1:i*s.nZeta,:);
    kappaCon((i-1)*s.nKappa+1:i*s.nKappa) = kappaK(:,1:s.nX,i)*xpCon((i-1)*s.nX+1:i*s.nX) ...
        + kappak(:,:,i)+ kappaZ(:,:,i)*zetaCon((i-1)*s.nZeta+1:i*s.nZeta);
    kappaD((i-1)*s.nKappa+1:i*s.nKappa,:) = kappaZ(:,:,i)*zetaD((i-1)*s.nZeta+1:i*s.nZeta,:);
    
    rhoVar((i-1)*s.nRho+1:i*s.nRho,:) = rhoK(:,:,i)*T + rhoZ(:,:,i)*zetaVar((i-1)*s.nZeta+1:i*s.nZeta,:);
    rhoCon((i-1)*s.nRho+1:i*s.nRho) = rhoK(:,1:s.nX,i)*xpCon((i-1)*s.nX+1:i*s.nX) ...
        + rhok(:,:,i)+ rhoZ(:,:,i)*zetaCon((i-1)*s.nZeta+1:i*s.nZeta);
    rhoD((i-1)*s.nRho+1:i*s.nRho,:) = rhoZ(:,:,i)*zetaD((i-1)*s.nZeta+1:i*s.nZeta,:);
    
    Th = [xpVar((i-1)*s.nX+1:i*s.nX,:);
          uVar((i-1)*s.nU+1:i*s.nU,:);
          zeros(s.nAlpha+s.nBeta,s.nX),eye(s.nAlpha+s.nBeta)];
      
    wVar((i-1)*s.nW+1:i*s.nW,:) = wM(:,:,i)*Th;
    wCon((i-1)*s.nW+1:i*s.nW) = wm(:,:,i)+wM(:,1:s.nX+s.nU,i)*[xpCon((i-1)*s.nX+1:i*s.nX);uCon((i-1)*s.nU+1:i*s.nU)];
    
    xpVar(i*s.nX+1:(i+1)*s.nX,:) = xpM(:,:,i)*Th;
    xpCon(i*s.nX+1:(i+1)*s.nX) = xpm(:,:,i)+xpM(:,1:s.nX+s.nU,i)*[xpCon((i-1)*s.nX+1:i*s.nX);uCon((i-1)*s.nU+1:i*s.nU)];
    
    lambdaVar((i-1)*s.nLambda+1:i*s.nLambda,:) = lambdaM(:,:,i)*Th ...
        + lambdaZ(:,:,i)*rhoVar((i-1)*s.nRho+1:i*s.nRho,:);
    lambdaCon((i-1)*s.nLambda+1:i*s.nLambda) = lambdam(:,:,i) ...
        + lambdaM(:,1:s.nX+s.nU,i)*[xpCon((i-1)*s.nX+1:i*s.nX);uCon((i-1)*s.nU+1:i*s.nU)] ...
        + lambdaZ(:,:,i)*rhoCon((i-1)*s.nRho+1:i*s.nRho);
    lambdaD((i-1)*s.nLambda+1:i*s.nLambda,:) = lambdaZ(:,:,i)*rhoD((i-1)*s.nRho+1:i*s.nRho,:);
    
    muVar((i-1)*s.nMu+1:i*s.nMu,:) = muM(:,:,i)*Th ...
        + muZ(:,:,i)*rhoVar((i-1)*s.nRho+1:i*s.nRho,:);
    muCon((i-1)*s.nMu+1:i*s.nMu) = mum(:,:,i) ...
        + muM(:,1:s.nX+s.nU,i)*[xpCon((i-1)*s.nX+1:i*s.nX);uCon((i-1)*s.nU+1:i*s.nU)] ...
        + muZ(:,:,i)*rhoCon((i-1)*s.nRho+1:i*s.nRho);
    muD((i-1)*s.nMu+1:i*s.nMu,:) = muZ(:,:,i)*rhoD((i-1)*s.nRho+1:i*s.nRho,:);
    
    zetaVar((i-1)*s.nZeta+1:i*s.nZeta,:) = zetaM(:,:,i)*Th ...
        + zetaZ(:,:,i)*rhoVar((i-1)*s.nRho+1:i*s.nRho,:);
    zetaCon((i-1)*s.nZeta+1:i*s.nZeta) = zetam(:,:,i) ...
        + zetaM(:,1:s.nX+s.nU,i)*[xpCon((i-1)*s.nX+1:i*s.nX);uCon((i-1)*s.nU+1:i*s.nU)] ...
        + zetaZ(:,:,i)*rhoCon((i-1)*s.nRho+1:i*s.nRho);
    zetaD((i-1)*s.nZeta+1:i*s.nZeta,:) = zetaZ(:,:,i)*rhoD((i-1)*s.nRho+1:i*s.nRho,:);
end