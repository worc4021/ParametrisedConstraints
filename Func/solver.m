function [wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,...
    zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ] = solver(actU,actW,s)

qc = zeros(s.ndelta,1);
Pc = blkdiag(s.P,zeros(s.nAlpha+s.nBeta));
Ph = zeros(s.nX+s.nU+s.nAlpha+s.nBeta);
qh = zeros(s.nX+s.nU+s.nAlpha+s.nBeta,1);
Ch = zeros(0,s.ndelta);
W = blkdiag(s.Q,s.R,zeros(s.nAlpha+s.nBeta));
phih = [];

initMats;

for n = s.N:-1:1    
    nlambda = numel(find(actW(n,:)));
    sol = maxWrapper(Pc,qc,Ch,phih,actW(n,:),s);
    
    wM(:,:,n) = sol.wM;
    wm(:,1,n) = sol.wm;
    
    xpM(:,:,n) = sol.xpM;
    xpm(:,1,n) = sol.xpm;
    
    defect = size(sol.lambdaZ,2);
    
    lambdaM(actW(n,:),:,n) = sol.lambdaM;
	lambdam(actW(n,:),1,n) = sol.lambdam;
    lambdaZ(actW(n,:),1:defect,n) = sol.lambdaZ;

    muM(:,:,n) = sol.muM;
    mum(:,1,n) = sol.mum;
    muZ(:,1:defect,n) = sol.muZ;
    
    nzeta = length(sol.zetam);
    zetaM(1:nzeta,:,n) = sol.zetaM;
    zetam(1:nzeta,1,n) = sol.zetam;
    zetaZ(1:nzeta,1:defect,n) = sol.zetaZ;
    
    T = [xpM(:,:,n);
         zeros(2,s.nX+s.nU),eye(2)];
	mh = [xpm(:,1,n);zeros(2,1)];
    
    Ph = T'*Pc*T-s.gammasq*wM(:,:,n)'*wM(:,:,n);
    qh = (qc'*T+mh'*Pc*T-s.gammasq*wm(:,1,n)'*wM(:,:,n))';
    
    Z = [sol.lambdaZ;
         sol.muZ;
         sol.zetaZ];
    
    if ~isempty(Z)
        C = -Z'*[zeros(nlambda,s.nX+s.nU+nAlpha),ones(nlambda,s.nBeta);
                       s.A,s.B,zeros(s.nAlpha+s.nBeta);
                       zeros(size(Ch,1),s.nX+s.nU),-Ch(s.nX+s.nU:end)];
        phi = Z'*[ones(nlambda,1);
                        zeros(s.nX,1);
                        phih];
    else
        C = zeros(0,s.nX+s.nU+s.nAlpha+s.nBeta);
        phi = [];
    end
               
	nkappa = numel(find(actU(n,:)));
    
    if nkappa>0
        0;
    end
    
    sol = minWrapper(Ph,qh,C,phi,actU(n,:),s);
    
    uK(:,:,n) = sol.uK;
    uk(:,:,n) = sol.uk;
    
    defect = max(size(sol.kappaZ,2),size(sol.rhoZ,2));
    
    kappaK(actU(n,:),:,n) = sol.kappaK;
    kappak(actU(n,:),:,n) = sol.kappak;
    kappaZ(actU(n,:),1:defect,n) = sol.kappaZ;
    
    nrho = size(sol.rhoK,1);
    rhoK(1:nrho,:,n) = sol.rhoK;
    rhok(1:nrho,:,n) = sol.rhok;
    rhoZ(1:nrho,1:defect,:) = sol.rhoZ;
    
    Z = [sol.kappaZ;
         sol.rhoZ];
	if ~isempty(Z)
        Ch = Z'*[zeros(nkappa,s.nX+s.nAlpha+s.nBeta),1;
                 -C(:,1:s.nX),-C(:,s.nX+s.nU+1),-C(:,s.nX+s.nU+s.nAlpha+1)];
        phih = Z'*[ones(nkappa,1);-phi];
	else
        Ch = zeros(0,s.nX+s.nAlpha+s.nBeta);
        phih = [];
	end
    
    V = [eye(s.nX),zeros(s.nX,s.nAlpha+s.nBeta);
         uK(:,:,n);
         zeros(s.nAlpha,s.nX),ones(s.nAlpha),zeros(s.nAlpha,s.nBeta);
         zeros(s.nBeta,s.nX),zeros(s.nBeta,s.nAlpha),ones(s.nBeta)];
    kh = [zeros(s.nX,1);uk(:,:,n);zeros(s.nAlpha+s.nBeta,1)];
    
    Pc = V'*(W+Ph)*V;
    qc = (qh'*V+kh'*(W+Ph)*V)';
end
    fprintf('\n')