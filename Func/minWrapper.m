function sol = minWrapper(P,q,C,phi,actU,s)

K = s.R+P(s.nX+1:s.nX+s.nU,s.nX+1:s.nX+s.nU);
L = [s.F(actU,:);
     C(:,s.nX+1:s.nX+s.nU)];

nkappa = numel(find(actU));
 
LHSconst = [-q(s.nX+1:s.nX+s.nU);ones(nkappa,1);-phi];
LHSvar = [-P(s.nX+1:s.nX+s.nU,1:s.nX),-P(s.nX+1:s.nX+s.nU,s.nX+s.nU+1:end);
          zeros(nkappa,s.nX),ones(nkappa,s.nAlpha),zeros(nkappa,s.nBeta);
          -C(:,1:s.nX),-C(:,s.nX+s.nU+1:end)];

      [Q,R] = qr(L');
        if and(~isempty(L),size(L,1)+1<=size(Q,2))
            EW = eig(Q(:,size(L,1)+1:end)'*K*Q(:,size(L,1)+1:end));
        elseif size(L,1)+1>size(Q,2)
            EW = zeros(size(K));
        else
            EW = eig(K);
        end
        
        fprintf(' ')
        for i = 1:length(EW)
            fprintf('%6.1f ',EW(i))
        end
        fprintf(':')
        
[RHSconst,RHSvar,RHSorth] = KKTSolver(K,L,LHSconst,LHSvar);


sol.uK = RHSvar(1:s.nX,:);
sol.uk = RHSconst(1:s.nX);


if nkappa>0
    sol.kappaK = RHSvar(s.nX+1:s.nX+nkappa,:);
    sol.kappak = RHSconst(s.nX+1:s.nX+nkappa);
    if ~isempty(RHSorth)
        sol.kappaZ = RHSorth(s.nX+1:s.nX+nkappa,:);
    else
        sol.kappaZ = [];
    end
    if ~isempty(C)
        sol.rhoK = RHSvar(s.nX+nkappa+1:end,:);
        sol.rhok = RHSconst(s.nX+nkappa+1:end);
        if ~isempty(RHSorth)
            sol.rhoZ = RHSorth(s.nX+nkappa+1:end,:); 
        else
            sol.rhoZ = [];
        end
    else
        sol.rhoK = zeros(0,s.nX+s.nAlpha+s.nBeta);
        sol.rhok = [];
        sol.rhoZ = []; 
    end
else
    sol.kappaK = zeros(0,s.nX+s.nAlpha+s.nBeta);
    sol.kappak = [];
    sol.kappaZ = [];
    
    if ~isempty(C)
        sol.rhoK = RHSvar(s.nX+1+nkappa:end,1:s.nX);
        sol.rhok = RHSconst(s.nX+1+nkappa:end);
        if ~isempty(RHSorth)
            sol.rhoZ = RHSorth(s.nX+1+nkappa:end,:); 
        else
            sol.rhoZ = [];
        end
    else
        sol.rhoK = zeros(0,s.nX+s.nAlpha+s.nBeta);
        sol.rhok = [];
        sol.rhoZ = []; 
    end
    
end
        