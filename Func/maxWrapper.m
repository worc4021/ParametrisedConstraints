function sol = maxWrapper(P,q,Ch,phih,actW,s)

% [K  L'][x]        =   [LHSconst]  + [LHSvar][theta]
% [L  0 ][lambda]   =   [LHSconst]  + [LHSvar][theta]

nlambda = length(find(actW));

K = blkdiag(-s.gammasq*eye(s.nW),P(1:s.nX,1:s.nX));
L = [s.G(actW,:),zeros(nlambda,s.nX);
    -s.D,eye(s.nX);
     zeros(size(Ch,1),s.nW),Ch(:,1:s.nX)];
 
LHSconst = [zeros(s.nW,1);
            -q(1:s.nX);
            ones(nlambda,1);
            zeros(s.nX,1);
            phih];
            
LHSvar = [zeros(s.nW,s.nX+s.nU+s.nAlpha+s.nBeta);
          zeros(s.nX,s.nX+s.nU),-P(1:s.nX,s.nX+1:end);
          zeros(nlambda,s.nX+s.nU+s.nAlpha),ones(nlambda,s.nBeta);
          s.A,s.B,zeros(s.nX,s.nAlpha+s.nBeta);
          zeros(size(Ch,1),s.nX+s.nU),-Ch(:,s.nX+1:end)];
          
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

sol.wM = RHSvar(1:s.nW,:);
sol.wm = RHSconst(1:s.nW);

sol.xpM = RHSvar(s.nW+1:s.nW+s.nX,:);
sol.xpm = RHSconst(s.nW+1:s.nW+s.nX);

if nlambda>0
    sol.lambdaM = RHSvar(s.nW+s.nX+1:s.nW+s.nX+nlambda,:);
    sol.lambdam = RHSconst(s.nW+s.nX+1:s.nW+s.nX+nlambda);
    
    
    sol.muM = RHSvar(s.nW+s.nX+nlambda+1:s.nW+s.nX+nlambda+s.nMu,:);
    sol.mum = RHSconst(s.nW+s.nX+1+nlambda:s.nW+s.nX+nlambda+s.nMu);
    
    if ~isempty(RHSorth)
        sol.lambdaZ = RHSorth(s.nW+s.nX+1:s.nW+s.nX+nlambda,:);
        sol.muZ = RHSorth(s.nW+s.nX+1+nlambda:s.nW+s.nX+nlambda+s.nMu,:);
    else
        sol.lambdaZ = [];
        sol.muZ = [];
    end
    
    if ~isempty(Ch)
    sol.zetaM = RHSvar(s.nW+s.nX+nlambda+s.nX+nlambda+1:end,:);
    sol.zetam = RHSconst(s.nW+s.nX+nlambda+s.nX+nlambda+1:end);
        if ~isempty(RHSorth)
            sol.zetaZ = RHSorth(s.nW+s.nX+nlambda+s.nX+nlambda+1:end,:);
        else
            sol.zetaZ = [];
        end
    else
        sol.zetaM = [];
        sol.zetam = [];
        sol.zetaZ = [];
    end
    
else
    sol.lambdaM = [];
    sol.lambdam = [];
    sol.lambdaZ = [];
    
    
    sol.muM = RHSvar(s.nW+s.nX+1:s.nW+s.nX+s.nX,:);
    sol.mum = RHSconst(s.nW+s.nX+1:s.nW+s.nX+s.nX);
    if ~isempty(RHSorth);
        sol.muZ = RHSorth(s.nW+s.nX+1:s.nW+s.nX+s.nX,:);
    else
        sol.muZ = [];
    end
    
    if ~isempty(Ch)
    sol.zetaM = RHSvar(s.nW+s.nX+s.nX+1:end,:);
    sol.zetam = RHSconst(s.nW+s.nX+s.nX:end);
        if ~isempty(RHSorth)
            sol.zetaZ = RHSorth(s.nW+s.nX+s.nX:end,:);
        else
            sol.zetaZ = [];
        end
    else
        sol.zetaM = [];
        sol.zetam = [];
        sol.zetaZ = [];
    end
end