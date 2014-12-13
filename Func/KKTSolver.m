function [RHSconst,RHSvar,RHSorth] = KKTSolver(K,L,LHSconst,LHSvar)

[m,n] = size(L);

if m == 0
    RHSconst = K\LHSconst;
    RHSvar = K\LHSvar;
    RHSorth = [];
else
    
    [Q,R] = qr(L);
    
    if and(m == 1,~prod(prod(L == zeros(size(L)))))
        r = 1;
    elseif prod(prod(L == zeros(size(L))))
        r = 0;
    else
        r = numel(find(abs(diag(R)) > 1e-9));
    end
    
    if r == m
        hel = [K,L';
               L,zeros(m)];
        RHSconst = hel\LHSconst;
        RHSvar = hel\LHSvar;
        RHSorth = [];
    elseif r == n
        R = R(1:r,:);
        H = R\(Q(:,1:r)');
        hel = [zeros(n),H;
               H',-H'*K*H];
        RHSconst = hel*LHSconst;
        RHSvar = hel*LHSvar;
        Q2 = Q(:,(r+1):end);
        RHSorth = [zeros(size(K,1),size(Q2,2));Q2/diag(sum(Q2))];
    else
        [cT,p] = chol(K);
        if ~p
          icTR = cT'\R(1:r,:)';
          RiTR = icTR'*icTR; % = R*inv(Theta)*R'
          iTR = cT\icTR; % = inv(Theta)*R'
          c2 = chol(RiTR);
          ic2RiT = c2'\iTR'; % = inv(c2')*R*inv(Theta)
          ic2Q1 = c2'\Q(:,1:r)'; % = inv(c2')*Q1'
          B12 = ic2RiT'*ic2Q1; % = inv(Theta)*R'*inv(R*inv(Theta)*R')*Q1'
          RHSvar = [cT\(cT'\eye(n))-ic2RiT'*ic2RiT, B12; B12', -ic2Q1'*ic2Q1]*LHSvar;
          RHSconst = [cT\(cT'\eye(n))-ic2RiT'*ic2RiT, B12; B12', -ic2Q1'*ic2Q1]*LHSconst;
        else
          warning('conx0:ThetaPD','Theta is not positive definite, min(eig(Theta)) = %d\n', min(eig(K)));
          [L,D,P] = ldl(K);
          iT = P*(L'\(D\(L\(P'*eye(n))))); % = inv(Theta)
          iTR = iT*R'; % = inv(Theta)*R'
          [L2,D2,P2] = ldl(R*iTR);
          iRiTR = P2*(L2'\(D2\(L2\(P2'*eye(r))))); % = inv(R*inv(Theta)*R')
          iRiTR_RiT = iRiTR*iTR'; % = inv(R*inv(Theta)*R')*R*inv(Theta)
          B11 = iT*(eye(n) - R'*iRiTR_RiT);
          B21 = Q(:,1:r)*iRiTR_RiT;
          B22 = -Q(:,1:r)*iRiTR*Q(:,1:r)';
          RHSvar = [B11, B21'; B21, B22]*LHSvar;
          RHSconst = [B11, B21'; B21, B22]*LHSconst;
        end
    Q2 = Q(:,(r+1):m);
    RHSorth = [zeros(size(K,1),size(Q,2));Q2/diag(sum(Q2))];
    end
end    

    
% 
%     
%     Q1 = Q(:,1:r);
%     Q2 = Q(:,r+1:end);
% 
%     if r == n
%         R = R(1:r,:);
%         H = R\(Q1');
%         hel = [zeros(n),H;
%                H',-H'*K*H];
%         RHSconst = hel*LHSconst;
%         RHSvar = hel*LHSvar;
%         RHSorth = [zeros(n,m-r);
%                    Q2];
%     elseif r == m
%         hel = [K,L';
%                L,zeros(m)];
%         RHSconst = hel\LHSconst;
%         RHSvar = LHSvar;
%         RHSorth = zeros(n,0);
%     else
%         R = R(1:r,:);
%         A = (R*(K\R))\eye(r);
%         C = K\eye(n);
%         MAT = zeros(n+m);
%         MAT(1:n,1:n) = C-C*R'*A*R*C;
%         MAT(1:n,n+1:end) = C*R'*A*Q1';
%         MAT(n+1:end,1:n) = (C*R'*A*Q1')';
%         MAT(n+1:end,n+1:end) = -Q1*A*Q1';
%         RHSconst = MAT*LHSconst;
%         RHSvar = MAT*LHSvar;
%         RHSorth = [zeros(n,n-r);Q2];
%     end
% end
% 
