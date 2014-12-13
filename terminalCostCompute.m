X = sdpvar(size(A,1));
Y = sdpvar(size(B,2),size(B,1));
PT = sdpvar(size(A,1));
gammasq = sdpvar(1);

BLK11 = blkdiag(X,gammasq*eye(size(PT)));
BLK12 = [X*A'+Y'*B',(sqrtm(Q)*X)',(sqrtm(R)*Y)';
         D',zeros(size(X)),zeros(size(Y'))];
BLK22 = blkdiag(X,eye(size(Q)),eye(size(R)));

SDP = [X > 0, [BLK11,BLK12;BLK12',BLK22] > 0, [PT,eye(size(PT));eye(size(PT)),X]>0];

opt = sdpsettings('solver','mosek-sdp','verbose',0);

info = solvesdp(SDP,trace(PT),opt);


P = double(PT);
K = double(Y)*P;
gammaSQ = double(gammasq);
KW = (gammaSQ*eye(size(P))-D'*P*D)\D'*P*(A+B*K);