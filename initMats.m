% w = wM*[x;u;alpha;beta]+wm

xpM = zeros([s.nX,s.nX + s.nU + s.nAlpha + s.nBeta,s.N]);
xpm = zeros(s.nX,1,s.N);


wM = zeros([s.nW,s.nX+s.nU+ s.nAlpha + s.nBeta,s.N]);
wm = zeros(s.nW,1,s.N);

uK = zeros(s.nU,s.nX+ s.nAlpha + s.nBeta,s.N);
uk = zeros(s.nU,1,s.N);


muM = zeros([s.nX,s.nX+s.nU+ s.nAlpha + s.nBeta,s.N]);
mum = zeros(s.nX,1,s.N);
muZ = zeros([s.nX,s.nX + s.nU,s.N]);


lambdaM = zeros(s.nLambda,s.nX+s.nU+ s.nAlpha + s.nBeta,s.N);
lambdam = zeros(s.nLambda,1,s.N);
lambdaZ = zeros(s.nLambda,s.nX + s.nW,s.N);

zetaM = zeros(s.nX + s.nU,s.nX + s.nU+ s.nAlpha + s.nBeta,s.N);
zetam = zeros(s.nX + s.nU,1,s.N);
zetaZ = zeros(s.nX + s.nU,s.nX + s.nU,s.N);

kappaK = zeros(s.nKappa,s.nX + s.nAlpha + s.nBeta,s.N);
kappak = zeros(s.nKappa,1,s.N);
kappaZ = zeros(s.nKappa,s.nU + s.nX,s.N);

rhoK = zeros(s.nX + s.nU,s.nX + s.nAlpha + s.nBeta,s.N);
rhok = zeros(s.nX + s.nU,1,s.N);
rhoZ = zeros(s.nX + s.nU,s.nU+s.nX,s.N);