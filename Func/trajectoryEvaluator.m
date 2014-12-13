function [x,u,w,kappa,lambda] = trajectoryEvaluator(theta,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s)

x = zeros(s.nX,s.N+1);
u = zeros(s.nU,s.N);
w = zeros(s.nW,s.N);
kappa = zeros(s.nKappa,s.N);
lambda = zeros(s.nLambda,s.N);
x(:,1) = theta(1:s.nX);
alpha = theta(s.nX+1:s.nX+s.nAlpha);
beta = theta(s.nX+s.nAlpha+1:end);

for i = 1:s.N
    u(:,i) = uK(:,:,i)*[x(:,i);alpha;beta]+uk(:,:,i);
    kappa(:,i) = kappaK(:,:,i)*[x(:,i);alpha;beta]+kappak(:,:,i);
    w(:,i) = wM(:,:,i)*[x(:,i);u(:,i);alpha;beta]+wm(:,:,i);
    lambda(:,i) = lambdaM(:,:,i)*[x(:,i);u(:,i);alpha;beta]+lambdam(:,:,i);
    x(:,i+1) = xpM(:,:,i)*[x(:,i);u(:,i);alpha;beta]+xpm(:,:,i);
end