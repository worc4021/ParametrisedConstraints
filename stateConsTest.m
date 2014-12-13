clear all
close all
clc

A = .5;
D = 1;

V = [1,0;
    -1,0;
    0,-1;
    .2,-1;
    -.1,-1;
    .1,-.5];
v = ones(size(V,1),1);

G = [1;-1];
g = [.2;.2];

Vn = zeros(size(V));
vn = zeros(size(v));

positivity = [zeros(1),-eye(1)];

Vpost = [];
vpost = [];

Vpre = [V;positivity];
vpre = [v;1];

iter = 1;
iterMax = 50;
while and(~isContained(Vpre,vpre,Vpost,vpost),iter<iterMax)

    if iter~=1
        Vpre = Vpost;
        vpre = vpost;
        Vn = zeros(size(Vpre));
        vn = zeros(size(vpre));
    end
    
    for i = 1:size(Vpre,1)
    [x,fval,exitflag,output,lambda] = cplexlp(-Vpre(i,1)*D,G,ones(size(G,1),1));
    Vn(i,:) = [Vpre(i,1)*(A),Vpre(i,2)+Vpre(i,2)*x];
    vn(i) = vpre(i)-Vpre(i,1)*x;
    end

[Vpost,vpost] = linReduce([Vpre;Vn;positivity],[vpre;vn;ones(1)]);
iter = iter+1;
end

% P = Polyhedron(V,v);
% P.minHRep;
% WN = P.A;
% wN = P.b;
% [VN,vN] = linReduce(V,v,1e-9);