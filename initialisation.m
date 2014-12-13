clear all
close all
clc

p = path();
if isempty(strfind(p,'Func'));
    path(p,'Func');
end


A = [1,2;0,1];
B = [1,1;1,-1];
D = [1,2;2,1];

uStatMax = 2;
F = [1,1;
     2,0;
     1,-1];
F = [F;-F]/uStatMax;


w1max = .1;
w2max = .15;
G = [diag([1/w1max,1/w2max]);-diag([1/w1max,1/w2max])];

% x^+ = Ax + Bu + Dw
% Fu \leq 1 + \alpha
% Gw \leq 1 + \beta

Q = diag([2,1]);
R = diag([1,1]);

terminalCostCompute;

N = 6;
wIdx = logical(zeros(N,size(G,1)));
uIdx = logical(zeros(N,size(F,1)));

sizes.N = N;
sizes.nX = size(A,2);
sizes.nU = size(B,2);
sizes.nW = size(D,2);
sizes.nMu = sizes.nX;
sizes.nAlpha = 1;
sizes.nBeta = 1;
sizes.ndelta = sizes.nX + sizes.nAlpha + sizes.nBeta;
sizes.nZeta = sizes.nU + sizes.nX;
sizes.nRho = sizes.nZeta;
sizes.nLambda = size(G,1);
sizes.nKappa = size(F,1);
sizes.A = A;
sizes.B = B;
sizes.D = D;
sizes.F = F;
sizes.G = G;
sizes.Q = Q;
sizes.R = R;
sizes.P = P;
sizes.gammasq = gammaSQ;
sizes.tol = 1e-6;