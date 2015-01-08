
% Computation with parametrised constraints. (Doesn't really work!)

V = [F*K,-ones(size(F,1),sizes.nAlpha),zeros(size(F,1),sizes.nBeta);
    G*KW,zeros(size(G,1),sizes.nAlpha),-ones(size(G,1),sizes.nBeta);
    zeros(sizes.nAlpha+sizes.nBeta,sizes.nX),-eye(sizes.nAlpha+sizes.nBeta)];
v = ones(size(F,1)+size(G,1)+sizes.nAlpha+sizes.nBeta,1);

[V,v] = PreConditioning(V,v);

Vn =zeros(size(V));
vn = zeros(size(v));

Vpost = [];
vpost = [];

Vpre = V;
vpre = v;

iter = 1;
iterMax = 50;
while and(~altIsContained(Vpre,vpre,Vpost,vpost),iter<iterMax)
    
    if iter~=1
       [Vpre,vpre] = PreConditioning(Vpost,vpost);
        Vn = zeros(size(Vpre));
        vn = zeros(size(vpre));
    end
    
    for i = 1:size(Vpre,1)
    [x,fval,exitflag,output,lambda] = cplexlp(-Vpre(i,1:sizes.nX)*D,...
        G,ones(size(G,1),1));
    Vn(i,:) = [Vpre(i,1:sizes.nX)*(A+B*K),Vpre(i,sizes.nX+1),...
        Vpre(i,sizes.nX+sizes.nAlpha+1)+Vpre(i,1:sizes.nX)*x];
    vn(i) = vpre(i)-Vpre(i,1:sizes.nX)*D*x;
    end

[Vpost,vpost] = rowReduce([Vpre;Vn],[vpre;vn]);
% Vpost = [Vpre;Vn];
% vpost = [vpre;vn];
fprintf('At iteration %2.0d the number of constraints was %3.0d.\n',iter,length(vpost))

iter = iter+1;
end

eqVec = [1;-1];

% [Vpost,vpost] = PreConditioning(Vpost,vpost);
[Vpost,vpost] = rowReduce(Vpost,vpost);
[Vout,vout] = polyProject([Vpost;[zeros(2),eqVec,zeros(2,1);zeros(2,3),eqVec]],...
    [vpost;zeros(4,1)],2);
[Vout,vout] = rowReduce(Vout,vout);

tInvSet = struct('A',Vpost,'b',vpost);

% [Vout,vout] = linReduce(Vout,vout);
[Vout,vout] = PreConditioning(Vout,vout);
pInvSet = struct('A',Vout,'b',vout);


fprintf('Is the terminal set in the correct halfspace (alpha,beta>-1): %d\n',...
    isContained(Vpost,vpost,[zeros(2),-eye(2)],ones(2)));



% Constant Terminal set computation. Terminal set computation is working.

V = [F*K;G*KW];
v = ones(size(V,1),1);

[V,v] = PreConditioning(V,v);

Vn =zeros(size(V));
vn = zeros(size(v));

Vpost = [];
vpost = [];

Vpre = V;
vpre = v;

iter = 1;
iterMax = 50;
while and(~altIsContained(Vpre,vpre,Vpost,vpost),iter<iterMax)
    
        if iter~=1
           [Vpre,vpre] = PreConditioning(Vpost,vpost);
            Vn = zeros(size(Vpre));
            vn = zeros(size(vpre));
        end

        for i = 1:size(Vpre,1)
        [x,fval,exitflag,output,lambda] = cplexlp(-Vpre(i,:)*D,G,ones(size(G,1),1));
        Vn(i,:) = Vpre(i,:)*(A+B*K);
        vn(i) = vpre(i)+fval;
        end
    
        if ~isempty(find(vn<0,1,'first'))
            0;
        end
        
%     [Vpost,vpost] = linReduce([Vpre;Vn],[vpre;vn]);
Vpost = [Vpre;Vn];
vpost = [vpre;vn];
    fprintf('At iteration %2.0d the number of constraints was %3.0d.\n',...
        iter,length(vpost))

    iter = iter+1;
end

[Vpost,vpost] = PreConditioning(Vpost,vpost);
[Vpost,vpost] = linReduce(Vpost,vpost);
cInvSet = struct('A',Vpost,'b',vpost);

figure(1)
plot(Polyhedron(cInvSet.A,cInvSet.b),'alpha',.2,...
    Polyhedron(pInvSet.A,pInvSet.b),'alpha',.2,...
    Polyhedron(tInvSet.A(:,1:2),tInvSet.b),'alpha',.2,'color','b')

figure(gcf)
nHor = 5;

Lam = cell(1,nHor);
lam = cell(1,nHor);

% For constant terminal constraint:

% Lam{1} = [Vpost,zeros(length(vpost),sizes.nAlpha+sizes.nBeta);
%           zeros(2*(sizes.nAlpha+sizes.nBeta),sizes.nX),...
% [eye(sizes.nAlpha+sizes.nBeta);
%           -eye(sizes.nAlpha+sizes.nBeta)]];
% lam{1} = [vpost;20*ones(sizes.nAlpha+sizes.nBeta,1);...
% ones(sizes.nAlpha+sizes.nBeta,1)];

% For parametrised terminal constraint:

% Lam{1} = Vpost;
% lam{1} = vpost;
% 
% 
% for  i = 2:nHor
%     [tLam,tlam] = stepBack(Lam{i-1},lam{i-1},sizes);
%     [tLam,tlam] = linReduce(tLam,tlam);
%     [Lam{i},lam{i}] = PreConditioning(tLam,tlam);
% end