function [Aout,Bout] = linReduce(varargin)

if nargin<2
    A = varargin{1};
    b = ones(size(A,1),1);
    tol = 1e-9;
elseif nargin<3
    A = varargin{1};
    b = varargin{2};
    tol = 1e-9;
elseif nargin==3
    A = varargin{1};
    b = varargin{2};
    tol = varargin{3};
else
    error('Wrong number of arguments for linReduce');
end


% n = size(A,2);

% tempA = unique([A,b],'rows');
% A = tempA(:,1:n);
% b = tempA(:,n+1);

redIdx = [];
for i = 1:length(b)
    hIdx = 1:length(b);
    hIdx([redIdx,i]) = [];
    [x, fval, exitflag, output, lambda] = cplexlp(-A(i,:), A(hIdx,:), b(hIdx));
    if -fval - tol  < b(i)
        redIdx = [redIdx,i];
    end
end

outIdx = 1:length(b);
outIdx(redIdx) = [];

Aout = A(outIdx,:);
Bout = b(outIdx);
