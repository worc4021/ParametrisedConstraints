function combos = combinator(vars,order)

combos = zeros(order^vars,vars);

for n = vars:-1:1
    combos(:,n) = repmat(kron((1:order)',ones(order^(vars-n),1)),order^(n-1),1);
end


IDX = [];
% parfor n = 1:(order^vars)
for n = 1:(order^vars)
    C = unique(combos(n,:));
    if length(C)<vars
        IDX = [IDX,n];
    else
        combos(n,:) = C;
    end
end
combos(IDX,:) = [];

[combos,dummy1,dummy2] = unique(combos,'rows');