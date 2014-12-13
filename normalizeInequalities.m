function An = normalizeInequalities(A,b)

An = zeros(size(A));

for i = 1:length(b)
    An(i,:) = A(i,:)./b(i);
end

An = unique(An,'rows');