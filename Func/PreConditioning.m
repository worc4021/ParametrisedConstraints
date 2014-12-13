function [Aout,bout] = PreConditioning(Ain,bin)

Aout = zeros(size(Ain));
bout = zeros(size(bin));

for i = 1:length(bin)
    Aout(i,:) = Ain(i,:)./norm(Ain(i,:));
    bout(i) = bin(i)/norm(Ain(i,:));
end