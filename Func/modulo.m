function o = modulo(x,y)
% returns shifted index o between 1 and y
r = mod(x,y);
if r == 0 
    o = y;
else
    o = r;
end