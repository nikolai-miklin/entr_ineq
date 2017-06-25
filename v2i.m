function I = v2i(V)
% V2I translates binary numbers from V to a column of decimals.
% If V = [0 0 ... 0], I is an empty set.
%
%   See also I2V.

S = size(V);
I = zeros(S(1),1);
for k=1:S(1)
    for l=1:S(2)
        I(k) = I(k) + 2^(S(2)-l)*V(k,l);
    end
end
if I==0
    I = zeros(0);
end
end