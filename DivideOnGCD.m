function A = DivideOnGCD(A)
% DIVIDEONGCD takes a system of inequalities A, finds if there is a 
% common integer divisor for each inequality, divides on it, and, 
% finally, removes identical inequalities from the system.
%
%    Example:
%    for the system:  2*x + 4*y <= 6
%                       x + 2*y <= 3
%    program returns only the latter as an output.
%
CheckSystem(A);
epsilon = 10^-9; % precision of this program
max_divisor = 50; % maximal disisor that is handeled   
AExt = [A.A,A.b]; % it's easier to work with the extended system
AExt = round(AExt/epsilon)*epsilon;
for l=1:size(AExt,1)
    [~,denoms_in_l] = rat(AExt(l,:)); % find denominators in the l'th row
    if max(VectorLCM(denoms_in_l))>max_divisor
        % denominator is too large -> skip this inequality
    else
        AExt(l,:) = AExt(l,:)*VectorLCM(denoms_in_l); % bring inequality
        % to the form with only integer coefficients 
        AExt(l,:) = AExt(l,:)./VectorGCD(AExt(l,AExt(l,:)~=0)); % divide on
        % gcd of all nonzero elements of the row
    end
end
% round it again to the precision and remove identical inequalities 
AExt = round(AExt/epsilon)*epsilon;
AExt = unique(AExt,'rows');
A.A = AExt(:,1:end-1);
A.b = AExt(:,end);
end

function g = VectorGCD(L)
% VECTORGCD finds greatest common divisor of numbers in the vector L
% 
g = min(L);
for k=1:numel(L)
    g = gcd(g,L(k));
    if g==1
        break
    end
end
end

function m = VectorLCM(L)
% VECTORLCM finds least common multiple of numbers in the vector L
%
%   L must not contain zeros.
%
m = max(L);
product_of_all = prod(L);
for k=1:numel(L)
    m = lcm(m,L(k));
    if m==product_of_all
        break
    end
end
end
        