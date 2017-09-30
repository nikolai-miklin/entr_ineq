function A = ShannonCone(N,n)
% SHANNONCONE gives all basic Shannon inequalities for a set of variables n.
%    A = ShannonCone(N) takes n = [1,...,N].
%    A = ShannonCone(N,n) generates inequalities for n, but writes them in 
%    the space of N variables.	
%    Each inequality is a row in A.A
%    A is stored as a sparse matrix.
%    A does not contain redundant inequalities.
%
%    See also VONNEUMANNCONE.
%
switch nargin
    case 1
        n = 1:N;
end
if any(n>N)
    error('Some of the elements in [n] is greater than N.');
elseif isequal(1:N,sort(n))==1
    A.A = [Monotonicity(N); SubModularity(N)];
else
    size_n = numel(n);
    powerset_n = zeros(2^size_n-1,N);
    powerset_n(:,n) = i2v(size_n); % powerset of [n] in [1,..,N] 
    index_n = v2i(powerset_n);
    monot_ineq = Monotonicity(size_n);
    subm_ineq = SubModularity(size_n);
    A.A = zeros(size(monot_ineq,1)+size(subm_ineq,1),2^N-1);
    A.A(:,index_n) = [monot_ineq; subm_ineq];
end
A.A = sparse(A.A);
A.b = zeros(size(A.A,1),1);
A.var = find(sum(abs(A.A),1)~=0); % indexing of variables
A.A = A.A(:,A.var);
end

function A = Monotonicity(N)
% MONOTONICITY generates monotonicity inequalities for set of N variables. 
%    Monotonicity inequalities are: 
%    H(S(k,:)) - H([N]) <=0 for each S(k,:) = [N]\k, k-single element in [N]
%
    A = zeros(N,2^N-1);
    A(:,2^N-1) = -1;
    S = ones(N)-eye(N);
    for k=1:N
        A(k,v2i(S(k,:)))=1;    
    end
end

function A = SubModularity(N)
% SUBMODULARITY generates submodularity inequalities for set of N variables
%   Submodularity inequalities are:
%   H(S)+H(S+I+J)-H(S+I)-H(S+J) <=0
%   S denotes here powerset_not_c(k,:).
%
if N==2 
    A = [-1 -1 1];
else
    C = nchoosek(1:N,2);
    A = zeros(size(C,1)*(2^(N-2)),2^N-1);
    for n=1:size(C,1) % for each unique pair I,J in [1,...,N].
        I = zeros(1,N);
        I(1,C(n,1)) = 1;
        J = zeros(1,N);
        J(1,C(n,2)) = 1;
        not_c = sdiff(1:N,C(n,:)); 
        powerset_not_c = zeros(2^(N-2),N);
        powerset_not_c(2:end,not_c) = i2v(N-2);
        for k=1:2^(N-2)
            A((n-1)*(2^(N-2))+k,[v2i(powerset_not_c(k,:)+I),v2i(powerset_not_c(k,:)+J)])=-1;
            A((n-1)*(2^(N-2))+k,[v2i(powerset_not_c(k,:)+I+J),v2i(powerset_not_c(k,:))])=1;
        end
    end
end
end
