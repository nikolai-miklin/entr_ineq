function A = vonNeumannCone(N,n,n_c)
% VONNEUMANNCONE gives all basic Shannon and von Neumann inequalities for a set of 
% variables n, where n_c are classical variables.
%    A = vonNeumannCone(N) takes n = [1,...,N] and n_c empty.
%    A = vonNeumannCone(N,n,n_c) generates inequalities for n, but writes them in
%    the space of N variables. 	
%    A contanins many redundant inequalities!
%
%    See also SHANNONCONE.
%
switch nargin
    case 1
        n = 1:N;
        n_c = [];
end
if any(n>N)
    error('One of the elements in [n] is greater than N.');
elseif isempty(setdiff(n_c,n))==0
    error('Some of the elements in [n_c] are not in [n].');
elseif isequal(1:N,sort(n))==1
    A.A = [WeakMonotonicity(N,n_c); SubModularity(N)];
else
    size_n = length(n);
    powerset_n = zeros(2^size_n-1,N);
    powerset_n(:,n) = i2v(size_n); % powerset of [1,..,ln] in N variables
    index_n = v2i(powerset_n);
    [~,index_n_c] = intersect(n,n_c);
    monot_ineq = WeakMonotonicity(size_n,index_n_c);
    subm_ineq = SubModularity(size_n);
    A.A = zeros(size([monot_ineq;subm_ineq],1),2^N-1);
    A.A(:,index_n) = [monot_ineq;subm_ineq];
end
A.A = unique(A.A,'rows');
Warning_message = ['Please note that the generated system may contain', ...
' redundant inequalities. Please consider removing them using RemoveRedundancy.'];
disp(Warning_message);
A.A = sparse(A.A);
A.b = zeros(size(A.A,1),1);
A.var = find(sum(abs(A.A),1)~=0); % indexing of variables
A.A = A.A(:,A.var);
end

function A = WeakMonotonicity(N,cl)
% WEAKMONOTONICITY generates monotonicity inequalities for set of N variables.
% cl is a set of classical variables
%    Weak monotonicity inequalities are:
%    H(X) +H(Y) -H(X+Z) -H(Y+Z) <=0
%    if Y is empty set and X is classical - use monotonicity:
%    H(X) -H(X+Z) <=0
% 
A = zeros(0,2^N-1);
for k=1:2^N-2
    Z = i2v(N,k); % take in terns Z - all subsets of 1:N
    not_Z = find(Z==0);
    size_not_Z = numel(not_Z);
    for l=1:2^size_not_Z-1
        x = i2v(size_not_Z,l);
        X = zeros(1,N);
        X(not_Z) = x; % take another subset X, not intersecting Z
        not_ZX = find(Z+X==0);
        size_not_ZX = numel(not_ZX);
        if size_not_ZX==0 % first the case when X+Z = 1:N
            if sum(X(cl))==sum(X) % if X is classical>use monotonicity
                A(end+1,v2i(X+Z)) = -1;
                A(end,v2i(X)) = 1;
            else % otherwise use weak monotinicity
                A(end+1,[v2i(X+Z),v2i(Z)]) = -1;
                A(end,v2i(X)) = 1;
            end
        else % now more general case of X+Z subset of 1:N
            for j=0:2^size_not_ZX-1
                y = i2v(size_not_ZX,j);
                Y = zeros(1,N);
                Y(not_ZX) = y; % yet another subset of 1:N
                if j==0 && sum(X(cl))==sum(X)
                    A(end+1,v2i(X+Z)) = -1;
                    A(end,v2i(X)) = 1;
                else
                    A(end+1,[v2i(X+Z),v2i(Y+Z)]) = -1;
                    A(end,[v2i(X),v2i(Y)]) = 1;
                end
            end
        end
    end
end
end

function A = SubModularity(N)
% SUBMODULARITY generates submodularity inequalities for set of N variables
%    Submodularity inequalities are:
%    H(Nsi) +H(Nsj) -H(Nsij) -H(Ns) <=0
%
%    Note that this function generates set of inequalities containing
%    redundancy.
%
A = zeros(0,2^N-1);
C = nchoosek(1:N,2);
for n=1:size(C,1)
    i = C(n,1);
    j = C(n,2);
    not_C = sdiff(1:N,C(n,:));
    for k=1:N-2 % first for the case of nonempty Ns (see description above)
        c = nchoosek(1:N-2,k);
        for l=1:size(c,1)
            s = zeros(1,N);
            s(not_C(c(l,:)))=1;
            Ns = v2i(s);
            si = s;
            si(i)=1;
            Nsi = v2i(si);
            sj = s;
            sj(j)=1;
            Nsj = v2i(sj);
            sij = s;
            sij([i,j])=1;
            Nsij = v2i(sij);
            A(end+1,:) = zeros(1,2^N-1);
            A(end,[Nsi,Nsj])=-1;
            A(end,[Nsij,Ns])=1;
        end 
    end
    % now for empty Ns (see description above)
    s = zeros(1,N);
    s(i) = 1;
    Ni = v2i(s);
    s = zeros(1,N);
    s(j) = 1;
    Nj = v2i(s);
    s = zeros(1,N);
    s([i,j]) = 1;
    Nij = v2i(s);
    A(end+1,:) = zeros(1,2^N-1);
    A(end,Nij)=1;
    A(end,[Ni,Nj])=-1;
end
end
