function [A,L] = ApplyConstraints(A,L,keep_list)
% APPLYCONSTRAINTS finds an intersection of cone given by A.A*x<=A.b with
% linear constraints given by L.A*x=L.b. Effectively this program 
% substitute variables in A using constraints L.
%
%    A -- structure of the initial system
%    L -- structures of linear constraints
%    keep_list -- list of variables which will not be substituted (often
%    marginal scenario).
%    keep_list can be as well given as a binary matrix of the marginal scenario.
%    priority_list -- list of variables to be subsituted first, e.g.
%    variables corresponding to entropies of more variables.
%    
%    See also D_SEPARATION.
%
switch nargin
    case 1
        error('Not enough input arguments.');
    case 2
        keep_list = [];
    case 3
        % do nothing
    otherwise
        error('Too many input arguments.');          
end
CheckSystem(A);
CheckSystem(L);
% Additionally we need to check that system of linear equation L.A*x=L.b
% is consitent.
if rank([L.A,L.b])>rank(L.A)
    error('System of linear constraints L is inconsistent.');
end
% check whether keep_list is given as a list of variables or as a matrix
if isempty(keep_list)==0
    if isempty(find(keep_list==0,1))==0 % if it contrains zeros
        if isempty(sdiff(unique(keep_list),[0,1])) % if it only contains 0 and 1
            keep_list = v2i(keep_list); % turn into indices
        else
            error('The list of variables to keep is given in the wrong form.');
        end
    end
end

% A and L can contain different lists of variables, hence we need to 
% change them such that A.var = L.var

var_common = sort(unique([A.var,L.var])); % common list of variables of A and L
A_new = zeros(size(A.A,1),numel(var_common));
[~,sub_index] = intersect(var_common,A.var);
A_new(:,sub_index) = A.A;
A.A = A_new;
L_new = zeros(size(L.A,1),numel(var_common));
[~,sub_index] = intersect(var_common,L.var);
L_new(:,sub_index) = L.A;
L.A = L_new;
A.var = var_common;
L.var = var_common;

% now proceed with new systems A and L

subst_list = sdiff(1:size(A.A,2),keep_list); % list of variables, which can be substituted.
subst_list = intersect(subst_list,find(sum(abs(L.A),1)~=0)); % remove variables not included in L.

% Next, use Gaussian elimination to find the best way to substitute variables.
% The point here is that some constraints might be linearly dependent in the
% subspace corresponding to the variables, which can be substituted. 

disp('Calculating the best way to apply the constraints ---  ');
[L_to_apply,eqs_list] = GaussElimination(L.A(:,subst_list));

kept_var = 1:numel(var_common); % keep track of variables that have been substituted

disp(' ');
disp('Applying constraints ---  ');
for k=1:size(L_to_apply,1)
    DisplayProgress(k,size(L_to_apply,1));
    var_to_subst = subst_list(find(L_to_apply(k,:)~=0,1)); % variable to substitute on step k
    kept_var = kept_var(kept_var~=var_to_subst);
    % substitution in A
    coeff = L.A(eqs_list(k),var_to_subst);
    vect = A.A(:,var_to_subst); % have to store as a separate variable
    A.A = A.A.*abs(coeff)-sign(coeff)*(vect*L.A(eqs_list(k),:));
    A.b = A.b.*abs(coeff)-sign(coeff)*(vect*L.b(eqs_list(k)));
    % substitution in L
    vect = L.A(:,var_to_subst); 
    L.A = L.A.*abs(coeff)-sign(coeff)*(vect*L.A(eqs_list(k),:));
    L.b = L.b.*abs(coeff)-sign(coeff)*(vect*L.b(eqs_list(k)));
end
% remove trivial constraints and inequalitites
eq_list = find(sum(abs(L.A),2)~=0);
L.A = L.A(eq_list,:);
L.b = L.b(eq_list);
eq_list = find(sum(abs(A.A),2)~=0);
A.A = A.A(eq_list,:);
A.b = A.b(eq_list);
% check if some constraints were left
if isempty(L.A)==0
    disp('Not all linear constraints were used!');
end
% now remove substrituted variables from both A and L
A.A = A.A(:,kept_var);
A.var = A.var(kept_var);
L.A = L.A(:,kept_var);
L.var = L.var(kept_var);
% last, but not least divide each inequality on gcd where possible
A = DivideOnGCD(A);
if isempty(L.A)==0
    L = DivideOnGCD(L);
end
end


function [A,eqs_list] = GaussElimination(A)
% GAUSSELIMINATION performs Gaussian elimination in the system A*x=0
% [PARTIALLY!] 
%
% For our puproses we only need the diagonal matrix which tells us
% which equation to use to substitute which variable.
%
eqs_list = zeros(size(A,1),1);
variable_list = 1:size(A,1);
disp('  ');
for l=1:size(A,2)
    DisplayProgress(l,size(A,2));
    nl = find(A(:,l)~=0);
    nl = sdiff(nl,eqs_list);
    if isempty(nl)==0
        eqs_list(l) = nl(1);
        variable_list = variable_list(variable_list~=eqs_list(l));
        A(variable_list,:) = A(variable_list,:)-A(variable_list,l)*A(eqs_list(l),:)./A(eqs_list(l),l);
    end
end
eqs_list = eqs_list(eqs_list~=0);
A = A(eqs_list,:);
end

