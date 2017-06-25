function [A,p] = RemoveRedundancy(A)
% REMOVEREDUNDANCY eliminates redundant inequalities in the system of 
% inequalities: A.A*x <= A.b
%    A must be a structure.
%    Output p is a list of indexes of irredundant inequalities.
%
%    See also ISREDUNDANT, IRREDUNDINEQS.
%
CheckSystem(A);
p = 1:size(A.A,1); % list of all inequalities
fprintf('\n');
str = 'Removing redundant inequalities ----- ';
fprintf(1,'%c',str);
for k=1:size(A.A,1)
    DisplayProgress(k,size(A.A,1));
    p = p(p~=k);
    exit_flag = isRedundant(A.A(k,:),A.b(k),A.A(p,:),A.b(p),[],[]);
    if exit_flag<0
        p = sort([p,k]);
    end
end
A.A = A.A(p,:);
A.b = A.b(p);
for k=1:length(str)+1;
    fprintf('\b');
end
end
