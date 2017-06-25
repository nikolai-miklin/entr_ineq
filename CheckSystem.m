function CheckSystem(A)
% CHECKSYSTEM checks whether system is composed correctly.
%
format_message = ['\nSystem of inequalities A*x<=b must be defined ',...
    'as a structure S, with S.A = A, S.b = b and S.var - a list of variables indexes.\n',...
    'For more details please see the package documentation.'];
if isstruct(A)
    names = fieldnames(A);
    if numel(names)==3
        if isempty(setdiff(names,{'A';'b';'var'}))
            if size(A.A,1)~=size(A.b,1)
                error_message = ['Sizes of the coefficient matrix .A and the vector of',...
                    ' constant .b are not consistent.'];
                error(sprintf(error_message));
            elseif size(A.A,2)~=numel(A.var)
                error_message = ['Number of variables in the coefficient matrix .A and',...
                    ' in the list of indexes .var must be the same.'];
                error(sprintf(error_message));
            end
        else
            error_message = ['Names of fields in the structure are not correct.',format_message];
            error(sprintf(error_message));
        end
    else
        error_message = ['Structure, which defines the system must contain',...
            ' three fields .A, .b, .var.',format_message];
        error(sprintf(error_message));
    end
else
    error_message = ['Input system must be of the type "structure". See MATLAB', ...
    ' documentation for structures.', format_message];
    error(sprintf(error_message));
end
end