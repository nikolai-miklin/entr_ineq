function c = sdiff(a,b)
% SDIFF does the same as setdiff for two lists a and b, but sdiff keeps
%    the order. 
%    Does not work for rows.
%
%    See also SETDIFF
%
[c,index_c] = setdiff(a,b);
[~,order_c] = sort(index_c);
c = c(order_c);
end