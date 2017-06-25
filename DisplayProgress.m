function DisplayProgress(k,N)
% DISPLAYPROGRESS is used to display progress in the precentage while k
%    runs from 1 to N.
%
if N>=2 % otherwise it does not make sense
    intervals_k = floor(N*(.01:.01:1));
    if isempty(find(intervals_k==k,1))==0 % if k is one of the intervals
        percent_to_diplay = floor(100*k/N);
        if percent_to_diplay==floor(100/N); 
            fprintf(1,' %i%%',percent_to_diplay);
        elseif percent_to_diplay<=10
            fprintf(1,'\b\b%i%%',percent_to_diplay);
        elseif percent_to_diplay==100
            fprintf(1,'\b\b\b\b\n');
        else
            fprintf(1,'\b\b\b%i%%',percent_to_diplay);
        end
    end
end
end