function [intervals_1, intervals_neg1] = find_intervals(x)
    % Initialize variables
    intervals_1 = [];
    intervals_neg1 = [];
    
    % Get the length of the input array
    n = length(x);
    
    % Start tracking intervals
    i = 1;
    while i < n
        if x(i) == 1
            % Start of interval containing 1
            start_1 = i;
            while i < n && x(i) == 1
                i = i + 1;
            end
            % End of interval containing 1
            intervals_1 = [intervals_1; start_1, i];
        elseif x(i) == -1
            % Start of interval containing -1
            start_neg1 = i;
            while i < n && x(i) == -1
                i = i + 1;
            end
            % End of interval containing -1
            intervals_neg1 = [intervals_neg1; start_neg1, i];
        else
            % Skip other values
            i = i + 1;
        end
    end
end