function idx = find_idx_po(bd)
    idx_LI = nan;
    idx_g  = nan;
    idx_stab = nan;
    idx_min_x    = nan;
    idx_max_x = nan;
    for i = 1:size(bd, 2)
        if strcmp(bd{1, i}, 'lambdaI')
            idx_LI = i;
        elseif strcmp(bd{1, i}, 'g')
            idx_g  = i;
        elseif strcmp(bd{1, i}, 'po.test.USTAB')
            idx_stab = i;
        elseif strcmp(bd{1, i}, 'MIN(x)')
            idx_min_x    = i;
        elseif strcmp(bd{1, i}, 'MAX(x)')
            idx_max_x    = i;
        end
    end
    idx = [idx_LI, idx_g, idx_stab, idx_min_x, idx_max_x];
end