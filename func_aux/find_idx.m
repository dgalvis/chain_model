function idx = find_idx(bd)

    idx_LI = nan;
    idx_g  = nan;
    idx_c  = nan;
    idx_stab = nan;
    idx_x    = nan;
    for i = 1:size(bd, 2)
        if strcmp(bd{1, i}, 'lambdaI')
            idx_LI = i;
        elseif strcmp(bd{1, i}, 'g')
            idx_g  = i;
        elseif strcmp(bd{1, i}, 'c')
            idx_c  = i;
        elseif strcmp(bd{1, i}, 'ep.test.USTAB')
            idx_stab = i;
        elseif strcmp(bd{1, i}, 'x')
            idx_x    = i;
        end
    end
    idx = [idx_LI, idx_g, idx_c, idx_stab, idx_x];
end