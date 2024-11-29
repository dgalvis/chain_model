function plot_curve(fin, col, uspec, flag)

    if nargin < 4
        flag = 'LI';
    end

    bd = coco_bd_read(fin);
    idx = find_idx(bd);
    [LI, gg, cc, x, stab] = find_curves(bd, idx);
    [int_stab, int_ustab] = find_intervals(stab);

    
    if strcmp(flag, 'LI')
        pars = LI;
    elseif strcmp(flag, 'gg')
        pars = gg;
    elseif strcmp(flag, 'cc')
        pars = cc;
    end

    for i = 1:size(int_stab, 1)
        plot(pars(int_stab(i,1):int_stab(i,2)),  x(1, int_stab(i,1):int_stab(i,2)), '-', 'color', col);
    end
    for i = 1:size(int_ustab, 1)
        plot(pars(int_ustab(i,1):int_ustab(i,2)),  x(1, int_ustab(i,1):int_ustab(i,2)), uspec, 'color', col);
    end
end