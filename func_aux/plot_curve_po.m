function plot_curve_po(fin, col, uspec, linewidth)

    if nargin < 4
        linewidth = 2;
    end

    bd = coco_bd_read(fin);
    idx = find_idx_po(bd);
    [LI, ~, MINx, MAXx, stab] = find_curves_po(bd, idx);
    [int_stab, int_ustab] = find_intervals(stab);
    for i = 1:size(int_stab, 1)
        plot(LI(int_stab(i,1):int_stab(i,2)),  MAXx(1, int_stab(i,1):int_stab(i,2)), '-', 'color', col, 'Linewidth', linewidth);
        plot(LI(int_stab(i,1):int_stab(i,2)),  MINx(1, int_stab(i,1):int_stab(i,2)), '-', 'color', col, 'Linewidth', linewidth);
    end
    for i = 1:size(int_ustab, 1)
        plot(LI(int_ustab(i,1):int_ustab(i,2)),  MAXx(1, int_ustab(i,1):int_ustab(i,2)), uspec, 'color', col, 'Linewidth', linewidth);
        plot(LI(int_ustab(i,1):int_ustab(i,2)),  MINx(1, int_ustab(i,1):int_ustab(i,2)), uspec, 'color', col, 'Linewidth', linewidth);
    end

end