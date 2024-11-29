function [LI, gg, cc, x] = find_special_points(fin, type)
    bd = coco_bd_read(fin);
    
    idx = find_idx(bd);
    [LI, gg, cc, x, ~] = find_curves(bd, idx);
    for i = 1:size(bd, 2)
        if strcmp(bd{1, i}, 'TYPE')
            idx_type = i;
        end
    end
    
    out = [];
    for i = 2:size(bd, 1)
        if strcmp(bd{i, idx_type}, type)
            out = [out, i];
        end
    end
    LI = LI(out);
    gg = gg(out);
    cc = cc(out);
    x  = x(out);
end