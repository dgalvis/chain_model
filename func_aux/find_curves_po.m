function [LI, gg, MINx, MAXx, stab] = find_curves_po(bd, idx)

    num_pts = size(bd, 1)- 1;

    LI = zeros(1, num_pts);
    gg = zeros(1, num_pts);
    MAXx  = zeros(1, num_pts);
    MINx = zeros(1, num_pts);
    stab = zeros(1, num_pts);

    for i = 2:size(bd, 1)
        LI(i-1) = bd{i, idx(1)};
        gg(i-1) = bd{i, idx(2)};
        MINx(i-1) = bd{i, idx(4)}(1);
        MAXx(i-1) = bd{i, idx(5)}(1);
        stab(i-1) = bd{i, idx(3)};%-sign(max(real(bd{i, idx(3)})));
    end

    stab_aux = nan(size(stab));
    stab_aux(stab == 0) = 1;
    stab_aux(stab ~= 0) = -1;
    stab = stab_aux;

end