function [LI, gg, cc, x, stab] = find_curves_all_x(bd, idx, N)

    num_pts = size(bd, 1)- 1;

    LI = zeros(1, num_pts);
    gg = zeros(1, num_pts);
    cc  = zeros(1, num_pts);
    x  = zeros(2*N-1, num_pts);
    stab = zeros(1, num_pts);

    for i = 2:size(bd, 1)
        LI(i-1) = bd{i, idx(1)};
        gg(i-1) = bd{i, idx(2)};
        cc(i-1) = bd{i, idx(3)};
        x(:, i-1) = bd{i, idx(5)};
        stab(i-1) = bd{i, idx(4)};%-sign(max(real(bd{i, idx(4)})));
    end

    stab_aux = nan(size(stab));
    stab_aux(stab == 0) = 1;
    stab_aux(stab ~= 0) = -1;
    stab = stab_aux;

end
