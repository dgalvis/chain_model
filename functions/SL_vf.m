%% Stuart-Landau (Cartesian) - ODE v1 (complex version)
% vectorized = true
% System:
% dA/dt = sigma A - L/2 * A * |A|^2
% sigma = lambda + 1j*(omega + lambda*c)
% L = 2*lambda * (1 + 1j*c)

% Coupling chain
% I_coup,j = g * (A_{j+1} - A_j) + g * (A_{j-1} - A_j)

function dydt = SL_vf(~, y, p, N)
    lambdaE = p(1,:);
    lambdaI = p(2,:);
    g       = p(3,:);
    c       = p(4,:);
    omega   = p(5,:);

    sigma = zeros(N,size(c,2));
    L = zeros(N,size(c,2));
    sigma([1,N],:) = repmat(lambdaE + (omega + lambdaE.*c)*1j, [2,1]);
    sigma(2:N-1,:) = repmat(lambdaI + (omega + lambdaI.*c)*1j, [N-2,1]);
    L([1, N],:)    = repmat(lambdaE .*  (1 + c*1j), [2,1]);
    L(2:N-1,:)     = repmat(lambdaI .*  (1 + c*1j), [N-2,1]);

    Ax = y(1:N, :);
    Ay = y(N+1:end, :);

    A = Ax + 1j*Ay;

    dydt = sigma .* A - L .* A .* abs(A).^2;

    dydt(2:end,:) = dydt(2:end,:) + g .* (A(1:end-1,:) - A(2:end,:));
    dydt(1:end-1,:) = dydt(1:end-1,:) + g .* (A(2:end,:) - A(1:end-1,:));


    dydt = [real(dydt); imag(dydt)];
end