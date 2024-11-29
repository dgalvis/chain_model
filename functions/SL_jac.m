%% Stuart-Landau (Cartesian) - Jacobian
% vectorized = true
% System:
% dA/dt = sigma A - L/2 * A * |A|^2
% sigma = lambda + 1j*(omega + lambda*c)
% L = 2*lambda * (1 + 1j*c)

% Coupling chain
% I_coup,j = g * (A_{j+1} - A_j) + g * (A_{j-1} - A_j)
function [J, DFDP] = SL_jac(~, y, p, N)
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


    J = zeros(size(y,1), size(y,1), size(y,2));

    gaux = zeros(N, size(lambdaE, 2));
    gaux([1,N], :) = -1*repmat(g, [2,1]);
    gaux(2:N-1, :) = -2*repmat(g, [N-2,1]);

    aux0 = real(sigma) - real(L).* (3*Ax.^2 + Ay.^2) + imag(L).* Ay .* 2 .* Ax + gaux;
    for i = 1:N
        J(i,i,:) = aux0(i,:);
    end

    aux1 = -imag(sigma) - real(L) .* (Ax .* 2 .* Ay) + imag(L) .* (3*Ay.^2 + Ax.^2);
    for i = 1:N
        J(i,i+N,:) = aux1(i,:);
    end

    aux2 = imag(sigma) - imag(L).* (3*Ax.^2 + Ay.^2) - real(L).* Ay .* 2 .* Ax;
    for i = 1:N
        J(i+N,i,:) = aux2(i,:);
    end

    aux3 = real(sigma) - imag(L) .* (Ax .* 2 .* Ay) - real(L) .* (3*Ay.^2 + Ax.^2) + gaux;
    for i = 1:N
        J(i+N,i+N,:) = aux3(i,:);
    end



    for i = 1:N
        if i ~= N
            J(i, i + 1, :) = g(1,:);
            J(i + N, i + N + 1, :) = g(1,:);
        end
        if i~=1
            J(i, i - 1, :) = g(1,:);
            J(i + N, i + N - 1, :) = g(1,:);
        end
    end

    if nargout == 2
        DFDP = SL_DFDP(0, y, p, N);
    end

end