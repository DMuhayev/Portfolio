function [valX, valY, valS, valSigma, Theta] = RegFind(err, k, eta, A, B_2, C_2, C_d, B_d, D_dd, D_2d, D_d2, r, f, k_r, m_f)
%Данная функция решает задачу поиска регулятора в системах с помехами 
%  Подробное описание алгоритма описано в книге
%   Д. В. Баландин, М. М. Коган "Синтез законов управления на основе
%   линейных матричных неравенств" (стр 116)
%   названия переменных сохранены такими же как в оригинале

    zeta = 1/eta;
    
    n_x = numel(A(1, :));
    n_u = numel(B_2(1, :));
    n_y = numel(C_2(:, 1));
    n_d = numel(B_d(1, :));
    
    %% some measurements for future!

    A_0 = [A, zeros(n_x, k); zeros(k, n_x), zeros(k, k)];
    B = [zeros(n_x, k), B_2; eye(k), zeros(k, n_u)];
    C = [zeros(k, n_x), eye(k); C_2, zeros(n_y, k)];
    B_0 = [B_d; zeros(k, n_d)];
    C_0 = [C_d, zeros(n_d, k)];
    D_1 = [zeros(k, n_d); D_2d];
    D_2 = [zeros(n_d, k), D_d2];

    P = [C, D_1, zeros(n_y + k, n_d)];
    R = [B', zeros(n_u + k, n_d), D_2'];

    W_p = null(P);
    W_r = null(R);

    %% initialization of X, Y and S, Sigma matrices 

    setlmis([])

    [X, ~, sX] = lmivar(1, [n_x + k, 1]); 
    [Y, ~, sY] = lmivar(1, [n_x + k, 1]);

    struct = zeros(r+f, 2);

    for i=1:f+r
        if (i <= r) 
            struct(i, :) = [k_r(i), 1];
        else
            struct(i, :) = [m_f(i), 0];
        end
    end

    [S, ~, sS] = lmivar(1, struct);
    [Sigma, ~, sSigma] = lmivar(1, struct);
    %added matrixes needed to evaluate X = Y^-1 and S = Sigma^-1
    %X_roof = diag(X, S), Y_roof = diag(Y, Sigma)
    [X_roof, ~, sX_roof] = lmivar(3, [sX, zeros(n_x + k, numel(sS(1, :)) ); zeros(numel(sS(1, :)), n_x + k), sS]);
    [Y_roof, ~, sY_roof] = lmivar(3, [sY, zeros(n_x + k, numel(sS(1, :)) ); zeros(numel(sS(1, :)), n_x + k), sSigma]);
    n_Xroof = numel(sX_roof(1, :) );
    [Lambda, ~, sLambda] = lmivar(2, [1, 1]);
    %% initialization of LMI sys (11.17)

    % first LMI of sys
    %W_p'*[matrix]*W_p < 0      [matrix] - given in lmiterm lower
    lmiterm([1 0 0 0], W_p);  %W_p'*[matrix]*W_p
    lmiterm([1 1 1 X], A_0', 1, 's');
    %lmiterm([1 1 2 X], 1, B_0);
    %lmiterm([1 1 3 0], C_0'); %first row  A_0'*X+X*A_0, X*B_0, C_0';
    lmiterm([1 2 1 X], B_0', 1);
    lmiterm([1 2 2 S], -zeta, 1);
    %lmiterm([1 2 3 0], D_dd'); %second row  B_0'*X, -zeta*S, D_dd';
    lmiterm([1 3 1 0], C_0);
    lmiterm([1 3 2 0], D_dd);
    lmiterm([1 3 3 Sigma], -zeta, 1); %third row C_0, D_dd, -zeta*Sigma;

    % second LMI of sys (likewise)
    %W_r'*[matrix]*W_r < 0      [matrix] - given in lmiterm lower
    lmiterm([2 0 0 0], W_r);  %W_r'*[matrix]*W_r
    lmiterm([2 1 1 Y], A_0, 1, 's');
    %lmiterm([2 1 2 0], B_0);
    %lmiterm([2 1 3 Y], 1, C_0'); %first row  Y*A_0'+A_0*Y, B_0, Y*C_0';
    lmiterm([2 2 1 0], B_0');
    lmiterm([2 2 2 S], -zeta, 1);
    %lmiterm([2 2 3 0], D_dd'); %second row  B_0', -zeta*S, D_dd';
    lmiterm([2 3 1 Y], C_0, 1);
    lmiterm([2 3 2 0], D_dd);
    lmiterm([2 3 3 Sigma], -zeta, 1); %third row C_0*Y, D_dd, -zeta*Sigma;

    % added LMIs needed to validate that X = Y^-1; S = Sigma^-1;
    %LMI3 is
    %   [X_roof   I; I   Y_roof] > 0

    lmiterm([3 1 1 X_roof], -1, 1);
    %lmiterm([3 1 2 0], -eye(n_Xroof) );
    lmiterm([3 2 1 0], eye(n_Xroof) );
    lmiterm([3 2 2 Y_roof], -1, 1);

    lmiterm([-4 1 1 X_roof], 1, 1);         %LMI2 : X > 0
    lmiterm([-5 1 1 Y_roof], 1, 1);         %LMI3 : Y > 0

    lmisys1 = getlmis;

    %% some additional preparations before solving
    %algorithm described in the book Д. В. Баландин, М. М. Коган 
    %  "Синтез законов управления на основе линейных матричных неравенств" on a page 153

    G1 = rand*sX_roof;
    G2 = rand*sY_roof;

    setlmis(lmisys1);
    %added LMI needed for algorithm

    addedlmi = newlmi; %this one is Gamma < lambda * I

    lmiterm([addedlmi 1 1 X_roof], 1, 1);
    lmiterm([addedlmi 1 1 0], 2*G1);
    lmiterm([addedlmi 1 1 Y_roof], G1, G1);
    lmiterm([addedlmi 1 1 Y_roof], 1, 1);
    lmiterm([addedlmi 1 1 0], 2*G2);
    lmiterm([addedlmi 1 1 X_roof], G2, G2);
    lmiterm([-addedlmi 1 1 Lambda], 1, 1);

    lmisys1 = getlmis;
    %% prepare x_init
    %{
    for j=1:n
        [Xj_roof, Yj_roof] = defcx(lmisys1, j, X_roof, Y_roof);

        c(j) = trace(Xj_roof) + trace(Yj_roof); 
    end

    options=zeros(1,5);

    [~, x_init] = mincx(lmisys1, c);
    %}
    %% find X_roof, Y_roof, G1, G2 such that X = Y^-1, S = Sigma^-1 
    n = decnbr(lmisys1);
    xinit = rand(1, n);
    copt = err +1;
    c = zeros(n, 1);
    c(sLambda) = 1;
    options=zeros(1,5);
    options(5)= 1;
    while (copt > err)

        [copt, xopt] = mincx(lmisys1, c, options, xinit, err);

        valX_roof = dec2mat(lmisys1, xopt, X_roof);
        valY_roof = dec2mat(lmisys1, xopt, Y_roof);

        G1 = -inv(valY_roof);
        G2 = -inv(valX_roof);
        %xinit = xopt/2;
        %update LMIsys with new G1 and G2
        lmisys1 = dellmi(lmisys1, addedlmi);
        setlmis(lmisys1);

        [Lambda, ~, ~] = lmivar(2, [1, 1]);

        addedlmi = newlmi;

        lmiterm([addedlmi 1 1 X_roof], 1, 1); %Gamma < lambda * I
        lmiterm([addedlmi 1 1 0], 2*G1);
        lmiterm([addedlmi 1 1 Y_roof], G1, G1);
        lmiterm([addedlmi 1 1 Y_roof], 1, 1);
        lmiterm([addedlmi 1 1 0], 2*G2);
        lmiterm([addedlmi 1 1 X_roof], G2, G2);
        lmiterm([-addedlmi 1 1 Lambda], 1, 1);
        %disp(copt);
        lmisys1 = getlmis;
    end


    valX = dec2mat(lmisys1, xopt, X);
    valY = dec2mat(lmisys1, xopt, Y);
    valS = dec2mat(lmisys1, xopt, S);
    valSigma = dec2mat(lmisys1, xopt, Sigma);

    Q = [B'*valX, zeros(n_u + k, n_d), D_2'];
    Psi = [A_0'*valX + valX*A_0, valX*B_0, C_0'; B_0'*valX, -zeta*valS, D_dd'];
    Psi = [Psi; C_0, D_dd, -zeta*valSigma];

    %finally find regulator
    Theta = basiclmi(Psi, Q, P);

end

