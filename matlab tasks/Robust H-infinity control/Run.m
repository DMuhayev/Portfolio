%% init data from example from book
%%{
err = 0.0001; %U can touch this
%example data

m_0 = 1;
b_0 = 1;
c_0 = 100;
w_d = 0.1;

%dim sizes

n_x = 2;
n_d = 3;
n_u = 1;
n_y = 1;

%matrices

A = [0, 1; -c_0/m_0, -b_0/m_0];
B_d = [0, 0, 0; w_d, w_d, w_d];
B_2 = [0; 1/m_0];
C_d = [c_0/m_0, b_0/m_0; -c_0/m_0, 0; 0 -b_0/m_0];
D_dd = [-w_d, -w_d, -w_d; 0, 0, 0; 0, 0, 0];
D_d2 = [-1/m_0; 0; 0];
D_2d = zeros(n_y, n_d);
C_2 = [1 0];

% /delta(t)
r = 3;
f = 0;

k_r = [1 1 1];
m_f = zeros(1, f);
%}

%% these can be changed by user (algorithm doesn't guarantee that good result exists)
eta = 0.4;
zeta = 1/eta;
k = 0;

[valX, valY, valS, valSigma, Theta]  = RegFind(err, k, eta, A, B_2, C_2, C_d, B_d, D_dd, D_2d, D_d2, r, f, k_r, m_f);

disp("X = ");
disp(valX);
disp("Y = ");
disp(valY);
disp("X*Y = ");
disp(valX*valY);
disp("S = ");
disp(valS);
disp("Sigma = ");
disp(valSigma);
disp("S*Sigma = ");
disp(valS*valSigma);
disp("Theta = ");
disp(Theta);