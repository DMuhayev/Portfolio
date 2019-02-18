function dd = Stepanova_sys(t, y, a, b, m1, m2, sigma)

% STEPANOVA represents Stepanova model of tumor-immune system
% a, b, m0, m1, sigma - constants from model (usually real numbers)

dd = zeros(2, 1);
dd(1) = y(1)*(a - y(2));
dd(2) = b*y(1)*y(2) - m1*y(2) - m2*(y(1)^2)*y(2) + sigma;