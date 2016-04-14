E = 1.3e10;
D = 480;
w = 0.3;
d = 0.03;
L = 2;
g = -9.81;

p = 100;

syms s(x);
s(x) = -p*g*sin(pi*x/L);


syms y(x);
I = (w*d^3)/12;
y(x) = (s(x)/(24*E*I))*x^2*(x^2 - 4*L*x + 6*L^2) - ((g*p*L)/(E*I*pi))*(L^3/pi^3 * sin(pi*x/L) - x^3/6 + L*x^2/2 - L^2*x/pi^2);

n = zeros(11, 1);
y_num_L = zeros(11, 1);
y_actual_L = zeros(11, 1);
error = zeros(11, 1);
cond_A = zeros(11, 1);
for i = (1:11)
    n(i) = 10*2^i;
    disp(n(i));
    b = eval((s(L/n(i):L/n(i):L)))';
    y_num = eulerbernoulli(E, D, w, d, L, n(i), b);
    y_num_L(i) = y_num(n(i));
    y_actual_L(i) = y(L);
    error(i) = abs(y(L) - max(y_num(n(i))));
    A = lagA(n(i));
    cond_A(i) = cond(A, inf);
end
T = table(n, y_num_L, y_actual_L, error, cond_A);
disp(T);

plot(log(n), log(error));
hold on;
plot(cond_A*eps);
hold on;
plot(L^2/n.^2);
