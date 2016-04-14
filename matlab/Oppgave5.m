% Oppgave 5

E = 1.3e10;
D = 480;
w = 0.3;
d = 0.03;
L = 2;

% faktisk løsning, kun egenvekt:
syms y(x);
I = (w*d^3)/12;
f = -9.81*D*w*d;
y(x) = (f/(24*E*I))*x^2*(x^2 - 4*L*x + 6*L^2);

i_max = 11;

n = zeros(i_max, 1);
y_num_L = zeros(i_max, 1);
y_actual_L = zeros(i_max, 1);
error = zeros(i_max, 1);
cond_A = zeros(i_max, 1);
for i = (1:i_max)
    n(i) = 10*2^i;
    disp(n(i));
    y_num = eulerbernoulli(E, D, w, d, L, n(i));
    y_num_L(i) = y_num(n(i));
    y_actual_L(i) = y(L);
    error(i) = abs(y(L) - y_num(n(i)));
    A = lagA(n(i));
    cond_A(i) = condest(A);       
end
T = table(n, y_num_L, y_actual_L, error, cond_A);
disp(T);
figure;
subplot(2, 1, 1);
plot(log(n), log(error));
title('Logaritmisk plot av n og feil');
subplot(2, 1, 2);
plot(log(cond_A));
title('Logaritmisk plot av kondisjonstallet til A av størrelse n'),