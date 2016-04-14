% % Oppgave 3
format long;


E = 1.3e10;
D = 480;
w = 0.3;
d = 0.03;
L = 2;



% numerisk løsning
n = 10;
y_num = eulerbernoulli(E, D, w, d, L, n);
disp(y_num);



% faktisk løsning, kun egenvekt:
syms y(x);
I = (w*d^3)/12;
f = -9.81*D*w*d;
y(x) = (f/(24*E*I))*x^2*(x^2 - 4*L*x + 6*L^2);

y_e = y(L/n:L/n:L)';

