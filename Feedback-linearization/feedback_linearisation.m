load('parameters')
load('linear_model')
%% Postavka problema
syms x1 x2 x1dot x2dot real
syms x xdot real
syms y h u g real
syms Lfh LgLfh Lgh LfLfh real

syms L C R E positive

%x1 - struja kalema iL
%x2 - napon kondenzatora uC
f = [1/L*x2; -1/C*x1 - x2/R/C];
g = [-x2/L + E/L; x1/C];
xdot = f + g*u;
h = 0.5*L*x1^2 + 0.5*C*(x2 - E)^2;

disp('f :');
disp(f);
disp('g :');
disp(g);
disp('h: ');
disp(h);


Lgh = simplify([diff(h,x1) diff(h,x2)]*g);
disp('Lgh :');
disp(Lgh);

Lfh = expand(simplify([diff(h,x1) diff(h,x2)]*f));
disp('Lfh :');
disp(Lfh);

LgLfh = expand(simplify([diff(Lfh,x1) diff(Lfh,x2)]*g));
disp('LgLfh :');
disp(LgLfh);

LfLfh = expand(simplify([diff(Lfh,x1) diff(Lfh,x2)]*f));
disp('LfLfh :');
disp(LfLfh);

%% eFL parametri
%f(s) = s^2 + K1*s + K0

syms s
Te = 0.025; %s

f = expand((Te/2*s + 1)^2/(Te/2)^2);

K0 = coeffs(f);
K1 = double(K0(2));
K0 = double(K0(1));

disp('f_desired :');
disp(f);
disp('K0 :');
disp(K0);
disp('K1 :');
disp(K1);

wp = K1;

disp('wp :');
disp(wp);

%% eFL + I parametri 
%f(s) = (s^2 + K1*s + K0)*s + Ki

syms s


f = expand((Te/3*s + 1)^3);
c = coeffs(f);
f = f/(c(end));

K0i = coeffs(f);
K1i = double(K0i(3));
Ki = double(K0i(1));
K0i = double(K0i(2));

disp('f_desired :');
disp(f);
disp('Ki :');
disp(Ki);
disp('K0i :');
disp(K0i);
disp('K1i :');
disp(K1i);

%%
save('feedback_linearisation_params', 'K0', 'K1', 'Ki', 'K0i', 'K1i', 'wp')