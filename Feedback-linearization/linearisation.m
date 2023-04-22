


%% Ravnotezno stanje
m = load('parameters');

syms iL uC u real
syms L C R E positive


diL = (1-u)*uC/L + u*E/L;
duC = -(1-u)*iL/C - uC/R/C;

uC = m.uC0;

u0 = solve(diL, u);
disp('u0 :');
disp(u0);
iL0 = solve(duC, iL);
disp('iL0 :');
disp(iL0);


R = m.R;
L = m.L;
C = m.C;
E = m.E;

u0 = vpa(subs(u0), 3);
u = u0;
disp('u0 :');
disp(u0);

iL0 = vpa(subs(iL0), 3);
disp('iL0 :');
disp(iL0);

%% Jakobijan Linearizacija

syms iL uC u real
syms L C R E positive

diL = (1-u)*uC/L + u*E/L;
duC = -(1-u)*iL/C - uC/R/C;

A = [diff(diL, iL) diff(diL, uC);
     diff(duC, iL) diff(duC, uC)];

B = [diff(diL, u); 
     diff(duC, u)];
 
Cm = [0 1];

disp('A:');
Ad = vpa(subs(A, u, u0), 2);
disp(Ad);

disp('B:');
Bd = vpa(subs(B, u, u0), 2);
Bd = vpa(subs(Bd, uC, m.uC0), 2);
Bd = vpa(subs(Bd, iL, iL0), 2);
disp(Bd);

disp('Cm:');
disp(Cm);
syms s 

G = Cm*(s*eye(2) - A)^-1*B;
G = simplify(G);
disp('G:');
disp(G);
Gd = subs(G, u, u0);
Gd = subs(Gd, uC, m.uC0);
Gd = vpa(subs(Gd, iL, iL0), 2);
disp(Gd);



%% Numeri?ki

uC = m.uC0;
iL = iL0;
u = u0;

R = m.R;
L = m.L;
C = m.C;
E = m.E;


A = vpa(eval(A), 3);
disp('A:');
disp(A);

B = vpa(eval(B), 3);
disp('B:');
disp(B);

Cm = [0 1];

disp('Cm:');
disp(Cm);

G = vpa(eval(G), 3);

disp('G:');
disp(G);

symExp(s) = G;
ExpFun = matlabFunction(symExp);
ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
G = tf(ExpFun(tf('s')));
%%
u0 = double(u0);
uC0 = m.uC0;
iL0 = double(iL0);
save('linear_model', 'G', 'u0', 'iL0', 'uC0');
