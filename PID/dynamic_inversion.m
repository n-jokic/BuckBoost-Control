

load('parameters')
load('linear_model')
%%

%resavamo problem nestabilne nule
A = zpk(G);
A.Z{1} = - A.Z{1};
a = 5;

w0 =abs(A.Z{1}/a); %rad/s
s = tf('s');
K = 1/s;

K = -minreal(w0/s*A^-1); %Zbog nestabilnih nula
disp('K :');
tf(K)
f = figure(002);
    f.Name = 'rlocus K_inv';
    h = rlocusplot(minreal(G*K));
    p = getoptions(h);
    p.XLabel.String = '\sigma';
    p.Ylabel.String = 'j\omega';
    p.Title.String = 'GMK';
    setoptions(h,p);
%podesavanje pomocu rlocus, za K ~ 1.34 dobijamo optimalno zeta


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%%
Kp = 0.835;
K = -Kp*minreal(w0/s*A^-1); 

%Kontroler je kauzalan

f = figure(003);
f.Name = 'Bode_K_inv';

bode(-K*G, 'black', G, 'r--');
    grid('on');
    legend({'KG', 'G'});
    
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(004);
f.Name = 'MsMt_K_inv';
find_MS_MT(K*G,1);

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
save('controler_inversion', 'K');

%% PID approximation
[num,den] = tfdata(K);
syms s
K_sym = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s)*s;

PID = vpa(taylor(K_sym, s, 0, 'Order', 3)/s,4);

symExp(s) = PID;
ExpFun = matlabFunction(symExp);
ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
K_PID = tf(ExpFun(tf('s')));

save('controler_inversion', 'K', 'K_PID');

  
f = figure(005);
f.Name = 'Bode K_PID';
bode(K_PID*G, 'black', G, 'r--');
    grid('on');
    legend({'KG', 'G'});
    
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(006);
f.Name = 'MsMt K_PID';
find_MS_MT(K_PID*G, 1);

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end