

load('parameters')
load('linear_model')
%% 
close all;
f = figure(001);
    f.Name = 'GMK_G';
    h = rlocusplot(minreal(-G));
    p = getoptions(h);
    p.XLabel.String = '\sigma';
    p.Ylabel.String = 'j\omega';
    p.Title.String = 'GMK';
    setoptions(h,p);
if(SAVE)
    saveas(f,[path '\' f.Name],'epsc');
end
%Gmk sistema, kao sto mozemo da vidimo izrazena oscilatornost, ocekujem los PF iz Bode-a 

f =figure(002);
f.Name = 'Bode G';


bode(-G, 'black');
    grid('on');
    [Gm,Pm,Wgm,Wpm] = margin(-G);
    text(1000, 200, {['$\Phi_{pf}$ = ' num2str(round(Pm, 3)) '$^{\circ}$' ...
    ' $\omega_{pf}$ = ' num2str(round(Wpm, 3)) ' $\mathrm{\frac{rad}{s}}$'], ...
    ['$d$ = ' num2str(round(Gm, 3)) ' $\omega_{0}$ = ' ...
    num2str(round(Wgm, 3)) ' $\mathrm{\frac{rad}{s}}$']}...
    , 'Interpreter', 'latex');
 
%Projektujemo lead-lag:
if(SAVE)
    saveas(f,[path '\' f.Name],'epsc');
end
%% 
% 
% Pf_0 = 100;%[deg]
% syms x;
% w0 = double(vpasolve(atan2(x,606.3) + atan2(40.92*x,(1.618e04 -x^2)) - (180-Pf_0)/180*pi));
w0 = 200;

disp('\omega_0 :');
disp(w0);


Kp = 1/abs(evalfr(-G,1j*w0));




disp('Kp: ');
disp(Kp);

s = tf('s');

f =figure(003);
f.Name = 'Bode_Kp';

K_p = Kp;
bode(-K_p*G, 'black', -G, 'r--');
    grid('on');
    legend({'KG', 'G'});
    [Gm,Pm,Wgm,Wpm] = margin(-K_p*G);
    text(10, 140, {['$\Phi_{pf}$ = ' num2str(round(Pm, 3)) '$^{\circ}$' ...
    ' $\omega_{pf}$ = ' num2str(round(Wpm, 3)) ' $\mathrm{\frac{rad}{s}}$'], ...
    ['$d$ = ' num2str(round(Gm, 3)) ' $\omega_{0}$ = ' ...
    num2str(round(Wgm, 3)) ' $\mathrm{\frac{rad}{s}}$']}...
    , 'Interpreter', 'latex');

w0 = Wpm;

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%%
wn = w0/10;

Klag = (s+wn)/s;
wi = wn;

disp('wn: ');
disp(wn);

f =figure(004);
f.Name = 'Bode_K_P+lag';
K  = Klag*K_p;
Kc = 1/abs(evalfr(-K*G,1j*w0));

K = Kc*Klag*K_p;

bode(-K*G, 'black', -G, 'r--');
    grid('on');
    legend({'KG', 'G'});
    [Gm,Pm,Wgm,Wpm] = margin(-K*G);
    text(1, 140, {['$\Phi_{pf}$ = ' num2str(round(Pm, 3)) '$^{\circ}$' ...
    ' $\omega_{pf}$ = ' num2str(round(Wpm, 3)) ' $\mathrm{\frac{rad}{s}}$'], ...
    ['$d$ = ' num2str(round(Gm, 3)) ' $\omega_{0}$ = ' ...
    num2str(round(Wgm, 3)) ' $\mathrm{\frac{rad}{s}}$']}...
    , 'Interpreter', 'latex');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%Hajde da vidimo kako izgleda rlocus...

f = figure(005);
    f.Name = 'GMK KG';
    h = rlocusplot(minreal(-K*G));
    p = getoptions(h);
    p.XLabel.String = '\sigma';
    p.Ylabel.String = 'j\omega';
    p.Title.String = 'GMK';
    p.XLim = [-250, 700];
    setoptions(h,p);
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
w0 = Wpm;

%%
% K = 2;
% dF = 10; %[deg]
% dF = dF/180*pi;
% a = (1+sin(dF))/(1-sin(dF));
% w1 = w0*a^(-1/(2*K));
% % w0 = w0*a^(1/(2*K));
% wn = w0/sqrt(a);
% wp = w0*sqrt(a);
% 
% Kp = 1/abs(evalfr(-G,1j*w1));

wn = w0/1.2;
% wn = 129;

wp = 10*wn;



Klead = (s/wn+1)/(s/wp+1);

disp('wn:');
disp(wn);
disp('wp:');
disp(wp);

K  = Klag*K_p*Klead;
Kc = 1/abs(evalfr(-K*G,1j*w0));
K  = Kc*Klag*K_p*Klead;

f =figure(006);
f.Name = 'Bode_K_P+lag+lead';

bode(-K*G, 'black', -G, 'r--');
    grid('on');
    legend({'KG', 'G'});
    [Gm,Pm,Wgm,Wpm] = margin(-K*G);
    text(1, 140, {['$\Phi_{pf}$ = ' num2str(round(Pm, 3)) '$^{\circ}$' ...
    ' $\omega_{pf}$ = ' num2str(round(Wpm, 3)) ' $\mathrm{\frac{rad}{s}}$'], ...
    ['$d$ = ' num2str(round(Gm, 3)) ' $\omega_{0}$ = ' ...
    num2str(round(Wgm, 3)) ' $\mathrm{\frac{rad}{s}}$']}...
    , 'Interpreter', 'latex');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%Hajde da vidimo kako izgleda rlocus...

f = figure(007);
    f.Name = 'GMK K_finalG';
    h = rlocusplot(minreal(-K*G));
    p = getoptions(h);
    p.XLabel.String = '\sigma';
    p.Ylabel.String = 'j\omega';
    p.Title.String = 'GMK';
    p.XLim = [-250, 700];
    setoptions(h,p);
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

%%
Kfinal = -K;

f = figure(9);
f.Name = 'MsMt_K_finalG';
    find_MS_MT(Kfinal*G, 1);
    

    
save('controler', 'Kfinal');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end   

%%

% wn = w0/0.9;
% wn = 129;



% wp = 10*wn;
% Kc = 1.6*2.2;


% 
% Klead = (s/wn+1)/(s/wp+1);


K1 = Klead;
% 
% Klag = (s+wi)/s;

K2 = -Kc*Kp*Klag;

f =  figure(12);
f.Name = 'MsMt_K_PID_exact';
find_MS_MT(K2*G, K1);
s = tf('s');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
save('controler', 'Kfinal', 'K1', 'K2');
