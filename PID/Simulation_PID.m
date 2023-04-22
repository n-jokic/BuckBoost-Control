%% Simulation

load('parameters')
load('linear_model')
load('controler_inversion')

t_step = 0.00001;

y0 = uC0;
iL0 = 0;
uC0 = 0;
r_on = 1;
d_delta = 0;

r_delta = 0;

[ num , ~ ] = tfdata(K_PID);
num = num{1};

Kp = num(2);
Ti = 1/(num(3)/Kp);
Td = num(1)/Kp;
[ num , den ] = tfdata(K);
num = num{1};
den = den{1};
%Parameter initialisation
d_on = 0;
u_manual = u0;
d_time = 0;
d_val = 0;

t_end = 1;
t_switch = 1;

n_on = 0;
sigma_n = 1;

r_time = t_end;
r_val = -1;

sel = 1;

% Kp = -  0.8139;
% Ti =  0.0011;
% Td = 4.0647e-04;

%% referenca i poremecaj
t_switch = 0.4;

r_val = 0.1*abs(y0);
r_on = 1;
r_time = t_switch;
r_delta = 1;
d_on = 1;
n_on = 0;
d_val = I0*0.2;
d_delta = r_delta;
d_time = t_switch + 20*r_delta;
sigma_n = abs(y0/100*5/3);

t_end = d_time + 10*d_delta;


sim('.\model\inverzija.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 4*r_delta, 1); 
f = figure(301);
f.Name = 'Reference Response K_PID';
figure(f);
    subplot(2, 1, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$');
        legend('boxoff');
        
    hold off;
    subplot(2, 1, 2);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);

N_low = find(t >=d_time*0.98, 1);
N_high = find(t >=d_time + 1.5*d_delta, 1); 
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
f = figure(302);
f.Name = 'Disturbance Response K_PID';
figure(f);
    subplot(2, 1, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$');
        legend('boxoff');
    hold off;
    subplot(2, 1, 2);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);
 
 if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
 %% sum
 n_on = 1;
 
 sim('.\model\inverzija.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 4*r_delta, 1); 
f = figure(303);
f.Name = 'Reference + Noise Response K_PID';
figure(f);
    subplot(2, 1, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$');
        legend('boxoff');
    hold off;
    subplot(2, 1, 2);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);

N_low = find(t >=d_time*0.98, 1);
N_high = find(t >=d_time + 1.5*d_delta, 1); 
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
f = figure(304);
f.Name = 'Disturbance + Noise Response K_PID';
figure(f);
    subplot(2, 1, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
   p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        legend(p, '$r(t)$');
        legend('boxoff');
        xlim([t(N_low) t(N_high)]);
    hold off;
    subplot(2, 1, 2);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);
 if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end       
%% Robustnost:

load('parameters')
L_r = [L*0.8, L*0.9, L, 1.1*L];
n_on = 0;



names = {'$0.8L$', '$0.9L$', '$L$', '$1.1L$', '$1.2L$'};

for k = 1 : length(L_r)
    L = L_r(k);


    sim('.\model\inverzija.slx');


    N_low = find(t >=t_switch*0.95, 1);
    N_high = find(t >=t_switch + 4*r_delta, 1); 
    f = figure(305);
    f.Name = 'Reference Response + Robust K_PID';
    figure(f);
    subplot(2, 1, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high));
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
    hold on;
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$');
        legend('boxoff');       
        
    subplot(2, 1, 2);
    hold on;
    plot(t(N_low : N_high), u_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);

    N_low = find(t >=d_time*0.95, 1);
    N_high = find(t >=d_time + 1.5*d_delta, 1); 


    f = figure(306);
    f.Name = 'Disturbance Response + Robust K_PID';
    figure(f);
        subplot(2, 1, 1);
        plot(t(N_low : N_high), y_out(N_low : N_high));   
        hold on;
        p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$');
        legend('boxoff');


    subplot(2, 1, 2);
    hold on;
    plot(t(N_low : N_high), u_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);
    
    
end

f = figure(305);
    legend(names);
    legend('boxoff');
    
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(306);
    legend(names);
    legend('boxoff');
    
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
    
%% Pocetni uslov poremecaj

d_on = 0;
r_on = 0;

t_switch = 0;

t_end = 0.3;
t_step = 0.0001;


L = L_r(3);

iL_r = -20:10:20;
uC_r = -80:10:20;

f = figure(307);
f.Name = 'State Space K_PID';
hold all;

for i = 1:length(iL_r)
    for j = 1:length(uC_r)
        iL0 = iL_r(i);
        uC0 = uC_r(j);
        sim('.\model\inverzija.slx');
        plot(iL_out, y_out);
    end
end

xlim([-24 26]);
ylim([-134 52]);
xlabel('$i_L(t)$');
ylabel('$u_C(t)$'); 
grid('on');

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end