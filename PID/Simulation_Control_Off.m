load('parameters')
load('linear_model')
load('controler_inversion')

y0 = uC0;
iL0 = 0;
uC0 = 0;
r_on = 1;
d_delta = 0;

r_delta = 0;


[ num , den ] = tfdata(K);

Kp = 1;
Ti = Kp;
Td = Kp;

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
r_on = 0;

r_time = t_end;
r_val = -1;

sel = 1;


%% Open loop step response:
t_end = 0.39;
t_step = 0.00001;
t_switch = t_end;

sim('.\model\inverzija.slx');


f = figure(231);
f.Name = 'Step Response';
figure(f);
    plot(t, y_out, 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$j(t)$ [V]');
        xlim([t(1) t(end)]);
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end       
%% Pocetni uslovi:
t_end = 1;
t_step = 0.0001;
t_switch = t_end;

iL_r = -6:6:6;
uC_r = -10:10:10;

f = figure(232);
f.Name = 'State Space';
hold all;

for i = 1:length(iL_r)
    for j = 1:length(uC_r)
        iL0 = iL_r(i);
        uC0 = uC_r(j);
        sim('.\model\inverzija.slx');
        plot(iL_out, y_out);
            xlabel('$i_L(t)$');
            ylabel('$u_C(t)$'); 
            grid('on');
    end
end
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end