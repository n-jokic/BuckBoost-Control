%% Simulation

load('parameters')
load('linear_model')
load('controler')

t_step = 0.00001;
I0 = iL0;
y0 = uC0;
iL0 = 0;
uC0 = 0;
r_on = 1;
d_delta = 0;

r_delta = 0;

[ num , ~ ] = tfdata(K2);
num = num{1};

Kp = num(1);
Ti = 1/(num(2)/Kp);

[ num , den ] = tfdata(K1);
num1 = num{1};
den1 = den{1};
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



[ num , den ] = tfdata(K1);
num = num{1};
den = den{1};
sigma_n = 0.4;

%% referenca i poremecaj
t_switch = 0.35;

r_val = 0.1*abs(y0);
r_on = 1;
r_time = t_switch;
r_delta = 0.6;
d_on = 1;
n_on = 0;
d_val = I0*0.2;
d_delta = r_delta;
d_time = t_switch + 6*r_delta;

t_end = d_time + 5*d_delta;
t_end = d_time + 10*d_delta;

sim('.\model\loop.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 1.98*r_delta, 1); 
N_mid = find(t >=t_switch + 0.98*r_delta, 1); 
N_s = find(t >=t_switch, 1);

[Ts,Tr] = findTsTr(y_out, t, N_s, N_mid);
disp('up step')
disp('Ts:');
disp(Ts);
disp('Tr:');
disp(Tr);

Ts_a = Ts;
Tr_a = Tr;

disp('----')
N_s = find(t >=t_switch + 2*r_delta, 1);
N_mid = find(t >=t_switch + 2.98*r_delta, 1); 
[Ts,Tr] = findTsTr(y_out, t, N_s, N_mid);
disp('down step');
disp('Ts:');
disp(Ts);
disp('Tr:');
disp(Tr);
disp('----')

Ts_a = Ts/2+Ts_a/2;
Tr_a = Tr/2+Tr_a/2;

disp('Ts_avg:');
disp(Ts_a);
disp('Tr_avg:');
disp(Tr_a);
disp('----')
%%
f = figure(401);
f.Renderer = 'painters';
f.Name = 'Reference_Response_K_PID';
figure(f);
    subplot(2, 2, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$', 'Location', 'best');
        legend('boxoff');
        
    hold off;
    subplot(2, 2, 3);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);
        ylim([0, 1]);
        
   subplot(2, 2, 2);  
        plot(t(N_low : N_high), iL_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([t(N_low) t(N_high)]);
    subplot(2, 2, 4);  
        plot(iL_out(N_low : N_high), y_out(N_low : N_high), 'black');
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
        
%     axes('position',[.135 .59 .135 .135]);
%         indexOfInterest = (t < r_time + 0.15) & (t > r_time); % range of t near perturbation
%         axis tight;
%         plot(t(indexOfInterest),y_out(indexOfInterest), 'black') % plot on new axes
%         set(gca,'xticklabel',[])
%         set(gca,'yticklabel',[])
%         grid on;
%         xlim([min(t(indexOfInterest)) max(t(indexOfInterest))]);
set(findall(gcf,'-property','FontSize'),'FontSize',font_size);

N_low = find(t >=d_time*0.98, 1);
N_high = find(t >=d_time + 1.98*d_delta, 1); 
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
f = figure(402);
f.Renderer = 'painters';
f.Name = 'Disturbance_Response_K_PID';
figure(f);
subplot(2, 2, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        legend(p, '$r(t)$', 'Location', 'best');
        legend('boxoff');
    hold off;
    subplot(2, 2, 3);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([t(N_low) t(N_high)]);
        ylim([0, 1]);
        
    subplot(2, 2, 2);  
        plot(t(N_low : N_high), iL_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([t(N_low) t(N_high)]);
    subplot(2, 2, 4);  
        plot(iL_out(N_low : N_high), y_out(N_low : N_high), 'black');
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
 
set(findall(gcf,'-property','FontSize'),'FontSize',font_size);
 
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end 
 %% sum
%  n_on = 1;
% 
%  sim('.\model\loop.slx');
% 
% 
% N_low = find(t >=t_switch*0.95, 1);
% N_high = find(t >=t_switch + 4*r_delta, 1); 
% f = figure(403);
% f.Name = 'Reference_+_Noise_Response_K_PID';
% figure(f);
%     subplot(2, 1, 1);
%     plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
%     hold on;
%     p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$y(t)$ [V]');
%         xlim([t(N_low) t(N_high)]);
%         legend(p, '$r(t)$', 'Location', 'best');
%         legend('boxoff');
%     hold off;
%     subplot(2, 1, 2);
%     plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$u(t)$');
%         xlim([t(N_low) t(N_high)]);
% 
% N_low = find(t >=d_time*0.98, 1);
% N_high = find(t >=d_time + 1.5*d_delta, 1); 
% if(SAVE)
%     saveas(f,[path '\' f.Name '.eps'],'epsc');
% end
% f = figure(404);
% f.Name = 'Disturbance_+_Noise_Response_K_PID';
% figure(f);
%     subplot(2, 1, 1);
%     plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
%     hold on;
%    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$y(t)$ [V]');
%         legend(p, '$r(t)$', 'Location', 'best');
%         legend('boxoff');
%         xlim([t(N_low) t(N_high)]);
%     hold off;
%     subplot(2, 1, 2);
%     plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$u(t)$');
%         xlim([t(N_low) t(N_high)]);
% if(SAVE)
%     saveas(f,[path '\' f.Name '.eps'],'epsc');
% end        
% %% Robustnost:
% 
% load('parameters')
% L_r = [L*0.8, L*0.9, L, L*1.1, L*1.2];
% n_on = 0;
% t_switch = 0.6;
% 
% r_on = 1;
% 
% d_on = 1;
% n_on = 0;
% 
% 
% names = {'$0.8L$', '$0.9L$', '$L$', '1.1L', '1.2L'};
% 
% for k = 1 : length(L_r)
%     L = L_r(k);
% 
% 
%     sim('.\model\loop.slx');
% 
% 
%     N_low = find(t >=t_switch*0.95, 1);
%     N_high = find(t >=t_switch + 4*r_delta, 1); 
%     f = figure(405);
%     f.Name = 'Reference_Response_+_Robust_K_PID';
%     figure(f);
%     subplot(2, 1, 1);
%     plot(t(N_low : N_high), y_out(N_low : N_high));
%     hold on;
%     p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
% 
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$y(t)$ [V]');
%         xlim([t(N_low) t(N_high)]);
%         legend(p, '$r(t)$', 'Location', 'best');
%         legend('boxoff');       
%         
%     subplot(2, 1, 2);
%     hold on;
%     plot(t(N_low : N_high), u_out(N_low : N_high));
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$u(t)$');
%         xlim([t(N_low) t(N_high)]);
% 
% N_low = find(t >=d_time*0.98, 1);
% N_high = find(t >=d_time + 1.5*d_delta, 1); 
% 
% 
%     f = figure(406);
%     f.Name = 'Disturbance_Response_+_Robust_K_PID';
%     figure(f);
%         subplot(2, 1, 1);
%         plot(t(N_low : N_high), y_out(N_low : N_high));
%         hold on;
%         p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
%         hold on;
%     
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$y(t)$ [V]');
%         xlim([t(N_low) t(N_high)]);
%         legend(p, '$r(t)$', 'Location', 'best');
%         legend('boxoff');
% 
% 
%     subplot(2, 1, 2);
%     hold on;
%     plot(t(N_low : N_high), u_out(N_low : N_high));
%         grid on;
%         xlabel('$t$ [s]');
%         ylabel('$u(t)$');
%         xlim([t(N_low) t(N_high)]);
%     
%     
% end
% f = figure(405);
%     legend(names);
%     legend('boxoff');
% if(SAVE)
%     saveas(f,[path '\' f.Name '.eps'],'epsc');
% end
% f = figure(406);
%     legend(names);
%     legend('boxoff');
% if(SAVE)
%     saveas(f,[path '\' f.Name '.eps'],'epsc');
% end
%     
% %% Pocetni uslov poremecaj
% 
% d_on = 0;
% r_on = 0;
% 
% t_switch = 0;
% 
% t_end = 0.3;
% 
% 
% 
% L = L_r(3);
% 
% iL_r = -20:10:20;
% uC_r = -80:10:20;
% 
% f = figure(408);
% f.Name = 'State_Space_K_PID';
% hold all;
% 
% for i = 1:length(iL_r)
%     for j = 1:length(uC_r)
%         iL0 = iL_r(i);
%         uC0 = uC_r(j);
%         sim('.\model\loop.slx');
%         plot(iL_out, y_out);
%     end
% end
% 
% xlim([-24 26]);
% ylim([-134 52]);
% xlabel('$i_L(t)$');
% ylabel('$u_C(t)$'); 
% grid('on');
% 
% if(SAVE)
%     saveas(f,[path '\' f.Name '.eps'],'epsc');
%end