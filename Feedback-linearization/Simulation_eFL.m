%% Simulation

load('parameters')
load('feedback_linearisation_params')
load('linear_model')
I0 = iL0;
y0 = uC0;
iL0 = 0;
uC0 = 0;
r_on = 1;
d_delta = 0;

r_delta = 0;

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


t_step = 0.00001;

AW = 0;

Kti = 1.2*Ki;
%% referenca i poremecaj AW off
t_switch = 0.4;


r_val = 0.1*abs(y0);
r_on = 1;
r_time = t_switch;
r_delta = 0.2;
d_on = 1;
n_on = 0;
d_val = I0*0.2;
d_delta = r_delta*1;
d_time = t_switch + 4*r_delta;
sigma_n = abs(y0/100*5/3);

t_end = d_time + 4*d_delta;



       

%% 


sim('.\model\feedback_lin.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 1.95*r_delta, 1); 
N_mid = find(t >=t_switch + 0.95*r_delta, 1); 
N_s = find(t >=t_switch, 1);

[Ts,Tr] = findTsTr(y_out, t, N_s, N_mid);
disp('up step')
disp('Ts:');
disp(Ts);
disp('Tr:');
disp(Tr);

disp('----')
N_s = find(t >=t_switch + 2*r_delta, 1);
N_mid = find(t >=t_switch + 2.95*r_delta, 1); 
[Ts,Tr] = findTsTr(y_out, t, N_s, N_mid);
disp('down step');
disp('Ts:');
disp(Ts);
disp('Tr:');
disp(Tr);
disp('----')


f = figure(101);
f.Renderer = 'painters';
f.Name = 'Reference_Response_K_eFL';
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
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
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

set(findall(gcf,'-property','FontSize'),'FontSize', 8)
N_low = find(t >=d_time*0.98, 1);
N_high = find(t >=d_time + 1.98*d_delta, 1); 

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
% 
f = figure(102);
f.Renderer = 'painters';
f.Name = 'Disturbance_Response_K_eFL';
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
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
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
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


 %% sum
 n_on = 1;
 
 
sim('.\model\feedback_lin.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 1.95*r_delta, 1); 
f = figure(103);
f.Name = 'Reference_+_Noise_Response_K_eFL';
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
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
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

    

N_low = find(t >=d_time*0.98, 1);
N_high = find(t >=d_time + 1.98*d_delta, 1); 

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(104);
f.Name = 'Disturbance_+_Noise_Response_K_eFL';
figure(f);
    subplot(2, 2, 1);
    plot(t(N_low : N_high), y_out(N_low : N_high), 'black');
    hold on;
    p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'r--');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        legend(p, '$r(t)$', 'Location', 'best');
        legend('boxoff');
        xlim([t(N_low) t(N_high)]);
    hold off;
    subplot(2, 2, 3);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
        
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

f = figure(104);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
        
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
        
%% Robustnost:
load('parameters')
L_r = [L*0.8, L, L*1.2];
n_on = 0;

d_on = 1;
n_on = 0;

zoom = [];
names = {'$0.8L$','$L$', '$1.2L$'};

for k = 1 : length(L_r)
    L = L_r(k);

    sim('.\model\feedback_lin.slx');


    N_low = find(t >=t_switch*0.95, 1);
    N_high = find(t >=t_switch + 1.98*r_delta, 1);  
    
    f = figure(105);
    f.Name = 'Reference_Response_+_Robust_K_eFL';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(N_low : N_high), y_out(N_low : N_high));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');
        
    subplot(2, 2, 3);
    hold on;
    plot(t(N_low : N_high), u_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
        
    subplot(2, 2, 2);  
    hold on;
        plot(t(N_low : N_high), iL_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([t(N_low) t(N_high)]);
    subplot(2, 2, 4);  
        hold on;
        plot(iL_out(N_low : N_high), y_out(N_low : N_high));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');

    N_low = find(t >=d_time*0.98, 1);
    N_high = find(t >=d_time + 1.98*d_delta, 1); 

    
    f = figure(106);
    f.Renderer = 'painters';
    f.Name = 'Disturbance_Response_+_Robust_K_eFL';
    figure(f);
        subplot(2, 2, 1);
        hold on;
        plot(t(N_low : N_high), y_out(N_low : N_high));
        hold on;
        %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([t(N_low) t(N_high)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');


    subplot(2, 2, 3);
    hold on;
    plot(t(N_low : N_high), u_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
     subplot(2, 2, 2);
     hold on;
        plot(t(N_low : N_high), iL_out(N_low : N_high));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([t(N_low) t(N_high)]);
    subplot(2, 2, 4);  
    hold on;
        plot(iL_out(N_low : N_high), y_out(N_low : N_high));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indexOfInterest = (t < r_time + 0.046) & (t > r_time-0.002); % range of t near perturbation

   f = figure(1050);
    time = t(indexOfInterest);
    f.Name = 'Reference_Response_+_Robust_K_eFL_zoom';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(indexOfInterest), total_ref(indexOfInterest), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([time(1) time(end)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');
        
    subplot(2, 2, 3);
    hold on;
    plot(t(indexOfInterest), u_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([time(1) time(end)]);
        
    subplot(2, 2, 2);  
    hold on;
        plot(t(indexOfInterest), iL_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([time(1) time(end)]);
    subplot(2, 2, 4);  
        hold on;
        plot(iL_out(N_low : N_high), y_out(N_low : N_high));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
        
    N_low = find(t >=d_time*0.98, 1);
    N_high = find(t >=d_time + 1.98*d_delta, 1); 

    indexOfInterest = (t < d_time + 0.09) & (t > d_time-0.002); % range of t near perturbation
    time = t(indexOfInterest);
    f = figure(1060);
    f.Renderer = 'painters';
    f.Name = 'Disturbance_Response_+_Robust_K_eFL_zoom';
    figure(f);
        subplot(2, 2, 1);
        hold on;
        plot(t(indexOfInterest), y_out(indexOfInterest));
        hold on;
        %p = plot(t(indexOfInterest), total_ref(indexOfInterest), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([time(1) time(end)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');


    subplot(2, 2, 3);
    hold on;
    plot(t(indexOfInterest), u_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([time(1) time(end)]);
     subplot(2, 2, 2);
     hold on;
        plot(t(indexOfInterest), iL_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([time(1) time(end)]);
    subplot(2, 2, 4);  
    hold on;
        plot(iL_out(indexOfInterest), y_out(indexOfInterest));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
    
    
end
f = figure(105);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
%       axes('position',[.135 .59 .135 .135]);
%         indexOfInterest = (t < r_time + 0.15) & (t > r_time); % range of t near perturbation
%         axis tight;
%         plot(t(indexOfInterest),zoom) % plot on new axes
%         set(gca,'xticklabel',[])
%         set(gca,'yticklabel',[])
%         grid on;
%         xlim([min(t(indexOfInterest)) max(t(indexOfInterest))]);
%         box on; % put box around new pair of axes
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    
f = figure(1050);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
%       axes('position',[.135 .59 .135 .135]);
%         indexOfInterest = (t < r_time + 0.15) & (t > r_time); % range of t near perturbation
%         axis tight;
%         plot(t(indexOfInterest),zoom) % plot on new axes
%         set(gca,'xticklabel',[])
%         set(gca,'yticklabel',[])
%         grid on;
%         xlim([min(t(indexOfInterest)) max(t(indexOfInterest))]);
%         box on; % put box around new pair of axes
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    
    
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
    f = figure(105);
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
f = figure(106);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    
f = figure(1060);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
    f = figure(106);
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
    
%% Pocetni uslov poremecaj

d_on = 0;
r_on = 0

t_switch = 0;

t_end = 0.3;
t_step = 0.0001;


L = L_r(3);

iL_r = -20:5:20;
uC_r = -80:5:20;

f = figure(107);
f.Name = 'State_Space_K_eFL';
hold all;

for i = 1:length(iL_r)
    for j = 1:length(uC_r)
        iL0 = iL_r(i);
        uC0 = uC_r(j);
        sim('.\model\feedback_lin.slx');
        plot(iL_out, y_out);

    end
end

xlim([-21 36]);
ylim([-170 30]);
xlabel('$i_L(t) [A]$');
ylabel('$u_C(t) [V]$'); 
grid('on');

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%


