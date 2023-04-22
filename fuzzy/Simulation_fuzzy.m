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

Kti = Ki;
Ki = 0.5*Ki;

K_r = 0;
Kif = 1;
Kn = 1;
Ku = 1;
Kp = 1;
Kd = 1;
K = 1;

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

t_end = d_time;
sigma_n = abs(y0/100*5/3);

%%
r_sel = 2;
sel = 1;
Kn = 1;

Ki = 0;

Ki = 0;
Kd = 1;
Kp = 1;
sim('.\model\feedback_lin_2018.slx');

t_start = find(t >= 0.08, 1);

er = total_ref-y_out;

error = (max(abs(er(t_start:end))));
Kn = 1/error;

Kp = 0.0089;
Kd = 0.000112;
Ki = Kp/0.0094;
Tf = 0.0011;


Ku = Kp/Kn;
Kp = -Kn;
Kd = -Kd/Ku;
Ki = -Ki/Ku;

Ki = 0.81*Ki;
Kd = 1.1*Kd;
Kp = 0.95*Kp;
r_sel = 2;
sel = 2;



sim('.\model\feedback_lin_2018.slx');



f = figure(100);
f.Renderer = 'painters';
f.Name = 'error_range_fuzzy';

figure(f);
    plot(t, e, 'k');
    hold on;
    xlim([t(1), t(end)]);
    plot(t, ones(length(t), 1), 'r--');
    plot(t, -ones(length(t), 1), 'r--');
    hold off;
        ylim([-2, 2]);
        xlabel('$t$ [s]');
        ylabel('$e(t)$');
        grid('on');
set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(101);
f.Renderer = 'painters';
f.Name = 'diff_error_range_fuzzy';
e = total_ref-y_out;
figure(f);
    plot(t, ed, 'k');
    hold on;
    xlim([t(1), t(end)]);
    plot(t, ones(length(t), 1), 'r--');
    plot(t, -ones(length(t), 1), 'r--');
    hold off;
        ylim([-2, 2]);
        xlabel('$t$ [s]');
        ylabel('$e_d(t)$');
        grid('on');
set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


%%
r_sel = 2;
sel = 1;
Kn = 1;

Ki = 0;
Kd = 1;
Kp = 1;

sim('.\model\feedback_lin_2018.slx');

t_start = find(t >= 0.08, 1);


er = total_ref-y_out;

error = (max(abs(er(t_start:end))));
Kn = 1/error;

Kp = 0.0089;
Kd = 0.000112;
Ki = Kp/0.0094;

Ku = Kp/Kn;
Kp = -Kn;
Kd = -Kd/Ku;
Ki = -Ki/Ku;
r_sel = 2;
sel = 3;



sim('.\model\feedback_lin_2018.slx');



f = figure(200);
f.Renderer = 'painters';
f.Name = 'error_range_fuzzy_modif';

figure(f);
    plot(t, e_1, 'k');
    hold on;
    plot(t, ones(length(t), 1), 'r--');
    plot(t, -ones(length(t), 1), 'r--');
    hold off;
        xlim([t(1), t(end)]);
        ylim([-2, 2]);
        xlabel('$t$ [s]');
        ylabel('$e(t)$');
        grid('on');
set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

f = figure(201);
f.Renderer = 'painters';
f.Name = 'diff_error_range_fuzzy_modif';
e = total_ref-y_out;
figure(f);
    plot(t, ed_1, 'k');
    xlim([t(1), t(end)]);
    hold on;
    plot(t, ones(length(t), 1), 'r--');
    plot(t, -ones(length(t), 1), 'r--');
    hold off;
        ylim([-2, 2]);
        xlabel('$t$ [s]');
        ylabel('$e_d(t)$');
        grid('on');
set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end




%%
standard_input_range = [-1 1];
number_of_inputs = 2;
number_of_input_mfs = 3;
min_ir = min(standard_input_range);
max_ir = max(standard_input_range);
input_triang_mfs = [-2 -1 0; -1 0 1; 0 1 2];
step = (max_ir - min_ir) / (number_of_input_mfs - 1);

min_out_range = number_of_inputs * min_ir;
max_out_range = number_of_inputs * max_ir;
standard_output_range = [min_out_range, max_out_range];
output_singleton_mfs = [min_out_range: step: max_out_range]';


lin_fpd = sugfis('Name','pravila');
lin_fpd.andMethod    = 'prod';
lin_fpd.orMethod     = 'probor';
lin_fpd.impMethod    = 'prod';
lin_fpd.aggMethod    = 'sum';
lin_fpd.defuzzMethod = 'wtaver';

lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
lin_fpd = addMF(lin_fpd, 'ê', 'trimf', input_triang_mfs(1,:), 'Name', 'NEG');
lin_fpd = addMF(lin_fpd, 'ê', 'trimf', input_triang_mfs(2,:), 'Name', 'NULA');
lin_fpd = addMF(lin_fpd, 'ê', 'trimf', input_triang_mfs(3,:), 'Name', 'POZ');

f = figure(1); plotmf(lin_fpd,'input',1); ylabel('Stepen pripadnosti');

f.Renderer = 'painters';
f.Name = 'error_membership_fPID';

set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
lin_fpd = addMF(lin_fpd, 'êd', 'trimf', input_triang_mfs(1,:), 'Name', 'NEG');
lin_fpd = addMF(lin_fpd, 'êd', 'trimf', input_triang_mfs(2,:), 'Name', 'NULA');
lin_fpd = addMF(lin_fpd, 'êd', 'trimf', input_triang_mfs(3,:), 'Name', 'POZ');

f = figure(2); plotmf(lin_fpd,'input',2); ylabel('Stepen pripadnosti');
f.Renderer = 'painters';
f.Name = 'dif_error_membership_fPID';
set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');


f = figure(3);
num_of_singletons = max(size(output_singleton_mfs));
stem(output_singleton_mfs, ones(num_of_singletons,1));
for i=1:num_of_singletons
    text(output_singleton_mfs(i),1.1,lin_fpd.output.MembershipFunctions(i).Name);
end

xlabel('û'); ylabel('Stepen pripadnosti'); a=axis; a(4)=1.2; axis(a);
xlim([-3, 3]);
f.Renderer = 'painters';
f.Name = 'out_fPID';set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end



ruleList = [];
for i = 1: number_of_input_mfs
  for j = 1: number_of_input_mfs
    output_val = input_triang_mfs(i, 2) + input_triang_mfs(j, 2);
    k = (output_val - min_out_range)/step + 1;
    current_rule = [i j k 1 1];
    ruleList = [ruleList; current_rule];
  end
end
lin_fpd = addRule(lin_fpd,ruleList);

f = figure(4); gensurf(lin_fpd);
writefis(lin_fpd,'pravila');

f.Renderer = 'painters';
f.Name = 'surf_fPID';set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end



%%

t_end = r_time + 0.1;
sigma_arr = [0.45, 0.5, 0.6];
sel = 3;
sim('.\model\feedback_lin_2018.slx');
indexOfInterest = (t < r_time + 0.046) & (t > r_time-0.002);
time = t(indexOfInterest);

names = {'linearan','$\sigma = 0.45$', '$\sigma = 0.50$', '$\sigma = 0.60$'};

f = figure(205);
    f.Name = 'Reference_Comparison_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
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
        




for i = 1 : length(sigma_arr)
    sigma = sigma_arr(i);


    standard_input_range = [-1 1];
    number_of_inputs = 2;
    number_of_input_mfs = 3;
    min_ir = min(standard_input_range);
    max_ir = max(standard_input_range);
    input_triang_mfs = [sigma -1 0; sigma 0 1; sigma 1 0];
    step = (max_ir - min_ir) / (number_of_input_mfs - 1);

    min_out_range = number_of_inputs * min_ir;
    max_out_range = number_of_inputs * max_ir;
    standard_output_range = [min_out_range, max_out_range];
    output_singleton_mfs = [min_out_range: step: max_out_range]';


    lin_fpd = sugfis('Name','pravila_modif');
    lin_fpd.andMethod    = 'prod';
    lin_fpd.orMethod     = 'probor';
    lin_fpd.impMethod    = 'prod';
    lin_fpd.aggMethod    = 'sum';
    lin_fpd.defuzzMethod = 'wtaver';

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');


    lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');




    lin_fpd = addRule(lin_fpd,ruleList);


    writeFIS(lin_fpd,'pravila_modif');











    sel = 4;
    sim('.\model\feedback_lin_2018.slx');

indexOfInterest = (t < r_time + 0.046) & (t > r_time-0.002);


f = figure(205);
    f.Name = 'sigma_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([time(1) time(end)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');
        hold off;
        
    subplot(2, 2, 3);
    hold on;
    plot(t(indexOfInterest), u_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([time(1) time(end)]);
        hold off;
        
    subplot(2, 2, 2);  
    hold on;
        plot(t(indexOfInterest), iL_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([time(1) time(end)]);
        hold off;
    subplot(2, 2, 4);  
        hold on;
        plot(iL_out(indexOfInterest), y_out(indexOfInterest));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
        hold off;




end

f = figure(205);
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');

    subplot(2, 2, 1);
    hold on;
    p = plot(t(indexOfInterest), total_ref(indexOfInterest), 'k--');
    
    legend(p, '$r(t)$', 'Location', 'best');
    legend('boxoff');
    set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%%

bounds_arr = [0.35, 0.4, 0.5];
t_end = r_time + 0.1;
sigma = 0.45;

    standard_input_range = [-1 1];
    number_of_inputs = 2;
    number_of_input_mfs = 3;
    min_ir = min(standard_input_range);
    max_ir = max(standard_input_range);
    input_triang_mfs = [sigma -1 0; sigma 0 1; sigma 1 0];
    step = (max_ir - min_ir) / (number_of_input_mfs - 1);

    min_out_range = number_of_inputs * min_ir;
    max_out_range = number_of_inputs * max_ir;
    standard_output_range = [min_out_range, max_out_range];
    output_singleton_mfs = [min_out_range: step: max_out_range]';


    lin_fpd = sugfis('Name','pravila_modif');
    lin_fpd.andMethod    = 'prod';
    lin_fpd.orMethod     = 'probor';
    lin_fpd.impMethod    = 'prod';
    lin_fpd.aggMethod    = 'sum';
    lin_fpd.defuzzMethod = 'wtaver';

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');


    lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');




    lin_fpd = addRule(lin_fpd,ruleList);


    writeFIS(lin_fpd,'pravila_modif');


sel = 4;
sim('.\model\feedback_lin_2018.slx');
indexOfInterest = (t < r_time + 0.046) & (t > r_time-0.002);
time = t(indexOfInterest);

names = {'$\sigma_z = 0.45$','$\sigma_z = 0.35$', '$\sigma_z = 0.4$', '$\sigma_z = 0.5$'};

f = figure(2051);
    f.Name = 'Reference_a_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
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
        




for i = 1 : length(sigma_arr)

    bounds = bounds_arr(i);
    standard_input_range = [-1 1];
    number_of_inputs = 2;
    number_of_input_mfs = 3;
    min_ir = min(standard_input_range);
    max_ir = max(standard_input_range);
    input_triang_mfs = [sigma -1 0; bounds 0 0; sigma 1 0];
    step = (max_ir - min_ir) / (number_of_input_mfs - 1);

    min_out_range = number_of_inputs * min_ir;
    max_out_range = number_of_inputs * max_ir;
    standard_output_range = [min_out_range, max_out_range];
    output_singleton_mfs = [min_out_range: step: max_out_range]';


    lin_fpd = sugfis('Name','pravila_modif');
    lin_fpd.andMethod    = 'prod';
    lin_fpd.orMethod     = 'probor';
    lin_fpd.impMethod    = 'prod';
    lin_fpd.aggMethod    = 'sum';
    lin_fpd.defuzzMethod = 'wtaver';

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');


    lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');




    lin_fpd = addRule(lin_fpd,ruleList);


    writeFIS(lin_fpd,'pravila_modif');











    sel = 4;
    sim('.\model\feedback_lin_2018.slx');

indexOfInterest = (t < r_time + 0.046) & (t > r_time-0.002);


f = figure(2051);
    f.Name = 'Reference_a_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
        grid on;
        xlabel('$t$ [s]');
        ylabel('$y(t)$ [V]');
        xlim([time(1) time(end)]);
        %legend(p, '$r(t)$', 'Location', 'best');
        %legend('boxoff');
        hold off;
        
    subplot(2, 2, 3);
    hold on;
    plot(t(indexOfInterest), u_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        xlim([time(1) time(end)]);
        hold off;
        
    subplot(2, 2, 2);  
    hold on;
        plot(t(indexOfInterest), iL_out(indexOfInterest));
        grid on;
        xlabel('$t$ [s]');
        ylabel('$i_L(t)$ [A]');
        xlim([time(1) time(end)]);
        hold off;
    subplot(2, 2, 4);  
        hold on;
        plot(iL_out(indexOfInterest), y_out(indexOfInterest));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
        hold off;




end

f = figure(2051);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    hold off;
    subplot(2, 2, 1);
    hold on;
    p = plot(t(indexOfInterest), total_ref(indexOfInterest), 'k--');
    
    legend(p, '$r(t)$', 'Location', 'best');
    legend('boxoff');
    set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

sigma_z = 0.35;
%%
sel = 4;
t_end = t_switch + 3.1*r_delta;


standard_input_range = [-1 1];
    number_of_inputs = 2;
    number_of_input_mfs = 3;
    min_ir = min(standard_input_range);
    max_ir = max(standard_input_range);
    input_triang_mfs = [sigma -1 0; sigma_z 0 0; sigma 1 0];
    step = (max_ir - min_ir) / (number_of_input_mfs - 1);

    min_out_range = number_of_inputs * min_ir;
    max_out_range = number_of_inputs * max_ir;
    standard_output_range = [min_out_range, max_out_range];
    output_singleton_mfs = [min_out_range: step: max_out_range]';


    lin_fpd = sugfis('Name','pravila_modif');
    lin_fpd.andMethod    = 'prod';
    lin_fpd.orMethod     = 'probor';
    lin_fpd.impMethod    = 'prod';
    lin_fpd.aggMethod    = 'sum';
    lin_fpd.defuzzMethod = 'wtaver';

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');


    lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
    lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');




    lin_fpd = addRule(lin_fpd,ruleList);


    writeFIS(lin_fpd,'pravila_modif');



sim('.\model\feedback_lin_2018.slx');



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




indexOfInterest = (t < r_time + 0.08) & (t > r_time-0.002);
time = t(indexOfInterest);

names = {'linearna pravila', 'nelinearna pravila'};

f = figure(206);
        f.Name = 'Reference_rules_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
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
        

    standard_input_range = [-1 1];
    number_of_inputs = 2;
    number_of_input_mfs = 3;
    min_ir = min(standard_input_range);
    max_ir = max(standard_input_range);
    input_triang_mfs = [sigma -1 0; sigma_z 0 0; sigma 1 0];
    step = (max_ir - min_ir) / (number_of_input_mfs - 1);

min_out_range = number_of_inputs * min_ir;
max_out_range = number_of_inputs * max_ir;
standard_output_range = [min_out_range, max_out_range];
output_singleton_mfs = [min_out_range: step: max_out_range]';


lin_fpd = sugfis('Name','pravila');
lin_fpd.andMethod    = 'min';
lin_fpd.orMethod     = 'max';
lin_fpd.impMethod    = 'prod';
lin_fpd.aggMethod    = 'sum';
lin_fpd.defuzzMethod = 'wtsum';

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'ê');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'ê', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');



f = figure(10); plotmf(lin_fpd,'input',1); ylabel('Stepen pripadnosti');

f.Renderer = 'painters';
f.Name = 'error_membership_fPID_modif';

set(findall(gcf,'-property','FontSize'),'FontSize', 8)


if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

    lin_fpd = addInput(lin_fpd, standard_input_range, 'Name', 'êd');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(1,:), 'Name', 'NEG');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(2,:), 'Name', 'NULA');
    lin_fpd = addMF(lin_fpd, 'êd', 'gaussmf', input_triang_mfs(3,:), 'Name', 'POZ');


f = figure(20); plotmf(lin_fpd,'input',2); ylabel('Stepen pripadnosti');
f.Renderer = 'painters';
f.Name = 'dif_error_membership_fPID_modif';
set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


lin_fpd = addOutput(lin_fpd, standard_output_range, 'Name', 'û');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(1,:), 'Name', 'NEG VEL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(2,:), 'Name', 'NEG MAL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(3,:), 'Name', 'NULA');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(4,:), 'Name', 'POZ MAL');
lin_fpd = addMF(lin_fpd, 'û', 'constant', output_singleton_mfs(5,:), 'Name', 'POZ VEL');


f = figure(30);
num_of_singletons = max(size(output_singleton_mfs));
stem(output_singleton_mfs, ones(num_of_singletons,1));
for i=1:num_of_singletons
    text(output_singleton_mfs(i),1.1,lin_fpd.output.MembershipFunctions(i).Name);
end

xlabel('û'); ylabel('Stepen pripadnosti'); a=axis; a(4)=1.2; axis(a);
xlim([-3, 3]);
f.Renderer = 'painters';
f.Name = 'out_fPID_modif';set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end



ruleList = [];
for i = 1: number_of_input_mfs
  for j = 1: number_of_input_mfs
    output_val = input_triang_mfs(i, 2) + input_triang_mfs(j, 2);
    k = (output_val - min_out_range)/step + 1;
    current_rule = [i j k 1 1];
    ruleList = [ruleList; current_rule];
  end
end
lin_fpd = addRule(lin_fpd,ruleList);

f = figure(40); gensurf(lin_fpd);
writefis(lin_fpd,'pravila_modif');

f.Renderer = 'painters';
f.Name = 'surf_fPID_modif';set(findall(gcf,'-property','FontSize'),'FontSize', 8)

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


sim('.\model\feedback_lin_2018.slx');



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





f = figure(206);
    f.Name = 'Reference_rules_K_fuzzy';
    f.Renderer = 'painters';
    figure(f);
    subplot(2, 2, 1);
    hold on;
    plot(t(indexOfInterest), y_out(indexOfInterest));
    hold on;
%   zoom = [zoom y_out(indexOfInterest)];
    %p = plot(t(N_low :  N_high), total_ref(N_low :  N_high), 'k--');
    
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

f = figure(206);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    hold off;

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end

%%

Ki = 0.81*Ki;
Kd = 1.1*Kd;
Kp = 0.95*Kp;


close all;
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

sel = 4;


%%

sim('.\model\feedback_lin_2018.slx');


N_low = find(t >=0, 1);
N_high = find(t >= t_switch - 0.01, 1); 

f = figure(201);
f.Renderer = 'painters';
f.Name = 'Transient_K_fuzzy';
figure(f);
    subplot(2, 1, 1);
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
    subplot(2, 1, 2);
    plot(t(N_low : N_high), u_out(N_low : N_high), 'black');
        grid on;
        xlabel('$t$ [s]');
        ylabel('$u(t)$');
        ylim([0, 1]);
        xlim([t(N_low) t(N_high)]);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 1.95*r_delta, 1); 


f = figure(201);
f.Renderer = 'painters';
f.Name = 'Reference_Response_K_fuzzy';
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
f = figure(202);
f.Renderer = 'painters';
f.Name = 'Disturbance_Response_K_fuzzy';
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
 
 
sim('.\model\feedback_lin_2018.slx');


N_low = find(t >=t_switch*0.95, 1);
N_high = find(t >=t_switch + 1.95*r_delta, 1); 
f = figure(203);
f.Name = 'Reference_+_Noise_Response_K_fuzzy';
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

f = figure(204);
f.Name = 'Disturbance_+_Noise_Response_K_fuzzy';
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

f = figure(204);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
        
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
        
%% Robustnost:
t_switch = 0.4;
%load('parameters')
L_r = [L*0.8, L, L*1.2];
n_on = 0;

d_on = 1;
n_on = 0;

names = {'$0.8L$','$L$', '$1.2L$'};

for i = 1 : length(L_r)
    L = L_r(i);

    sim('.\model\feedback_lin_2018.slx');


    N_low = find(t >=t_switch*0.95, 1);
    N_high = find(t >=t_switch + 1.98*r_delta, 1);  
    
    f = figure(205);
    f.Name = 'Reference_Response_+_Robust_K_fuzzy';
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

    
    f = figure(206);
    f.Renderer = 'painters';
    f.Name = 'Disturbance_Response_+_Robust_K_fuzzy';
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

   f = figure(2050);
    time = t(indexOfInterest);
    f.Name = 'Reference_Response_+_Robust_K_fuzzy_zoom';
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
        plot(iL_out(indexOfInterest), y_out(indexOfInterest));
        grid on;
        xlabel('$i_L(t)$ [A]');
        ylabel('$u_C(t)$ [V]');
        
    N_low = find(t >=d_time*0.98, 1);
    N_high = find(t >=d_time + 1.98*d_delta, 1); 

    indexOfInterest = (t < d_time + 0.09) & (t > d_time-0.002); % range of t near perturbation
    time = t(indexOfInterest);
    f = figure(2060);
    f.Renderer = 'painters';
    f.Name = 'Disturbance_Response_+_Robust_K_fuzzy_zoom';
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
f = figure(205);
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
    
f = figure(2050);
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
    f = figure(205);
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
f = figure(206);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
    
f = figure(2060);
set(findall(gcf,'-property','FontSize'),'FontSize', 8)
f.Renderer = 'painters';
    subplot(2, 2, 3); 
    legend(names);
    legend('boxoff');
    legend('Location', 'best');
if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
    f = figure(206);
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
    
%% Pocetni uslov poremecaj

d_on = 0;
r_on = 0;

t_switch = 0;

t_end = 0.3;
t_step = 0.0001;


L = L_r(3);

iL_r = -10:5:10;
uC_r = -10:5:10;

f = figure(207);
f.Name = 'State_Space_K_fuzzy';
hold all;

for i = 1:length(iL_r)
    for j = 1:length(uC_r)
        iL0 = iL_r(i);
        uC0 = uC_r(j);
        sim('.\model\feedback_lin_2018.slx');
        f = figure(207);
        plot(iL_out, y_out);

    end
end

xlim([-21 36]);
ylim([-220 110]);
xlabel('$i_L(t)$');
ylabel('$u_C(t)$'); 
grid('on');

if(SAVE)
    saveas(f,[path '\' f.Name '.eps'],'epsc');
end
%
