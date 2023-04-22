%%
clear all;
close all;


% This script changes all interpreters from tex to latex. 

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));

for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%%

SAVE = 0; 
path = 'C:\Users\milos\OneDrive\VII semestar\NSU2\projekat\projekat1\izvestaj\slike\zad3';


uC0 = -22.5; %V
L = 15.91e-3; %H
C = 470e-6; %F
R = 52; %Omh
E = 12; %V

save('parameters', 'R', 'L', 'E', 'C', 'uC0');

load('linear_model')
y0 = uC0;
load('controler')

%Parameter initialisation

font_size = 9;

%%

loop_shaping;

%%
sel = 1;
Simulation_loop;
%%
% sel = 1;
% Simulation_PID_mcl_loop;
%%
sel = 2;
Simulation_PID_loop;
