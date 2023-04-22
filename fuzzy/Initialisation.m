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
path = 'C:\Users\milos\OneDrive\VII semestar\NSU2\projekat\izvestaj\slike\projekat3\zad1';

uC0 = -22.5; %V
L = 15.91e-3; %H
C = 470e-6; %F
R = 52; %Omh
E = 12; %V


save('parameters', 'R', 'L', 'E', 'C', 'uC0');

load('linear_model')
load('feedback_linearisation_params')
K_r = 0;

y0 = uC0;

%Parameter initialisation

%%
sel = 1;

Simulation_fuzzy;






