%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to load the model parameters                                 %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,d,u,M_logmean,M_logsigma,M,d_chemotherapy_mean,d_chemotherapy_std,d_chemotherapy,alpha1,alpha2,beta,t_chemo] = loadParameters()
%== growth dynamics parameters
b = 0.667;      % birth rate
d = 0.661;      % death rate
u  = 10^-5;     % transition rate
%== tumor burden at diagnosis 
M_logmean  = 26.5934; % parameter of log-normal prob. dist. fun.: Log mean
M_logsigma = 0.4709;  % parameter of log-normal prob. dist. fun.: Log standard deviation
M = lognstat(M_logmean,M_logsigma);
%== treatment parameters
d_chemotherapy_mean = -0.9; % (chemo induced death rate) par. of log-normal prob. dist. fun.: Log mean
d_chemotherapy_std = 1.2208;% (chemo induced death rate) par. of log-normal prob. dist. fun.: Log mean
d_chemotherapy = lognstat(d_chemotherapy_mean,d_chemotherapy_std);
alpha1 = 0.02;              % weight of drug effect on partially resistant cells (one res. mechanism)
alpha2 = 0.01;              % weight of drug effect on partially resistant cells (two res. mechanism)
beta = .99;                 % reduction of tumor burden by IDS
t_chemo = 3*21;             % time of three cycles of chemotherapy

end