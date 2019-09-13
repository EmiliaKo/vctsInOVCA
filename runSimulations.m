%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main script to simulate Virtual clinical trials simulations           %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load model parameters
[b,d,u,M_logmean,M_logsigma,M,d_chemotherapy_mean,d_chemotherapy_std,d_chemotherapy,alpha1,alpha2,beta,t_chemo] = loadParameters();
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('none');
% compile mex file
t_stop = 20*12*30; % max simulation time
Cells_init = [1 0 0 0 0 0 0 0];
codegen Drugresmodel.m -args {Cells_init,t_stop,b,d,u,M,d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a}

%% run Standard-of-care simulation and show the model fit to the data
out = SimulateVirtualPatientsCohort('soc.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
plotKaplanMeier('SOC_run.mat')

%% perform virtual randomized controlled trials simulations
% simulate SOC + trientine
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('trientine');
t = SimulateVirtualPatientsCohort('trientine.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + wee1i
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('wee1i');
w = SimulateVirtualPatientsCohort('wee1i.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + birinapant
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('birinapant');
b = SimulateVirtualPatientsCohort('birinapant.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + trientine + wee1i
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('t_w');
tw = SimulateVirtualPatientsCohort('trientine_wee1i.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + trientine + birinapant
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('t_b');
tb = SimulateVirtualPatientsCohort('trientine_birinapant.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + wee1i + trientine
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('w_b');
wb = SimulateVirtualPatientsCohort('wee1i_birinapant.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);
% simulate SOC + wee1i + trientine + birinapant
[alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('t_w_b');
twb = SimulateVirtualPatientsCohort('trientine_wee1i_birinapant.mat',b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a);

%% perform stratified virtual clinical trials simulations
nPatients = 1000;
[PFI_trientine, Comp_trientine]          = SimulateStratifiedCTS(nPatients,1,'positive');
[PFI_trientine_negative, Comp_trientine_negative] = SimulateStratifiedCTS(nPatients,1,'negative');

[PFI_wee1i, Comp_wee1i]          = SimulateStratifiedCTS(nPatients,2,'positive');
[PFI_wee1i_negative, Comp_wee1i_negative] = SimulateStratifiedCTS(nPatients,2,'negative');

[PFI_birinapant, Comp_birinapant]          = SimulateStratifiedCTS(nPatients,3,'positive');
[PFI_birinapant_negative, Comp_birinapant_negative] = SimulateStratifiedCTS(nPatients,3,'negative');


