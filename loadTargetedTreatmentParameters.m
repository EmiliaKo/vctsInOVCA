%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function load parameters related to targeted treatment                %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters(drug)
switch drug
    case 'trientine'
        alpha1_pre = 0.02+0.13;
        alpha2_pre = 0.01+0.13;
        u_ddr = 0;
        d_a = 0;
    case 'wee1i'
        alpha1_pre = 0.02;
        alpha2_pre = 0.01;
        u_ddr = 0.07;
        d_a = 0;        
    case 'birinapant'
        alpha1_pre = 0.02;
        alpha2_pre = 0.01;
        u_ddr = 0;
        d_a = 0.07;
    case 'none'
        alpha1_pre = 0.02;
        alpha2_pre = 0.01;
        u_ddr = 0;
        d_a = 0;      
    case 't_w'
        alpha1_pre = 0.02+0.13;
        alpha2_pre = 0.01+0.13;
        u_ddr = 0.07;
        d_a = 0;
    case 't_b'
        alpha1_pre = 0.02+0.13;
        alpha2_pre = 0.01+0.13;
        u_ddr = 0;
        d_a = 0.07;       
    case 'w_b'
        alpha1_pre = 0.02;
        alpha2_pre = 0.01;
        u_ddr = 0.07;
        d_a = 0;  
    case 't_w_b'
        alpha1_pre = 0.02;
        alpha2_pre = 0.01;
        u_ddr = 0.07;
        d_a = 0.07;      
end