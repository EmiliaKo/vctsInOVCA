%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to simulate stratified virtual clinical trials simulations   %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PFI, TumorCompositionAfterTreatment] = SimulateStratifiedCTS(nPatients,drug,arm)
[b,d,u,M_logmean,M_logsigma,d_chemotherapy_mean,d_chemotherapy_std,alpha1,alpha2,beta] = loadParameters();

% == set outputs 
PFI = zeros(nPatients,1);
TumorCompositionAfterTreatment = zeros(nPatients,8);

% == simulate virtual patient cohort
i = 1;
while true 
    disp(i)
    % sample from lognormal distribution two parameters: M & d_chemotherapy
    M              = lognrnd(M_logmean,M_logsigma);  
    d_chemotherapy = lognrnd(d_chemotherapy_mean,d_chemotherapy_std); 
    
    % simulate pre-treatment phase
    [CellsDiagnosis]= Drugresmodel_mex([1 0 0 0 0 0 0 0],10e9,b,d,u,M, 0,0,0,0,0,0,0);
    
    % identify dominant platinum resistance mechanism
    PRE  = sum([CellsDiagnosis(2) CellsDiagnosis(5) CellsDiagnosis(6) CellsDiagnosis(8)]);
    ON   = sum([CellsDiagnosis(3) CellsDiagnosis(5) CellsDiagnosis(7) CellsDiagnosis(8)]);
    POST = sum([CellsDiagnosis(4) CellsDiagnosis(6) CellsDiagnosis(7) CellsDiagnosis(8)]);
    
    [~,idx] = sort([PRE ON POST]);
    dominantResMech = idx(end);
    
    % set treatment parameters to a given stratified VCTS
    if strcmp(arm,'positive') % biomarker positive clinical trial arm
        if dominantResMech ~= drug
            continue     
        end
        switch drug
            case 1
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('trientine');
            case 2 
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('wee1i');
            case 3
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('birinapant');
        end
    else % biomarker negative clinical trial arm
        if dominantResMech == drug
            continue     
        end
        switch drug
            case 1
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('trientine');
            case 2 
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('wee1i');
            case 3
                [alpha1_pre,alpha2_pre,u_ddr,d_a] = loadTargetedTreatmentParameters('birinapant');
        end        
    end       
    %primary treatment
    [CellsNACT]= Drugresmodel_mex(CellsDiagnosis,3*7,b,d,u,M, d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a);
    CellsSurgery = (1-beta).*sum(CellsNACT).*CellsNACT./sum(CellsNACT);    
    [CellsADJ]= Drugresmodel_mex(CellsSurgery,3*7,b,d,u,M, d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a);
    %post-treatment
    if sum(CellsADJ) > 0 % if residual tumor after treatment
        [~, pfi]= Drugresmodel_mex(CellsADJ,10e9,b,d,u,10^9, 0,0,0,0,0,0,0);
    else
        pfi = 10*365;
    end
    
    if pfi > 10*365
        pfi = 10*365;
    end
    
    %outputs
    PFI(i) = pfi/30;
    TumorCompositionAfterTreatment(i,:) = CellsADJ;
     
    if i >= nPatients
       break; 
    end
    
    i = i +1;
end

end