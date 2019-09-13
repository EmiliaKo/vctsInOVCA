%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to simulate virtual HGSOC patient cohort                     %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = SimulateVirtualPatientsCohort(fileOut,b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% == inputs
% fileOut <- MAT file name where the results will be saved
% b,d,u,beta,M_logmean,M_logsigma,alpha1,alpha2,d_chemotherapy_mean,d_chemotherapy_std,alpha1_pre,alpha2_pre,u_ddr,d_a -> model parameters (see loadParameters function)
%% outputs
% out <- struct with results from simulations: PFI, tumor compositionat
% diagnosis/after primary treatment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare rng
rng('shuffle','twister');
%% set cohort parameters & initial conditions (start from a single sensitive cell)
nPatients = 1000; % cohort size
Cells_init = [1   % 000
              0   % 100
              0   % 010
              0   % 001
              0   % 110
              0   % 101
              0   % 011
              0]';% 111
                    
%% set outputs          
CellsDiag = zeros(nPatients,8);  % number of cells at diagnosis 
PFI = zeros(nPatients,1);        % platinum-free interval [months]
CellsAfter = zeros(nPatients,8); % number of cells after treatment

%% == iterate over patients
for i = 1:nPatients
    disp(i)
    % sample M & d_chemotherapy from log-normal distribution
    M              = lognrnd(M_logmean,M_logsigma);  
    d_chemotherapy = lognrnd(d_chemotherapy_mean,d_chemotherapy_std); 
%% == simulations start here
    %1)  pre-treatment phase
    t_stop = 10e10;
    [CellsDiagnosis]= Drugresmodel_mex(Cells_init,t_stop,b,d,u,M, 0,0,0,0,0,0,0);
    %2) primary treatment
    t_chemotherapy = 21; % time interval between two chemotherapy cycles
    cycles_chemotherapy = 3;% # of chemotherapy cycles
    CellsNACT    = Drugresmodel_mex(CellsDiagnosis,t_chemotherapy*cycles_chemotherapy,b,d,u,M, d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a);
    CellsSurgery = (1-beta).*sum(CellsNACT).*CellsNACT./sum(CellsNACT);
    CellsADJ     = Drugresmodel_mex(CellsSurgery,t_chemotherapy*cycles_chemotherapy,b,d,u,M, d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a);
    %3) post-treatment 
    if sum(CellsADJ) > 0
        M_relapse = 10^9; % tumor burden at relapse
        [~, pfi]= Drugresmodel_mex(CellsADJ,t_stop,b,d,u,M_relapse, 0,0,0,0,0,0,0);
    else
        pfi = 10*365; % PFI equals 10 years
    end
%% == simulations ends here    
    % outputs
    PFI(i,:)        = pfi/30;
    CellsDiag(i,:)  = CellsDiagnosis;
    CellsAfter(i,:) = CellsADJ;
end
% save outputs
out.pfi = PFI;
out.cellsAtDiagnosis = CellsDiag;
out.cellAfterTreatment = CellsAfter;
save(fileOut)
end