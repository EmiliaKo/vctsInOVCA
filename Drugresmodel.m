%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Core function which is responsible for the model simulations          %%
%% author: Emilia Kozlowska                                              %%
%% the last update 13/09/2019                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cells, time]= Drugresmodel(Cells_init,t_stop,b,d,u,M,d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% == inputs
% Cells_init <- vector of size 8 indicating # of cells for each subclone at
% start of simulation
% t_stop <- maximum simulation time
% b,d,u,M,d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a -> model parameters (see loadParameters function)
%% outputs
% Cells <- vector of size 8 indicating # of cells for each subclone at the
% end of simulation
% time <- time at the end of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% == set transition matrix
[transition_slow, transition_fast] = transitions();
transition = [transition_fast; transition_slow];    
                   
%% == set parameters
Cells = Cells_init;
time = 0;

%% simulation start here
while(true)
    %== calculate propensity
    [propensity_slow,propensity_fast] = propensity_calculator(Cells,b,d,u,d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a); 
    propensity= [propensity_fast; propensity_slow];      
    
    if(sum(Cells) < 10^7) % if tumor size is below 10^7 cells
        % time to next slow event
        Tau = exprnd(1/(sum(propensity_slow))); 
        % which next slow event will occur?
        normProb = cumsum(propensity_slow./sum(propensity_slow));
        next_reaction_ind = find(normProb>rand(),1,'first');
        
        % update slow reaction
        if next_reaction_ind ~= 0
        Cells = Cells + transition_slow(next_reaction_ind,:);
        else
            Tau = .1;
        end
        % update fast reactions
        lambda=propensity_fast.*Tau;
        R = transition_fast;
        for idx = 1:size(transition_fast,1)
            R(idx,:) = transition_fast(idx,:).*poissrnd(lambda(idx));
        end
        Cells = Cells + sum(R,1);
        Cells(Cells<0) = 0; % to avoid negative # of cells  
        time = time + Tau;  % update time
        
    else % if tumor size is above 10^5 cells (use tau-leaping)
        Tau = .1;
        lambda=propensity.*Tau;
        R = transition;
        for idx = 1:size(transition,1)
            R(idx,:) = transition(idx,:).*poissrnd(lambda(idx));
        end
        Cells = Cells + sum(R,1);
        Cells(Cells<0) = 0; % to avoid negative # of cells  
        time = time + Tau; % update time
    end

    if (sum(Cells,2)) <= 0 || isnan(sum(Cells,2)) % start from the beginning
            Cells = Cells_init;
            time = 0; 
    end
      
    if  sum(Cells,2) > M || time >t_stop % when to stop simulations?
        break;    
    end

end
%% simulation ends here

end

%% helper function: function to calcule events propensity
function [propensity_slow,propensity_fast] = propensity_calculator(Cells, b,d,u,d_chemotherapy,alpha1,alpha2,alpha1_pre,alpha2_pre,u_ddr,d_a)
    propensity_slow = [b*u*Cells(1) 
                       b*u*Cells(1) 
                       b*u*Cells(1) 
                       b*u*Cells(2)
                       b*u*Cells(2)
                       b*u*Cells(3)
                       b*u*Cells(3)                          
                       b*u*Cells(4)
                       b*u*Cells(4)                     
                       b*u*Cells(5)                          
                       b*u*Cells(6)
                       b*u*Cells(7)];  
      
    propensity_fast = [b*(1-3*u)*Cells(1)
                       b*(1-2*u)*Cells(2)
                       b*(1-2*u)*Cells(3)
                       b*(1-2*u)*Cells(4)
                       b*(1-u)*Cells(5)
                       b*(1-u)*Cells(6)
                       b*(1-u)*Cells(7)
                       b*Cells(8)                  
                       d*Cells(1)
                       d*Cells(2)
                       d*Cells(3)
                       d*Cells(4)
                       d*Cells(5)
                       d*Cells(6)
                       d*Cells(7)
                       d*Cells(8)                  
                       d_chemotherapy*Cells(1)
                       d_chemotherapy*alpha1_pre*Cells(2)
                       d_chemotherapy*alpha1*Cells(3)
                       d_chemotherapy*alpha1*Cells(4)
                       d_chemotherapy*alpha2_pre*Cells(5)
                       d_chemotherapy*alpha2_pre*Cells(6)
                       d_chemotherapy*alpha2*Cells(7)
                       d_a*Cells(8)
                       d_a*Cells(7)  
                       d_a*Cells(6)                       
                       d_a*Cells(4)
                       u_ddr*Cells(8)
                       u_ddr*Cells(8)
                       u_ddr*Cells(8)
                       u_ddr*Cells(7)
                       u_ddr*Cells(7)
                       u_ddr*Cells(6)
                       u_ddr*Cells(6)                      
                       u_ddr*Cells(5)
                       u_ddr*Cells(5)                        
                       u_ddr*Cells(3)                      
                       u_ddr*Cells(2)
                       u_ddr*Cells(4)                        
                       ];
end

%% helper function: define transitions in all events
function [transition_slow, transition_fast] = transitions()
% division with mutation is a slow process
    transition_slow = [0 1 0 0 0 0 0 0  % 000 -> 000 + 100
                       0 0 1 0 0 0 0 0  % 000 -> 000 + 010
                       0 0 0 1 0 0 0 0  % 000 -> 000 + 001
                       0 0 0 0 1 0 0 0  % 100 -> 100 + 110
                       0 0 0 0 0 1 0 0  % 100 -> 100 + 101
                       0 0 0 0 1 0 0 0  % 010 -> 000 + 100
                       0 0 0 0 0 0 1 0  % 010 -> 000 + 100
                       0 0 0 0 0 1 0 0  % 001 -> 000 + 100
                       0 0 0 0 0 0 1 0  % 001 -> 000 + 100
                       0 0 0 0 0 0 0 1  % 110 -> 110 + 111
                       0 0 0 0 0 0 0 1  % 011 -> 011 + 111
                       0 0 0 0 0 0 0 1];% 110 -> 110 + 111    
% other reactions are fast                       
 transition_fast =    [1 0 0 0 0 0 0 0   % 000 -> 000 + 000
                       0 1 0 0 0 0 0 0   % 100 -> 100 + 100
                       0 0 1 0 0 0 0 0   % 010 -> 010 + 010
                       0 0 0 1 0 0 0 0   % 001 -> 001 + 001
                       0 0 0 0 1 0 0 0   % 110 -> 110 + 110 
                       0 0 0 0 0 1 0 0   % 101 -> 101 + 101
                       0 0 0 0 0 0 1 0   % 110 -> 110 + 110
                       0 0 0 0 0 0 0 1   % 111 -> 111 + 111
                       -1 0 0 0 0 0 0 0  % 000 -> 
                       0 -1 0 0 0 0 0 0  % 100 -> 
                       0 0 -1 0 0 0 0 0  % 010 -> 
                       0 0 0 -1 0 0 0 0  % 001 -> 
                       0 0 0 0 -1 0 0 0  % 110 -> 
                       0 0 0 0 0 -1 0 0  % 101 -> 
                       0 0 0 0 0 0 -1 0  % 011 -> 
                       0 0 0 0 0 0 0 -1  % 111 ->                
                       -1 0 0 0 0 0 0 0  % 000 -> 
                       0 -1 0 0 0 0 0 0  % 100 -> 
                       0 0 -1 0 0 0 0 0  % 010 -> 
                       0 0 0 -1 0 0 0 0  % 001 -> 
                       0 0 0 0 -1 0 0 0  % 110 -> 
                       0 0 0 0 0 -1 0 0  % 101 -> 
                       0 0 0 0 0 0 -1 0  % 011 -> 
                       0 0 0 0 0 0 0 -1  % 111 ->   
                       0 0 0 0 0 0 -1 0  % 011 ->
                       0 0 0 0 0 -1 0 0  % 101 ->
                       0 0 0 -1 0 0 0 0  % 001 ->
                       0 0 0 0 0 1 0 -1  % 111 -> 101
                       0 0 0 0 0 0 1 -1  % 111 -> 011
                       0 0 0 0 1 0 0 -1  % 111 -> 110
                       0 0 0 1 0 -1 0 0  % 011 -> 001 
                       0 0 1 0 0 -1 0 0  % 011 -> 010 
                       0 0 0 1 0 -1 0 0  % 101 -> 001 
                       0 1 0 0 0 -1 0 0  % 101 -> 100 
                       0 1 0 0 -1 0 0 0  % 110 -> 100 
                       0 0 1 0 -1 0 0 0  % 110 -> 010 
                       1 0 -1 0 0 0 0 0  % 010 -> 000                    
                       1 -1 0 0 0 0 0 0  % 100 -> 000                    
                       1 0 0 -1 0 0 0 0];% 001 -> 000                    
                       
end

