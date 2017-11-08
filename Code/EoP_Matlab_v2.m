%% EOP

% The following code simulates the game from the paper,
% and produces the figures used.



clear all; clc; close all;

cd 'C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab' % Set current directory

%To keep the code manageable, all variable names will follow a three-part
%structure: 
% 1) The first part of each will be either pop (population, agt (agent),
% prn (principal) or oth (other)
%) 2) The second part will be a common denominater, e.g. mean, var
%(variance), est (estimate) or similar. 
% 3) Finally a fitting name that allows one to understand what it is
%    (e.g. the letters used in the text)
% 4) If relevant, there can be a fourth part with a _t. 

%%
% It has the following section
% GLOBAL PARAMETERS
% Preallocate matrices for storage 
% SIMULATE GAME
% FIGURES
% CLEANUP


%% GLOBAL PARAMETERS
%Randomization
seed = 52; 
    rng(seed);

%Number of iterations of simulations
T = 100; %Game runs for T periods 
N = 1000; %We simulate a population of N. 

%Population means
pop_mean_thetabar_0         = unifrnd(5,10);
    pop_mean_etabar         = unifrnd(0.01,0.02);
pop_mean_vbar               = unifrnd(2,3);
pop_mean_mubar              = unifrnd(1,2); %Large enough that a high share has positive motivation

%Population variances 
pop_var_sigmatheta_0          = unifrnd(5,8);
    pop_var_sigmaeta         = unifrnd(0.1, 0.2);
pop_var_s2v               = unifrnd(2,3);
pop_var_sigmamu              = unifrnd(2,3);
pop_var_s2mu                 = unifrnd(2,3);
pop_var_s2theta              = unifrnd(15,20);

%Population constants
pop_value_w                  = unifrnd(1,2); %Must be lower than pop_mean_vbar
pr_value_lambda              = 0.5;
pr_value_alpha               = 0.8; 
pr_value_kp                  = 1; 
pr_value_b                   = 1;
prn_value_x_constant         = 1; 

%Figures?
fig = 1; 

%Shock in period?
shock = 0; %True or false
shock_period = 50; 
shock_multiplier = 2; %Multiplier on theta_bar_0

%Pick theta_0?
pick_theta = 0; 
picked_theta = 10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE ALL CHANGES IN THE SCRIPT, AND COPY. 
% FROM HERE ON AND DOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preallocate matrices for storage
    %We will generally shift rows between periods, and columns between
    %agents.
pop_value_theta_t    = zeros(T,1); %The value of the social good in each period
pop_mean_mu_t        = zeros(T,1);
%pop_mean_thetabar_t  = zeros(T,1);
pop_var_sigmatheta_t = zeros(T,1);
pop_value_rho_t      = zeros(T,1);
prn_value_x_t        = zeros(T,1);   
pop_value_chix_t     = zeros(T,1);
pop_value_abar_t     = zeros(T,1);

agt_value_v_it       = zeros(T,N); 
agt_value_theta_it   = zeros(T,N);    
agt_value_mu_it      = zeros(T,N);   
agt_value_a_it       = zeros(T,N);   
agt_value_abaractual_t = zeros(T,1);
agt_exp_theta_t     = zeros(T,1); 
agt_value_utility   = zeros(T,N);
agt_value_utilitynorm   = zeros(T,N);
agt_value_reputil     =zeros(T,N); 
agt_value_contFB    = zeros(T,1);
agt_value_utilityFB = zeros(T,1);

%pr_value_thetahat_t             = zeros(T,1);
pr_value_thetahataactual_t      = zeros(T,1);
%pr_matrix_sigmatheta_t          = zeros(T,1);
pr_value_gammax                  = zeros(T,1);
%pr_matrix_covariance            = zeros(T,T); 
%pr_matrix_invcovariance         = zeros(T,T); %We meed the matrices to be the smaller size at intermediate calculations
pr_value_aP_t                    = zeros(T,1);
pr_value_aPFB_t                  = zeros(T,1); 
pr_value_utility                  = zeros(T,1);
pr_value_utilitynorm                  = zeros(T,1);
pr_value_utilityFB              = zeros(T,1);

%% SIMULATE GAME
    %We reset variables, since they need proper dimensions
    %clearvars pr_matrix_sigmatheta_t pr_matrix_covariance 
for t=1:T

% PRELIMNARY DRAWS AND THE LIKE
%The value of the public good evolves as a stochastic process from second
%period. 
    %We also update the variance expressions 
if t==1
    if pick_theta == 1
    pop_value_theta_t(t,1)      = picked_theta;    
    else
    pop_value_theta_t(t,1)      = normrnd(pop_mean_thetabar_0, sqrt(pop_var_sigmatheta_0));
    end 
    pop_mean_thetabar_t(t,1)    = pop_mean_thetabar_0; 
    pop_var_sigmatheta_t(t,1)   = pop_var_sigmatheta_0;
    elseif (shock==1 && t==shock_period)
        shock_size = shock_multiplier*pop_mean_thetabar_0; 
    pop_value_theta_t(t,1)      = pop_value_theta_t(t-1,1)+shock_size;
    pop_mean_thetabar_t(t,1)    = pop_mean_thetabar_t(t-1,1) + pop_mean_etabar;
    pop_var_sigmatheta_t(t,1)   = pop_var_sigmatheta_t(t-1,1)+pop_var_sigmaeta;       
    else 
    pop_value_theta_t(t,1)      = pop_value_theta_t(t-1,1)+normrnd(pop_mean_etabar, sqrt(pop_var_sigmaeta));
    pop_mean_thetabar_t(t,1)    = pop_mean_thetabar_t(t-1,1) + pop_mean_etabar;
    pop_var_sigmatheta_t(t,1)   = pop_var_sigmatheta_t(t-1,1)+pop_var_sigmaeta;
end

% PRINCIPALE SETS x_t
    %ALTERNATIVE: Exogenous x_t
    prn_value_x_t(t,1) = prn_value_x_constant; 

%Draw average reputational motivation
pop_mean_mu_t(t,1) = normrnd(pop_mean_mubar,sqrt(pop_var_sigmamu));

%Compute the rho used in expectations
pop_value_rho_t(t,1) = pop_var_sigmatheta_t(t,1)/(pop_var_sigmatheta_t(t,1)+pop_var_s2theta);

%Solve chi(x)
    syms y  
    eqn = y ==  pop_var_s2theta/...
                 (prn_value_x_t(t,1)^2*y^2*pop_var_s2mu ...
                 + pop_var_s2v ...
                 + pop_value_rho_t(t,1)*pop_var_s2theta);
    sol = solve(eqn, y,'Real',true);
    %chi_x_0 = feval(symengine,'numeric::solve','1-0.99-(u+1)*exp(-u)','AllRealRoots')
pop_value_chix_t(t,1) = double(sol);
    clearvars y eqn sol;


% AGENTS
for i=1:N
% Draw type 
agt_value_v_it(t,i) = normrnd(pop_mean_vbar, sqrt(pop_var_s2v));
agt_value_theta_it(t,i)  = pop_value_theta_t(t,1) + normrnd(0, sqrt(pop_var_s2theta));
agt_value_mu_it(t,i) = normrnd(pop_mean_mu_t(t,1), sqrt(pop_var_s2mu));

%Form expectation on theta
agt_exp_theta_t(t,i)  = (1-pop_value_rho_t(t,1))    *pop_mean_thetabar_t(t,1) ...
                        + pop_value_rho_t(t,1)       *agt_value_theta_it(t,i);

%Choose contribution level
agt_value_a_it(t,i) =   agt_value_v_it(t,i) ...
                        +agt_exp_theta_t(t,i) ...
                        +prn_value_x_t(t,1)*agt_value_mu_it(t,i)*pop_value_chix_t(t,1);

%Average contribution levels
pop_value_abar_t(t,1)   =   pop_mean_vbar ...
                            + pop_value_rho_t(t,1)*pop_value_theta_t(t,1) ...
                            +(1-pop_value_rho_t(t,1))*pop_mean_thetabar_t(t,1) ...
                            + prn_value_x_t(t,1)*pop_mean_mu_t(t,1)*pop_value_chix_t(t,1); 
agt_value_abaractual_t(t,1)   = mean(agt_value_a_it(t,:));

end %End of loop over agents
    

% PRINCIPAL SETS CONTRIBUTION

%First, she computes the newest signal
%Define gamma as a helper value
pr_value_gammax(t,1) = ((prn_value_x_t(t,1)*pop_value_chix_t(t,1))/pop_value_rho_t(t,1)); 

pr_value_thetahat_t(t,1) =  pop_value_theta_t(t,1) ...
                           + pr_value_gammax(t,1)*(pop_mean_mu_t(t,1)-pop_mean_mubar);

%TODO: Compute based on actual observed values.                        
                       
%privalue_thetahat_actual 

pr_matrix_sigmatheta_t(t,1) =  pop_var_sigmatheta_t(t,1);

%Fill out covariance matrix iteratively with a few clever codes 
pr_matrix_covariance(t,t) =     pop_var_sigmatheta_t(t,1) ...
                                + pr_value_gammax(t,1)^2*pop_var_sigmamu;
    for k=1:t-1
    pr_matrix_covariance(t,k) = pop_var_sigmatheta_t(k,1);
    pr_matrix_covariance(k,t) = pop_var_sigmatheta_t(k,1);
    end
    clearvars k

    
%pr_value_phi =     
pr_exp_theta_t(t,1) =  pop_mean_thetabar_t(t,1) ...
                     + (transpose(pr_matrix_sigmatheta_t)/pr_matrix_covariance)...
                     *(pr_value_thetahat_t...
                         - pop_mean_thetabar_t);
                     
%Then, she computes her contribution
pr_value_varphi          = pr_value_lambda+(1-pr_value_lambda)*pr_value_b;
pr_value_aP_t(t,1)       =  ((pop_value_w+pr_exp_theta_t(t,1))*pr_value_varphi)/...
                            (1-pr_value_lambda)*pr_value_kp;

pr_value_aPFB_t(t,1)    =  ((pop_value_w+pop_value_theta_t(t,1))*pr_value_varphi)/...
                            (1-pr_value_lambda)*pr_value_kp;                      
        
% CALCULATE UTILITIES
%Agents (need principals contribution)
%We only calculate direct utilities, as reputational ones cancel out over
%population
for i=1:N
agt_value_utility(t,i)  = (agt_value_v_it(t,i)+pop_value_theta_t(t,1))*agt_value_a_it(t,i)...
                            +(pop_value_w+pop_value_theta_t(t,1))*...
                            (pop_value_abar_t(t,1)+pr_value_aP_t(t,1))...
                            - 0.5*agt_value_a_it(t,i)^2;
 


%Normalize utility by value of normal good 
agt_value_utilitynorm(t,i) =  agt_value_utility(t,i)/pop_value_theta_t(t,1);     

end %End short utility loop

%Check that reputation cancels out 
%mean(mean(agt_value_reputil,2)-pop_mean_vbar,1);

%Principal
pr_value_utility(t,1) = pr_value_lambda*(pr_value_alpha*((pop_mean_vbar+pop_value_theta_t(t,1))*pop_value_abar_t(t,1)) ...
                        +(pop_value_w+pop_value_theta_t(t,1))*...
                            (pop_value_abar_t(t,1)+pr_value_aP_t(t,1)) ...
                            - 0.5*pop_value_abar_t(t,1)^2) ...
                        +(1-pr_value_lambda)*(pr_value_b*((pop_value_w+pop_value_theta_t(t,1))*(pop_value_abar_t(t,1)+pr_value_aP_t(t,1)))...
                        - pr_value_kp*0.5*pr_value_aP_t(t,1)^2);

pr_value_utilitynorm(t,1) =  pr_value_utility(t,1)/(pop_value_theta_t(t,1));

%     if t == 1
%         pr_value_utilitydetrend = pr_value_utility(t,1); 
%     else
%         Delta_Theta = (pop_value_theta_t(t,1)-pop_value_theta_t(t-1,1));
%     pr_value_utilitydetrend(t,1) = pr_value_utility(t,1) - ...
%                                   pr_value_lambda*( ...
%                                   pr_value_alpha*(Delta_Theta*pop_value_abar_t(t-1,1)) ...
%                                   +Delta_Theta*(pop_value_abar_t(t-1,1)+pr_value_aP_t(t-1,1))) ...
%                                   + (1-pr_value_lambda)*pr_value_b*Delta_Theta*(pop_value_abar_t(t-1,1)+pr_value_aP_t(t-1,1));
%                                 
%         clearvars Delta_Theta
%     
%     end

%Agent aggregate first-best
agt_value_contFB(t,1) = pop_mean_vbar + pop_value_w+ pop_value_theta_t(t,1);

%First-best utilities 
agt_value_utilityFB(t,1)   = (pop_mean_vbar+pop_value_theta_t(t,1))*agt_value_contFB(t,1) ...
                            +(pop_value_w+pop_value_theta_t(t,1))*...
                            (agt_value_contFB(t,1)+pr_value_aPFB_t(t,1)) ...
                            - 0.5*agt_value_contFB(t,1)^2;

pr_value_utilityFB(t,1) = pr_value_lambda*(pr_value_alpha*((pop_mean_vbar+pop_value_theta_t(t,1))*agt_value_contFB(t,1)) ...
                        +(pop_value_w+pop_value_theta_t(t,1))*...
                            (agt_value_contFB(t,1)+pr_value_aPFB_t(t,1)) ...
                            - 0.5*agt_value_contFB(t,1)^2) ...
                        +(1-pr_value_lambda)*(pr_value_b*((pop_value_w+pop_value_theta_t(t,1))*(agt_value_contFB(t,1)+pr_value_aPFB_t(t,1)))...
                        - pr_value_kp*0.5*pr_value_aPFB_t(t,1)^2);                        


end %End of loop over time periods

%Compute a few averages we use in the the figures for ease
agt_exp_avgexptheta_t      = mean(agt_exp_theta_t,2);
agt_value_avgutility       = mean(agt_value_utility,2);
 

%% FIGURES
   if fig == 1
%set(groot,'DefaultAxesLineStyleOrder','-|-|--|--')
%set(groot,'FontName','CMU Serif')

%%%% FIGURE 1. How do expectations evolve over time?
figure;  set(gcf, 'Position', [80,80, 800, 800])
        subplot(2,2,[1 2])
        fig1 = plot([pop_value_theta_t, pop_mean_thetabar_t, pr_exp_theta_t,agt_exp_avgexptheta_t]);
        lgd = legend('True value','Prior belief', 'Principal', 'Average agent',...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Evolution of expectations','FontName','CMU Serif');
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Value of public good','FontName','CMU Serif');
        %fig1.linestyle = '-|--|-|-';
        set(fig1, {'LineStyle'}, {'-';'--';'-';'-'});
        set(fig1, {'Marker'}, {'none';'none';'o';'square'}); 
        set(fig1, 'Markersize', 3); 
        set(fig1, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        set(fig1, 'Linewidth', 1);
        set(fig1, {'Color'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        
        
        
        subplot(2,2,3)
        fig2 = plot([pr_value_aPFB_t,pr_value_aP_t,agt_value_contFB,agt_value_abaractual_t]);
        lgd = legend('Principal (FB)','Principal', 'Avg. agent (FB)', 'Avg. agent',...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Contributions','FontName','CMU Serif');
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Contributions','FontName','CMU Serif');
        set(fig2, {'Color'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[  0.4660    0.6740    0.1880];[0.8500    0.3250    0.0980]});
        set(fig2, {'Marker'}, {'none';'o';'x';'square'}); 
        set(fig2, 'Markersize', 2); 
        set(fig2, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[  0.4660    0.6740    0.1880];[0.8500    0.3250    0.0980]});
         set(fig1, 'Linewidth', 1);
   
        subplot(2,2,4)
        fig3 = plot([pr_value_utilityFB, pr_value_utility, agt_value_utilityFB, agt_value_avgutility]);
        lgd = legend('Principal (FB)','Principal', 'Avg. agent (FB)', 'Avg. agent', ...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Utility','FontName','CMU Serif'); 
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Utility','FontName','CMU Serif');
        set(fig3, {'Color'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[  0.4660    0.6740    0.1880];[0.8500    0.3250    0.0980]});
        set(fig3, {'Marker'}, {'none';'o';'x';'square'}); 
        set(fig3, 'Markersize', 2); 
        set(fig3, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[  0.4660    0.6740    0.1880];[0.8500    0.3250    0.0980]});
        
        print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\1_WorkingFile', '-depsc'); 
   end %End of figures

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% END OF FUNCTION