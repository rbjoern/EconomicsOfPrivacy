clear all; clc; close all;  

%% GLOBAL PARAMETER DEFAULTS
%Randomization
seed = 13; 
    rng(seed);

%Number of iterations of simulations
T = 100; %Game runs for T periods 
N = 1000; %We simulate a population of N. 

%Population means
pop_mean_thetabar_0         = 5;        %unifrnd(5,10);
pop_mean_etabar             = 0.03;     %unifrnd(0.01,0.02);
pop_mean_vbar               = 10;          %unifrnd(2,3);
pop_mean_mubar              = 10;           %unifrnd(1,2); %Large enough that a high share has positive motivation

%Population variances
pop_var_sigmatheta_0         = 5;       %unifrnd(5,8);
    pop_var_sigmaeta         = 0.5;     %unifrnd(0.1, 0.2);
pop_var_s2v                  = 5;         %unifrnd(2,3);
pop_var_sigmamu              = 5;         %unifrnd(2,3);
pop_var_s2mu                 = 5;
pop_var_s2theta              = 15;

%Population constants
pop_value_w                  = 5;             %unifrnd(1,2); %Must be lower than pop_mean_vbar
pr_value_lambda              = 0.5;
pr_value_alpha               = 0.5; 
pr_value_kp                  = 1; 
pr_value_b                   = 0.5;
prn_value_x_constant         = 2; 

%Figures?
fig =0; 

%Shock in period?
shock = 0; %True or false
shock_period = 50; 
shock_multiplier = 4; %Multiplier on theta_bar_0

%Pick theta_0?
pick_theta = 1; 
picked_theta = 10; 

%% RUN POSSIBLE LOOP

%% CHANGE SPECIFIC PARAMETERS

%% RUN GAME FROM FUNCTION
[   pop_value_theta_t ...
                ,pop_value_abar_t ...
                , agt_exp_theta_t ...
                , pr_value_thetahat_t  ...
                , pr_exp_theta_t ...
                , pop_mean_thetabar_t ...
                , pr_value_aPFB_t ...
                , pr_value_aP_t ...
                , agt_value_abaractual_t ...
                , agt_value_utility ...
                , pr_value_utility ...
                , agt_exp_avgexptheta_t ... 
                , agt_value_avgutility ... 
            ] = ...
EoP.SimulateGame( seed ...
                 ,T ...
                 ,N ....
                 ,pop_mean_thetabar_0 ...
                 ,pop_mean_etabar ...
                 ,pop_mean_vbar ...
                 ,pop_mean_mubar ...
                 ,pop_var_sigmatheta_0 ...
                 ,pop_var_sigmaeta ...
                 ,pop_var_s2v ...
                 ,pop_var_sigmamu ...
                 ,pop_var_s2mu ...
                 ,pop_var_s2theta ...
                 ,pop_value_w ...
                 ,pr_value_lambda ...
                 ,pr_value_alpha ...
                , pr_value_kp ...
                , pr_value_b ... 
                , prn_value_x_constant ...
                , fig ...
                , shock ...
                , shock_period ...
                , shock_multiplier ...
                , pick_theta...
                , picked_theta ...
                );

%% DEFINE OUTPUT OF LOOP IF RELEVANT

%%FIGURES IF RELEVANT
fig_here = 1;

if fig_here == 1;

%%%% FIGURE 1. How do expectations evolve over time?
figure;  set(gcf, 'Position', [80,80, 800, 800])
        subplot(2,2,[1 2])
        fig1 = plot([pop_value_theta_t, pop_mean_thetabar_t, pr_exp_theta_t,agt_exp_avgexptheta_t]);
        lgd = legend('True value','Prior belief', 'Principal', 'Average agent',...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        title('Evolution of expectations','FontName','CMU Serif');
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Value of public good','FontName','CMU Serif');
        %fig1.linestyle = '-|--|-|-';
        set(fig1, {'LineStyle'}, {'-';'--';'-';'-'});
        set(fig1, {'Marker'}, {'none';'none';'o';'square'}); 
        set(fig1, 'Markersize', 3); 
        set(fig1, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        %set(fig1, 'Linewidth', 1.2);
        set(fig1, {'Color'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        
        
        
        subplot(2,2,3)
        fig2 = plot([pr_value_aPFB_t,pr_value_aP_t,agt_value_abaractual_t]);
        lgd = legend('Principals first-best','Principal', 'Average agent',...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        title('Contributions','FontName','CMU Serif');
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Contributions','FontName','CMU Serif');
        set(fig2, {'Color'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[0.8500    0.3250    0.0980]});
        set(fig2, {'Marker'}, {'none';'o';'square'}); 
        set(fig2, 'Markersize', 3); 
        set(fig2, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                    ;[0.8500    0.3250    0.0980]});
        %set(fig1, {'LineStyle'}, {'-';'-';'-';'-'});
   
        subplot(2,2,4)
        fig3 = plot([pr_value_utility, agt_value_avgutility]);
        lgd = legend('Principal', 'Average agent', ...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        title('Utility','FontName','CMU Serif'); 
        xlabel('Time period','FontName','CMU Serif');
        ylabel('Utility','FontName','CMU Serif');
        set(fig3, {'Color'}, {[0   0.4470    0.7410]...
                    ;[0.8500    0.3250    0.0980]});
        set(fig3, {'Marker'}, {'o';'square'}); 
        set(fig3, 'Markersize', 3); 
        set(fig3, {'MarkerFaceColor'}, {[0   0.4470    0.7410]...
                    ;[0.8500    0.3250    0.0980]});
        
        print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\1_WorkingFile', '-depsc'); 
   end %End of figures
 
     
  

print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\Figure_Default', '-depsc'); 
