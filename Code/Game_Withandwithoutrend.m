clear all; clc; close all;

%% GLOBAL PARAMETER DEFAULTS
%Randomization
seed = 66; 
    rng(seed);

%Number of iterations of simulations
T = 100; %Game runs for T periods 
N = 10000; %We simulate a population of N. 

%Population means
pop_mean_thetabar_0         = unifrnd(5,10);
pop_mean_etabar             = unifrnd(0.05,0.1);
pop_mean_vbar               = unifrnd(5,10);
pop_mean_mubar              = unifrnd(5,10); %Large enough that a high share has positive motivation

%Population variances 
pop_var_sigmatheta_0          = unifrnd(5,8);
    pop_var_sigmaeta         = unifrnd(0.1, 0.5);
pop_var_s2v                  = unifrnd(2,8);
pop_var_sigmamu              = unifrnd(2,8);
pop_var_s2mu                 = unifrnd(2,8);
pop_var_s2theta              = unifrnd(15,20);

%Population constants
pop_value_w                  = unifrnd(1,2); %Must be lower than pop_mean_vbar
pr_value_lambda              = 0.5;
pr_value_alpha               = 0.8; 
pr_value_kp                  = 1; 
pr_value_b                   = 1;
prn_value_x_constant         = 1; 

%Figures?
fig = 0; 

%Shock in period?
shock = 0; %True or false
shock_period = 50; 
shock_multiplier = 4; %Multiplier on theta_bar_0

%Pick theta_0?
pick_theta = 0; 
picked_theta = 10; 

%% RUN POSSIBLE LOOP
figure; set(gcf, 'Position', [80,80, 800, 800]) 

seed = 1336;
%T = 300;
%N=100;
for n=1:2 

    
if n==2 
   pop_mean_etabar = 0; 
   pop_var_sigmaeta = 0;
end
 

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
%%%% FIGURE 1. How do expectations evolve over time?
        
        subplot(2,1,n)
       fig1 = plot([pop_value_theta_t, pop_mean_thetabar_t, pr_exp_theta_t, mean(agt_exp_theta_t,2)]);
       if n ==1  
       lgd = legend('True value','Prior belief', 'Principal', 'Avg. agent',...
            'Location', 'northwest');
        lgd.FontSize = 10;
        lgd.FontName = 'CMU Serif';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Evolution of expectations');
        ylabel('Random walk','FontName','CMU Serif'); 
        xlabel('Time period','FontName','CMU Serif');
       end
       if n ==2 
       title('Evolution of expectations');    
       xlabel('Time period','FontName','CMU Serif');
       ylabel('Fixed value','FontName','CMU Serif');
       end 
        set(fig1, {'LineStyle'}, {'-';'--';'-';'-'});
        set(fig1, {'Marker'}, {'none';'none';'o';'square'}); 
        set(fig1, 'Markersize', 2); 
        set(fig1, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        %set(fig1, 'Linewidth', 1);
        set(fig1, {'Color'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.8500    0.3250    0.0980]});
        
       


%%FIGURES IF RELEVANT

end

print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\FigurWithandWithoutTrend', '-depsc'); 