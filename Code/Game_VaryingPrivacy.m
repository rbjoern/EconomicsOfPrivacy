clear all; clc; close all;  

%% GLOBAL PARAMETER DEFAULTS
%Randomization
seed = 11; 
    rng(seed);

%Number of iterations of simulations
T = 100; %Game runs for T periods 
N = 1000; %We simulate a population of N. 

%Population means
pop_mean_thetabar_0         = 5;        %unifrnd(5,10);
pop_mean_etabar             = 0.03;     %unifrnd(0.01,0.02);
pop_mean_vbar               = 5;          %unifrnd(2,3);
pop_mean_mubar              = 10;           %unifrnd(1,2); %Large enough that a high share has positive motivation

%Population variances
pop_var_sigmatheta_0         = 5;       %unifrnd(5,8);
    pop_var_sigmaeta         = 0.5;     %unifrnd(0.1, 0.2);
pop_var_s2v                  = 5;         %unifrnd(2,3);
pop_var_sigmamu              = 5;         %unifrnd(2,3);
pop_var_s2mu                 = 5;
pop_var_s2theta              = 15;

%Population constants
pop_value_w                  = 10;             %unifrnd(1,2); %Must be greater than pop_mean_vbar
pr_value_lambda              = 0.5;
pr_value_alpha               = 0.5; 
pr_value_kp                  = 1; 
pr_value_b                   = 0.5;
prn_value_x_constant         = 0.3; 

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
for n = 1:2 
   

%% CHANGE SPECIFIC PARAMETERS
    if n==2
    prn_value_x_constant = 10*prn_value_x_constant;
    end 

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
                , agt_value_contFB... 
                , agt_value_utilityFB ...
                , pr_value_utilityFB ...
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
%pop_value_theta_t(:,n)             = pop_value_theta_t; %Check 
x_check(:,n) = prn_value_x_constant;

pr_exp_theta_t_x(:,n)                 = pr_exp_theta_t;
agt_exp_avgexptheta_t_x(:,n)          = agt_exp_avgexptheta_t;
pr_value_aPFB_tx(:,n)                       = pr_value_aPFB_t;  
pr_value_aP_tx(:,n)                         = pr_value_aP_t;      
agt_value_abaractual_tx(:,n)                = agt_value_abaractual_t;  
pr_value_utilityx(:,n)                            =    pr_value_utility;
agt_value_avgutilityx(:,n)                        = agt_value_avgutility; 

end % LOOP 

%%FIGURES IF RELEVANT
fig_here = 1;

if fig_here == 1;

%%%% FIGURE 1. How do expectations evolve over time?
figure;  set(gcf, 'Position', [80,80, 800, 800])
        subplot(2,2,[1 2])
        fig1 = plot([pop_value_theta_t(:,1), pop_mean_thetabar_t(:,1), pr_exp_theta_t_x,agt_exp_avgexptheta_t(:,1)]);
        lgd = legend('True value','Prior belief', 'Principal (low x)','Principal (high x)', 'Average agent',...
            'Location', 'northwest');
        lgd.FontSize = 8;
        lgd.FontName = 'Century';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Evolution of expectations','FontName','Century');
        xlabel('Time period','FontName','Century');
        ylabel('Value of public good','FontName','Century');
        %fig1.linestyle = '-|--|-|-';
        set(fig1, {'LineStyle'}, {'-';'--';'-';'-';'-'});
        set(fig1, {'Marker'}, {'none';'none';'o';'^';'square'}); 
        set(fig1, 'Markersize', 2); 
        set(fig1, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.3010    0.7450    0.9330];[0.8500    0.3250    0.0980]});
        %set(fig1, 'Linewidth', 1.2);
        set(fig1, {'Color'}, {[0.4940    0.1840    0.5560];[0.6350    0.0780    0.1840]...
                                ;[0    0.4470    0.7410];[0.3010    0.7450    0.9330];[0.8500    0.3250    0.0980]});
        %LineStyleOrder('-','-','-','--');
        
        
        subplot(2,2,[3 4])
        fig2 = plot([pr_value_aPFB_t,pr_value_aP_tx,agt_value_contFB, agt_value_abaractual_tx]);
        lgd = legend('Principal (FB)','Principal (low x)','Principal (high x)','Average agent (FB)', 'Average agent (low x)','Average agent (high x)',...
            'Location', 'northwest');
        lgd.FontSize = 8;
        lgd.FontName = 'Century';
        lgd_pos = get(lgd,'position'); 
        set(lgd,'position',[lgd_pos(1), lgd_pos(2), 1.4*lgd_pos(3), lgd_pos(4)])
        title('Contributions','FontName','Century');
        xlabel('Time period','FontName','Century');
        ylabel('Contributions','FontName','Century');
        set(fig2, {'Color'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                                ;[0.3010    0.7450    0.9330];[  0.4660    0.6740    0.1880] ...
                                ;[0.8500    0.3250    0.0980];[0.9290    0.6940    0.1250]});
        set(fig2, {'Marker'}, {'none';'o';'^';'x';'square';'diamond'}); 
        set(fig2, 'Markersize', 2); 
        set(fig2, {'MarkerFaceColor'}, {[0.4940    0.1840    0.5560];[0   0.4470    0.7410]...
                                ;[0.3010    0.7450    0.9330];[  0.4660    0.6740    0.1880] ...
                                ;[0.8500    0.3250    0.0980];[0.9290    0.6940    0.1250]});
        %set(fig1, {'LineStyle'}, {'-';'-';'-';'-'});
   
%         subplot(2,2,4)
%         fig3 = plot([pr_value_utilityx, agt_value_avgutilityx]);
%         lgd = legend('Principal (low x)','Principal (high x)', 'Average agent (low x)','Average agent (high x)', ...
%             'Location', 'northwest');
%         lgd.FontSize = 8;
%         lgd.FontName = 'Century';
%         title('Utility','FontName','Century'); 
%         xlabel('Time period','FontName','Century');
%         ylabel('Utility','FontName','Century');
%         set(fig3, {'Color'}, {[  0    0.4470    0.7410]...
%                                 ;[0.3010    0.7450    0.9330];[0.8500    0.3250    0.0980];[0.9290    0.6940    0.1250]});
%         set(fig3, {'Marker'}, {'o';'^';'square';'diamond'}); 
%         set(fig3, 'Markersize', 3); 
%         set(fig3, {'MarkerFaceColor'}, {[  0    0.4470    0.7410]...
%                                 ;[0.3010    0.7450    0.9330];[0.8500    0.3250    0.0980];[0.9290    0.6940    0.1250]});
%         
        %print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\1_Expectations', '-depsc'); 
   end %End of figures    
     
 

print('C:\Users\rbjoe\Dropbox\Kugejl\9. semester\Economics of Privacy\Matlab\Images\Figure_VaryingPrivacy', '-depsc'); 
