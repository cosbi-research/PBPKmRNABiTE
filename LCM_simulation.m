% Clear the workspace
clear all 
% Close all the figures
close all
% Clear the command window
clc 

%% ************* Source code for the model simulation for both the recombinant and mRNA-encoded protein administration (LCM) *********************

%% *************** Brief descriptio of the code structure: *********************************

% ***************** Specific quantities to reproduce the experimental setting presented in Huang_2023 ******************
% Dose injected and time-window of observation

% ****************** Experimental data collected from Huang_2023 **************************************************
% Mean and confidence intervals for both the recombinant and the
% mRNA-encoded BiTE administration

% ********************** General parameters for the Protein and the mRNA Transportation Modules **************************
% PBPK model parameters that do not change based on protein size (Li_2019)

% ************************ Fitted parameters ***********************************************
% In the case of the Long Chain Model (LCM) we fit 8 parameters
% the fluid phase pinocytosis rate CL_UP
% the transfer rate of the mRNA injected from blood to the liver: k_BL2LV
% the degradation rate of the mRNA: k_d_mRNA
% the rates at which BiTEs are produced via translation in each mRNA phase:
% k_tr_mRNA_LV_1, k_tr_mRNA_LV_2, k_tr_mRNA_LV_3
% the phenomenological transition rates between mRNA phases:
% k_s_12 (from phase 1 to phase 2), k_s_23 (from phase 2 to phase 3)

% ************************ Specific species parameters for the mRNA and the Protein Transportation modules **************************
% Parameters collection related to the mouse physiology
% Parameters collection related to the trasport of 55kDa proteins

% ***************** Setting the model simulation *********************************
% Time interval for the simulation -> time window of the observed data
% Options about the model computational resolution

% ************* Model simulation in the case of Recombinant proteins administration ************
% 1) Setting the initial condition of the simulation 
% 2) Model simulation
% 3) Model results

% ************* Model simulation in the case of mRNA-encoded proteins administration ************
% 1) Setting the initial condition of the simulation 
% 2) Model simulation
% 3) Model results


%% ***************** Specific quantities to reproduce the experimental setting presented in Huang_2023 ******************
% Recombinat dose injected - 6 mg/Kg -> 6 mg*0.02 (mouse weight 0.02 Kg)
dose_BiTE = 6*0.02; % Unit measure: mug (Huang_2023)
% mRNA dose injected - 1.5 mg/Kg -> 1.5 mg*0.02 (mouse weight 0.02 Kg)
dose_mRNA = 1.5*0.02; % Unit measure: mug (Huang_2023)
% Time-window of observation 
T_initial = 0; % Unit measure: h (Huang_2023)
T_final = 168; % Unit measure: h (Huang_2023)

%% ****************** Experimental data collected from Huang_2023 **************************************************
% Time points
t_data = [1, 4, 6, 12, 24, 48, 72, 144, 168]; %Unit measure: hours

% --------------------------------------------------------------------------------------
% Mean concentration of BiTEs in PLASMA in the case of Recombinant proteins
% administration
mean_data = [6227.272727272725, 1772.727272727272, 1090.90909090909, 295.454524545454413,  204.54545454545405, 90.9090909090919, 90.9090909090919, 45.45454545454595, 45.45454545454595]; %Unit measure: ng/mL
% Confidence intervals
plus_data = [7772.72727272727, 2181.818181818181, 1454.545454545454, 409.0909090909081, 318.181818181818, 204.54545454545405, 181.81818181818198, 136.36363636363603, 90.9090909090919]; %Unit measure: ng/mL
minus_data = [4681.818181818181, 1272.727272727273, 727.2727272727261, 159.09090909090992, 0, 0, 0, 0, 0]; %Unit measure: ng/mL

% ------------------------------------------------------------------------------------------
% Mean concentration of BiTEs in PLASMA in the case of mRNA-encoded proteins
% administration
 mean_data_LNP = [0, 3221.0526315789484, 6463.1578947368425, 5284.21052631579,  4547.368421052633, 3957.894736842106, 3136.8421052631584, 547.3684210526335, 84.21052631579005]; %Unit measure: ng/mL
% Confidence intervals
 plus_data_LNP = [378.9473684210534, 3894.7368421052633, 7557.894736842106, 5957.894736842107, 5157.894736842107, 4421.052631578948, 4168.421052631579, 842.1052631578968, 84.21052631579005]; %Unit measure: ng/mL
 minus_data_LNP = [0, 2547.3684210526326, 5368.42105263158, 4589.473684210528, 3957.894736842106, 3536.8421052631584, 2147.3684210526326, 231.57894736842172, 0]; %Unit measure: ng/mL

%% ********************** General parameters for the Protein and the mRNA Transportation Modules **************************
% PBPK model parameters that do not change based on protein size

% Proportionality constant between the rate at which antibody transfers 
% from the lymph compartment to the blood compartment and the plasma flow of the given species
C_LNLF = 9.1;  % Unit measure: pure number (Shah_2012)
% Large pore radius
r_L = 22.85; % Unit measure: nm
% Small pore radius
r_S = 4.44;  % Unit measure: nm
% Fractional hydraulic conductance of large pores 
alpha_L = 0.042; % Unit measure: pure number (Li_2019)
% Fractional hydraulic conductance of small pores 
alpha_S = 0.958; % Unit measure: pure number (Li_2019)
% Lysosome degradation rate for all the proteins (Model for recombinant
% proteins)
K_deg = 32.2; % Unit measure: 1/h (Li_2019)
% Lysosome degradation rate for all the proteins (Model for mRNA-encoded
% proteins)
K_deg_LNP = 32.2; 
% Constant on the relative pore abundance (Supp. Mat. Li_2019)
X_J = 0.38; % Unit measure: pure number (Supp. Mat. Li_2019)
% Constant dependent on pore size, relative hydraulic conduntance of large
% and small pores and protein size
Xp = 13197; % Unit measure: nm^3 (Supp. Mat. Li_2019)

%% ************************ Fitted parameters ***********************************************
% In the case of the Long Chain Model (LCM) we fit 8 parameters
% the fluid phase pinocytosis rate CL_UP
% the transfer rate of the mRNA injected from blood to the liver: k_BL2LV
% the degradation rate of the mRNA: k_d_mRNA
% the rates at which BiTEs are produced via translation in each mRNA phase:
% k_tr_mRNA_LV_1, k_tr_mRNA_LV_2, k_tr_mRNA_LV_3
% the phenomenological transition rates between mRNA phases:
% k_s_12 (from phase 1 to phase 2), k_s_23 (from phase 2 to phase 3)

load('LCM_parameters_estimated.mat')
parameters = parameters_estimated;


%% ************************ Specific species parameters for the mRNA and the Protein Transportation modules **************************
% Parameters collection related to the mouse physiology
pars = LCM_parameters_collection_mouse(parameters,C_LNLF,K_deg,K_deg_LNP);
% Parameters collection related to the trasport of 55kDa proteins
pars_dim = Parameters_collection_BiTEs_55kDa(X_J, Xp, alpha_L, alpha_S, pars);

%% -------------------------------------------------------------------------------------------------------------------------------

%% ***************** Setting the model simulation *********************************
% Time interval for the simulation -> time window of the observed data
tspan=[T_initial T_final]; %time interval of resolution
% Options about the model computational resolution
opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% ************* Model simulation in the case of Recombinant proteins administration ************

%% Setting the initial condition of the simulation (the model computes the BiTE concentration in Kg/L)
% In Huang_2023: intra-venous injection 
% The concentration of BiTEs in all the compartments at exception of the
% blood is setted to 0
C0 = zeros(50,1);
% Recobinant BiTE dose -> from mg to Kg we multiple by 10^(-6)
C0_plasma = (dose_BiTE*10^(-6))/(pars.V_TOT_PLASMA); % Unit measure: Kg/L
% BiTE concentration in plasma is computed in Eq.46 in the model
C0(46)=C0_plasma; 

%% Model simulation
% T is the vector that contains the time points at which the solution is computed
% SV is the matrix of the ODE system solution: each column corresponds to
% the time-evolution of the corresponding quantity
[T,SV] = ode15s(@(t,Y) recombinant_model(t,Y,pars,pars_dim), tspan, C0, opt);
% Computed BiTE concentration in PLASMA
C_PLASMA = SV(:,46); % Unit measure: Kg/L
% Computed BiTE concentration in PLASMA comparable with experimental data
C_TOT_ngmL_PLASMA = ((C_PLASMA/(10^3)))*10^12; % Unit measure: ng/mL

%% ************* Model results **********
% Figure 1: Time-evolution of BiTE concentration in PLASMA in the case of
% Recombinant proteins administration
figure(1)
p1=plot(T,C_TOT_ngmL_PLASMA);
p1.LineStyle='-';
p1.Marker='none';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0 0.1 0.4];
p1.Color=[0 0.1 0.4];
p1.LineWidth=2;
hold on
for i=1:length(t_data)
p1=line([t_data(i) t_data(i)],[minus_data(i) plus_data(i)]);
p1.LineStyle='-';
p1.Marker='none';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=1;
hold on
end
hold on
p1=plot(t_data(1:end),plus_data(1:end));
p1.LineStyle='none';
p1.Marker='*';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=3;
hold on
p1=plot(t_data(1:end),minus_data(1:end));
p1.LineStyle='none';
p1.Marker='*';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=3;
hold on
p1=plot(t_data(1:end),mean_data(1:end));
p1.LineStyle='none';
p1.Marker='s';
p1.MarkerSize=4;
p1.MarkerFaceColor=[0.5 0.7 0.9];
p1.Color=[0.5 0.7 0.9];
p1.LineWidth=3;
hold on
xlim([0 T_final])
ylim([0.1 max(plus_data)+0.1*max(plus_data)])
axis square
xlabel('time [hours]')
ylabel('concentration [ng/ml]')
title('BiTE concentration in PLASMA')
set(gca,'FontSize',20)

%% ----------------------------------------------------------------------------------------------------------------------------
%% ************* Model simulation in the case of mRNA-encoded proteins administration ************

%% Setting the initial condition of the simulation (the model computes the BiTE concentration in Kg/L)
% In Huang_2023: intra-venous injection 
% The concentration of BiTEs in all the compartments is setted to 0
C0 = zeros(54,1);
% mRNA-encoded BiTE dose -> from mg to Kg we multiple by 10^(-6)
C0_mRNA_plasma = (dose_mRNA*10^(-6))/(pars.V_TOT_PLASMA); % Unit measure: Kg/L
% mRNA concentration in plasma is computed in Eq.51 in the model
C0(51)=C0_mRNA_plasma;

%% Model simulation
% T_LNP is the vector that contains the time points at which the solution is computed
% SV_LNP is the matrix of the ODE system solution: each column corresponds to
% the time-evolution of the corresponding quantity
[T_LNP,SV_LNP] = ode15s(@(t,Y) LCM(t,Y,pars,pars_dim), tspan, C0, opt);

% Computed BiTE concentration in PLASMA
C_PLASMA_LNP = SV_LNP(:,46); % Unit measure: Kg/L
% Computed BiTE concentration in PLASMA comparable with experimental data
C_TOT_ngmL_PLASMA_LNP = ((C_PLASMA_LNP/(10^3)))*10^12; % Unit measure: ng/mL

%% ************* Model results **********
% Figure 2: Time-evolution of BiTE concentration in PLASMA in the case of
% mRNA-encoded proteins administration

figure(2)
p1=plot(T_LNP,C_TOT_ngmL_PLASMA_LNP);
p1.LineStyle='-';
p1.Marker='none';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0 0.1 0.4];
p1.Color=[0 0.1 0.4];
p1.LineWidth=2;
hold on
for i=1:length(t_data)
p1=line([t_data(i) t_data(i)],[minus_data_LNP(i) plus_data_LNP(i)]);
p1.LineStyle='-';
p1.Marker='none';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=1;
hold on
end
hold on
p1=plot(t_data(1:end),plus_data_LNP(1:end));
p1.LineStyle='none';
p1.Marker='*';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=3;
hold on
p1=plot(t_data(1:end),minus_data_LNP(1:end));
p1.LineStyle='none';
p1.Marker='*';
p1.MarkerSize=2.5;
p1.MarkerFaceColor=[0.8 0.9 1];
p1.Color=[0.8 0.9 1];
p1.LineWidth=3;
hold on
p1=plot(t_data(1:end),mean_data_LNP(1:end));
p1.LineStyle='none';
p1.Marker='s';
p1.MarkerSize=4;
p1.MarkerFaceColor=[0.5 0.7 0.9];
p1.Color=[0.5 0.7 0.9];
p1.LineWidth=3;
xlim([0 T_final])
ylim([0.1 max(plus_data_LNP)+0.1*max(plus_data_LNP)])
axis square
xlabel('time [hours]')
ylabel('concentration [ng/ml]')
title('BiTE concentration in PLASMA')
set(gca,'FontSize',20)
