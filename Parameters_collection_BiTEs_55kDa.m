function pars_dim = Parameters_collection_BiTEs_55kDa(X_J, Xp, alpha_L, alpha_S, pars)

%% Table 4 of Li_2019 for proteins of 55kDa 

%% ******************** General parameters for 55kDa proteins ********************************
% Stokes–Einstein radius
pars_dim.a_e = 3.26; % Unit measure: nm
% Small pore vascular reflection coefficient
pars_dim.sigma_S = 0.906; % Unit measure: pure number
% Large pore vascular reflection coefficient
pars_dim.sigma_L = 0.0877; % Unit measure: pure number
% Peclet number of small pore
pars_dim.Pe_S = 0.113; % Unit measure: pure number
% Peclet number of large pore
pars_dim.Pe_L =2.254425255333927; % Unit measure: pure number
% Glomerular sieving coefficient
pars_dim.theta = 0.0986; % Unit measure: pure number
% Fractional accessible pore size of small pore
pars_dim.AonA0_S = 2.45e-03; % Unit measure: pure number
% Fractional accessible pore size of large pore
pars_dim.AonA0_L = 0.522; % Unit measure: pure number
% Relative ratio between protein size and small pore size
pars_dim.r_S = pars_dim.a_e / 0.735; % Unit measure: pure number
% Relative ratio between protein size and large pore size
pars_dim.r_L = pars_dim.a_e / 0.143; % Unit measure: pure number
% Fractional hydraulic conductance of large pores 
pars_dim.alpha_L = alpha_L;
% Fractional hydraulic conductance of small pores 
pars_dim.alpha_S = alpha_S;
% Renal clearence - Eq.6 Li_2019
pars_dim.CL_R = pars.GFR * pars_dim.theta; % Unit measure: L/h
% Constant related to two pores dynamics for small pores Eq.23 in Li_2019
% Supplementary materials
pars_dim.Xp_S = Xp * (1/pars_dim.a_e) * pars_dim.AonA0_S * (alpha_S/pars_dim.r_S^2); % Unit measure: pure number
% Constant related to two pores dynamics for large pores Eq.23 in Li_2019
% Supplementary materials
pars_dim.Xp_L = Xp * (1/pars_dim.a_e) * pars_dim.AonA0_L * (alpha_L/pars_dim.r_L^2); % Unit measure: pure number

%% ******************** Parameters that regulate the proteins transport for 55kDa proteins ********************************

%% ++++++++++++++++  Permeability-surface area product of large pore ++++++++++++++++++++++++++++++++++
pars_dim.PS_L_HEART = pars_dim.Xp_L * pars.L_HEART;
pars_dim.PS_L_LUNG = pars_dim.Xp_L * pars.L_LUNG;
pars_dim.PS_L_MUSCLE = pars_dim.Xp_L * pars.L_MUSCLE;
pars_dim.PS_L_SKIN = pars_dim.Xp_L * pars.L_SKIN;
pars_dim.PS_L_ADIPOSE = pars_dim.Xp_L * pars.L_ADIPOSE;
pars_dim.PS_L_BONE = pars_dim.Xp_L * pars.L_BONE;
pars_dim.PS_L_BRAIN = pars_dim.Xp_L * pars.L_BRAIN;
pars_dim.PS_L_KIDNEY = pars_dim.Xp_L * pars.L_KIDNEY;
pars_dim.PS_L_LIVER = pars_dim.Xp_L * pars.L_LIVER;
pars_dim.PS_L_SM_INT = pars_dim.Xp_L * pars.L_SM_INT;
pars_dim.PS_L_LG_INT = pars_dim.Xp_L * pars.L_LG_INT;
pars_dim.PS_L_PANCREAS = pars_dim.Xp_L * pars.L_PANCREAS;
pars_dim.PS_L_THYMUS = pars_dim.Xp_L * pars.L_THYMUS;
pars_dim.PS_L_SPLEEN = pars_dim.Xp_L * pars.L_SPLEEN;
pars_dim.PS_L_OTHER = pars_dim.Xp_L * pars.L_OTHER;
pars_dim.PS_L_TUMOUR = pars_dim.Xp_L * pars.L_TUMOUR;

%% +++++++++++++++++++++++++++ Permeability-surface area product of small pore ++++++++++++++++++++++++++++++++++++++++++++++++++++
pars_dim.PS_S_HEART = pars_dim.Xp_S * pars.L_HEART;
pars_dim.PS_S_LUNG = pars_dim.Xp_S * pars.L_LUNG;
pars_dim.PS_S_MUSCLE = pars_dim.Xp_S * pars.L_MUSCLE;
pars_dim.PS_S_SKIN = pars_dim.Xp_S * pars.L_SKIN;
pars_dim.PS_S_ADIPOSE = pars_dim.Xp_S * pars.L_ADIPOSE;
pars_dim.PS_S_BONE = pars_dim.Xp_S * pars.L_BONE;
pars_dim.PS_S_BRAIN = pars_dim.Xp_S * pars.L_BRAIN;
pars_dim.PS_S_KIDNEY = pars_dim.Xp_S * pars.L_KIDNEY;
pars_dim.PS_S_LIVER = pars_dim.Xp_S * pars.L_LIVER;
pars_dim.PS_S_SM_INT = pars_dim.Xp_S * pars.L_SM_INT;
pars_dim.PS_S_LG_INT = pars_dim.Xp_S * pars.L_LG_INT;
pars_dim.PS_S_PANCREAS = pars_dim.Xp_S * pars.L_PANCREAS;
pars_dim.PS_S_THYMUS = pars_dim.Xp_S * pars.L_THYMUS;
pars_dim.PS_S_SPLEEN = pars_dim.Xp_S * pars.L_SPLEEN;
pars_dim.PS_S_OTHER = pars_dim.Xp_S * pars.L_OTHER;
pars_dim.PS_S_TUMOUR = pars_dim.Xp_S * pars.L_TUMOUR;

%% ++++++++++++++++++++++++++++ Circular isogravimetric flow ++++++++++++++++++++++++++++++++++++++++++++++++++++
% J_ISO for large pore. For small pores we use -J_ISO
J_ISO_HEART = X_J * pars.L_HEART;
J_ISO_LUNG = X_J * pars.L_LUNG;
J_ISO_MUSCLE = X_J * pars.L_MUSCLE;
J_ISO_SKIN = X_J * pars.L_SKIN;
J_ISO_ADIPOSE = X_J * pars.L_ADIPOSE;
J_ISO_BONE = X_J * pars.L_BONE;
J_ISO_BRAIN = X_J * pars.L_BRAIN;
J_ISO_KIDNEY = X_J * pars.L_KIDNEY;
J_ISO_LIVER = X_J * pars.L_LIVER;
J_ISO_SM_INT = X_J * pars.L_SM_INT;
J_ISO_LG_INT = X_J * pars.L_LG_INT;
J_ISO_PANCREAS = X_J * pars.L_PANCREAS;
J_ISO_THYMUS = X_J * pars.L_THYMUS;
J_ISO_SPLEEN = X_J * pars.L_SPLEEN;
J_ISO_OTHER = X_J * pars.L_OTHER;
J_ISO_TUMOUR = X_J * pars.L_TUMOUR;

%% ++++++++++++++++++ Lymph flow through large pores ++++++++++++++++++++++++
% General law: J_L = J_ISO + alpha_L*J where J is the lymphatic flux -> Eq.
% 8 in Li_2019 Supplementary materials
pars_dim.J_L_HEART = J_ISO_HEART + alpha_L*pars.L_HEART;
pars_dim.J_L_LUNG = J_ISO_LUNG + alpha_L*pars.L_LUNG;
pars_dim.J_L_MUSCLE = J_ISO_MUSCLE + alpha_L*pars.L_MUSCLE;
pars_dim.J_L_SKIN = J_ISO_SKIN + alpha_L*pars.L_SKIN;
pars_dim.J_L_ADIPOSE = J_ISO_ADIPOSE + alpha_L*pars.L_ADIPOSE;
pars_dim.J_L_BONE = J_ISO_BONE + alpha_L*pars.L_BONE;
pars_dim.J_L_BRAIN = J_ISO_BRAIN + alpha_L*pars.L_BRAIN;
pars_dim.J_L_KIDNEY = J_ISO_KIDNEY + alpha_L*pars.L_KIDNEY;
pars_dim.J_L_LIVER = J_ISO_LIVER + alpha_L*pars.L_LIVER;
pars_dim.J_L_SM_INT = J_ISO_SM_INT + alpha_L*pars.L_SM_INT;
pars_dim.J_L_LG_INT = J_ISO_LG_INT + alpha_L*pars.L_LG_INT;
pars_dim.J_L_PANCREAS = J_ISO_PANCREAS + alpha_L*pars.L_PANCREAS;
pars_dim.J_L_THYMUS = J_ISO_THYMUS + alpha_L*pars.L_THYMUS;
pars_dim.J_L_SPLEEN = J_ISO_SPLEEN + alpha_L*pars.L_SPLEEN;
pars_dim.J_L_OTHER = J_ISO_OTHER + alpha_L*pars.L_OTHER;
pars_dim.J_L_TUMOUR = J_ISO_TUMOUR + alpha_L*pars.L_TUMOUR;

%% ++++++++++++++++++ Lymph flow through small pores ++++++++++++++++++++++++
% General law: J_L = - J_ISO + alpha_s*J where J is the lymphatic flux -> Eq.
% 9 in Li_2019 Supplementary materials
pars_dim.J_S_HEART = -J_ISO_HEART + alpha_S*pars.L_HEART;
pars_dim.J_S_LUNG = -J_ISO_LUNG + alpha_S*pars.L_LUNG;
pars_dim.J_S_MUSCLE = -J_ISO_MUSCLE + alpha_S*pars.L_MUSCLE;
pars_dim.J_S_SKIN = -J_ISO_SKIN + alpha_S*pars.L_SKIN;
pars_dim.J_S_ADIPOSE = -J_ISO_ADIPOSE + alpha_S*pars.L_ADIPOSE;
pars_dim.J_S_BONE = -J_ISO_BONE + alpha_S*pars.L_BONE;
pars_dim.J_S_BRAIN = -J_ISO_BRAIN + alpha_S*pars.L_BRAIN;
pars_dim.J_S_KIDNEY = -J_ISO_KIDNEY + alpha_S*pars.L_KIDNEY;
pars_dim.J_S_LIVER = -J_ISO_LIVER + alpha_S*pars.L_LIVER;
pars_dim.J_S_SM_INT = -J_ISO_SM_INT + alpha_S*pars.L_SM_INT;
pars_dim.J_S_LG_INT = -J_ISO_LG_INT + alpha_S*pars.L_LG_INT;
pars_dim.J_S_PANCREAS = -J_ISO_PANCREAS + alpha_S*pars.L_PANCREAS;
pars_dim.J_S_THYMUS = -J_ISO_THYMUS + alpha_S*pars.L_THYMUS;
pars_dim.J_S_SPLEEN = -J_ISO_SPLEEN + alpha_S*pars.L_SPLEEN;
pars_dim.J_S_OTHER = -J_ISO_OTHER + alpha_S*pars.L_OTHER;
pars_dim.J_S_TUMOUR = -J_ISO_TUMOUR + alpha_S*pars.L_TUMOUR;

%% **************** Useful fractions for implementation - ratios in Eq. 10 and Eq. 9 in Li_2019
pars_dim.Pe_S_ratio = pars_dim.Pe_S/(exp(pars_dim.Pe_S) - 1); %frazione in eq 10
pars_dim.Pe_L_ratio = pars_dim.Pe_L/(exp(pars_dim.Pe_L) - 1); %frazione in eq 9
end