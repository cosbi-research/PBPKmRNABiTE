function dydt = recombinant_model(t,y,pars,pars_dim)

dydt=zeros(50,1);

%% Concentration of mAb in vascular compartment for each tissue 
C_V_HEART = y(1);
C_V_LUNG = y(2);
C_V_MUSCLE = y(3);
C_V_SKIN = y(4);
C_V_ADIPOSE = y(5);
C_V_BONE = y(6);
C_V_BRAIN = y(7);
C_V_KIDNEY = y(8);
C_V_LIVER = y(9);
C_V_SM_INT = y(10);
C_V_LG_INT = y(11);
C_V_PANCREAS = y(12);
C_V_THYMUS = y(13);
C_V_SPLEEN = y(14);
C_V_OTHER = y(15);

%% Concentration of mAb in endosomal compartment for each tissue 
C_E_HEART = y(16);
C_E_LUNG = y(17);
C_E_MUSCLE = y(18);
C_E_SKIN = y(19);
C_E_ADIPOSE = y(20);
C_E_BONE = y(21);
C_E_BRAIN = y(22);
C_E_KIDNEY = y(23);
C_E_LIVER = y(24);
C_E_SM_INT = y(25);
C_E_LG_INT = y(26);
C_E_PANCREAS = y(27);
C_E_THYMUS = y(28);
C_E_SPLEEN = y(29);
C_E_OTHER = y(30);

%% Concentration of mAb in interstitial compartment for each tissue 
C_I_HEART = y(31);
C_I_LUNG = y(32);
C_I_MUSCLE = y(33);
C_I_SKIN = y(34);
C_I_ADIPOSE = y(35);
C_I_BONE = y(36);
C_I_BRAIN = y(37);
C_I_KIDNEY = y(38);
C_I_LIVER = y(39);
C_I_SM_INT = y(40);
C_I_LG_INT = y(41);
C_I_PANCREAS = y(42);
C_I_THYMUS = y(43);
C_I_SPLEEN = y(44);
C_I_OTHER = y(45);

%% Compartments without absorption
C_PLASMA = y(46);
C_LYMPH = y(47);

%% Tumour compartment
C_V_TUMOUR = y(48);
C_E_TUMOUR = y(49);
C_I_TUMOUR = y(50);

%% Fluxes in entrance in each organ
% Notation PLQ, L are respectevely Q, J in the paper
% Fluxes in vascular space plasma
M_VI_HEART     = C_V_LUNG*pars.PLQ_HEART;   
M_VI_MUSCLE    = C_V_LUNG*pars.PLQ_MUSCLE;  
M_VI_SKIN      = C_V_LUNG*pars.PLQ_SKIN ;   
M_VI_ADIPOSE   = C_V_LUNG*pars.PLQ_ADIPOSE ;
M_VI_BONE      = C_V_LUNG*pars.PLQ_BONE  ;  
M_VI_BRAIN     = C_V_LUNG*pars.PLQ_BRAIN ;  
M_VI_KIDNEY    = C_V_LUNG*pars.PLQ_KIDNEY;  
M_VI_SM_INT    = C_V_LUNG*pars.PLQ_SM_INT ; 
M_VI_LG_INT    = C_V_LUNG*pars.PLQ_LG_INT ; 
M_VI_PANCREAS  = C_V_LUNG*pars.PLQ_PANCREAS;
M_VI_THYMUS    = C_V_LUNG*pars.PLQ_THYMUS  ;
M_VI_SPLEEN    = C_V_LUNG*pars.PLQ_SPLEEN  ;
M_VI_OTHER     = C_V_LUNG*pars.PLQ_OTHER  ;         
M_VI_LIVER     = C_V_LUNG*pars.PLQ_LIVER + (pars.PLQ_PANCREAS-pars.L_PANCREAS)*C_V_PANCREAS + (pars.PLQ_SM_INT-pars.L_SM_INT)*C_V_SM_INT + (pars.PLQ_LG_INT-pars.L_LG_INT)*C_V_LG_INT + (pars.PLQ_SPLEEN-pars.L_SPLEEN)*C_V_SPLEEN;
M_VI_LUNG      = C_PLASMA*(pars.PLQ_LUNG+pars.L_LUNG); 
M_VI_TUMOUR    = C_V_LUNG*pars.PLQ_TUMOUR  ;


%% Fluxes in exit in each organ

% Fluxes in vascular space plasma
Q_VOUT_HEART       = (pars.PLQ_HEART     - pars.L_HEART    );
Q_VOUT_MUSCLE      = (pars.PLQ_MUSCLE    - pars.L_MUSCLE   );
Q_VOUT_SKIN        = (pars.PLQ_SKIN      - pars.L_SKIN     );
Q_VOUT_ADIPOSE     = (pars.PLQ_ADIPOSE   - pars.L_ADIPOSE  );
Q_VOUT_BONE        = (pars.PLQ_BONE      - pars.L_BONE     );
Q_VOUT_BRAIN       = (pars.PLQ_BRAIN     - pars.L_BRAIN    );
Q_VOUT_KIDNEY      = (pars.PLQ_KIDNEY    - pars.L_KIDNEY   );
Q_VOUT_SM_INT      = (pars.PLQ_SM_INT    - pars.L_SM_INT   );
Q_VOUT_LG_INT      = (pars.PLQ_LG_INT    - pars.L_LG_INT   );
Q_VOUT_PANCREAS    = (pars.PLQ_PANCREAS  - pars.L_PANCREAS );
Q_VOUT_THYMUS      = (pars.PLQ_THYMUS    - pars.L_THYMUS   );
Q_VOUT_SPLEEN      = (pars.PLQ_SPLEEN    - pars.L_SPLEEN   );
Q_VOUT_OTHER       = (pars.PLQ_OTHER     - pars.L_OTHER    );
PLQ_LIVER_UPSTREAM = (pars.PLQ_PANCREAS-pars.L_PANCREAS) + (pars.PLQ_SM_INT-pars.L_SM_INT) + (pars.PLQ_LG_INT-pars.L_LG_INT) + (pars.PLQ_SPLEEN-pars.L_SPLEEN);
Q_VOUT_LIVER       = (pars.PLQ_LIVER     - pars.L_LIVER    )  + PLQ_LIVER_UPSTREAM;
Q_VOUT_LUNG        = pars.PLQ_LUNG;
Q_VOUT_TUMOUR      = (pars.PLQ_TUMOUR     - pars.L_TUMOUR  );

%%scaling CL_UP accordingly to the original Shah 2012
CL_UP_HEART = pars.CL_UP*pars.V_E_HEART;
CL_UP_LUNG = pars.CL_UP*pars.V_E_LUNG;
CL_UP_MUSCLE = pars.CL_UP*pars.V_E_MUSCLE;
CL_UP_SKIN = pars.CL_UP*pars.V_E_SKIN;
CL_UP_ADIPOSE = pars.CL_UP*pars.V_E_ADIPOSE;
CL_UP_BONE = pars.CL_UP*pars.V_E_BONE;
CL_UP_BRAIN = pars.CL_UP*pars.V_E_BRAIN;
CL_UP_KIDNEY = pars.CL_UP*pars.V_E_KIDNEY;
CL_UP_LIVER = pars.CL_UP*pars.V_E_LIVER;
CL_UP_SM_INT = pars.CL_UP*pars.V_E_SM_INT;
CL_UP_LG_INT = pars.CL_UP*pars.V_E_LG_INT;
CL_UP_PANCREAS = pars.CL_UP*pars.V_E_PANCREAS;
CL_UP_THYMUS = pars.CL_UP*pars.V_E_THYMUS;
CL_UP_SPLEEN = pars.CL_UP*pars.V_E_SPLEEN;
CL_UP_OTHER = pars.CL_UP*pars.V_E_OTHER;
CL_UP_TUMOUR = pars.CL_UP*pars.V_E_TUMOUR;



%% Equations for two-pore clearance
%Two-pore clearance for large pores - eq 9 +13
CL_TP_L_HEART = pars_dim.PS_L_HEART * (C_V_HEART - C_I_HEART) * pars_dim.Pe_L_ratio + pars_dim.J_L_HEART * (1 - pars_dim.alpha_L)*C_V_HEART;
CL_TP_L_LUNG = pars_dim.PS_L_LUNG * (C_V_LUNG - C_I_LUNG) * pars_dim.Pe_L_ratio + pars_dim.J_L_LUNG * (1 - pars_dim.alpha_L)*C_V_LUNG;
CL_TP_L_MUSCLE = pars_dim.PS_L_MUSCLE * (C_V_MUSCLE - C_I_MUSCLE) * pars_dim.Pe_L_ratio + pars_dim.J_L_MUSCLE * (1 - pars_dim.alpha_L)*C_V_MUSCLE;
CL_TP_L_SKIN = pars_dim.PS_L_SKIN * (C_V_SKIN  - C_I_SKIN ) * pars_dim.Pe_L_ratio + pars_dim.J_L_SKIN * (1 - pars_dim.alpha_L)*C_V_SKIN ;
CL_TP_L_ADIPOSE = pars_dim.PS_L_ADIPOSE * (C_V_ADIPOSE - C_I_ADIPOSE) * pars_dim.Pe_L_ratio + pars_dim.J_L_ADIPOSE * (1 - pars_dim.alpha_L)*C_V_ADIPOSE;
CL_TP_L_BONE = pars_dim.PS_L_BONE * (C_V_BONE - C_I_BONE) * pars_dim.Pe_L_ratio + pars_dim.J_L_BONE * (1 - pars_dim.alpha_L)*C_V_BONE;
CL_TP_L_BRAIN = pars_dim.PS_L_BRAIN * (C_V_BRAIN - C_I_BRAIN) * pars_dim.Pe_L_ratio + pars_dim.J_L_BRAIN * (1 - pars_dim.alpha_L)*C_V_BRAIN;
CL_TP_L_KIDNEY = pars_dim.PS_L_KIDNEY * (C_V_KIDNEY - C_I_KIDNEY) * pars_dim.Pe_L_ratio + pars_dim.J_L_KIDNEY * (1 - pars_dim.alpha_L)*C_V_KIDNEY;
CL_TP_L_LIVER = pars_dim.PS_L_LIVER * (C_V_LIVER - C_I_LIVER) * pars_dim.Pe_L_ratio + pars_dim.J_L_LIVER * (1 - pars_dim.alpha_L)*C_V_LIVER;
CL_TP_L_SM_INT = pars_dim.PS_L_SM_INT *(C_V_SM_INT - C_I_SM_INT) * pars_dim.Pe_L_ratio + pars_dim.J_L_SM_INT * (1 - pars_dim.alpha_L)*C_V_SM_INT;
CL_TP_L_LG_INT = pars_dim.PS_L_LG_INT *(C_V_HEART - C_I_LG_INT) * pars_dim.Pe_L_ratio + pars_dim.J_L_LG_INT * (1 - pars_dim.alpha_L)*C_V_LG_INT;
CL_TP_L_PANCREAS = pars_dim.PS_L_PANCREAS * (C_V_PANCREAS - C_I_PANCREAS) * pars_dim.Pe_L_ratio + pars_dim.J_L_PANCREAS * (1 - pars_dim.alpha_L)*C_V_PANCREAS;
CL_TP_L_THYMUS = pars_dim.PS_L_THYMUS * (C_V_THYMUS - C_I_THYMUS) * pars_dim.Pe_L_ratio + pars_dim.J_L_THYMUS * (1 - pars_dim.alpha_L)*C_V_THYMUS;
CL_TP_L_SPLEEN = pars_dim.PS_L_SPLEEN * (C_V_SPLEEN - C_I_SPLEEN) * pars_dim.Pe_L_ratio + pars_dim.J_L_SPLEEN * (1 - pars_dim.alpha_L)*C_V_SPLEEN;
CL_TP_L_OTHER = pars_dim.PS_L_OTHER * (C_V_OTHER - C_I_OTHER) * pars_dim.Pe_L_ratio + pars_dim.J_L_OTHER * (1 - pars_dim.alpha_L)*C_V_OTHER;
CL_TP_L_TUMOUR = pars_dim.PS_L_TUMOUR * (C_V_TUMOUR - C_I_TUMOUR) * pars_dim.Pe_L_ratio + pars_dim.J_L_TUMOUR * (1 - pars_dim.alpha_L)*C_V_TUMOUR;

%Two-pore clearance for small pores - eq 10 + 14
CL_TP_S_HEART = pars_dim.PS_S_HEART * (C_V_HEART - C_I_HEART) * pars_dim.Pe_S_ratio + pars_dim.J_S_HEART * (1 - pars_dim.alpha_S)*C_V_HEART;
CL_TP_S_LUNG = pars_dim.PS_S_LUNG * (C_V_LUNG - C_I_LUNG) * pars_dim.Pe_S_ratio + pars_dim.J_S_LUNG * (1 - pars_dim.alpha_S)*C_V_LUNG;
CL_TP_S_MUSCLE = pars_dim.PS_S_MUSCLE * (C_V_MUSCLE - C_I_MUSCLE) * pars_dim.Pe_S_ratio + pars_dim.J_S_MUSCLE * (1 - pars_dim.alpha_S)*C_V_MUSCLE;
CL_TP_S_SKIN = pars_dim.PS_S_SKIN * (C_V_SKIN  - C_I_SKIN ) * pars_dim.Pe_S_ratio + pars_dim.J_S_SKIN * (1 - pars_dim.alpha_S)*C_V_SKIN ;
CL_TP_S_ADIPOSE = pars_dim.PS_S_ADIPOSE * (C_V_ADIPOSE - C_I_ADIPOSE) * pars_dim.Pe_S_ratio + pars_dim.J_S_ADIPOSE * (1 - pars_dim.alpha_S)*C_V_ADIPOSE;
CL_TP_S_BONE = pars_dim.PS_S_BONE * (C_V_BONE - C_I_BONE) * pars_dim.Pe_S_ratio + pars_dim.J_S_BONE * (1 - pars_dim.alpha_S)*C_V_BONE;
CL_TP_S_BRAIN = pars_dim.PS_S_BRAIN * (C_V_BRAIN - C_I_BRAIN) * pars_dim.Pe_S_ratio + pars_dim.J_S_BRAIN * (1 - pars_dim.alpha_S)*C_V_BRAIN;
CL_TP_S_KIDNEY = pars_dim.PS_S_KIDNEY * (C_V_KIDNEY - C_I_KIDNEY) * pars_dim.Pe_S_ratio + pars_dim.J_S_KIDNEY * (1 - pars_dim.alpha_S)*C_V_KIDNEY;
CL_TP_S_LIVER = pars_dim.PS_S_LIVER * (C_V_LIVER - C_I_LIVER) * pars_dim.Pe_S_ratio + pars_dim.J_S_LIVER * (1 - pars_dim.alpha_S)*C_V_LIVER;
CL_TP_S_SM_INT = pars_dim.PS_S_SM_INT *(C_V_SM_INT - C_I_SM_INT) * pars_dim.Pe_S_ratio + pars_dim.J_S_SM_INT * (1 - pars_dim.alpha_S)*C_V_SM_INT;
CL_TP_S_LG_INT = pars_dim.PS_S_LG_INT *(C_V_HEART - C_I_LG_INT) * pars_dim.Pe_S_ratio + pars_dim.J_S_LG_INT * (1 - pars_dim.alpha_S)*C_V_LG_INT;
CL_TP_S_PANCREAS = pars_dim.PS_S_PANCREAS * (C_V_PANCREAS - C_I_PANCREAS) * pars_dim.Pe_S_ratio + pars_dim.J_S_PANCREAS * (1 - pars_dim.alpha_S)*C_V_PANCREAS;
CL_TP_S_THYMUS = pars_dim.PS_S_THYMUS * (C_V_THYMUS - C_I_THYMUS) * pars_dim.Pe_S_ratio + pars_dim.J_S_THYMUS * (1 - pars_dim.alpha_S)*C_V_THYMUS;
CL_TP_S_SPLEEN = pars_dim.PS_S_SPLEEN * (C_V_SPLEEN - C_I_SPLEEN) * pars_dim.Pe_S_ratio + pars_dim.J_S_SPLEEN * (1 - pars_dim.alpha_S)*C_V_SPLEEN;
CL_TP_S_OTHER = pars_dim.PS_S_OTHER * (C_V_OTHER - C_I_OTHER) * pars_dim.Pe_S_ratio + pars_dim.J_S_OTHER * (1 - pars_dim.alpha_S)*C_V_OTHER;
CL_TP_S_TUMOUR = pars_dim.PS_S_TUMOUR * (C_V_TUMOUR - C_I_TUMOUR) * pars_dim.Pe_S_ratio + pars_dim.J_S_TUMOUR * (1 - pars_dim.alpha_S)*C_V_TUMOUR;

%total two-pore clearance
CL_TP_HEART = CL_TP_L_HEART + CL_TP_S_HEART;
CL_TP_LUNG = CL_TP_L_LUNG + CL_TP_S_LUNG;
CL_TP_MUSCLE = CL_TP_L_MUSCLE + CL_TP_S_MUSCLE;
CL_TP_SKIN = CL_TP_L_SKIN + CL_TP_S_SKIN;
CL_TP_ADIPOSE = CL_TP_L_ADIPOSE + CL_TP_S_ADIPOSE;
CL_TP_BONE = CL_TP_L_BONE + CL_TP_S_BONE;
CL_TP_BRAIN = CL_TP_L_BRAIN + CL_TP_S_BRAIN;
CL_TP_KIDNEY = CL_TP_L_KIDNEY + CL_TP_S_KIDNEY;
CL_TP_LIVER = CL_TP_L_LIVER + CL_TP_S_LIVER;
CL_TP_SM_INT = CL_TP_L_SM_INT + CL_TP_S_SM_INT;
CL_TP_LG_INT = CL_TP_L_LG_INT + CL_TP_S_LG_INT;
CL_TP_PANCREAS = CL_TP_L_PANCREAS + CL_TP_S_PANCREAS;
CL_TP_THYMUS = CL_TP_L_THYMUS + CL_TP_S_THYMUS;
CL_TP_SPLEEN = CL_TP_L_SPLEEN + CL_TP_S_SPLEEN;
CL_TP_OTHER = CL_TP_L_OTHER + CL_TP_S_OTHER;
CL_TP_TUMOUR = CL_TP_L_TUMOUR + CL_TP_S_TUMOUR;


%% ODEs for organs
%Vascular space - eq 3 in the paper
dydt(1) = (M_VI_HEART - Q_VOUT_HEART*C_V_HEART - CL_TP_HEART - CL_UP_HEART*C_V_HEART)/pars.V_V_LUNG;
dydt(2) = (M_VI_LUNG - Q_VOUT_LUNG*C_V_LUNG - CL_TP_LUNG - CL_UP_LUNG*C_V_LUNG)/pars.V_V_LUNG;
dydt(3) = (M_VI_MUSCLE - Q_VOUT_MUSCLE*C_V_MUSCLE - CL_TP_MUSCLE - CL_UP_MUSCLE*C_V_MUSCLE)/pars.V_V_MUSCLE;
dydt(4) = (M_VI_SKIN - Q_VOUT_SKIN*C_V_SKIN - CL_TP_SKIN - CL_UP_SKIN*C_V_SKIN)/pars.V_V_SKIN;
dydt(5) = (M_VI_ADIPOSE - Q_VOUT_ADIPOSE*C_V_ADIPOSE - CL_TP_ADIPOSE - CL_UP_ADIPOSE*C_V_ADIPOSE)/pars.V_V_ADIPOSE;
dydt(6) = (M_VI_BONE - Q_VOUT_BONE*C_V_BONE - CL_TP_BONE - CL_UP_BONE*C_V_BONE)/pars.V_V_BONE;
dydt(7) = (M_VI_BRAIN - Q_VOUT_BRAIN*C_V_BRAIN - CL_TP_BRAIN - CL_UP_BRAIN*C_V_BRAIN)/pars.V_V_BRAIN;
%renal clearance - eq 7
dydt(8) = (M_VI_KIDNEY - Q_VOUT_KIDNEY*C_V_KIDNEY - CL_TP_KIDNEY - (CL_UP_KIDNEY + pars_dim.CL_R)*C_V_KIDNEY)/pars.V_V_KIDNEY;
dydt(9) = (M_VI_LIVER - Q_VOUT_LIVER*C_V_LIVER - CL_TP_LIVER - CL_UP_LIVER*C_V_LIVER)/pars.V_V_LIVER;
dydt(10) = (M_VI_SM_INT - Q_VOUT_SM_INT*C_V_SM_INT - CL_TP_SM_INT - CL_UP_SM_INT*C_V_SM_INT)/pars.V_V_SM_INT;
dydt(11) = (M_VI_LG_INT - Q_VOUT_LG_INT*C_V_LG_INT - CL_TP_LG_INT - CL_UP_LG_INT*C_V_LG_INT)/pars.V_V_LG_INT;
dydt(12) = (M_VI_PANCREAS - Q_VOUT_PANCREAS*C_V_PANCREAS - CL_TP_PANCREAS - CL_UP_PANCREAS*C_V_PANCREAS)/pars.V_V_PANCREAS;
dydt(13) = (M_VI_THYMUS - Q_VOUT_THYMUS*C_V_THYMUS - CL_TP_THYMUS - CL_UP_THYMUS*C_V_THYMUS)/pars.V_V_THYMUS;
dydt(14) = (M_VI_SPLEEN - Q_VOUT_SPLEEN*C_V_SPLEEN - CL_TP_SPLEEN - CL_UP_SPLEEN*C_V_SPLEEN)/pars.V_V_SPLEEN;
dydt(15) = (M_VI_OTHER - Q_VOUT_OTHER*C_V_OTHER - CL_TP_OTHER - CL_UP_OTHER*C_V_OTHER)/pars.V_V_OTHER;


%endosomal space - eq 4 in the paper
dydt(16) = (C_V_HEART + C_I_HEART)*CL_UP_HEART/pars.V_E_HEART - pars.K_deg*C_E_HEART;
dydt(17) = (C_V_LUNG + C_I_LUNG)*CL_UP_LUNG/pars.V_E_LUNG - pars.K_deg*C_E_LUNG;
dydt(18) = (C_V_MUSCLE + C_I_MUSCLE)*CL_UP_MUSCLE/pars.V_E_MUSCLE - pars.K_deg*C_E_MUSCLE;
dydt(19) = (C_V_SKIN + C_I_SKIN)*CL_UP_SKIN/pars.V_E_SKIN - pars.K_deg*C_E_SKIN;
dydt(20) = (C_V_ADIPOSE + C_I_ADIPOSE)*CL_UP_ADIPOSE/pars.V_E_ADIPOSE - pars.K_deg*C_E_ADIPOSE;
dydt(21) = (C_V_BONE + C_I_BONE)*CL_UP_BONE/pars.V_E_BONE - pars.K_deg*C_E_BONE;
dydt(22) = (C_V_BRAIN + C_I_BRAIN)*CL_UP_BRAIN/pars.V_E_BRAIN - pars.K_deg*C_E_BRAIN;
dydt(23) = (C_V_KIDNEY + C_I_KIDNEY)*CL_UP_KIDNEY/pars.V_E_KIDNEY - pars.K_deg*C_E_KIDNEY;
dydt(24) = (C_V_LIVER + C_I_LIVER)*CL_UP_LIVER/pars.V_E_LIVER - pars.K_deg*C_E_LIVER;
dydt(25) = (C_V_SM_INT + C_I_SM_INT)*CL_UP_SM_INT/pars.V_E_SM_INT - pars.K_deg*C_E_SM_INT;
dydt(26) = (C_V_LG_INT + C_I_LG_INT)*CL_UP_LG_INT/pars.V_E_LG_INT - pars.K_deg*C_E_LG_INT;
dydt(27) = (C_V_PANCREAS + C_I_PANCREAS)*CL_UP_PANCREAS/pars.V_E_PANCREAS - pars.K_deg*C_E_PANCREAS;
dydt(28) = (C_V_THYMUS + C_I_THYMUS)*CL_UP_THYMUS/pars.V_E_THYMUS - pars.K_deg*C_E_THYMUS;
dydt(29) = (C_V_SPLEEN + C_I_SPLEEN)*CL_UP_SPLEEN/pars.V_E_SPLEEN - pars.K_deg*C_E_SPLEEN;
dydt(30) = (C_V_OTHER + C_I_OTHER)*CL_UP_OTHER/pars.V_E_OTHER - pars.K_deg*C_E_OTHER;

% %interstitial space - eq 5 in the paper
% dydt(31) =   CL_TP_HEART*C_V_HEART/pars.V_I_HEART - (1.0-pars.SIGMA_I_HEART)*pars.L_HEART*C_I_HEART/pars.V_I_HEART - CL_UP_HEART*(C_I_HEART + C_V_HEART)/pars.V_I_HEART ;
% dydt(32) =   CL_TP_LUNG*C_V_LUNG/pars.V_I_LUNG - (1.0-pars.SIGMA_I_LUNG)*pars.L_LUNG*C_I_LUNG/pars.V_I_LUNG - CL_UP_LUNG*(C_I_LUNG + C_V_LUNG)/pars.V_I_LUNG ;
% dydt(33) =   CL_TP_MUSCLE*C_V_MUSCLE/pars.V_I_MUSCLE - (1.0-pars.SIGMA_I_MUSCLE)*pars.L_MUSCLE*C_I_MUSCLE/pars.V_I_MUSCLE - CL_UP_MUSCLE*(C_I_MUSCLE + C_V_MUSCLE)/pars.V_I_MUSCLE ;
% dydt(34) =   CL_TP_SKIN*C_V_SKIN/pars.V_I_SKIN - (1.0-pars.SIGMA_I_SKIN)*pars.L_SKIN*C_I_SKIN/pars.V_I_SKIN - CL_UP_SKIN*(C_I_SKIN + C_V_SKIN)/pars.V_I_SKIN ;
% dydt(35) =   CL_TP_ADIPOSE*C_V_ADIPOSE/pars.V_I_ADIPOSE - (1.0-pars.SIGMA_I_ADIPOSE)*pars.L_ADIPOSE*C_I_ADIPOSE/pars.V_I_ADIPOSE - CL_UP_ADIPOSE*(C_I_ADIPOSE + C_V_ADIPOSE)/pars.V_I_ADIPOSE ;
% dydt(36) =   CL_TP_BONE*C_V_BONE/pars.V_I_BONE - (1.0-pars.SIGMA_I_BONE)*pars.L_BONE*C_I_BONE/pars.V_I_BONE - CL_UP_BONE*(C_I_BONE + C_V_BONE)/pars.V_I_BONE ;
% dydt(37) =   CL_TP_BRAIN*C_V_BRAIN/pars.V_I_BRAIN - (1.0-pars.SIGMA_I_BRAIN)*pars.L_BRAIN*C_I_BRAIN/pars.V_I_BRAIN - CL_UP_BRAIN*(C_I_BRAIN + C_V_BRAIN)/pars.V_I_BRAIN ;
% dydt(38) =   CL_TP_KIDNEY*C_V_KIDNEY/pars.V_I_KIDNEY - (1.0-pars.SIGMA_I_KIDNEY)*pars.L_KIDNEY*C_I_KIDNEY/pars.V_I_KIDNEY - CL_UP_KIDNEY*(C_I_KIDNEY + C_V_KIDNEY)/pars.V_I_KIDNEY ;
% dydt(39) =   CL_TP_LIVER*C_V_LIVER/pars.V_I_LIVER - (1.0-pars.SIGMA_I_LIVER)*pars.L_LIVER*C_I_LIVER/pars.V_I_LIVER - CL_UP_LIVER*(C_I_LIVER + C_V_LIVER)/pars.V_I_LIVER ;
% dydt(40) =   CL_TP_SM_INT*C_V_SM_INT/pars.V_I_SM_INT - (1.0-pars.SIGMA_I_SM_INT)*pars.L_SM_INT*C_I_SM_INT/pars.V_I_SM_INT - CL_UP_SM_INT*(C_I_SM_INT + C_V_SM_INT)/pars.V_I_SM_INT ;
% dydt(41) =   CL_TP_LG_INT*C_V_LG_INT/pars.V_I_LG_INT - (1.0-pars.SIGMA_I_LG_INT)*pars.L_LG_INT*C_I_LG_INT/pars.V_I_LG_INT - CL_UP_LG_INT*(C_I_LG_INT + C_V_LG_INT)/pars.V_I_LG_INT ;
% dydt(42) =   CL_TP_PANCREAS*C_V_PANCREAS/pars.V_I_PANCREAS - (1.0-pars.SIGMA_I_PANCREAS)*pars.L_PANCREAS*C_I_PANCREAS/pars.V_I_PANCREAS - CL_UP_PANCREAS*(C_I_PANCREAS + C_V_PANCREAS)/pars.V_I_PANCREAS ;
% dydt(43) =   CL_TP_THYMUS*C_V_THYMUS/pars.V_I_THYMUS - (1.0-pars.SIGMA_I_THYMUS)*pars.L_THYMUS*C_I_THYMUS/pars.V_I_THYMUS - CL_UP_THYMUS*(C_I_THYMUS + C_V_THYMUS)/pars.V_I_THYMUS ;
% dydt(44) =   CL_TP_SPLEEN*C_V_SPLEEN/pars.V_I_SPLEEN - (1.0-pars.SIGMA_I_SPLEEN)*pars.L_SPLEEN*C_I_SPLEEN/pars.V_I_SPLEEN - CL_UP_SPLEEN*(C_I_SPLEEN + C_V_SPLEEN)/pars.V_I_SPLEEN ;
% dydt(45) =   CL_TP_OTHER*C_V_OTHER/pars.V_I_OTHER - (1.0-pars.SIGMA_I_OTHER)*pars.L_OTHER*C_I_OTHER/pars.V_I_OTHER - CL_UP_OTHER*(C_I_OTHER + C_V_OTHER)/pars.V_I_OTHER ;

%interstitial space - eq 5 in the paper
dydt(31) =   CL_TP_HEART/pars.V_I_HEART - (1.0-pars.SIGMA_I_HEART)*pars.L_HEART*C_I_HEART/pars.V_I_HEART - CL_UP_HEART*(C_I_HEART)/pars.V_I_HEART ;
dydt(32) =   CL_TP_LUNG/pars.V_I_LUNG - (1.0-pars.SIGMA_I_LUNG)*pars.L_LUNG*C_I_LUNG/pars.V_I_LUNG - CL_UP_LUNG*(C_I_LUNG)/pars.V_I_LUNG ;
dydt(33) =   CL_TP_MUSCLE/pars.V_I_MUSCLE - (1.0-pars.SIGMA_I_MUSCLE)*pars.L_MUSCLE*C_I_MUSCLE/pars.V_I_MUSCLE - CL_UP_MUSCLE*(C_I_MUSCLE)/pars.V_I_MUSCLE ;
dydt(34) =   CL_TP_SKIN/pars.V_I_SKIN - (1.0-pars.SIGMA_I_SKIN)*pars.L_SKIN*C_I_SKIN/pars.V_I_SKIN - CL_UP_SKIN*(C_I_SKIN)/pars.V_I_SKIN ;
dydt(35) =   CL_TP_ADIPOSE/pars.V_I_ADIPOSE - (1.0-pars.SIGMA_I_ADIPOSE)*pars.L_ADIPOSE*C_I_ADIPOSE/pars.V_I_ADIPOSE - CL_UP_ADIPOSE*(C_I_ADIPOSE)/pars.V_I_ADIPOSE ;
dydt(36) =   CL_TP_BONE/pars.V_I_BONE - (1.0-pars.SIGMA_I_BONE)*pars.L_BONE*C_I_BONE/pars.V_I_BONE - CL_UP_BONE*(C_I_BONE)/pars.V_I_BONE ;
dydt(37) =   CL_TP_BRAIN/pars.V_I_BRAIN - (1.0-pars.SIGMA_I_BRAIN)*pars.L_BRAIN*C_I_BRAIN/pars.V_I_BRAIN - CL_UP_BRAIN*(C_I_BRAIN)/pars.V_I_BRAIN ;
dydt(38) =   CL_TP_KIDNEY/pars.V_I_KIDNEY - (1.0-pars.SIGMA_I_KIDNEY)*pars.L_KIDNEY*C_I_KIDNEY/pars.V_I_KIDNEY - CL_UP_KIDNEY*(C_I_KIDNEY)/pars.V_I_KIDNEY ;
dydt(39) =   CL_TP_LIVER/pars.V_I_LIVER - (1.0-pars.SIGMA_I_LIVER)*pars.L_LIVER*C_I_LIVER/pars.V_I_LIVER - CL_UP_LIVER*(C_I_LIVER)/pars.V_I_LIVER ;
dydt(40) =   CL_TP_SM_INT/pars.V_I_SM_INT - (1.0-pars.SIGMA_I_SM_INT)*pars.L_SM_INT*C_I_SM_INT/pars.V_I_SM_INT - CL_UP_SM_INT*(C_I_SM_INT)/pars.V_I_SM_INT ;
dydt(41) =   CL_TP_LG_INT/pars.V_I_LG_INT - (1.0-pars.SIGMA_I_LG_INT)*pars.L_LG_INT*C_I_LG_INT/pars.V_I_LG_INT - CL_UP_LG_INT*(C_I_LG_INT)/pars.V_I_LG_INT ;
dydt(42) =   CL_TP_PANCREAS/pars.V_I_PANCREAS - (1.0-pars.SIGMA_I_PANCREAS)*pars.L_PANCREAS*C_I_PANCREAS/pars.V_I_PANCREAS - CL_UP_PANCREAS*(C_I_PANCREAS)/pars.V_I_PANCREAS ;
dydt(43) =   CL_TP_THYMUS/pars.V_I_THYMUS - (1.0-pars.SIGMA_I_THYMUS)*pars.L_THYMUS*C_I_THYMUS/pars.V_I_THYMUS - CL_UP_THYMUS*(C_I_THYMUS)/pars.V_I_THYMUS ;
dydt(44) =   CL_TP_SPLEEN/pars.V_I_SPLEEN - (1.0-pars.SIGMA_I_SPLEEN)*pars.L_SPLEEN*C_I_SPLEEN/pars.V_I_SPLEEN - CL_UP_SPLEEN*(C_I_SPLEEN)/pars.V_I_SPLEEN ;
dydt(45) =   CL_TP_OTHER/pars.V_I_OTHER - (1.0-pars.SIGMA_I_OTHER)*pars.L_OTHER*C_I_OTHER/pars.V_I_OTHER - CL_UP_OTHER*(C_I_OTHER)/pars.V_I_OTHER ;

%% ODEs for fluids

%Plasma - eq 1
dydt(46) =  (Q_VOUT_HEART*C_V_HEART/pars.V_TOT_PLASMA + Q_VOUT_KIDNEY*C_V_KIDNEY/pars.V_TOT_PLASMA + Q_VOUT_MUSCLE*C_V_MUSCLE/pars.V_TOT_PLASMA + Q_VOUT_SKIN*C_V_SKIN/pars.V_TOT_PLASMA + Q_VOUT_BRAIN*C_V_BRAIN/pars.V_TOT_PLASMA + Q_VOUT_ADIPOSE*C_V_ADIPOSE/pars.V_TOT_PLASMA + Q_VOUT_THYMUS*C_V_THYMUS/pars.V_TOT_PLASMA + Q_VOUT_LIVER*C_V_LIVER/pars.V_TOT_PLASMA + Q_VOUT_BONE*C_V_BONE/pars.V_TOT_PLASMA + Q_VOUT_OTHER*C_V_OTHER/pars.V_TOT_PLASMA + Q_VOUT_TUMOUR*C_V_TUMOUR/pars.V_TOT_PLASMA + pars.L_LYMPH*C_LYMPH/pars.V_TOT_PLASMA - (pars.PLQ_LUNG+pars.L_LUNG)*C_PLASMA/pars.V_TOT_PLASMA );

%Lymphnode - eq 2
dydt(47) =   (1.0-pars.SIGMA_I_HEART)*pars.L_HEART*C_I_HEART/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_KIDNEY)*pars.L_KIDNEY*C_I_KIDNEY/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_MUSCLE)*pars.L_MUSCLE*C_I_MUSCLE/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_SKIN)*pars.L_SKIN*C_I_SKIN/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_BRAIN)*pars.L_BRAIN*C_I_BRAIN/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_ADIPOSE)*pars.L_ADIPOSE*C_I_ADIPOSE/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_THYMUS)*pars.L_THYMUS*C_I_THYMUS/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_LIVER)*pars.L_LIVER*C_I_LIVER/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_SPLEEN)*pars.L_SPLEEN*C_I_SPLEEN/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_PANCREAS)*pars.L_PANCREAS*C_I_PANCREAS/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_SM_INT)*pars.L_SM_INT*C_I_SM_INT/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_LG_INT)*pars.L_LG_INT*C_I_LG_INT/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_BONE)*pars.L_BONE*C_I_BONE/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_OTHER)*pars.L_OTHER*C_I_OTHER/pars.V_TOT_LYMPH + (1.0-pars.SIGMA_I_LUNG)*pars.L_LUNG*C_I_LUNG/pars.V_TOT_LYMPH - pars.L_LYMPH*C_LYMPH/pars.V_TOT_LYMPH;


%% tumour compartment
dydt(48) = (M_VI_TUMOUR - Q_VOUT_TUMOUR*C_V_TUMOUR - CL_TP_TUMOUR - CL_UP_TUMOUR*C_V_TUMOUR)/pars.V_V_TUMOUR;
dydt(49) = (C_V_TUMOUR + C_I_TUMOUR)*CL_UP_TUMOUR/pars.V_E_TUMOUR - pars.K_deg*C_E_TUMOUR;
dydt(50) =   CL_TP_TUMOUR/pars.V_I_TUMOUR - (1.0-pars.SIGMA_I_TUMOUR)*pars.L_TUMOUR*C_I_TUMOUR/pars.V_I_TUMOUR - CL_UP_TUMOUR*(C_I_TUMOUR)/pars.V_I_TUMOUR ;
end