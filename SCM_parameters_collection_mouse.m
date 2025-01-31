function pars = SCM_parameters_collection_mouse(par_fit,C_LNLF,K_deg,K_deg_LNP)
    
    %% ***************** Generic species specific parameters (Li_2019) *************************************
    % Body weight animal model (mouse)
    pars.BW = 0.028; % Unit measure: Kg
    % Lymph reflection coefficient
    pars.SIGMA_I = 0.2; % Unit measure: pure number
    % Glomerular filtration rate in mouse
    pars.GFR = 0.0167; %0.278 mL/min ---> 0.0167 L/h  % Unit measure: L/h

    %% ***************** Product specific parameters ***********************************
    % Rate of pinocytosis per unit endosomal space
    pars.CL_UP = par_fit(1); % Unit measure: 1/h
    % Proportionality constant between the rate at which antibody transfers 
    % from the lymph compartment to the blood compartment and the plasma flow of the given species
    pars.C_LNLF = C_LNLF; % Unit measure: pure number (Shah_2012)
    % Lysosome degradation rate for all the proteins (Model for recombinant
    % proteins)
    pars.K_deg = K_deg; % Unit measure: 1/h (Li_2019)
    % Lysosome degradation rate for all the proteins (Model for mRNA-encoded
    % proteins)
    pars.K_deg_LNP = K_deg_LNP; % Unit measure: 1/h (Li_2019)
    

    %% ++++++++++++++++ Body parameters ++++++++++++++++++++++++++++++++++++++++++++++++

    %% ----------------- VOLUMES -----------------------------------------

    %% ***************** Species specific total volumes of the compartments (Unit measure: L) ***************************
    % Heart
    pars.V_TOT_HEART = 0.152*10^(-3); % Unit measure: L
    % Lungs
    pars.V_TOT_LUNG = 0.204*10^(-3);  % Unit measure: L
    % Muscles
    pars.V_TOT_MUSCLE = 11.3*10^(-3); % Unit measure: L
    % Skin
    pars.V_TOT_SKIN = 5.02*10^(-3); % Unit measure: L
    % Adipose tissue
    pars.V_TOT_ADIPOSE = 1.98*10^(-3); % Unit measure: L
    % Bones
    pars.V_TOT_BONE = 2.82*10^(-3); % Unit measure: L 
    % Brain
    pars.V_TOT_BRAIN = 0.485*10^(-3); % Unit measure: L 
    % Kidneys
    pars.V_TOT_KIDNEY = 0.525*10^(-3); % Unit measure: L
    % Liver
    pars.V_TOT_LIVER = 1.93*10^(-3);  % Unit measure: L  
    % Small intestine
    pars.V_TOT_SM_INT = 0.728*10^(-3); % Unit measure: L 
    % Large intestine
    pars.V_TOT_LG_INT = 0.314*10^(-3); % Unit measure: L
    % Pancreas
    pars.V_TOT_PANCREAS = 0.0970*10^(-3);  % Unit measure: L
    % Thymus
    pars.V_TOT_THYMUS = 0.009*10^(-3);  % Unit measure: L   
    % Spleen 
    pars.V_TOT_SPLEEN = 0.127*10^(-3); % Unit measure: L
    % All the other organs not represented
    pars.V_TOT_OTHER = 0.465*10^(-3); % Unit measure: L
    % Lymphatuc system
    pars.V_TOT_LYMPH = 0.000113; % Unit measure: L
    % Plasma
    pars.V_TOT_PLASMA = 0.00094435; % Unit measure: L

    %% *************** Species specific vascular volumes of the compartments (Unit measure: L) *********************************
    pars.V_V_HEART = 0.000005852;   % Unit measure: L           
    pars.V_V_LUNG = 0.00002945509193;    % Unit measure: L        
    pars.V_V_MUSCLE = 0.000249018; % Unit measure: L             
    pars.V_V_SKIN = 0.0001877106;  % Unit measure: L              
    pars.V_V_ADIPOSE = 0.000021802; % Unit measure: L            
    pars.V_V_BONE = 0.000062128;   % Unit measure: L              
    pars.V_V_BRAIN = 0.00001067;  % Unit measure: L                
    pars.V_V_KIDNEY = 0.000028875; % Unit measure: L         
    pars.V_V_LIVER = 0.00016410625;  % Unit measure: L            
    pars.V_V_SM_INT = 0.0000116116;  % Unit measure: L
    pars.V_V_LG_INT = 0.000005000325;   % Unit measure: L  
    pars.V_V_PANCREAS = 0.000005335;   % Unit measure: L           
    pars.V_V_THYMUS = 0.000000495;   % Unit measure: L            
    pars.V_V_SPLEEN = 0.000015367;   % Unit measure: L          
    pars.V_V_OTHER = 0.000019535609097;  % Unit measure: L   

    %% **************** Species specific interstitial volumes of the compartments (Unit measure: L) ************************************                                     
    pars.V_I_HEART = 0.000021736;     % Unit measure: L           
    pars.V_I_LUNG = 0.0000384285724; % Unit measure: L          
    pars.V_I_MUSCLE = 0.00147147;    % Unit measure: L            
    pars.V_I_SKIN = 0.00165627;     % Unit measure: L             
    pars.V_I_ADIPOSE = 0.00033694;   % Unit measure: L           
    pars.V_I_BONE = 0.000525264;     % Unit measure: L          
    pars.V_I_BRAIN = 0.0000873;      % Unit measure: L           
    pars.V_I_KIDNEY = 0.00007875;    % Unit measure: L           
    pars.V_I_LIVER = 0.000384595163;   % Unit measure: L         
    pars.V_I_SM_INT = 0.000126672;     % Unit measure: L        
    pars.V_I_LG_INT = 0.000054549;     % Unit measure: L          
    pars.V_I_PANCREAS = 0.000016878;    % Unit measure: L         
    pars.V_I_THYMUS = 0.00000153;      % Unit measure: L        
    pars.V_I_SPLEEN = 0.0000254;      % Unit measure: L          
    pars.V_I_OTHER = 0.0000797245183;  % Unit measure: L         
                                                                   
    %% ************************ Species specific endosomal volumes of the compartments (Unit measure: L)*************************                                         
    pars.V_E_HEART = 0.00000076;  % Unit measure: L             
    pars.V_E_LUNG = 0.0000010220365;    % Unit measure: L        
    pars.V_E_MUSCLE = 0.000056595;    % Unit measure: L        
    pars.V_E_SKIN = 0.000025095;     % Unit measure: L          
    pars.V_E_ADIPOSE = 0.00000991;   % Unit measure: L           
    pars.V_E_BONE = 0.00001412;     % Unit measure: L            
    pars.V_E_BRAIN = 0.000002425;   % Unit measure: L           
    pars.V_E_KIDNEY = 0.000002625;   % Unit measure: L          
    pars.V_E_LIVER = 0.000009625;    % Unit measure: L          
    pars.V_E_SM_INT = 0.00000364;    % Unit measure: L         
    pars.V_E_LG_INT = 0.0000015675;  % Unit measure: L           
    pars.V_E_PANCREAS = 0.000000485;  % Unit measure: L         
    pars.V_E_THYMUS = 0.000000045;    % Unit measure: L       
    pars.V_E_SPLEEN = 0.000000635;    % Unit measure: L       
    pars.V_E_OTHER = 0.000002326591;   % Unit measure: L
          
     %% ************** Tumor volumes (Unit measure: L) **************************************************
    pars.V_TUMOUR = 0.472*10^(-3); % Unit measure: L
    pars.V_V_TUMOUR = 0.07*pars.V_TUMOUR; % Unit measure: L
    pars.V_E_TUMOUR = 0.005*pars.V_TUMOUR; % Unit measure: L
    pars.V_I_TUMOUR = 0.55*pars.V_TUMOUR;% Unit measure: L
                                                       
  
    %% ++++++++++++++++ Tissue absorpition dynamics parameters (Protein Transportation Module) ++++++++++++++++++++++++++++++++++++++++++++++++

    % Arterial blood to each tissue compartment is delivered by the efferent
    % blood supply from the lung and venous return from most tissues
    % except small intestine, large intestine, spleen, pancreas, is
    % delivered to the blood compartment, which represents a venous pool. 
    % Venous return from small intestine, large intestine, spleen and pancreas is delivered to the liver.
    % The flow circuitry is completed by delivering the blood from the
    % blood (venous pool) compartment to the lung compartment.

    %% ********************** PLQ = Plasma Flower Rate to the tissues *****************************************                   
    pars.PLQ_HEART        =0.036498;  % Unit measure: L/h  
    pars.PLQ_MUSCLE       =0.08613;   % Unit measure: L/h  
    pars.PLQ_SKIN         =0.027819;  % Unit measure: L/h    
    pars.PLQ_ADIPOSE      =0.013431;  % Unit measure: L/h  
    pars.PLQ_BONE         =0.01518;   % Unit measure: L/h    
    pars.PLQ_BRAIN        =0.011781;  % Unit measure: L/h    
    pars.PLQ_KIDNEY       =0.068508;  % Unit measure: L/h   
    pars.PLQ_LIVER        =0.010263;  % Unit measure: L/h    
    pars.PLQ_SM_INT       =0.05808;   % Unit measure: L/h  
    pars.PLQ_LG_INT       =0.017292;  % Unit measure: L/h  
    pars.PLQ_PANCREAS     =0.006237;  % Unit measure: L/h 
    pars.PLQ_THYMUS       =0.001188;  % Unit measure: L/h   
    pars.PLQ_SPLEEN       =0.008184;  % Unit measure: L/h  
    pars.PLQ_OTHER        =0.01254;   % Unit measure: L/h 
    pars.PLQ_TUMOUR       =12.7*10^(-3); % Unit measure: L/h 

    % Conservation law for fluxes to compute the lung flux according with
    % the physiology
    pars.PLQ_LUNG =(pars.PLQ_HEART+ pars.PLQ_MUSCLE+pars.PLQ_SKIN+pars.PLQ_ADIPOSE+pars.PLQ_BONE+pars.PLQ_BRAIN+pars.PLQ_KIDNEY+pars.PLQ_LIVER+pars.PLQ_SM_INT+pars.PLQ_LG_INT+pars.PLQ_PANCREAS+pars.PLQ_THYMUS+pars.PLQ_SPLEEN+pars.PLQ_OTHER+pars.PLQ_TUMOUR);
   
    %% ********************** L = Lymph flow from the tissue *****************************************    
    % Lymph flow from all the tissues is delivered to the blood compartment via a lymph node compartment.
    % Lympathic flow is 500 times slower than the plasma one i.e. 1/500=0.002)
    pars.L_HEART     = 0.002*pars.PLQ_HEART;   % Unit measure: L/h 
    pars.L_MUSCLE    = 0.002*pars.PLQ_MUSCLE;  % Unit measure: L/h 
    pars.L_SKIN      = 0.002*pars.PLQ_SKIN;    % Unit measure: L/h 
    pars.L_ADIPOSE   = 0.002*pars.PLQ_ADIPOSE ; % Unit measure: L/h 
    pars.L_BONE      = 0.002*pars.PLQ_BONE;   % Unit measure: L/h  
    pars.L_BRAIN     = 0.002*pars.PLQ_BRAIN ;  % Unit measure: L/h 
    pars.L_KIDNEY    = 0.002*pars.PLQ_KIDNEY ; % Unit measure: L/h 
    pars.L_SM_INT    = 0.002*pars.PLQ_SM_INT;  % Unit measure: L/h 
    pars.L_LG_INT    = 0.002*pars.PLQ_LG_INT;  % Unit measure: L/h 
    pars.L_PANCREAS  = 0.002*pars.PLQ_PANCREAS; % Unit measure: L/h 
    pars.L_THYMUS    = 0.002*pars.PLQ_THYMUS;  % Unit measure: L/h 
    pars.L_SPLEEN    = 0.002*pars.PLQ_SPLEEN ; % Unit measure: L/h 
    pars.L_OTHER     = 0.002*pars.PLQ_OTHER ; % Unit measure: L/h 
    pars.L_LUNG      = 0.002*pars.PLQ_LUNG;   % Unit measure: L/h 
    pars.L_TUMOUR    = 0.002*pars.PLQ_TUMOUR; % Unit measure: L/h

    % Special computation for liver according with the physiology
    % represented: the venous return from small intestine, large intestine,
    % spleen and pancreas is delivered to the liver contributing to the
    % dynamics
    pars.PLQ_LIVER_UPSTREAM = (pars.PLQ_PANCREAS-pars.L_PANCREAS) + (pars.PLQ_SM_INT-pars.L_SM_INT) + (pars.PLQ_LG_INT-pars.L_LG_INT) + (pars.PLQ_SPLEEN-pars.L_SPLEEN);
    pars.L_LIVER     = 0.002*(pars.PLQ_LIVER + pars.PLQ_LIVER_UPSTREAM);
    
    %% ***************** Lymphatic system dynamics **********************************************
    pars.L_LYMPH = pars.PLQ_LUNG*pars.C_LNLF;% Unit measure: L/h 

   %% ****************  Pinocytosis rate for each tissue *********************************************
    % General formula: CL_UP_{ORG} = CL_UP*V_E_{ORG};
    pars.CL_UP_HEART = pars.CL_UP*pars.V_E_HEART; % Unit measure: L/h
    pars.CL_UP_LUNG = pars.CL_UP*pars.V_E_LUNG; % Unit measure: L/h
    pars.CL_UP_MUSCLE = pars.CL_UP*pars.V_E_MUSCLE; % Unit measure: L/h
    pars.CL_UP_SKIN = pars.CL_UP*pars.V_E_SKIN; % Unit measure: L/h
    pars.CL_UP_ADIPOSE = pars.CL_UP*pars.V_E_ADIPOSE; % Unit measure: L/h
    pars.CL_UP_BONE = pars.CL_UP*pars.V_E_BONE; % Unit measure: L/h
    pars.CL_UP_BRAIN = pars.CL_UP*pars.V_E_BRAIN; % Unit measure: L/h
    pars.CL_UP_KIDNEY = pars.CL_UP*pars.V_E_KIDNEY; % Unit measure: L/h
    pars.CL_UP_LIVER = pars.CL_UP*pars.V_E_LIVER; % Unit measure: L/h
    pars.CL_UP_SM_INT = pars.CL_UP*pars.V_E_SM_INT; % Unit measure: L/h
    pars.CL_UP_LG_INT = pars.CL_UP*pars.V_E_LG_INT; % Unit measure: L/h
    pars.CL_UP_PANCREAS = pars.CL_UP*pars.V_E_PANCREAS; % Unit measure: L/h
    pars.CL_UP_THYMUS = pars.CL_UP*pars.V_E_THYMUS; % Unit measure: L/h
    pars.CL_UP_SPLEEN = pars.CL_UP*pars.V_E_SPLEEN; % Unit measure: L/h
    pars.CL_UP_OTHER = pars.CL_UP*pars.V_E_OTHER; % Unit measure: L/h
    pars.CL_UP_TUMOUR = pars.CL_UP*pars.V_E_TUMOUR; % Unit measure: L/h

    %% ******************** Reflection interstitial coefficient for each organ *****************************************
    % General formula: SIGMA_I_{ORG} = SIGMA_I;
    pars.SIGMA_I_HEART = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_LUNG = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_MUSCLE = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_SKIN = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_ADIPOSE = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_BONE = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_BRAIN = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_KIDNEY = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_LIVER = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_SM_INT = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_LG_INT = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_PANCREAS = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_THYMUS = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_SPLEEN = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_OTHER = pars.SIGMA_I; % Unit measure: pure number
    pars.SIGMA_I_TUMOUR = pars.SIGMA_I; % Unit measure: pure number

    %% ++++++++++++++++ mRNA transport dynamics parameters (mRNA Transportation Module) ++++++++++++++++++++++++++++++++++++++++++++++++
    % Ratio of the mRNA injected that reach all the other organs at
    % exception of the liver - estimated using the data showed in Figure 5c
    % in Huang_2023
    pars.ratio_in_other =0.0110; % Unit measure: pure number
    % Transfer rate of the mRNA injected from blood to the liver:
    pars.k_BL2LV = par_fit(2); % Unit measure: 1/h
    % Degradation rate of the mRNA: k_d_mRNA
    pars.k_d_mRNA = par_fit(3); % Unit measure: 1/h
    % Rates at which BiTEs are produced via translation in each mRNA phase
    pars.k_tr_mRNA_LV_1 = par_fit(4); % Unit measure: 1/h
    pars.k_tr_mRNA_LV_2 = par_fit(5); % Unit measure: 1/h
    % Phenomenological transition rates between mRNA phases
    pars.k_s_12 = par_fit(6); % from phase 1 to phase 2  % Unit measure: 1/h
end