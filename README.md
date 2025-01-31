# A Multi-Scale Physiologically Based Pharmacokinetics Model to Support mRNA-Encoded BiTE Therapy in Cancer Treatment

**Authors:** Giada Fiandaca, Elio Campanile, Lorena Leonardelli, Elisa Pettinà, Stefano Giampiccolo, Elizabeth J Carstens, Lorenzo Dasti, Natascia Zangani, Luca Marchetti

This repository contains the complete MATLAB source codes for reproducing the dynamics of mRNA trafficking, translation, and the subsequent production of functional bispecific antibodies. It also includes their trafficking in both of the two modeling approaches proposed in the paper: the **Short Chain Model (SCM)** and the **Long Chain Model (LCM)**.

---

## 📂 CONTENTS

### **Short Chain Model (SCM)**

- `Parameters_collection_BiTEs_55kDa.m`: Function containing the parameters regulating the transport of 55 kDa proteins (BiTE dimension is 55 kDa).
- `recombinant_model.m`: Function containing the equations of the model for recombinant protein transport.

- `SCM_simulation.m`: Source code for the model simulation for both recombinant and mRNA-encoded protein administration in the Short Chain Model. Fixing the dose injected and the observation time window allows computing the time evolution of BiTEs' concentration in blood.
- `SCM.m`: Function containing the equations of the Short Chain Model.
- `SCM_parameters_collection_mouse.m`: Function containing parameters related to mouse physiology, including species-specific protein transport, body parameters, tissue absorption dynamics, and mRNA transport dynamics for the Short Chain Model.
- `SCM_parameters_estimated.mat`: `.mat` file containing estimated parameter values (e.g., rate of pinocytosis, mRNA degradation, and mRNA transfer rates) to fit data presented in *Huang\_2023* for the Short Chain Model.

- `LCM_simulation.m`: Source code for model simulation for both recombinant and mRNA-encoded protein administration in the Long Chain Model. Fixing the dose injected and observation time window allows computing the time evolution of BiTEs' concentration in blood.
- `LCM.m`: Function containing the equations of the Long Chain Model.
- `LCM_parameters_collection_mouse.m`: Function containing parameters related to mouse physiology, including species-specific protein transport, body parameters, tissue absorption dynamics, and mRNA transport dynamics for the Long Chain Model.
- `LCM_parameters_estimated.mat`: `.mat` file containing estimated parameter values (e.g., rate of pinocytosis, mRNA degradation, mRNA transfer rates) to fit data presented in *Huang\_2023* for the Long Chain Model.

---

## 🚀 HOW TO RUN A SIMULATION

To perform a simulation using the **Short Chain Model** or **Long Chain Model**:

1. Navigate to the directory containing the model source files.
2. In MATLAB, run: SCM_simulation or LCM_simulation

---

For any questions or issues, feel free to open an [Issue](https://github.com/) on this repository! 🎯

