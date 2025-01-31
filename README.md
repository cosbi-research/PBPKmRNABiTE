# A Multi-Scale Physiologically Based Pharmacokinetics Model to Support mRNA-Encoded BiTE Therapy in Cancer Treatment

**Authors:** Giada Fiandaca, Elio Campanile, Lorena Leonardelli, Elisa Pettinà, Stefano Giampiccolo, Elizabeth J Carstens, Lorenzo Dasti, Natascia Zangani, Luca Marchetti

This repository contains the complete MATLAB source codes for reproducing the dynamics of mRNA trafficking, translation, and the subsequent production of functional bispecific antibodies. It also includes their trafficking in both of the two modeling approaches proposed in the paper: the **Short Chain Model (SCM)** and the **Long Chain Model (LCM)**.

---

## 📂 CONTENTS

### **Short Chain Model (SCM)**

- ``: Function containing the parameters regulating the transport of 55 kDa proteins (BiTE dimension is 55 kDa).
- ``: Function containing the equations of the model for recombinant protein transport.
- ``: Source code for the model simulation for both recombinant and mRNA-encoded protein administration in the Short Chain Model. Fixing the dose injected and the observation time window allows computing the time evolution of BiTEs' concentration in blood.
- ``: Function containing the equations of the Short Chain Model.
- ``: Function containing parameters related to mouse physiology, including species-specific protein transport, body parameters, tissue absorption dynamics, and mRNA transport dynamics for the Short Chain Model.
- ``: `.mat` file containing estimated parameter values (e.g., rate of pinocytosis, mRNA degradation, and mRNA transfer rates) to fit data presented in *Huang\_2023* for the Short Chain Model.

### **Long Chain Model (LCM)**

- ``: Source code for model simulation for both recombinant and mRNA-encoded protein administration in the Long Chain Model. Fixing the dose injected and observation time window allows computing the time evolution of BiTEs' concentration in blood.
- ``: Function containing the equations of the Long Chain Model.
- ``: Function containing parameters related to mouse physiology, including species-specific protein transport, body parameters, tissue absorption dynamics, and mRNA transport dynamics for the Long Chain Model.
- ``: `.mat` file containing estimated parameter values (e.g., rate of pinocytosis, mRNA degradation, mRNA transfer rates) to fit data presented in *Huang\_2023* for the Long Chain Model.

---

## 🚀 HOW TO RUN A SIMULATION

To perform a simulation using the **Short Chain Model** or **Long Chain Model**:

1. Navigate to the directory containing the model source files.
2. In MATLAB, run:
   ```matlab
   SCM_simulation
   ```
   or
   ```matlab
   LCM_simulation
   ```

This will execute the respective model and simulate the BiTEs' concentration dynamics in the bloodstream.

---

## 📌 Citation

If you use this model in your research, please cite our paper: **A Multi-Scale Physiologically Based Pharmacokinetics Model to Support mRNA-Encoded BiTE Therapy in Cancer Treatment**.

---

For any questions or issues, feel free to open an [Issue](https://github.com/) on this repository! 🎯

