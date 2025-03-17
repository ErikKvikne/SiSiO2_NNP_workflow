# **Si/SiO2 Neural Network Potential Training Workflow**

This repository contains the **automated workflow** for generating **Si/SiO₂ interfaces**, running **Molecular Dynamics (MD) simulations**, and performing **Density Functional Theory (DFT) calculations**. The goal is to create a **high-quality dataset** for training **Neural Network Potentials (NNPs)**, enabling accurate and efficient large-scale simulations of silicon/silica interfaces.

## **Project Overview**  
Silicon/silica (Si/SiO₂) interfaces are crucial in semiconductor devices, influencing properties such as **charge distribution**, **defect states**, and **electronic structure**. Traditional methods like **DFT** are computationally expensive for large systems, while empirical force fields lack accuracy. **NNPs** offer a promising alternative, but their effectiveness depends on the quality of the training dataset.

This workflow automates:  
- **Si/SiO₂ interface generation** using **Pymatgen**.  
- **MD simulations** using **LAMMPS** and **MACE-MP-0**.  
- **DFT calculations** using **CP2K**.  
- **Data collection** for training an **NNP**.  

## **Repository Structure**  
```
📂 SiSiO2_NNP_workflow/  
│── src/ (Main source code)  
│   │── structures/ (Interface generation scripts)  
│   │   ├── si_sio2_generator.py (Class for Si/SiO2 interface generation)  
│   │── md_simulations/ (Molecular Dynamics (MD) scripts)  
│   │   ├── md_runner.py (Runs LAMMPS simulations)  
│   │   ├── analysis.py (Post-processing tools)  
│   │── dft_calculations/ (DFT simulation scripts)  
│   │   ├── cp2k_runner.py (Runs CP2K calculations)  
│   │   ├── input_generator.py (Generates CP2K input files)  
│   │── utils/ (Helper functions)  
│   │   ├── file_io.py (File reading/writing utilities)  
│   │   ├── visualization.py (Visualization tools)  
│── scripts/ (Standalone execution scripts)
│   │── md/ (Interface generation scripts) 
│   │   ├── md_runner.py (Runs LAMMPS simulations)  
│   │   ├── analysis.py (Post-processing tools)  
│── README.md (Project documentation)  
│── .gitignore (Files to exclude from version control)  
```
