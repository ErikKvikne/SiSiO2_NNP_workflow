# **Si/SiO2 Neural Network Potential Training Workflow**

This repository contains the **automated workflow** for generating **Si/SiO₂ interfaces**, running **Molecular Dynamics (MD) simulations**, and performing **Density Functional Theory (DFT) calculations**. The goal is to create a **high-quality dataset** for training **Neural Network Potentials (NNPs)**, enabling accurate and efficient large-scale simulations of silicon/silica interfaces.

## **Project Overview**  
Silicon/silica (Si/SiO₂) interfaces are crucial in semiconductor devices, influencing properties such as **charge distribution**, **defect states**, and **electronic structure**. Traditional methods like **DFT** are computationally expensive for large systems, while empirical force fields lack accuracy. **NNPs** offer a promising alternative, but their effectiveness depends on the quality of the training dataset.

This workflow automates:  
- **Si/SiO₂ interface generation** using **Pymatgen**.  
- **MD simulations** using **LAMMPS** and **MACE-MP-0**.  
- **DFT calculations** using **CP2K**.  
- **Data collection** for training an **NNP**.  
