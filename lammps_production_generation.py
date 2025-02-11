import os
import shutil

# List of configurations
configurations = [
    {'nodes': 1, 'tasks': 1, 'cpus_per_task': 8}
]

temperature_initial = 300

# Define target temperature
temperature = 1500 

# Define the folder containing the .data files
data_folder = '/path/to/structure/folder'

# Populate the atomic_configs list with all .data files in the folder
atomic_configs = [
    file_name for file_name in os.listdir(data_folder) if file_name.endswith('.data')
]

# Print the list to verify
print("Atomic Configurations:")
print(atomic_configs)

# LAMMPS input file template with placeholders
in_file_template = """# LAMMPS Input Script for Si/SiO2 Interface Simulation using MACE

# Initialization
units         metal
atom_style    atomic
atom_modify   map yes
newton        on

boundary      p p p

# Atom Definition
read_data     {data_file}

# Masses (if not defined in data file)
mass          1 15.9994    # O
mass          2 28.0855    # Si

# Potential Definition
pair_style    mace no_domain_decomposition
pair_coeff    * * /projappl/jeakola/ekvikne/force_fields/mace-small-density-agnesi-stress.model-lammps.pt O Si

# Neighbor List Settings
neighbor      2.0 bin
neigh_modify  delay 0 every 1 check yes

# Thermodynamic Output
thermo        100
thermo_style  custom step temp etotal pe ke press vol

# Dump Snapshots
dump           dmp all custom 100 snapshot.lammpstrj id type x y z

# -----------------------------
# Energy Minimization
# -----------------------------

minimize      1e-10 1e-10 1000 10000

# ---------------------------------
# Initial Equilibration at lower timestep
# ---------------------------------
# Use a smaller timestep for initial equilibration at the starting temperature
timestep      0.00025     # smaller timestep for initial stability

# Velocity Initialization
velocity      all create {temp_init} 12345 dist gaussian

# Equilibrate at {temp_init} with a smaller timestep
fix           eq all nvt temp {temp_init} {temp_init} 0.1
run           5000
unfix         eq

# ---------------------------------------------------------
# Increase the timestep for heating and main simulation
# ---------------------------------------------------------
timestep      0.001      # now we increase the timestep

# Gradual Heating from {temp_init} to {temp}
fix           mynpt all npt temp {temp_init} {temp} 0.2 aniso 0.0 0.0 1.0
run           5000    # gradual heating phase
unfix         mynpt

# Equilibration at {temp} with increased timestep
fix           mynpt all npt temp {temp} {temp} 0.2 aniso 0.0 0.0 1.0
run           10000  # main simulation
unfix         mynpt

write_data    si_sio2_snapshot_{data_name}.data

"""

# Submission script template
submit_script_template = """#!/bin/bash
#SBATCH --account=jeakola
#SBATCH --job-name={data_name}
#SBATCH --output=my_job_%j.out
#SBATCH --error=my_job_%j.err
#SBATCH --partition=large
#SBATCH --time=2-00:00:00
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem=16G

module purge
module load gcc/13.2.0 openmpi/5.0.5 intel-oneapi-mkl/2024.0.0

# Update path to your LAMMPS installation if needed
export PATH=/projappl/jeakola/ekvikne/lammps/bin:$PATH

# Set the number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Set LD_LIBRARY_PATH to include the directory with libtorch.so
export LD_LIBRARY_PATH=/projappl/jeakola/ekvikne/libtorch/lib:$LD_LIBRARY_PATH

# Run LAMMPS with the appropriate input file
srun lmp -sf omp -in in.lammps

# Capture the job ID
JOBID=$SLURM_JOB_ID

# Wait a moment to ensure the job data is available
sleep 5

# Write `seff` output to a file
seff $JOBID > seff_output_${{JOBID}}.txt
"""

# Main directory name
main_dir = "main_dir_name"

# Ensure the main directory exists
os.makedirs(main_dir, exist_ok=True)

# Create the master launcher script
launcher_script_path = os.path.join(main_dir, "launch_all_jobs.sh")
with open(launcher_script_path, 'w') as launcher_script:
    launcher_script.write("#!/bin/bash\n\n")

    # Iterate over configurations and atomic starting files
    for config in configurations:
        for atomic_config in atomic_configs:
            # Get the atomic configuration base name
            data_name = os.path.splitext(atomic_config)[0]

            # Create a unique folder name for each atomic configuration
            folder_name = f"{data_name}"
            folder_path = os.path.join(main_dir, folder_name)
            os.makedirs(folder_path, exist_ok=True)

            # Write the LAMMPS input file with placeholders replaced
            in_file_content = in_file_template.format(
                temp_init=temperature_initial,
                temp=temperature,
                data_file=atomic_config,
                data_name=data_name
            )
            in_file_path = os.path.join(folder_path, "in.lammps")
            with open(in_file_path, 'w') as in_file:
                in_file.write(in_file_content)

            # Customize the submission script
            submit_script_content = submit_script_template.format(
                nodes=config['nodes'],
                tasks=config['tasks'],
                cpus_per_task=config['cpus_per_task'],
                data_name=data_name
            )
            submit_script_path = os.path.join(folder_path, "submit.sh")
            with open(submit_script_path, 'w') as submit_script:
                submit_script.write(submit_script_content)

            # Make the submission script executable
            os.chmod(submit_script_path, 0o755)

            # Copy the atomic configuration file from the data folder
            source_data_file = os.path.join(data_folder, atomic_config)
            dest_data_file = os.path.join(folder_path, atomic_config)
            if os.path.exists(source_data_file):
                shutil.copy(source_data_file, dest_data_file)
            else:
                print(f"Warning: {source_data_file} not found in the data folder.")

            # Add commands to the launcher script
            launcher_script.write(f"cd {folder_name}\n")
            launcher_script.write("sbatch submit.sh\n")
            launcher_script.write("cd ..\n\n")

# Make the launcher script executable
os.chmod(launcher_script_path, 0o755)

print(f"Folder structure created in '{main_dir}'.")
print(f"Run 'bash launch_all_jobs.sh' inside '{main_dir}' to submit all jobs.")
