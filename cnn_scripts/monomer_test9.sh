#!/bin/bash

#SBATCH --job-name=monomer_cnag_csplus2_100G_1gpu_16xCPU_xla_memf_5.0_nocontainer       #Job Name

#SBATCH --output=%j_monomer_cnag_csplus2_100G-1gpu_16xCPU_xla_memf_5.0_nocontainer.out

#SBATCH --error=%j_monomer_cnag_csplus2_100G-1gpu_16xCPU_xla_memf_5.0_nocontainer.err   #Error file will allow you to follow progress

#SBATCH -p compsci-gpu            # Need the GPU partition

#SBATCH --mem=200G          #set memory usage

#SBATCH --gres=gpu:v100:1              # request a GPU

#SBATCH -c 16                      #alphafold uses 8 cpus

#SBATCH --mail-type=END

#SBATCH --mail-user=bmp40@duke.edu

##SBATCH --exclusive

 

export SINGULARITY_TMPDIR=/usr/xtmp/derf/tmp             #You need to set the tmp and cache files for singularity

export SINGULARITY_CACHEDIR=/usr/xtmp/derf/cache

export PATH=/usr/project/csplus2/dietrich/pkg/hhsuite/bin:/usr/project/csplus2/dietrich/pkg/conda/bin:$PATH

export TF_FORCE_UNIFIED_MEMORY=1

export XLA_PYTHON_CLIENT_MEM_FRACTION=9.0

export OPENMM_CPU_THREADS=8

export XLA_PYTHON_CLIENT_ALLOCATOR=platform

unset PYTHONPATH

 

proteinName=AER304C                          #Set protein name, this will be folder title with output

faFile=/usr/xtmp/derf/alphafold/data/test/${proteinName}    #Set path to fasta file

outputPath=/usr/xtmp/derf/alphafold/output               #Set path to output

ALPHAFOLD_DATA_PATH=/usr/project/csplus2/dietrich/datadir/    #This should not change

 

date

hostname

nvidia-smi

echo "The job is getting started..."

/usr/project/csplus2/dietrich/pkg/conda/bin/python3 /usr/project/csplus2/dietrich/alphafold/run_alphafold.py \

            --fasta_paths=$faFile \

            --data_dir=$ALPHAFOLD_DATA_PATH \

            --run_relax=True \

            --use_gpu_relax=True\

            --model_preset=monomer \

            --db_preset=full_dbs \

          --max_template_date=2021-11-01 \

            --uniref90_database_path=$ALPHAFOLD_DATA_PATH/uniref90/uniref90.fasta \

            --mgnify_database_path=$ALPHAFOLD_DATA_PATH/mgnify/mgy_clusters_2018_12.fa \

            --uniclust30_database_path=$ALPHAFOLD_DATA_PATH/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \

            --bfd_database_path=$ALPHAFOLD_DATA_PATH/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \

            --template_mmcif_dir=$ALPHAFOLD_DATA_PATH/pdb_mmcif/mmcif_files \

            --obsolete_pdbs_path=$ALPHAFOLD_DATA_PATH/pdb_mmcif/obsolete.dat \

            --pdb70_database_path=$ALPHAFOLD_DATA_PATH/pdb70/pdb70 \

            --output_dir=$outputPath

 

echo "The analysis is done."

date