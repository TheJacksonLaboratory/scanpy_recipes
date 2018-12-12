import subprocess


SAMPLE_SCRIPT = """
#!/usr/bin/env bash
#PBS -N {sampleid}_rds_creation
#PBS -o {wd}/{sampleid}_rds_creation.log
#PBS -j oe
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1
module load R

Rscript - <<eog
setwd("{wd}")
log2cpm <- read.csv("{sampleid}_counts.csv")
tsne.data <- read.csv("{sampleid}_umap3d.csv")
featuredata <- read.csv("{sampleid}_features.csv")
save(log2cpm, featuredata, tsne.data, file="{rdsfile}")
print("Done")
eog
"""

def submit_rds_job(sampleid, output_dir, rds_filename):
    sample_script = SAMPLE_SCRIPT.format(
        sampleid=sampleid,
        wd=output_dir,
        rdsfile=rds_filename
    )
    cmd = f"""
    qsub <<eof
    {sample_script}
    eof"""
    output = subprocess.check_output(cmd, shell=True).decode("ascii").strip()
    print(f"Rds creation submitted as helix job {output}.")
    print(f"Output will be located in [{output_dir}].")
