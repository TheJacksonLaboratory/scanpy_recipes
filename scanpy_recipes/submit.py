import os
import subprocess

SAMPLE_SCRIPT = """#!/usr/bin/env bash
{header}

Rscript - <<eog
setwd("{wd}")
options(stringsAsFactors = FALSE, row.names = 1, as.is = T)
log2cpm <- read.csv("{sampleid}_counts.csv", row.names = 1, stringsAsFactors = FALSE)
colnames(log2cpm) <- gsub(".", "-", colnames(log2cpm), fixed = T)
tsne.data <- read.csv("{sampleid}_umap3d.csv", row.names = 1, stringsAsFactors = FALSE)
featuredata <- read.csv("{sampleid}_features.csv", row.names = 1, stringsAsFactors = FALSE)
save(log2cpm, featuredata, tsne.data, file="{rdsfile}")
print("Done")
eog
"""

class Submitter(object):
    header = ""
    sample_script = SAMPLE_SCRIPT

    def __init__(self):
        self.scheduler = None
        self.submit_command = None
        self.header_kwargs = {"sampleid", "wd", "walltime"}
        self.script_kwargs = {"sampleid", "wd", "rdsfile"}

    def format(self, **kwargs):
        missing_header_kwargs = self.header_kwargs - set(kwargs.keys())
        assert missing_header_kwargs == set(), \
            f"Need to supply header kwargs: [{missing_header_kwargs}]"
        self.header = self.header.format(**kwargs)
        self.output_dir = kwargs.get("wd")

        missing_script_kwargs = self.script_kwargs - set(kwargs.keys())
        assert missing_script_kwargs == set(), \
            f"Need to supply script kwargs: [{missing_script_kwargs}]"
        self.sample_script = self.sample_script.format(header=self.header, **kwargs)

    def submit(self, additional_script="", **kwargs):
        save_file = os.path.join(self.output_dir, ".submit_rds_creation.sh")
        with open(save_file, "w") as fout:
            fout.write(self.sample_script)

        cmd = f"{self.submit_command} {save_file}"
        output = subprocess.check_output(cmd, shell=True).decode("ascii").strip()
        print(f"Rds creation submitted as {self.scheduler} job {output}.")
        print(f"Output will be located in [{self.output_dir}].")


class PBSSubmitter(Submitter):
    header = (
        "#PBS -N {sampleid}_rds_creation\n"
        "#PBS -o {wd}/{sampleid}_rds_creation.log\n"
        "#PBS -j oe\n"
        "#PBS -l walltime={walltime}\n"
        "#PBS -l vmem={mem}mb\n"
        "#PBS -l nodes=1:ppn=1\n"
        "module load R"
    )
    def __init__(self, ):
        super().__init__()
        self.scheduler = "pbs"
        self.submit_command = "qsub"


class SlurmSubmitter(Submitter):
    header = (
        "#SBATCH --job-name={sampleid}_rds_creation\n"
        "#SBATCH --output={wd}/{sampleid}_rds_creation.log\n"
        "#SBATCH --mail-type=FAIL\n\n"
        "#SBATCH --partition=batch\n"
        "#SBATCH --time={walltime}\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks-per-node=1\n"
        "#SBATCH --mem={mem}\n"
        "#SBATCH --export=NONE\n"
        "module load R"
    )
    def __init__(self, ):
        super().__init__()
        self.scheduler = "slurm"
        self.submit_command = "sbatch"


class LocalSubmitter(Submitter):
    # string formatting will only barf if keys present in string are not
    # present in format call, but not the other way around.  We can exploit
    # that here.
    header = ""
    def __init__(self, ):
        super().__init__()
        self.scheduler = "local"
        self.submit_command = "$SHELL"


SUBMITTERS = {
    "pbs": PBSSubmitter,
    "slurm": SlurmSubmitter,
    "local": LocalSubmitter,
}


def submit_rds_job(
        sampleid,
        output_dir,
        rds_filename,
        walltime="00:15:00",
        mem=64000,
        scheduler="pbs"
    ):
    if scheduler not in SUBMITTERS.keys():
        raise ValueError(f"Cannot handle scheduler [{scheduler}]!")

    submitter = SUBMITTERS[scheduler]()
    submitter.format(
        sampleid=sampleid, wd=output_dir, walltime=walltime, mem=mem,
        rdsfile=rds_filename
    )
    submitter.submit()
