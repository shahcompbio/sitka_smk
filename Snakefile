import pandas as pd
import numpy as np
configfile: "config.yaml"
df = pd.read_csv(config["samplesheet"])
intypes = ["tcn", "ascn",  "tcnchrom", "tcnploidygrow"]

rule all:
    input: 
        expand("results/output/{sample}-{inputtype}-heatmap.png", sample = PATIENTS, inputtype = intypes),
        expand("results/output/{sample}-heatmaps.pdf", sample = PATIENTS)

include: "rules/sitka_inference.smk"
include: "rules/sitka_inference_wgd.smk"