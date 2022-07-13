import pandas as pd
configfile: "config.yaml"
df = pd.read_csv(config["samplesheet"])

rule all:
    input: 
        expand("results/output/{sample}-heatmap.png", sample = df["patient"])

include: "rules/sitka_inference.smk"