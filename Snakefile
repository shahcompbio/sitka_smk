import pandas as pd
configfile: "config.yaml"
df = pd.read_csv(config["samplesheet"])
intypes = ["tcn", "ascn",  "tcnchrom"]

rule all:
    input: 
        expand("results/output/{sample}-{inputtype}-heatmap.png", sample = df["patient"], inputtype = intypes),
        expand("results/output/{sample}-summary.html", sample = df["patient"])

include: "rules/sitka_inference.smk"