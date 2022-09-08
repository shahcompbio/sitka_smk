import pandas as pd
configfile: "config.yaml"
df = pd.read_csv(config["samplesheet"])
intypes = ["tcn", "ascn"]#,  "tcnpadded"]

rule all:
    input: 
        expand("results/output/{sample}-{inputtype}-heatmap.png", sample = df["patient"], inputtype = intypes)

include: "rules/sitka_inference.smk"