import pandas as pd

configfile: "config.yaml"
df = pd.read_csv(config["patientlist"])
print(df)

rule all:
    input: 
        expand("results/plots/{sample}-heatmap.png", sample = df["patient"])

include: "rules/sitka_inference.smk"