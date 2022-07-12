
def _get_hscnpath(wildcards):
    df_filt = df[df["patient"] == wildcards.sample]
    path = df_filt["path"].to_list()[0]
    return path

rule hmmcopy_to_sitka_tree:
    input:
        cn_df = _get_hscnpath,
    output:
        sitka_input = "results/{sample}/inputdata/{sample}_sitka.csv",
        sitka_segs = "results/{sample}/inputdata/{sample}_sitka_segs.csv.gz",
        sitka_transitions = "results/{sample}/inputdata/{sample}_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliam1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/phylo/cnbins_to_sitka.R"

rule sitka_tree_inference:
    input:
        "results/{sample}/inputdata/{sample}_sitka.csv",
    output:
        posterior='results/{sample}/results/{sample}-phylo.csv',
        fnr='results/{sample}/results/{sample}-fnr.csv',
        fpr='results/{sample}/results/{sample}-fpr.csv',
        logdensity='results/{sample}/results/{sample}-sitka_tree_output/logDensity.csv'
    threads: 10
    resources:
        mem_mb=1024*2
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    params:
        nscans = config["nscans"],
    shell:
        '''
        sitka-infer-with-noisy-params \
            -Xmx200G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --model.globalParameterization true \
            --model.binaryMatrix {input} \
            --model.fprBound 0.1 \
            --model.fnrBound 0.5 \
            --engine PT \
            --engine.initialization FORWARD \
            --engine.ladder Polynomial \
            --engine.nScans {params.nscans} \
            --engine.nPassesPerScan 1 \
            --engine.nChains 10 \
            --engine.nThreads Fixed \
            --engine.nThreads.number 50;
        mv samples/phylo.csv {output.posterior};
        mv samples/fnr.csv {output.fnr};
        mv samples/fpr.csv {output.fpr};
        mv samples/logDensity.csv {output.logdensity};
        '''

rule sitka_tree_consensus:
    input:
        'results/{sample}/results/{sample}-phylo.csv',
    output:
        'results/{sample}/results/{sample}-consensus.newick'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        sitka-l1-decode \
            -Xmx30G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --samples {input};
        mv consensus.newick {output};
        '''

rule sitka_tree_average_tip_indicators:
    input:
        'results/{sample}/results/{sample}-phylo.csv'
    output:
        'results/{sample}/results/{sample}-average.csv'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        sitka-average \
            -Xmx30G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --csvFile {input} \
            --logisticTransform false;
        mv average.csv {output};
        '''

rule sitka_tree_decode:
    input:
        'results/{sample}/results/{sample}-average.csv'
    output:
        'results/{sample}/results/{sample}-tree.newick'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        sitka-greedy \
            -Xmx30G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --tipInclusionProbabilities ReadOnlyCLMatrix {input};
        mv tree.newick {output};
        '''

rule formatsitka:
    input:
        tree = 'results/{sample}/results/{sample}-tree.newick',
    output:
        tree = 'results/processed/{sample}-tree.newick'
    threads: 1
    resources:
        mem_mb=1024 * 10
    singularity: "docker://marcjwilliam1/signals:v0.7.6"
    script: "../scripts/phylo/processtree.R"

rule plotsitka:
    input:
        tree = 'results/processed/{sample}-tree.newick',
        hscn = _get_hscnpath
    output:
        plot = report("results/plots/{sample}-heatmap.pdf", category = "sitka"),
        plotpng = "results/plots/{sample}-heatmap.png"
    threads: 1
    resources:
        mem_mb=1024 * 50
    singularity: "docker://marcjwilliam1/signals:v0.7.6"
    script: "../scripts/phylo/plotheatmap.R"