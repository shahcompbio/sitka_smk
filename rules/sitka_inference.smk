def _get_hscnpath(wildcards):
    df_filt = df[df["patient"] == wildcards.sample]
    path = df_filt["path"].to_list()[0]
    return path

rule hmmcopy_to_sitka_tree:
    input:
        hscn = _get_hscnpath,
    output:
        sitka_input = "results/input/{sample}_sitka.csv",
        sitka_segs = "results/input/{sample}_sitka_segs.csv.gz",
        sitka_transitions = "results/input/{sample}_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/cnbins_to_sitka.R"

rule sitka_tree_inference:
    input:
        "results/input/{sample}_sitka.csv",
    output:
        posterior='results/output/{sample}-phylo.csv',
        fnr='results/output/{sample}-fnr.csv',
        fpr='results/output/{sample}-fpr.csv',
        logdensity='results/output/{sample}-sitka_tree_output/logDensity.csv'
    threads: 25
    resources:
        mem_mb=1024*2
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    params:
        nscans = config["nscans"],
        fpr = config["fpr"],
        fnr = config["fnr"]
    shell:
        '''
        corrupt-infer-with-noisy-params \
            -Xmx200G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --model.globalParameterization true \
            --model.binaryMatrix {input} \
            --model.fprBound {params.fpr} \
            --model.fnrBound {params.fnr} \
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
        'results/output/{sample}-phylo.csv',
    output:
        'results/output/{sample}-consensus.newick'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        corrupt-l1-decode \
            -Xmx30G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --samples {input};
        mv consensus.newick {output};
        '''

rule sitka_tree_average_tip_indicators:
    input:
        'results/output/{sample}-phylo.csv'
    output:
        'results/output/{sample}-average.csv'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        corrupt-average \
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
        'results/output/{sample}-average.csv'
    output:
        'results/output/{sample}-tree.newick'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*50
    shell:
        '''
        corrupt-greedy \
            -Xmx30G \
            --experimentConfigs.saveStandardStreams false \
            --experimentConfigs.recordGitInfo false \
            --experimentConfigs.managedExecutionFolder false \
            --tipInclusionProbabilities ReadOnlyCLMatrix {input};
        mv tree.newick {output};
        '''

rule formatsitka:
    input:
        tree = 'results/output/{sample}-tree.newick',
    output:
        tree = 'results/output/{sample}-tree-processed.newick'
    threads: 1
    resources:
        mem_mb=1024 * 10
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    script: "../scripts/processtree.R"

rule plotsitka:
    input:
        tree = 'results/output/{sample}-tree-processed.newick',
        hscn = _get_hscnpath
    output:
        plotpng = "results/output/{sample}-heatmap.png"
    threads: 1
    resources:
        mem_mb=1024 * 50
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    script: "../scripts/plotheatmap.R"