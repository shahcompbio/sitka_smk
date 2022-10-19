def _get_hscnpath(wildcards):
    df_filt = df[df["patient"] == wildcards.sample]
    path = df_filt["path"].to_list()[0]
    return path

rule hmmcopy_to_sitka_tree_tcn:
    input:
        hscn = _get_hscnpath,
    output:
        sitka_input = "results/input/{sample}-tcn_sitka.csv",
        sitka_segs = "results/input/{sample}-tcn_sitka_segs.csv.gz",
        sitka_transitions = "results/input/{sample}-tcn_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/cnbins_to_sitka_tcn.R"
    
rule hmmcopy_to_sitka_tree_ascn:
    input:
        hscn = _get_hscnpath,
    output:
        sitka_input = "results/input/{sample}-ascn_sitka.csv",
        sitka_segs = "results/input/{sample}-ascn_sitka_segs.csv.gz",
        sitka_transitions = "results/input/{sample}-ascn_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/cnbins_to_sitka_ascn.R"

rule hmmcopy_to_sitka_tree_tcnchrom:
    input:
        hscn = _get_hscnpath,
    output:
        sitka_input = "results/input/{sample}-tcnchrom_sitka.csv",
        sitka_segs = "results/input/{sample}-tcnchrom_sitka_segs.csv.gz",
        sitka_transitions = "results/input/{sample}-tcnchrom_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/cnbins_to_sitka_tcnchrom.R"

rule sitka_tree_inference:
    input:
        "results/input/{sample}-{inputtype}_sitka.csv",
    output:
        posterior='results/output/{sample}-{inputtype}-phylo.csv',
        fnr='results/output/{sample}-{inputtype}-fnr.csv',
        fpr='results/output/{sample}-{inputtype}-fpr.csv',
        logdensity='results/output/{sample}-{inputtype}-sitka_tree_output/logDensity.csv'
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
        'results/output/{sample}-{inputtype}-phylo.csv',
    output:
        'results/output/{sample}-{inputtype}-consensus.newick'
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
        'results/output/{sample}-{inputtype}-phylo.csv'
    output:
        'results/output/{sample}-{inputtype}-average.csv'
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
        'results/output/{sample}-{inputtype}-average.csv'
    output:
        'results/output/{sample}-{inputtype}-tree.newick'
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
        tree = 'results/output/{sample}-{inputtype}-tree.newick',
    output:
        tree = 'results/output/{sample}-{inputtype}-tree-processed.newick'
    threads: 1
    resources:
        mem_mb=1024 * 10
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    script: "../scripts/processtree.R"

rule plotsitka:
    input:
        tree = 'results/output/{sample}-{inputtype}-tree-processed.newick',
        hscn = _get_hscnpath
    output:
        plotpng = "results/output/{sample}-{inputtype}-heatmap.png"
    threads: 1
    resources:
        mem_mb=1024 * 50
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    script: "../scripts/plotheatmap.R"

rule plotheatmaps:
    input:
        trees_filt = expand("results/output/{{sample}}-{inputtype}-tree-processed.newick", inputtype = intypes),
        trees_unfilt = expand("results/output/{{sample}}-{inputtype}-tree.newick", inputtype = intypes),
        hscn = _get_hscnpath,
    output:
        pdf = "results/output/{sample}-heatmaps.pdf"
    threads: 2
    resources: mem_mb=1024 * 25
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    script: "../scripts/plotallheatmaps.R"

