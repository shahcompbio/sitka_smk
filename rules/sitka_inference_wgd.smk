rule hmmcopy_to_sitka_tree_tcn_wgd:
    input:
        hscn = _get_hscnpath,
    output:
        sitka_input_modeploidy = "results/input/{sample}-tcnploidy_sitka.csv",
        sitka_input_other = "results/input/{sample}-tcnploidy_sitka_other.csv",
        sitka_segs = "results/input/{sample}-tcnploidy_sitka_segs.csv.gz",
        sitka_transitions = "results/input/{sample}-tcnploidy_sitka_transitions.csv.gz",
    singularity: "docker://marcjwilliams1/signals:v0.7.6"
    threads: 15
    resources:
        mem_mb=1024*10
    script:
        "../scripts/cnbins_to_sitka_tcnwgd.R"

#get fn/fp rate by taking last 100 values of posterior
def _get_fnrate(wildcards):
    file = "results/output/" + wildcards.sample + "-tcnploidy-fnr.csv"
    df = pd.read_csv(file)
    return np.mean(df.tail(100).value.tolist())

def _get_fprate(wildcards):
    file = "results/output/" + wildcards.sample + "-tcnploidy-fpr.csv"
    df = pd.read_csv(file)
    return np.mean(df.tail(100).value.tolist())

rule sitka_grow:
    input:
        sitka_input = "results/input/{sample}-tcnploidy_sitka_other.csv",
        tree = 'results/output/{sample}-tcnploidy-tree.newick'
    output:
        tree = 'results/output/{sample}-tcnploidygrow-tree.newick',
    params:
        fprate = _get_fprate,
        fnrate = _get_fnrate,
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*25
    shell:
        '''
        corrupt-grow -Xmx20G \
            --matrix NoisyBinaryCLMatrix \
            --matrix.binaryMatrix {input.sitka_input} \
            --matrix.fpRate {params.fprate} \
            --matrix.fnRate {params.fnrate} \
            --phylo file {input.tree}
        mv results/all/*/grown.newick {output.tree}
        rm -r results/all/
        '''