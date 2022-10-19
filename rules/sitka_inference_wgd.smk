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

rule sitka_grow:
    input:
        sitka_input = "results/input/{sample}-tcnploidy_sitka_other.csv",
        tree = 'results/output/{sample}-tcnploidy-tree.newick'
    output:
        tree = 'results/output/{sample}-tcnploidygrow-tree.newick',
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*25
    shell:
        '''
        corrupt-grow -Xmx20G \
            --matrix NoisyBinaryCLMatrix \
            --matrix.binaryMatrix {input.sitka_input} \
            --matrix.fpRate 0.1 \
            --matrix.fnRate 0.9 \
            --phylo file {input.tree}
        pwd
        ls -lth results/
        mv results/all/*/grown.newick {output.tree}
        rm -r results/all/
        '''

rule sitka_grow2:
    input:
        sitka_input = "results/input/{sample}-tcnploidy_sitka_other.csv",
        tree = 'results/output/{sample}-tcnploidy-tree.newick'
    output:
        tree2 = 'results/output/{sample}-tcnploidygrow-tree2.newick'
    shadow: 'shallow'
    singularity: 'shub://funnell/nowellpack_singularity'
    resources:
        mem_mb=1024*25
    shell:
        '''
        corrupt-grow -Xmx20G \
            --matrix ReadOnlyCLMatrix {input.sitka_input} \
            --phylo file {input.tree}
        pwd
        ls -lth results/
        mv results/all/*/grown.newick {output.tree2}
        rm -r results/all/
        '''