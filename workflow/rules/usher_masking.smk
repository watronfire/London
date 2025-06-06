import pandas as pd
import numpy as np

candidates = pd.read_csv( "data/candidates.csv" )
CANDIDATES = candidates.set_index( "id" )["sequence_id"].to_dict()
REFERENCE = "data/cholera_reference.fasta"
TREE = "/PHShome/nm104/analysis/2025.01.15_class/pruned.pb"
FREQUENCIES = range(5,100,5)

rule extract_sequence:
    input:
        vcf = rules.prepare_VCF.output.compressed_sequences,
        reference = REFERENCE
    params:
        seq_id = lambda wildcards: CANDIDATES[wildcards.sample]
    output:
        sequence = temp( "intermediates/sequences/{sample}.fasta" )
    shell:
        """
        bcftools consensus --missing N --mark-del N -f {input.reference} -s '{params.seq_id}' {input.vcf} |\
        sed -e 's/*/N/g' > {output.sequence}
        """

rule randomize_sample:
    input:
        sequence = rules.extract_sequence.output.sequence
    params:
        trials = 100,
        frequency = lambda wildcards: int( wildcards.freq ) / 100,
        method = "smooth"
    output:
        randomized_sequences = temp( "intermediates/random_mask/{sample}_{freq}.fasta" )
    threads: 2
    shell:
        """
        RandomSeqMask \
            --sequence {input.sequence} \
            --fraction {params.frequency} \
            --method {params.method} \
            --count {params.trials} \
            --outfile {output.randomized_sequences}
        """

rule run_vibecheck_usher_masked:
    input:
        sequences = rules.randomize_sample.output.randomized_sequences,
    params:
        outdir = "intermediates/usher_masked/{sample}_{freq}/"
    output:
        parsed_output = "intermediates/usher_masked/{sample}_{freq}/{sample}_{freq}.csv"
    threads: 8
    shell:
        """
        vibecheck \
            --outdir {params.outdir} \
            --outfile {wildcards.sample}_{wildcards.freq}.csv \
            --threads {threads} \
            --max-ambiguity 1 \
            {input.sequences} 
        """


rule mask_variants_usher:
    input:
        sequences = TEST_SEQUENCES,
        canidates = "data/candidates.csv"
    params:
        frequencies = FREQUENCIES,
        trials = 100
    output:
        masked_variants = "results/usher-variant-masking.csv"
    log:
        notebook = "results/notebooks/mask-variants-usher.ipynb"
    threads: 16
    notebook: "../notebooks/mask-variants-usher.py.ipynb"


rule combine_usher_results:
    input:
        usher = expand( "intermediates/usher_masked/{sample}_{freq}/{sample}_{freq}.csv", sample=CANDIDATES, freq=FREQUENCIES )
    output:
        results = "results/usher-masking.csv"
    run:
        import pandas as pd
        import os

        results = list()
        for f in input.usher:
            name = os.path.splitext( os.path.basename( f ) )[0]
            lineage = name.split( "_" )[-2]
            frac = name.split( "_" )[-1]
            df = pd.read_csv( f )
            df["trial"] = df["sequence_id"].apply( lambda x: x.split( "_" )[1] )
            df["file"] = name
            df["actual"] = lineage
            df["frac_masked"] = frac
            results.append( df )
        results = pd.concat( results, ignore_index=True )
        results.to_csv( output.results, index=False )


rule plot_masking_results:
    input:
        results = rules.combine_usher_results.output.results
    output:
        accuracy_plot = "results/plots/accuracy-vs-masking.pdf",
        parsimony_plot = "results/plots/parsimony-vs-masking.pdf"
    log:
        notebook = "results/notebooks/plot_usher_masking.ipynb"
    notebook: "../notebooks/plot_usher_masking.py.ipynb"

