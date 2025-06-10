
REFERENCE = "/PHShome/nm104/scripts/vibecheck/vibecheck/resources/reference.fasta"
BARCODES = "/PHShome/nm104/scripts/vibecheck/vibecheck/resources/o1_barcodes.feather"
READ1 = "/data/wohllab/2025.01.28_class/input/CHS3677_R1.fastq.gz"
READ2 = "/data/wohllab/2025.01.28_class/input/CHS3677_R1.fastq.gz"
SAMPLES = {"CHS3677" : [READ1, READ2]}

rule align_reads:
    input:
        read1 = lambda wildcards: SAMPLES[wildcards.sample][0],
        read2 = lambda wildcards: SAMPLES[wildcards.sample][1],
        reference = REFERENCE
    output:
        alignment = "intermediates/alignments/{sample}.sorted-filtered.bam",
        alignment_index = "intermediates/alignments/{sample}.sorted-filtered.bam.bai"
    threads: 8
    shell:
        """
        minimap2 -ax sr -t {threads} {input.reference} {input.read1} {input.read2} |\
        samtools view -b -F 4 -F 2048 -F 8 -F 256 - |\
        samtools sort -o {output.alignment} - && \
        samtools index {output.alignment}
        """


rule subsample_reads:
    input:
        alignment = rules.align_reads.output.alignment
    params:
        reads = lambda wildcards: int( wildcards.reads )
    output:
        sampled_alignment = temp( "intermediates/subsamples/{sample}.{reads}.{trial}.bam" )
    threads: 4
    benchmark: "intermediates/benchmarks/subsample_reads/{sample}.{reads}.{trial}.txt"
    shell:
        """
        cat <(samtools view -H {input.alignment}) <(samtools view {input.alignment} | shuf -n {params.reads}) |\
        samtools view -b - |\
        samtools sort -o {output.sampled_alignment} - &&\
        samtools index {output.sampled_alignment}
        """

rule calculate_depth:
    input:
        alignment = rules.subsample_reads.output.sampled_alignment,
        reference = REFERENCE
    output:
        depth = temp( "intermediates/depth/{sample}.{reads}.{trial}.txt" )
    shell:
        """
        samtools mpileup -aa -A -d 2000 -Q 20 -q 0 -B -f {input.reference} {input.alignment} |\
        cut -f1-4 > {output.depth}
        """

rule call_variants:
    input:
        alignment = rules.subsample_reads.output.sampled_alignment,
        reference = REFERENCE
    output:
        variants = temp( "intermediates/variants/{sample}.{reads}.{trial}.vcf" ),
        variants_filled = "intermediates/variants/{sample}.{reads}.{trial}_filled.vcf"
    threads: 4
    benchmark: "intermediates/benchmarks/call_variants/{sample}.{reads}.{trial}.txt"
    shell:
        """
        bcftools mpileup --count-orphans --threads {threads} -d 600000 -Q 20 -q 0 -B -a INFO/AD,INFO/ADF,INFO/ADR -Ou -f {input.reference} {input.alignment} |\
        bcftools call --threads {threads} -mv -Ov --ploidy 1 -o {output.variants} && \
        bcftools +fill-tags {output.variants} -Ou -o {output.variants_filled} -- -t AF
        """

rule freyja_demix:
    input:
        variants = rules.call_variants.output.variants_filled,
        depth = rules.calculate_depth.output.depth,
        barcodes = BARCODES
    output:
        results = temp( "intermediates/freyja/{sample}.{reads}.{trial}.txt" )
    benchmark: "intermediates/benchmarks/freyja_demix/{sample}.{reads}.{trial}.txt"
    shell:
        """
        freyja demix {input.variants} {input.depth} --output {output.results} --barcodes {input.barcodes} ||\
        echo -e "\\t{input.variants}\\nsummarized\\t[('Other', 0.00)]\\nlineages\\tOther\\nabundances\\t1.0\\nresid\\t0.00\\ncoverage\\t0.0" > {output.results}
        """

rule parse_results:
    input:
        results = rules.freyja_demix.output.results,
        depths = rules.calculate_depth.output.depth
    output:
        parsed_results = "intermediates/parsed_freyja/{sample}.{reads}.{trial}.csv"
    benchmark: "intermediates/benchmarks/parse_results/{sample}.{reads}.{trial}.txt"
    run:
        import numpy as np

        field_parsers = {
            "lineages": lambda x: x[1:],  # Take everything after 'lineages'
            "abundances": lambda x: list(map(float, x[1:])),  # Convert to float
            "resid": lambda x: float(x[-1]),  # Take last value as float
            "coverage": lambda x: float(x[-1]),
        }

        results = {}

        with open(input.results) as f:
            for line in f:
                # Split line into tokens and find matching parser
                tokens = line.strip().split()
                field = next((k for k in field_parsers if line.startswith(k)), None)

                if field:
                    results[field] = field_parsers[field](tokens)

        # Verify all required fields were found
        missing_fields = set(field_parsers) - set(results)
        if missing_fields:
            sys.exit(-67)

        # Calculate conflict score, top lineage, and construct a summary
        confidence = np.exp(sum(a * np.log(a) for a in results["abundances"] if a > 0))
        top_lineage = results["lineages"][np.argmax(results["abundances"])]
        summary = "Freyja results: " + " ".join(
            f"{lin}({abun:.1%})"
            for lin, abun in zip(results["lineages"], results["abundances"])
        )
        
        depth = 0
        count = 0
        with open( input.depths, "rt" ) as depth_file:
            for line in depth_file:
                if int( line.split()[3] ) > 0:
                    depth += 1
                count += 1
        depth /= count

        with open(output.parsed_results, "wt") as out_file:
            out_file.write("sequence_id,lineage,confidence,classification_notes,depth\n")
            out_file.write(f"{input.results},{top_lineage},{confidence:.3f},{summary},{depth:.2%}\n")

rule combine_freyja_results:
    input:
        results = expand( "intermediates/parsed_freyja/{sample}.{reads}.{trial}.csv", sample=SAMPLES, reads=[100, 1000, 10000, 100000, 1000000], trial=range(1,11 ) )
    output:
        combined_results = "results/freyja_masking.csv"
    run:
        import pandas as pd
        df = pd.concat( [pd.read_csv(f) for f in input.results], ignore_index=True )
        df.to_csv( output.combined_results, index=False )


rule combine_freyja_benchmarks:
    input:
        subsamples_reads = expand( "intermediates/benchmarks/subsample_reads/{sample}.{reads}.{trial}.txt", sample=SAMPLES, reads=[
        100, 1000, 10000, 100000, 1000000], trial=range(1,11 ) ),
        call_variants = expand( "intermediates/benchmarks/call_variants/{sample}.{reads}.{trial}.txt",sample=SAMPLES,reads=[
        100, 1000, 10000, 100000, 1000000],trial=range( 1,11 ) ),
        freyja_demix = expand( "intermediates/benchmarks/freyja_demix/{sample}.{reads}.{trial}.txt",sample=SAMPLES,reads=[
        100, 1000, 10000, 100000, 1000000],trial=range( 1,11 ) ),
        parse_results = expand( "intermediates/benchmarks/parse_results/{sample}.{reads}.{trial}.txt",sample=SAMPLES,reads=[
        100, 1000, 10000, 100000, 1000000],trial=range( 1,11 ) )
    output:
        benchmarks = "results/freyja_benchmarks.csv"
    run:
        import pandas as pd
        from pathlib import Path

        bm = list()
        for step in [input.subsamples_reads, input.call_variants, input.freyja_demix, input.parse_results]:
            for file in step:
                path = Path(file)
                df = pd.read_csv( path, sep="\t" )
                df["file"] = path.name
                df["step"] = path.parent.name
                bm.append( df )
        bm = pd.concat( bm, ignore_index=True )
        bm[["reads", "trial"]] = bm["file"].str.extract( r"CHS3677.(\d+).(\d+).txt" )
        for col in ["reads", "trial"]:
            bm[col] = pd.to_numeric( bm[col] )
        bm.to_csv( output.benchmarks, index=False )


rule plot_freyja_benchmarks:
    input:
        results = rules.combine_freyja_benchmarks.output.benchmarks
    output:
        benchmark_plots = "results/plots/freyja_benchmarks.pdf"
    log:
        notebook = "results/notebooks/plot_freyja_benchmarks.ipynb"
    notebook: "../notebooks/plot_freyja_benchmarks.py.ipynb"


rule plot_freyja_results:
    input:
        results = rules.combine_freyja_results.output.combined_results
    output:
        accuracy_coverage_plot = "results/plots/accuracy-coverage-vs-reads.pdf"
    log:
        notebook = "results/notebooks/plot_freyja_masking.ipynb"
    notebook: "../notebooks/plot_freyja_masking.py.ipynb"