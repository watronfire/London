LINEAGE_CALLS = "data/test.lineages.tsv"
TEST_SEQUENCES = "data/test.sequences.vcf"
REFERENCE = "data/cholera_reference.fasta"
TERRA_RESULTS = "data/taxon_results.tsv"

rule prepare_VCF:
    input:
        sequences = TEST_SEQUENCES
    output:
        compressed_sequences = "intermediates/sequences.bcf.gz",
        sequences_index = "intermediates/sequences.bcf.gz.csi"
    shell:
        """
        bcftools view -Ob {input.sequences} > {output.compressed_sequences} &&
        bcftools index {output.compressed_sequences}
        """

rule prepare_raw_sequences:
    input:
        sequences = rules.prepare_VCF.output.compressed_sequences,
        reference = REFERENCE
    params:
        script = "workflow/scripts/vcf_to_fasta.py"
    output:
        raw_sequences = "intermediates/raw_sequences.fasta"
    shell:
        """
        python {params.script} \
            --reference={input.reference} \
            --vcf {input.sequences} \
            --output={output.raw_sequences}
        """

rule subsample_sequences:
    input:
        sequences = rules.prepare_raw_sequences.output.raw_sequences
    output:
        subsampled_sequences = "intermediates/subsamples/sample{sample}-trail{trial}.fasta"
    group: "benchmark"
    shell:
        """
        seqtk sample {input.sequences} {wildcards.sample} > {output.subsampled_sequences}
        """

rule benchmark_vibecheck_speed:
    input:
        subsampled_sequences = rules.subsample_sequences.output.subsampled_sequences
    params:
        outdir = "intermediates/speed_runs/sample{sample}-trail{trial}/"
    output:
        lineage_call = "intermediates/speed_runs/sample{sample}-trail{trial}/lineage_report.csv"
    benchmark: "intermediates/benchmarks/sample{sample}-trail{trial}.tsv"
    group: "benchmark"
    shell:
        """
        vibecheck --outdir {params.outdir} {input.subsampled_sequences}
        """


rule combine_speed_benchmark:
    input:
        benchmarks = expand( "intermediates/benchmarks/sample{sample}-trail{trial}.tsv", sample=SPEED_BENCHMARKS, trial=range(3) )
    output:
        benchmark_summary = "results/benchmark-summary.csv",
        benchmark_plot = "results/plots/speed-memory-benchmark.pdf"
    run:
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit

        res = list()
        for file in input.benchmarks:
            df = pd.read_csv( file,sep='\t' )
            df["file"] = os.path.basename( file )
            res.append( df )

        res = pd.concat( res )
        res[["sequences", "trial"]] = res["file"].str.extract( r"sample(\d+)-trail(\d+).tsv" )
        res["sequences"] = pd.to_numeric( res["sequences"] )
        res["trial"] = pd.to_numeric( res["trial"] )
        res = res.drop( columns=["file"] )
        res.to_csv( output.benchmark_summary, index=False )


        def basic_exp( x, a, b ):
            return a * np.power( x,b )


        fig, ax = plt.subplots( dpi=200,figsize=(12, 3.5),ncols=3 )

        fit, _ = curve_fit( basic_exp,xdata=res["sequences"],ydata=res["cpu_time"] )
        xs = np.linspace( 1,103,100 )
        ys = [basic_exp( x,*fit ) for x in xs]

        points = ax[0].scatter( res["sequences"],res["cpu_time"],s=10,zorder=100,c="#009E73" )
        points.set_clip_on( False )
        ax[0].plot( xs,ys,color="black",linewidth=1,linestyle="--",zorder=50 )
        ax[0].set_ylim( 0 )
        ax[0].set_xlim( -3,103 )
        ax[0].set_xlabel( "Sequences",fontweight="bold" )
        ax[0].set_ylabel( "CPU time (seconds)",fontweight="bold" )
        ax[0].grid( color="#EFEFEF" )

        fit2, _ = curve_fit( basic_exp,xdata=res["sequences"],ydata=res["cpu_time"] / res["sequences"] )
        ys2 = [basic_exp( x,*fit2 ) for x in xs]
        points = ax[1].scatter( res["sequences"],res["cpu_time"] / res["sequences"],s=10,zorder=100,c="#009E73" )
        points.set_clip_on( False )
        ax[1].plot( xs,ys2,color="black",linewidth=1,linestyle="--",zorder=50 )
        ax[1].set_ylim( 0 )
        ax[1].set_xlim( -3,103 )
        ax[1].set_xlabel( "Sequences",fontweight="bold" )
        ax[1].set_ylabel( "CPU time per sequence (seconds)",fontweight="bold" )

        ax[1].grid( color="#EFEFEF" )

        fit3, _ = curve_fit( basic_exp,xdata=res["sequences"],ydata=res["max_rss"] )
        ys3 = [basic_exp( x,*fit3 ) for x in xs]

        points = ax[2].scatter( res["sequences"],res["max_rss"],s=10,zorder=100,c="#009E73" )
        points.set_clip_on( False )
        ax[2].plot( xs,ys3,color="black",linewidth=1,linestyle="--",zorder=50 )
        ax[2].set_ylim( 0 )
        ax[2].set_xlim( -3,103 )
        ax[2].set_xlabel( "Sequences",fontweight="bold" )
        ax[2].set_ylabel( "Max memory (MB)",fontweight="bold" )

        ax[2].grid( color="#EFEFEF" )

        plt.tight_layout()
        plt.savefig( output.benchmark_plot )


rule classify_all_sequences:
    input:
        sequences = rules.prepare_raw_sequences.output.raw_sequences
    params:
        outdir = "intermediates/accuracy_runs/"
    output:
        lineage_call = "intermediates/accuracy_runs/lineage_report.csv"
    threads: 4
    shell:
        """
        vibecheck --threads {threads} --outdir {params.outdir} {input.sequences}
        """


rule plot_classification_stats:
    input:
        lineage_calls = rules.classify_all_sequences.output.lineage_call,
        actual_lineages = LINEAGE_CALLS
    output:
        accuracy_plot = "results/plots/accuracy-confusion-matrix.pdf",
        parsimony_plot = "results/plots/parsimoneous-placements.pdf"
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker
        import numpy as np
        from statsmodels.stats.proportion import proportion_confint

        actual = pd.read_csv( input.actual_lineages,sep="\t" )
        actual = actual.rename( columns={ "taxa": "sequence_id", "te": "lineage_actual" } )
        actual["lineage_actual"] = actual["lineage_actual"].replace( { "UNK": "UNDEFINED" } )

        res = pd.read_csv( input.lineage_calls )
        res = res.merge( actual,on="sequence_id",how="left" )

        assert res.loc[res["lineage_actual"].isna()].shape[0] == 0

        res["correct"] = res["lineage"] == res["lineage_actual"]
        res["parsimony_placements"] = res["classification_notes"].str.extract( r"[A-Z0-9.]+\([0-9]+\/([0-9]+)\)" )
        res["parsimony_placements"] = pd.to_numeric( res["parsimony_placements"] )

        successes = res["correct"].sum()
        observations = res.shape[0]
        accuracy = successes / observations
        ci = proportion_confint( successes,observations,alpha=0.05,method="jeffries" )

        accuracy_sum = res.groupby( "lineage_actual" )["correct"].agg( ["sum", "count"] )
        accuracy_sum.columns = ["successes", "observations"]
        accuracy_sum = accuracy_sum.reset_index()
        accuracy_sum["accuracy"] = accuracy_sum["successes"] / accuracy_sum["observations"]
        accuracy_sum[["accuracy_low", "accuracy_high"]] = accuracy_sum.apply(
            lambda x: pd.Series(
                proportion_confint(
                    x["successes"],x["observations"],alpha=0.05,method="jeffreys"
                )
            ),axis=1
        )

        accuracy_sum["numeric_lineage"] = accuracy_sum[
            "lineage_actual"].str.extract( r"(\d+)" ).fillna( 99 ).astype( int )

        accuracy_sum = accuracy_sum.sort_values( by="numeric_lineage" ).reset_index( drop=True )

        confusion_matrix = res.pivot_table( columns="lineage_actual",index="lineage",values="correct",aggfunc="count",fill_value=0 )
        confusion_matrix = confusion_matrix.reindex(
            columns=accuracy_sum["lineage_actual"],index=accuracy_sum["lineage_actual"]
        )

        fig, ax = plt.subplots( dpi=200,figsize=(10, 4),ncols=2 )

        lineages = accuracy_sum.shape[0]

        points = ax[0].scatter( accuracy_sum.index,accuracy_sum["accuracy"],color="#009E73",zorder=100 )
        points.set_clip_on( False )
        points = ax[0].scatter( lineages,accuracy,color="#009E73",zorder=100 )
        points.set_clip_on( False )

        ax[0].vlines(
            accuracy_sum.index,
            accuracy_sum["accuracy_low"],accuracy_sum["accuracy_high"],color="#009E73",capstyle="round"
        )
        ax[0].vlines( lineages,ci[0],ci[1],color="#009E73",zorder=100 )

        ax[0].set_xticks(
            range( lineages + 1 ),
            accuracy_sum["lineage_actual"].replace( { "UNDEFINED": "Other" } ).to_list() + ["All"],fontsize=7
            )
        ax[0].yaxis.set_major_formatter( mticker.PercentFormatter( 1 ) )

        ax[0].set_ylim( 0,1 )
        ax[0].set_xlim( -0.5,lineages + 0.5 )

        ax[0].axvline( lineages - 0.5,color="black",linewidth=0.75 )

        ax[0].set_ylabel( "Accuracy",fontweight="bold" )

        ax[0].grid( color="#EFEFEF" )

        ax[1].imshow( confusion_matrix,cmap="Blues",vmax=100 )

        ax[1].set_xticks(
            range( confusion_matrix.shape[0] ),[i.replace( "UNDEFINED","Other" ) for i in
                                                confusion_matrix.columns],fontsize=7
        )
        ax[1].set_yticks(
            range( confusion_matrix.shape[1] ),[i.replace( "UNDEFINED","Other" ) for i in
                                                confusion_matrix.index],fontsize=7
        )

        for i in range( confusion_matrix.shape[0] ):
            for j in range( confusion_matrix.shape[1] ):
                value = confusion_matrix.iloc[i, j]
                if value > 0:
                    ax[1].text(
                        j,i,confusion_matrix.iloc[
                            i, j],ha="center",va="center",color="black" if value < 60 else "white",fontsize=8
                    )

        ax[1].set_xticks( np.arange( 0,confusion_matrix.shape[0] ) + 0.5,minor=True )
        ax[1].set_yticks( np.arange( 0,confusion_matrix.shape[0] ) + 0.5,minor=True )

        ax[1].grid( which="minor",color="white",linewidth=1,zorder=100 )
        ax[1].tick_params( axis="both",which="minor",size=0 )

        ax[1].set_xlabel( "Actual lineage",fontweight="bold" )
        ax[1].set_ylabel( "Assigned lineage",fontweight="bold" )

        plt.tight_layout()
        plt.savefig( output.accuracy_plot )

        pp = res.groupby( "lineage_actual" )["parsimony_placements"].describe().reset_index()
        pp["numeric_lineage"] = pp["lineage_actual"].str.extract( r"(\d+)" ).fillna( 99 ).astype( int )
        pp = pp.sort_values( "numeric_lineage" ).reset_index( drop=True )

        fig, ax = plt.subplots( dpi=200,figsize=(5, 4) )

        points = ax.scatter( pp["50%"],pp.index,zorder=100 )
        points.set_clip_on( False )
        ax.hlines( pp.index,pp["min"],pp["max"],capstyle="round",zorder=99 )

        ax.set_yticks( range( pp.shape[0] ),pp["lineage_actual"].str.replace( "UNDEFINED","Other" ) )
        ax.set_xticks( range( 15 ),minor=True )

        ax.set_xlim( 0,15 )

        ax.set_xlabel( "Parsimoneous placements",fontweight="bold" )

        ax.grid( which="both",color="#EFEFEF" )

        plt.tight_layout()
        plt.savefig( output.parsimony_plot )


rule plot_freyja_accuracy:
    input:
        terra_results = TERRA_RESULTS
    output:
        accuracy_plot = "results/plots/freyja-accuracy-confusion-matrix.pdf"
    log:
        notebook = "results/notebooks/freyja-accuracy.ipynb"
    notebook: "../notebooks/plot_freyja_accuracy.py.ipynb"