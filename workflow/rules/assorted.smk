rule plot_lineage_defining_mutations:
    input:
        barcodes = "/PHShome/nm104/scripts/vibecheck/vibecheck/resources/o1_barcodes.feather"
    output:
        barcode_plot = "results/plots/barcode_plots.pdf"
    log:
        notebook = "results/notebooks/plot_lineage_defining_mutations.ipynb"
    notebook: "../notebooks/plot_lineage_defining_mutations.py.ipynb"