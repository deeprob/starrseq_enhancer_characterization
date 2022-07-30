import argparse
import os
import pandas as pd
import utils as ut


def main(
    lib_short,
    peak_dir,
    enhancer_gene_dir,
    rnaseq_dir,
    chip_dir,
    rpp_file,
    tf_bind_file,
    diff_activity_type,
    store_dir
    ):

    # get enhancer activity df
    peak_file = ut.get_lib_dapeak_filepath(peak_dir, lib_short, diff_activity_type)
    enhancer_activity_df = ut.link_enhancers_to_libwise_rpp(rpp_file, peak_file, lib_short, "CC").set_index(["chrom", "start", "end"])
    # enhancer gene mapped df
    target_gene_file = ut.get_enhancer_gene_mapped_file(enhancer_gene_dir, "diff_peak", lib_short, "great", diff_activity_type)
    rna_file = os.path.join(rnaseq_dir, f"{lib_short}_vs_CC.tsv")
    enhancer_gene_df = ut.link_enhancers_to_genes(target_gene_file, rna_file)
    # binding motif df
    enhancer_motif_df = ut.link_enhancers_to_binding_motifs(tf_bind_file, peak_file)
    # histone df
    histones = sorted([f.name for f in os.scandir(os.path.join(chip_dir, "histone"))])
    enhancer_histone_df = ut.link_enhancers_to_multiple_chips(chip_dir, "histone", histones, peak_file)
    # tf df
    tfs = sorted([f.name for f in os.scandir(os.path.join(chip_dir, "tf"))])
    enhancer_tf_df = ut.link_enhancers_to_multiple_chips(chip_dir, "tf", tfs, peak_file)
    # concatenated meta df
    meta_df = pd.concat((enhancer_activity_df, enhancer_gene_df, enhancer_motif_df, enhancer_tf_df, enhancer_histone_df), axis=1)
    store_file = os.path.join(store_dir, lib_short, diff_activity_type, "enhancer_meta.csv")
    os.makedirs(os.path.dirname(store_file))
    meta_df.to_csv(store_file, index=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq MEA analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("peak_dir", type=str, help="Dir where the library peaks are stored")
    parser.add_argument("rnaseq_dir", type=str, help="Dir where the rnaseq de results are stored")
    parser.add_argument("enhancer_gene_dir", type=str, help="Dir where the target genes mapped to enhancers are stored")
    parser.add_argument("chip_dir", type=str, help="Dir where theencode chipseq files are stored")
    parser.add_argument("rpp_file", type=str, help="The meta reads per plasmid matrix file where rpp info for all enhancers across all libraries is stored")
    parser.add_argument("tf_bind_file", type=str, help="The meta file that contains tf binding information for all enhancers across all libraries")
    parser.add_argument("store_dir", type=str, help="Dir to store great results")
    parser.add_argument("diff_activity_type", type=str, help="type of differential enhancer activity, use this argument to compare induced,repressed and constitutive peaks between lib1 and control")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)

    main(
        lib_args.library_short,
        cli_args.peak_dir,
        cli_args.enhancer_gene_dir,
        cli_args.rnaseq_dir,
        cli_args.chip_dir,
        cli_args.rpp_file,
        cli_args.tf_bind_file,
        cli_args.diff_activity_type,
        cli_args.store_dir
    )
