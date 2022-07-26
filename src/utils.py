import os
import json
from argparse import Namespace
import subprocess
import pandas as pd


###############################
# read meta file; create args #
###############################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)
        
    args = Namespace(
        # from metadata file
        library_prefix = meta_dict[lib_name]["prefix"],
        library_reps = meta_dict[lib_name]["replicates"],
        library_pair= meta_dict[lib_name]["read_pairs"],
        library_umi = meta_dict[lib_name]["umi"],
        library_suffix = meta_dict[lib_name]["suffix"],
        library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["filtered"]
    )

    return args


###################
# filepath parser #
###################

def get_lib_peak_parsed_filepath(store_dir, lib_short):
    peak_filepath = os.path.join(
        store_dir, "lib_peak",lib_short, "peaks.bed"
        )
    return peak_filepath

def get_lib_dapeak_parsed_filepath(store_dir, lib_short, da_type):
    peak_filepath = os.path.join(
        store_dir, "diff_peak", lib_short, f"{da_type}.bed"
        )
    return peak_filepath

def get_enhancer_gene_mapped_file(store_dir, peak_desc, lib_short, method, da_type):
    return os.path.join(store_dir, peak_desc, lib_short, method, f"{da_type}.tsv")

def get_rnaseq_de_file(store_dir, lib_short):
    return os.path.join(store_dir, f"{lib_short}_vs_CC.tsv")


#####################
# peak to genes lfc #
#####################

def parse_rnaseq_de_file(rnaseq_de_file):
    df = pd.read_csv(rnaseq_de_file, sep="\t")
    df.columns = [c.strip('"') for c in df.columns]
    return df.loc[:, ["gene_symbol", "logFC", "FDR"]]

def parse_great_output(great_outfile):
    """
    Returns enhancer mapped to their genes
    """
    df = pd.read_csv(great_outfile, skiprows=3, sep="\t")
    df = df.loc[df.iloc[:, 0]=="Ensembl Genes"]
    df.Regions = df.Regions.str.split(",")
    df = df.explode("Regions", ignore_index=True)
    df = df.loc[:, ["Regions", "Genes"]]
    df[["chrm", "start"]] = df.Regions.str.split("-", expand=True)[0].str.split(":", expand=True)
    df["end"] = df.Regions.str.split("-", expand=True)[1]
    df = df.astype({"start": int, "end": int})
    return df.loc[:, ["chrm", "start", "end", "Genes"]]


def link_enhancers_to_genes(great_outfile, rnaseq_de_file):
    df_great = parse_great_output(great_outfile)
    df_rnaseq = parse_rnaseq_de_file(rnaseq_de_file)
    df = df_great.merge(df_rnaseq, left_on="Genes", right_on="gene_symbol").drop(columns=["gene_symbol"])
    return df.sort_values(["chrm", "start", "end"]).reset_index(drop=True)


#####################
# peak to peaks lfc #
#####################

def parse_peak_file(peak_file):
    df = pd.read_csv(peak_file, sep="\t", header=None, usecols=[0,1,2], names=["chrm", "start", "end"])
    return df

def parse_meta_rpp_file(rpp_file, lib_short, cc_short):
    df_rpp = pd.read_csv(rpp_file, index_col=[0,1,2])
    df_rpp.index = df_rpp.index.rename(["chrm", "start", "end"])
    required_columns = [c for c in df_rpp.columns if ((c.startswith(lib_short)) or (c.startswith(cc_short)))]
    df_rpp = df_rpp.loc[:, required_columns]
    return df_rpp.reset_index()
    
def link_enhancers_to_libwise_rpp(meta_rpp_file, peak_file, lib_short, cc_short):
    df_peak = parse_peak_file(peak_file)
    df_rpp = parse_meta_rpp_file(meta_rpp_file, lib_short, cc_short)
    df = df_peak.merge(df_rpp, on=["chrm", "start", "end"])
    return df


##################
# binding motifs #
##################

