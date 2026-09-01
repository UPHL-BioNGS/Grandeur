import pandas as pd

# quast : combining both files
q_df  = pd.DataFrame()
qc_df = pd.DataFrame()
def summarize_quast_from_reads(summary_df, quast):
    print("Adding results for " + quast)
    file = quast
    analysis = str(file).split("_")[0]
    q_df = pd.read_table(file, dtype = str, index_col= False)
    q_df = q_df.add_prefix(analysis + "_")
    q_df.columns = [x.lower() for x in q_df.columns]

def summarize_quast_contig(summary_df, quast_contig):
    print("Adding results for " + quast_contig)
    file = quast_contig
    analysis = str(file).split("_")[0]
    qc_df = pd.read_table(file, dtype = str, index_col= False)
    qc_df = qc_df.add_prefix(analysis + "_")
    qc_df.columns = [x.lower() for x in qc_df.columns]

def summarize_quast(summary_df, quast_reads, quast_contigs):
    q_df  = pd.DataFrame()
    qc_df = pd.DataFrame()
    analysis = "quast"
    new_df = pd.concat([q_df, qc_df])
    new_df[analysis + "_sample"] = new_df[analysis + "_sample"].astype(str)

    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how='left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)

    return summary_df