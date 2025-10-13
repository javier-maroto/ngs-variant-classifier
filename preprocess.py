import pandas as pd
import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt


class Dataloader:

    MERGE_COLS = ['CHROM', 'POS', 'REF', 'ALT']
    VARL_COLS = []

    def __init__(self):
        pass

    def load_varlociraptor(self, inpdir):
        df_varlociraptor = pd.read_csv(f'{inpdir}/calls_filt.tsv', sep='\t')
        df_varlociraptor = self.create_pred_columns(df_varlociraptor)
        df_varlociraptor = df_varlociraptor.rename(columns={
            'seqnames': 'CHROM', 'position': 'POS', 'ref': 'REF', 'alt': 'ALT',
            'AF': 'AF_VARL', 'DP': 'DP_VARL'
        })
        df_varlociraptor['artifact_reason'] = df_varlociraptor.apply(
            Dataloader.create_artifact_reason_column, axis=1)
        return df_varlociraptor
    
    def load_varlociraptor_normal(self, inpdir):
        df_varlociraptor = pd.read_csv(f'{inpdir}/calls_filt_normal.tsv', sep='\t')
        df_varlociraptor = df_varlociraptor.rename(columns={
            'seqnames': 'CHROM', 'position': 'POS', 'ref': 'REF', 'alt': 'ALT',
            'AF': 'AF_NORMAL_VARL'
        })
        # df_varlociraptor = df_varlociraptor.drop_duplicates()
        return df_varlociraptor
    
    def load_vcf(self, inpdir):
        df_vcf = pd.read_csv(f'{inpdir}/candidates.tsv', sep='\t')
        df_vcf = df_vcf.rename(columns={
            'seqnames': 'CHROM', 'position': 'POS', 'ref': 'REF', 'alt': 'ALT'
        })
        # df_vcf['AF'] = df_vcf['AO'] / df_vcf['DP']
        return df_vcf

    def create_pred_columns(self, df):
        prob_cols = [col for col in df.columns if col.startswith('PROB_')]
        for col in prob_cols:
            df[col] = 10 ** -(df[col] / 10)
        class_labels = [col[5:] for col in prob_cols] + ['GERMLINE/SOMATIC']
        varlociraptor_label_mapping = {i: v for i, v in enumerate(class_labels)}
        preds = np.argmax(df[prob_cols], axis=1)
        df['pred'] = pd.Categorical(
            [varlociraptor_label_mapping[pred] for pred in preds],
            categories=class_labels
        )
        df['prob_variant'] = df['PROB_SOMATIC'] + df['PROB_GERMLINE']
        return df
    
    def create_artifact_reason_column(row):
        active_biases = (
            Dataloader.sb_bias(row) +
            Dataloader.rob_bias(row) +
            Dataloader.rpb_bias(row) +
            Dataloader.scb_bias(row) +
            Dataloader.he_bias(row) +
            Dataloader.alb_bias(row)
        )
        return ';'.join(active_biases)
        
    def sb_bias(row):
        if row['SB'] == '+':
            return ['Forward strand bias']
        elif row['SB'] == '-':
            return ['Reverse strand bias']
        return []
    
    def rob_bias(row):
        if row['ROB'] == '>':
            return ['F1R2 read orientation bias']
        elif row['ROB'] == '<':
            return ['F2R1 read orientation bias']
        return []
    
    def rpb_bias(row):
        if row['RPB'] == '^':
            return ['Read position bias (systematic sequencing errors)']
        return []
    
    def scb_bias(row):
        if row['SCB'] == '$':
            return ['Softclip bias (systematic alignment errors)']
        return []
        
    def he_bias(row):
        if row['HE'] == '*':
            return ['Homopolymer error (systematic PCR amplification errors)']
        return []
    
    def alb_bias(row):
        if row['ALB'] == '*':
            return ['Low MAPQ or major alternative alignment (XA tags)']
        return []


def create_scatter_plot_af(df):
    """
    Creates a scatter plot of AF_OBS_VARL (Estimated) vs AF (Observed)
    with a dashed reference line y = x.

    Parameters:
        df (pandas.DataFrame): Must contain columns 'AF_OBS_VARL' and 'AF'.

    Returns:
        matplotlib.figure.Figure: The figure object for saving or further editing.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(df['AF_OBS_VARL'], df['AF'], alpha=0.7, edgecolor='k')
    ax.plot([0, 1], [0, 1], '--', color='red', linewidth=1)

    ax.set_xlabel('Estimated')
    ax.set_ylabel('Observed')
    ax.set_title('Estimated vs Observed AF')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.6)

    return fig
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Load data with Dataloader")
    parser.add_argument(
        "--id",
        required=True,
        help="Identifier for varlociraptor data to load"
    )
    parser.add_argument(
        "--output_folder",
        required=True,
        help="Folder where all the outputs are saved"
    )
    parser.add_argument(
        "--tumor_cell_quantity",
        type=float,
        required=True,
        help="Tumor cell quantity (ranges from 0.0 to 1.0)"
    )
    params = parser.parse_args()

    inpdir = os.path.join(params.output_folder, params.id, "data")
    dl = Dataloader()
    df = dl.load_varlociraptor(inpdir)
    df_vcf = dl.load_vcf(inpdir)
    df = df[dl.MERGE_COLS + ['pred', 'prob_variant', 'artifact_reason'] + [col for col in df.columns if col.startswith('PROB_')] + ['AF_VARL', 'DP_VARL', 'SAOBS', 'SROBS']]
    df = df.merge(df_vcf[dl.MERGE_COLS + ['AF', 'AO', 'FAO', 'DP', 'FDP']], on=dl.MERGE_COLS)
    if os.path.exists(os.path.join(inpdir, "calls_filt_normal.tsv")):
        df = df.merge(dl.load_varlociraptor_normal(inpdir), on=dl.MERGE_COLS)
    else:
        df['AF_NORMAL_VARL'] = 0.0
    df['AF_OBS_VARL'] = df['AF_VARL'] * params.tumor_cell_quantity + df['AF_NORMAL_VARL'] * (1 - params.tumor_cell_quantity)
    df.loc[(df['pred'] == "SOMATIC") & (df['PROB_GERMLINE'] + 0.05 > df['PROB_SOMATIC']), 'pred'] = 'GERMLINE/SOMATIC'
    df.loc[(df['pred'] == "GERMLINE") & (df['PROB_SOMATIC'] + 0.05 > df['PROB_GERMLINE']), 'pred'] = 'GERMLINE/SOMATIC'
    outdir = os.path.join(params.output_folder, params.id, "results")
    df.to_csv(os.path.join(outdir, "varlociraptor_probs.tsv"), sep="\t", index=False)

    fig = create_scatter_plot_af(df)
    fig.savefig(os.path.join(outdir, 'estimated_vs_observed_AF.png'), dpi=300, bbox_inches='tight')

    print(df.to_string(index=False))
    n_higher_than_1_estimated_AF = len(df.loc[(df['AF_VARL'] > 0.95) & (df['AF'] > df['AF_OBS_VARL'] * 1.05) & (df['PROB_GERMLINE'] < 0.1)])
    if n_higher_than_1_estimated_AF > 0:
        print(f'Warning: tumor cell content prior might be too low (based on {n_higher_than_1_estimated_AF} calls)', file=sys.stderr)