import csv
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

def write_log(*args):
    with open(args[0], "a") as log_file:
        line = ' '.join([str(a) for a in args[1:]])
        log_file.write(line+'\n')
        print(line)


def parse_parameters_comma_separate(parameter): # in future we can use any delimiter
    new_list = []
    parameter_split = parameter.split(",")
    for para in parameter_split:
        new_list.append(para.strip())
    return new_list

# make a function to compute the de up or down regulated
def degs_up_down(df, pval_col="P.Value", fc_col= "logFC", p_val_cutoff = 0.05, logFC_cutoff=0.58, gene_col="geneSymbol"): #0.58 = 1.5 for log2FC
  up_deg = (df[(df[fc_col] > logFC_cutoff) & (df[pval_col] < p_val_cutoff)])[gene_col].tolist()
  down_deg = (df[(df[fc_col] < logFC_cutoff) & (df[pval_col] < p_val_cutoff)])[gene_col].tolist()
  return up_deg, down_deg



def make_plots(gene_sets,enr_obj, outdir, cmap="viridis",column="Adjusted P-value",title = "Downregulated in Rubicon KO vs WT",figsize=(4,12), size=3, top_term=10): #enr_obj = the result from analyze_plot_gsea ; for top 10 term default should be good #use column = P-value for pvalues
  # Choose a colormap
  cmap = plt.get_cmap(cmap)

  # Generate colors from the colormap
  colors = cmap(np.linspace(0, 1, len(gene_sets)))

  # Create a color dictionary mapping each element to a color
  color_dict = {element: color for element, color in zip(gene_sets, colors)}

  ax = dotplot(enr_obj.results,
              column=column,
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=size,
              top_term=top_term,
              figsize=figsize,
              title = title,
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
              ofname=f'{outdir}/categorical_scatter_dotplot_{title}.png'
             )
  ax = barplot(enr_obj.results,
              column=column,
              group='Gene_set', # set group, so you could do a multi-sample/library comparsion
              size=size,
              top_term=top_term,
              figsize=figsize,
              #color=['darkred', 'darkblue'] # set colors for group
              color = color_dict,
              ofname=f'{outdir}/horizontal_barplot_{title}.png'
             )
  # to save your figure, make sure that ``ofname`` is not None
  ax = dotplot(enr_obj.res2d, title=f'{title}/n{"_".join(gene_sets)}',cmap=cmap,size=size, column=column, figsize=figsize,ofname=f'{outdir}/dotplot_combined_score_{title}.png')

  ax = barplot(enr_obj.res2d,title=f'{title}/n{"_".join(gene_sets)}', figsize=figsize,ofname=f'{outdir}/barplot_res2d_{title}.png', color='darkred')

#make_plots(gene_sets,enr_obj, outdir, cmap="viridis",column="Adjusted P-value",title = "Downregulated in Rubicon KO vs WT",figsize=(4,12), size=3, top_term=10)
#analyze_plot_gsea(gene_sets,background="background.txt", save_dir='Young_SKO_Rubicon_up')

def analyze_plot_gsea(gene_list, gene_set,cutoff=0.05,background="background.txt", save_dir='Young_SKO_Rubicon_up'):
  # Assuming 'variable' is the variable you want to check
  is_list = isinstance(gene_set, list)
  if is_list:
    folder_result = f'{save_dir}_{"_".join(gene_set)}'
  else:
    folder_result = f'{save_dir}_{gene_set}'
  # check if the gene_set can be changed to either one
  #DEGs_up_1d
  # backgound only reconigized a gene list input.
  try:
    enr_bg = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                    gene_sets=gene_set,
                    # organism='human', # organism argment is ignored because user input a background
                    background=background,
                    outdir=folder_result, #None =  don't write to disk
                    cutoff=cutoff,
                    )
    print(f'...Displaying the top 10 enriched results for {gset}')
    print(enr_bg.results.head(10))

    return enr_bg,folder_result


  except:
    print(f'ValueError: Warning: No enrich terms when cutoff = 0.05 for {gene_set}')
    return None, None


"""### Prepare input file"""

def prepare_prerank_input(df,pval_col="P.Value", fc_col= "logFC"):
  df2=df.copy()
  # make gene symbol upper
  df2['rank'] = -np.log10(df2[pval_col]) * np.sign(df2[fc_col])

  # df2["rank"] = df2["P.Value"].apply(lambda x: -1*np.log10(x))
  df2["gene_Upper"] = df2.geneSymbol.apply(lambda x: str.upper(x))
  rnk_up=pd.DataFrame(list(df2["rank"]), index=df2["gene_Upper"])

  rnk_up=rnk_up.reset_index()
  #rename columns
  rnk_up.columns=[0,1]
  rnk_up=rnk_up.set_index(0)
  #avoid duplicate genes
  rnk_up = rnk_up.groupby(rnk_up.index).mean()
  rnk_up = rnk_up.sort_values(by=[1], ascending=False)
  print(f'... the input dataframe or file should be only 2 columns index is 0 and rank score is 1')
  print(rnk_up.head())
  return rnk_up

"""Generate a dataframe using the function"""


"""### Analyze using prerank"""

from sys import maxsize
def analyze_prerank(rank_df,gene_set, save_dir='Young_SKO_Rubicon_prerank',threads=4, minsize=5, maxsize=1000):

  # # run prerank
  # # enrichr libraries are supported by prerank module. Just provide the name
  # # use 4 process to acceralate the permutation speed
  try:
    pre_res = gp.prerank(rnk=rank_df, # or rnk = rnk,
                     gene_sets=gene_set,
                     threads=threads,
                     min_size=minsize,
                     max_size=maxsize,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=f'{save_dir}', # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )

    print(f'...Displaying the top 10 enriched results for {save_dir}')
    print(pre_res.res2d.head(10))

    return pre_res,f'{save_dir}'


  except:
    print(f'ValueError: Warning: No enrich terms when cutoff = 0.05 for {gene_set}')
    return None, None



