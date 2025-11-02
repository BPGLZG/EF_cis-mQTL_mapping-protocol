
### Mapping with TensorQTL
## 6.2. Compute mQTLs with RNT values

# To run this script, the following conda environment was activated:
# /data/genotools/conda_envs/tensorqtl/

#Working directory
#"/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_perm")

# Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans

print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))
print('TensorQTL {}'.format(tensorqtl.__version__))

# Set number of threads for tensorqtl
num_threads = 16
torch.set_num_threads(num_threads)

 #Define paths to data
plink_prefix_path = '/data/projects/endometriosis/tensorQTL_final/final_20250210/results/imp_130_final/imp_130_final'
expression_bed = '/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_130_met.bed.gz'
covariates_file = '/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_covs.tsv'
prefix = '/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_perm/cis_mqtls_130_perm_20250210'

#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Sort phenotype sample names
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

#Run TensorQTL - permuted
cis_df = cis.map_cis(genotype_df,
                     variant_df,
                     phenotype_df,
                     phenotype_pos_df,
                     covariates_df = covariates_df,
                     seed = 2024,
                     window = 500000)
                  
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)

cis_df.to_csv('/data/projects/endometriosis/tensorQTL_final/final_20250210/results/tensorqtl_perm/tensorQTL_130_20250210_rnt_perm.txt.gz', header = True, index = True, sep = '\t')

 
