
#Author: David J. Bradshaw, Ph.D. 
#Emails: dbradshaw2015@fau.edu, dbradshaw3366@gmail.com

conda create -n rescript -c conda-forge -c bioconda -c qiime2 -c defaults \
qiime2=2021.2.0=py36_0 q2cli=2021.2.0=py36_0 q2templates=2021.2.0=py36_0 q2-types=2021.2.0=py36_0 \
q2-longitudinal=2021.2.0=py36_0 q2-feature-classifier=2021.2.0=py36_0 "pandas>=0.25.3" xmltodict=0.12.0=py_0
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
conda install mamba=0.7.14=py36h05d92e0_0
mamba install -c conda-forge -c bioconda snakemake


conda create -n rescript \
  -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-longitudinal q2-feature-classifier "pandas>=0.25.3" xmltodict

#Get SILVA files
qiime rescript get-silva-data \
--p-version '138.1' \
--p-target 'SSURef_NR99' \
--p-include-species-labels \
--o-silva-sequences silva-138.1-ssu-nr99-seqs.qza \
--o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

qiime rescript cull-seqs \
--p-n-jobs 4 \
--i-sequences silva-138.1-ssu-nr99-seqs.qza \
--o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
--i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
--i-taxonomy silva-138.1-ssu-nr99-tax.qza \
--p-labels Archaea Bacteria Eukaryota \
--p-min-lens 900 1200 1400 \
--o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
--o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

qiime rescript dereplicate \
--i-sequences silva-138.1-ssu-nr99-seqs-filt.qza \
--i-taxa silva-138.1-ssu-nr99-tax.qza \
--p-rank-handles 'silva' \
--p-mode 'uniq' \
--o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
--o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

qiime feature-classifier extract-reads \
--i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza  \
--p-f-primer CCTACGGGNGGCWGCAG \
--p-r-primer GACTACHVGGGTATCTAATCC \
--p-n-jobs 4 \
--p-read-orientation 'forward' \
--o-reads silva-138.1-ssu-nr99-341f-805r-seqs.qza

qiime rescript dereplicate \
--i-sequences silva-138.1-ssu-nr99-341f-805r-seqs.qza \
--i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
--p-rank-handles 'silva' \
--p-mode 'super' \
--o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-341f-805r-derep-super.qza \
--o-dereplicated-taxa silva-138.1-ssu-nr99-tax-341f-805r-derep-super.qza

qiime rescript evaluate-fit-classifier \
--i-sequences silva-138.1-ssu-nr99-seqs-341f-805r-derep-super.qza \
--i-taxonomy silva-138.1-ssu-nr99-tax-341f-805r-derep-super.qza \
--o-classifier silva-138.1-99-341f-805r-2021.8-classifier.qza \
--o-observed-taxonomy silva-138-99-341f-805r--derep-super-taxonomy-predicted-taxonomy.qza \
--o-evaluation silva-138-99-341f-805r--derep-super-taxonomy-fit-classifier-evaluation.qzv \
--p-reads-per-batch 10000
