import sys
import pandas as pd

sampname = sys.argv[1]

def comb_anno(sampname, path='./data/annotated/', vcf_path='./data/vcf/', out_path='./data/annotated/'):
    ### parse raw data
    #samp = sampname.split('_')[0]
    tokb = pd.read_csv(path+sampname+'_oncokb_maf.txt', sep='\t', skiprows=[1], header=0)
    tav = pd.read_csv(path+sampname+'.hg19_multianno.txt', sep='\t', usecols=list(range(41)), index_col=False)
    tvcf = pd.read_csv(vcf_path+sampname+'_gatkfilt_vanfilt.vcf', sep='\t', comment='#', header=None)
    tvcf.columns=['CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'samp']

    ### extract VAF AD DP from VCF
    tad = [int(i.split(':')[1].split(',')[1]) for i in tvcf.samp]
    tdp = [int(i.split(':')[3]) for i in tvcf.samp]
    tvaf = [round(a/d, 3) for a,d in zip(tad,tdp)]
    tvcf = tvcf.assign(VAF=tvaf)
    tvcf = tvcf.assign(AD=tad)
    tvcf = tvcf.assign(DP=tdp)

    ### merge dfs
    tokb_vid_cols = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
    tav_vid_cols = ['Chr', 'Start', 'Ref', 'Alt']
    tokb_vids = gen_vid(tokb, tokb_vid_cols)
    tav_vids = gen_vid(tav, tav_vid_cols)
    tvcf_vids = gen_vcf_vid(tvcf)
    tokb = tokb.assign(vid=tokb_vids)
    tav = tav.assign(vid=tav_vids)
    tvcf = tvcf.assign(vid=tvcf_vids)
    mdf = pd.merge(tokb, tav, on='vid', how='outer')
    mdf = pd.merge(mdf, tvcf, on='vid', how='outer')

    ### reduced columns
    cols = ['Chromosome', 'Start_Position', 'End_Position',  'Reference_Allele', 'Tumor_Seq_Allele2', 'VAF', 'AD', 'DP', 'Hugo_Symbol', 'Variant_Type', 'Consequence', 'HGVSc', 'HGVSp_Short', 'RefSeq', 'dbSNP_RS', 'CLIN_SIG', 'InterVar_automated', 'IMPACT', 'mutation_effect', 'oncogenic', 'Highest_level']
    tmdf = mdf[cols]
    tmdf.columns = ['Chr', 'Start', 'End',  'Ref', 'Alt', 'VAF', 'AD', 'DP', 'Gene', 'Variant_Type', 'Consequence', 'HGVSc', 'HGVSp', 'RefSeq', 'dbSNP_RS', 'ClinVar_Sig', 'InterVar_automated', 'VEP_Impact', 'VEP_mutation_effect', 'VEP_oncogenic', 'OncoKB_Level']

    ### filter and sort
    tmdf = tmdf[tmdf.Variant_Type=='SNP']
    oncokb_dict = {'LEVEL_1':1, 'LEVEL_2':2, 'LEVEL_3A':3, 'LEVEL_3B':3, 'LEVEL_4':4, 'LEVEL_R1':1, 'LEVEL_R2':2, '':4}
    vep_impact_dict = {'LOW':3, 'MODIFER':2, 'MODERATE':2, 'HIGH':1, '':4}
    vep_oncogenic_dict = {'Unknown':3, 'Inconclusive':3, 'Predicted Oncogenic':2, 'Likely Oncogenic':2, 'Oncogenic':1, '':4}
    VEP_mutation_effect_dict = {'Unknown':3, 'Inconclusive':3, 'Likely Neutral':3, 'Likely Gain-of-function':2, 'Likely Loss-of-function':2, 'Gain-of-function':1, 'Loss-of-function':1, '':4}

    tmdf = tmdf.assign(oncokb_cat = [oncokb_dict[i] for i in tmdf['OncoKB_Level']])
    tmdf = tmdf.assign(vep_impact_cat = [vep_impact_dict[i] for i in tmdf['VEP_Impact']])
    tmdf = tmdf.assign(vep_oncogenic_cat = [vep_oncogenic_dict[i] for i in tmdf['VEP_oncogenic']])
    tmdf = tmdf.assign(vep_mutation_effect_cat = [vep_mutation_effect_dict[i] for i in tmdf['VEP_mutation_effect']])

    tmdf = tmdf.sort_values(by=['VAF', 'oncokb_cat', 'vep_impact_cat', 'vep_oncogenic_cat', 'vep_mutation_effect_cat'], ascending=[False, True, True, True, True])
    tmdf = tmdf.drop(columns=['oncokb_cat', 'vep_impact_cat', 'vep_oncogenic_cat', 'vep_mutation_effect_cat'])
    tmdf.to_csv(out_path+sampname+'_annotated.txt', sep='\t', index=False)

    return

comb_anno(sampname)
