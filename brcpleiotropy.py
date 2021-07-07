#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 13:04:15 2021

@author: meredithcarpenter
"""

import glob, os   
import pandas as pd
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
from collections import Counter
import collections
import math
from Bio.PDB.Polypeptide import *
from pyrosetta import *
from pyrosetta.toolbox import *
init()

#Removing any value from list containing val and returning a list exclusing those values
def remove_values(the_list, val):
    return[value for value in the_list if val.casefold() not in value.casefold()]

def remove(list): 
    pattern = '[A-Z,a-z,*]'
    list = [re.sub(pattern, '', i) for i in list] 
    return list

def shannon_entropy(sequence):
    m = len(sequence)
    bases = collections.Counter([tmp_base for tmp_base in sequence])
 
    shannon_entropy_value = 0
    for base in bases:
        # number of residues
        n_i = bases[base]
        # n_i (# residues type i) / M (# residues in column)
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i
 
    return shannon_entropy_value * -1

def eta(data, unit='natural'):
    base = {
        'shannon' : 2.,
        'natural' : math.exp(1),
        'hartley' : 10.
    }

    if len(data) <= 1:
        return 0

    counts = Counter()

    for d in data:
        counts[d] += 1

    ent = 0

    probs = [float(c) / len(data) for c in counts.values()]
    for p in probs:
        if p > 0.:
            ent -= p * math.log(p, base[unit])

    return ent

pd.set_option('display.max_colwidth',None)

#Dataframe from clinvar snps associated with breast cancer
brdf = pd.read_csv('./diseases/clinvar_breastcancer.csv')
brdf.columns = ['name', 'genes', 'proteins', 'conditions', 'clinsig', 'review', 'accession', 'grch37chr', 'grch37loc', 'grch38chr', 'grch38loc','varid', 'alleleid', 'dbsnp', 'spdi']
brdf = brdf[(brdf['conditions'].str.contains(r'(?=.*breast)', regex=True, case=False, na=False))]
len(brdf)
brcgenes = brdf['genes'].to_list()
brcgenes = set(brcgenes)
len(brcgenes)

#only considering variants with dbsnp identifiers
brdf_dbsnp = brdf.dropna(subset=['dbsnp'])
brcgenes_dbsnp = brdf_dbsnp['genes'].to_list()
brcgenes_dbsnp = set(brcgenes_dbsnp)
len(brcgenes_dbsnp)

#Complete list of genes associated with breast cancer (split on multiple listings)
splitgs = set()
for gene in brcgenes_dbsnp:
    gene = str(gene)
    if "|" in gene:
        gsplit = gene.split("|")
        splitgs.update(gsplit)
    else:
        splitgs.add(gene)
len(splitgs)

#Secondary diseases list:
dis_set = ['Fanconi anemia', 'Aplastic anemia', 'Thrombocytopenia', 'Dyskeratosis congenita',
              'Nephrolithiasis', 'Cystinuria', 'Renal hypodysplasia', 'Renal aplasia', 'Microcephaly',
              'Macrocephaly', 'Costello', 'Hemimegalencephaly', 'Lissencephaly', 'Nijmegen', 'Ataxia',
              'Cockayne', 'Xeroderma pigmentosum', 'Blepharocheilodontic', 'Smith-Kingsmore', 'Crouzon',
              'Cortical dysplasia', 'Polymicrogyria', 'Hirschsprung', 'Hutchinson-Gilford',
              'Megalencephaly-capillary', 'Cardiofaciocutaneous', 'Seizure', 'Tracheoesophageal Fistula',
              'VACTERL', 'Diaphyseal dysplasia', 'Dactyl', 'Noonan', 'Tuberous sclerosis', 'Thyroidism',
              'Immunodeficiency', 'Bloom', 'Cystic Fibrosis', 'Pulmonary Fibrosis']

#Identify SNPs associated with both breast cancer and disease
dis = '|'.join(dis_set)
pleios = brdf_dbsnp[brdf_dbsnp['conditions'].str.contains(dis, case=False)]
len(pleios)

#Keep only the first gene name assocaited with the variants
firstgene = pleios['genes'].str.split("|", n=1, expand = True)[0]
pleios['genes'] = firstgene
pleios = pleios.replace('C11orf65','ATM')
pleios = pleios.replace('AOPEP','FANCC')

#Identify genes associted with pleiotropic snps
pgenes = pleios['genes'].unique()

###list of genes for heatmap = pgenes###

####GENERATE DENSITY PLOTS FOR ALL SNPS####
#Complete disease data frame from all clinvar snps associated with the listed secondary diseases
diseases = pd.concat(map(pd.read_csv, glob.glob(os.path.join('./diseases', "*.csv"))), ignore_index=True)
diseases.columns = ['name', 'genes', 'proteins', 'conditions', 'clinsig', 'review', 'accession', 'grch37chr', 'grch37loc', 'grch38chr', 'grch38loc','varid', 'alleleid', 'dbsnp', 'spdi']
diseases = diseases.drop_duplicates()

#Create a new dataframe with only the genes in 'genes' 
genes_re = re.compile("|".join(pgenes))
genedf = diseases[diseases['genes'].str.contains(genes_re, na=False)]

for i in genedf:
    loc_list = genedf['grch38loc'].tolist()
    modlocs = list()
    for location in loc_list:
        if isinstance(location, int):
            modlocs.append(location)
        elif isinstance(location, float):
            modlocs.append(location)
        elif "-" in location:
            low, high = location.split("-",1)
            low.strip()
            high.strip()
            low = int(low)
            high = int(high)
            location = (high + low) / 2
            modlocs.append(location)
        else:
            modlocs.append(location)

#new dataframe with just the subset of information needed for the density plots
gene_list = genedf['genes'].tolist()
conditions_list = genedf['conditions'].tolist()
densitydf = pd.DataFrame(list(zip(modlocs,gene_list,conditions_list)), columns = ['location', 'gene', 'conditions'])

#Search density data frame by gene, determine number of snps for each condition at each location, build new df
#New dataframe: location, condition, count
brc = 'breast'
hcp = 'hereditary cancer-predisposing'

cancer_set = ['Breast', 'Other Cancer']
dis_cancer_set = cancer_set + dis_set


for i in pgenes:
    bygene = densitydf[densitydf['gene'].str.contains(i)]
    bygene = bygene.dropna()
    countdf = pd.DataFrame(columns = ['location', 'condition'])
    for index,row in bygene.iterrows():
        for dis in dis_cancer_set:
            if dis == 'Other Cancer':
                matches = ["cancer", "Cancer", "carcinoma", "Carcinoma", "tumor", "Tumor", "oma"]
                if any(x in row['conditions'] for x in matches) and not (brc.casefold() in row['conditions']) and not (hcp.casefold() in row['conditions']):
                    newdf = pd.DataFrame({"location":[row['location']], "condition":'Other Cancer'})
                    countdf = countdf.append(newdf, ignore_index=True)
            elif dis == 'Breast':
                if dis.casefold() in row['conditions'].casefold():
                    newdf = pd.DataFrame({"location":[row['location']], "condition":'Breast Cancer'})
                    countdf = countdf.append(newdf, ignore_index=True)
            elif dis.casefold() in row['conditions'].casefold():
                newdf = pd.DataFrame({"location":[row['location']], "condition":dis})
                countdf = countdf.append(newdf, ignore_index=True)
    countdf = countdf.dropna()
    countdf['location'] = pd.to_numeric(countdf['location'])
    print("complete", i)
    if countdf.empty:
        continue
    else:
        cond_list = countdf['condition'].unique() #generate list of the conditions in the dataframe to be plotted
        #iterate through the conditions
        for condition in cond_list:
            #Subsent to the condition
            subset = countdf[countdf['condition'] == condition]
            
            #Draw the density plot
            agd = sns.distplot(subset['location'], hist = False, kde = True, kde_kws = {'linewidth':2}, label=condition)
        
        #Plot formatting
        plt.legend(prop={'size':16}, title = 'Condition')
        plt.title(i)
        plt.xlabel('Chromosome location (GRCH38)')
        plt.ylabel('Density')
        plt.show()
        figname = i + '.png'
        agd.figure.savefig(figname, dpi=100)



#####GENERATE HEATMAP FOR COLOCALZIED SNPs####

#Iterate over list, search for matches and count number of times each disease appears on each gene, save to dataframe merged by gene

newallgenecounts = pd.DataFrame(pgenes, columns=['Gene'])
for i in dis_cancer_set:
    if i == "Breast":
        genecount = pleios['genes'].value_counts().rename_axis('Gene').reset_index(name='Breast Cancer')
        newallgenecounts = newallgenecounts.merge(genecount, how='left', on=['Gene'])
    elif i == "Other Cancer":
        modcons = list()
        for entry in pleios['conditions']:
            entry = str(entry)
            lst = entry.split("|")
            nob = [x for x in lst if "breast" not in x]
            nob = [x for x in nob if "Breast" not in x]
            nohcp = [x for x in nob if "Hereditary cancer-predisposing" not in x]
            nohcp_join = ''.join(nohcp)
            modcons.append(nohcp_join)
        pleios['modcons'] = modcons
        subset = pleios[(pleios['modcons'].str.contains(r'(?=.*cancer)', regex=True, case=False, na=False)) | (pleios['modcons'].str.contains(r'(?=.*carcinoma)', regex=True, case=False, na=False)) | (pleios['modcons'].str.contains(r'(?=.*tumor)', regex=True, case=False, na=False)) | (pleios['modcons'].str.contains(r'(?=.*oma)', regex=True, case=False, na=False))]
        genecount = subset['genes'].value_counts().rename_axis('Gene').reset_index(name='Other Cancer')
        newallgenecounts = newallgenecounts.merge(genecount, how='left', on=['Gene'])
    else:
        subset = pleios[pleios['conditions'].str.contains(i, case=False, na=False)]
        genecount = subset['genes'].value_counts().rename_axis('Gene').reset_index(name=i)
        newallgenecounts = newallgenecounts.merge(genecount, how='left', on=['Gene'])

#Remove genes with zero counts
newallgenecounts = newallgenecounts.loc[(newallgenecounts.sum(axis=1) != 0)]
#Set index as gene names
newallgenecounts = newallgenecounts.set_index('Gene')
#Replace NaN with 0s
newallgenecounts = newallgenecounts.fillna(0)
#Sort gene rows by sum of row (total counts)
newallgenecounts = (newallgenecounts.assign(sum=newallgenecounts.sum(axis=1)).sort_values(by='sum', ascending=False).iloc[:,:-1])

#Generate heatmap
sns.set_theme()
sns.set(rc={'figure.figsize':(20,15)})
agc = sns.heatmap(newallgenecounts, cmap="coolwarm", annot=False, robust=True, fmt="1g")
#plt.title('Heatmap of All Disease SNPs on Key Genes')
agc.figure.savefig('agcall.png', dpi=1000) ##To save heatmap file

###List of SNPs and locations for MendelVar
###Mendelvar - coassociated SNPs, location interval
#pleios df = snps assocaited with brc and additional condition
#pleiosnps = pleios.drop_duplicates(subset="dbsnp")

##Generate location interval
loc_list = pleios['grch38loc'].tolist()
loc1 = list()
loc2 = list()
for location in loc_list:
    if isinstance(location, int):
        loc1.append(location)
        loc2.append(location)
    elif isinstance(location, float):
        loc1.append(location)
        loc2.append(location)
    elif "-" in location:
        low, high = location.split("-",1)
        low.strip()
        high.strip()
        low = int(low)
        high = int(high)
        loc1.append(low)
        loc2.append(high)
    else:
        loc1.append(location)
        loc2.append(location)

chr_list = pleios['grch38chr'].tolist()
string = "chr"
chrs = [string + x for x in chr_list]
db_list = pleios['dbsnp'].tolist()
snpdf = pd.DataFrame(list(zip(chr_list,loc1,loc2,db_list)), columns = ['chromosome', 'location1', 'location2', 'dbsnp'])
snpdf.to_csv('All SNPs', header=False, index=False, sep=" ")

###PROTEIN SEQUENCE ALIGNMENTS

### GENERATION OF DATA SET FOR SNPS WITH BREAST CANCER + DISEASE : CO-LOCALIZED SNPS
#Data set = pleiosnps

sfp = pleios[['genes','grch38loc','proteins','dbsnp']]
sfp.dropna(subset = ['proteins'], inplace=True) #drop variants missing protein change data
sfp = sfp[~sfp['proteins'].str.contains('[a-z]')] #drop variants with fs and del mutations
sfp = sfp.reset_index(drop=True) #reset index
sfp = sfp.sort_values('genes') #sort by gene
sfp_genes = sfp['genes'].unique() #genes with in-frame mutations

#splitting protein change into start AA, location, end AA
AA_old = list()
AA_pos = list()
AA_new = list()

for i in sfp['proteins']:
    proteins = re.split(',', i)
    if len(proteins) == 1:
        res = re.split('(\d+)',i)
        pos = int(res[1])
        AA_old.append(res[0])
        AA_pos.append(pos)
        AA_new.append(res[2])
    else:    
        residues = remove(proteins)
        residues = [x.strip(' ') for x in residues]
        residues.sort(key=int)
        for j in proteins:
            if residues[-1] in j:
                res = re.split('(\d+)',j)
                res = [x.strip(' ') for x in res]
                pos = int(res[1])
                AA_old.append(res[0])
                AA_pos.append(pos)
                AA_new.append(res[2])
                break

sfp_g = sfp['genes'].tolist()
sfp_d = sfp['dbsnp'].tolist()
sfp_c = sfp['grch38loc'].tolist()

pro_sfp = pd.DataFrame(list(zip(sfp_g, sfp_d, sfp_c, AA_pos, AA_old, AA_new)), columns = ['gene', 'dbsnp', 'grch38loc','AA_pos', 'wt_AA', 'mut_AA'])
##Make necessary corrections
#BRCA1 rs587782026 AA_pos to 1824
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs587782026', 'AA_pos'] = 1824
#BRCA1 rs1555579648 AA_pos to 1686
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1555579648', 'AA_pos'] = 1686
#BRCA1 rs1265352633 AA_pos to 1552
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1265352633', 'AA_pos'] = 1552
#BRCA1 rs45553935 AA_pos to 1736
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs45553935', 'AA_pos'] = 1736
#BRCA1 rs80357123 AA_pos to 1751
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs80357123', 'AA_pos'] = 1751
#BRCA1 rs80356885 AA_pos to 1508
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs80356885', 'AA_pos'] = 1508
#BRCA1 rs80356962 AA_pos to 1815
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs80356962', 'AA_pos'] = 1815
#BRCA1 rs55770810 AA_pos to 1699
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs55770810', 'AA_pos'] = 1699
#PTEN rs121909219 AA_pos to 233
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs121909219', 'AA_pos'] = 233
#CASR rs1801725 AA_pos to 986
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1801725', 'AA_pos'] = 986
#TGFB1 rs1800470 AA_old to L and AA_new to P
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1800470', 'wt_AA'] = 'L'
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1800470', 'mut_AA'] = 'P'
#FGFR2 rs1057519045 AA_pos to 549
pro_sfp.loc[pro_sfp['dbsnp'] == 'rs1057519045', 'AA_pos'] = 549


snpres = pro_sfp.drop(columns=['grch38loc'])
snpres.to_csv('snp_residues.csv', index=False) #save SNP list file

for filename in os.listdir('./fastas'):
    if filename.endswith(".fasta"):
        genename = filename.rsplit('.',1)[0]
        filepath = os.path.join('./fastas',filename)
        snp_proteins = []
        for seq_record in SeqIO.parse(filepath, 'fasta'):
            snp_proteins.append(seq_record)
        for index, row in pro_sfp.iterrows():
            if row['gene'] == genename:
                pos = int(row['AA_pos']) - 1
                for seq_record in SeqIO.parse(filepath, 'fasta'):
                    if pos < len(seq_record) and seq_record[pos] == row['wt_AA']:
                        if row['mut_AA'] == 'fs':
                            continue
                        else:
                            seq_new = seq_record
                            seq_new = seq_new[:pos] + row['mut_AA'] + seq_new[pos+1:]
                            seq_new.id = str(row['gene']) + '|' + str(row['dbsnp'])
                            seq_new.name = str(row['gene']) + '|' + str(row['dbsnp'])
                            seq_new.description = str(row['gene']) + '|' + str(row['dbsnp'])
                            snp_proteins.append(seq_new)
                    else:
                        print('WT AA does not match isoform entry for SNP:') ##warning for mis-aligned SNPS
                        print(row['gene'],row['dbsnp'])
                        
        saveas = "coassoc_" + filename
        SeqIO.write(snp_proteins, filename, "fasta")

### RESIDUE CONSERVATION CALCULATIONS

# Generate list of all genes represented in file
genes = snpres.gene.unique()

#Drop duplicate dbSNP IDs when SNP associated with multiple protein changes
snpres_nr = snpres.drop_duplicates(subset='dbsnp')


# Single letter amino acid data for dataframe
aaslc = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T', '-']

finalsum = pd.DataFrame(columns=['gene','dbSNP','Residue','Isoform 1 AA', 'SNP AA', 'Highest AA:% at Res','2nd Highest AA:% at Res', '3rd Highest AA:% as Res', 'Shannon Entropy'])
# Loop through each gene to make subset dataframe from only that gene's SNPs
for gene in genes:
    a = {'AA': aaslc}
    filename = gene + '.aln'
    filepath = os.path.join('./alignments',filename) #filepath to search for gene alignment file in alignment folder
    alldf = pd.DataFrame(data=a) #initialize alldf with the aaslc column
    shanen = list()
    for index, row in snpres_nr.iterrows(): #loop through rows in snpres_nr df
        snp_dict = {}
        if row['gene'] == gene: #for any row containing the gene name
            snp = []
            pos = int(row['AA_pos']) - 1 #save that residue poistion as pos
            for record in SeqIO.parse(filepath, 'fasta'): #for each sequence in the alignment file 
                if pos < len(record):
                    seq = record.seq[pos] # save the pos residue AA as seq
                    snp.append(seq) # add the seq AA to the snp list
            snp_dict[row['dbsnp']] = snp
            for key, value in snp_dict.items(): #for each item in the snp_dict
                separator = ""
                res_join = separator.join(value)
                shan = eta(res_join, unit='natural')
                shan = round(shan, 4)
                shanen.append(shan)
                count= Counter(res_join) #count the number of times each AA appears
                df = pd.DataFrame.from_dict(count,orient='index') #create a df from the dictionary
                df = df.reset_index()
                df.columns = ['AA', key] # rename the columns
                alldf = alldf.merge(df, how='left', on='AA').replace(np.nan,0) #merge the new df onto the alldf by AA position
    
    ##TO SAVE AND VIEW COUNTS            
    #saveas = gene + ' Conservation Counts.csv'
    #alldf.to_csv(saveas, index=False) # save to .csv file
    
    # generate percentage counts table
    for column in alldf: #for each column in the alldf
        if column != 'AA': #not the AA column
            alldf[column] = alldf[column].transform(lambda x: (x / x.sum())*100) #transform counts into percentages of total
    alldf = alldf.round(2) #round to 2 decimal places
    
    ##TO SAVE AND VIEW PERCENTAGES
    #saveas = gene + ' Conservation Percentages.csv'
    #alldf.to_csv(saveas, index=False) # save to .csv file
    
    # create a summary table from the results
    sumdf = alldf.set_index('AA')
    columns = list(sumdf)
    aadict = {}
    for i in columns:
        newdf = sumdf[i].nlargest(3)
        newdf = newdf.reset_index()
        newdf[i] = newdf[i].astype(str) 
        newdf['aap'] = newdf['AA'].str.cat(newdf[i], sep=" : ")
        topaas = newdf['aap'].to_list()
        aadict[i] = topaas
    
    topaadf = pd.DataFrame.from_dict(aadict, orient='index', columns = ['Highest AA:% at Res', '2nd Highest AA:% at Res', '3rd Highest AA:% as Res'])
    topaadf = topaadf.reset_index()
    topaadf = topaadf.rename(columns={'index':'dbSNP'})
    topaadf['Shannon Entropy'] = shanen
    genedf = snpres[snpres['gene'].str.contains(gene)]
    genedf = genedf.rename(columns={"dbsnp": "dbSNP", "AA_pos":"Residue", "wt_AA": "Isoform 1 AA", "mut_AA": "SNP AA"})
    
    summarydf = genedf.merge(topaadf, how='left', on="dbSNP")
    summarydf.loc[(summarydf['3rd Highest AA:% as Res'] == 'G : 0.0'), '3rd Highest AA:% as Res'] = 'NaN'
    summarydf.loc[(summarydf['3rd Highest AA:% as Res'] == 'P : 0.0'), '3rd Highest AA:% as Res'] = 'NaN'
            
    ##TO SAVE AND VIEW PER GENE SUMMARIES
    #saveas = gene + ' Conservation Summary.csv'
    #summarydf.to_csv(saveas, index=False) # save to .csv file
    
    finalsum = pd.concat([finalsum,summarydf], ignore_index=True)

##TO SAVE SUMMARY FILE TO CSV
finalsum.to_csv('All Genes Summary Table.csv', index=False)


###PyRosetta Score Generation
from rosetta.protocols.relax import FastRelax
sfxn = get_score_function(True)
fr = FastRelax()
fr.set_scorefxn(sfxn)

##Dataframe of conserved, missense mutations (SE <= 1.0 and excluding nonsense * mutations)
prsnps = finalsum[['gene', 'dbSNP', 'Residue', 'Isoform 1 AA', 'SNP AA', 'Shannon Entropy']]
prsnps = prsnps.loc[prsnps['Shannon Entropy'] <= 1]
prsnps = prsnps.loc[prsnps['SNP AA'] != '*']
prsnps = prsnps.drop(columns=['Shannon Entropy'])
prsnps = prsnps.rename(columns={'dbSNP': 'dbsnp', 'Residue': 'res', 'Isoform 1 AA': 'i1aa', 'SNP AA': 'snpaa'})


scoredf = pd.DataFrame(columns = ['dbsnp','score'])
genescores = pd.DataFrame(columns = ['gene', 'score'])

#For each pdb file named by gene
for filename in os.listdir('./pdbs'):
    if filename.endswith(".pdb"):
        genename = filename.rsplit('.',1)[0]
        filepath = os.path.join('./pdbs',filename)
        
        #read in the pdb file
        pose = pose_from_pdb(filepath)
        fr.apply(pose)
        rname = genename + "_relaxed.pdb"
        pose.dump_pdb(rname)
        ascore = sfxn(pose)
        genescores = genescores.append({'gene': genename, 'score': ascore}, ignore_index=True)
        
        for index,row in prsnps.iterrows():
            original_pose = pose.clone()
            if row['gene'] == genename:
                fres = pose.pdb_info().pdb2pose('A', row['res'])
                if fres == 0:
                    print(row['dbsnp'] + ' out of range')
                    continue
                else:
                    tres = pose.residue(fres).name()
                    if tres == 'HIS_D':
                        t1 = 'H'
                    elif tres == 'MET:NtermProteinFull':
                        t1 = 'M'
                    elif tres == 'LYS:NtermProteinFull':
                        t1 = 'L'
                    else:
                        t1 = three_to_one(tres)
                if row['i1aa'] == t1:                
                    mutate_residue(original_pose, fres, row['snpaa'],pack_radius=4.5)
                    ##TO SAVE MUTATED PROTEINS
                    #mname = genename + "_mutated.pdb" 
                    #original_pose.dump_pdb(mname)
                    mscore = sfxn(original_pose)
                    netscore = mscore - ascore
                    scoredf = scoredf.append({'dbsnp': row['dbsnp'], 'score': netscore}, ignore_index=True)
                else:
                    print(row['dbsnp'] + ' no match')

finalscores = prsnps.merge(scoredf, on="dbsnp")
finalscores.to_csv('final_scores.csv', index=False)
genescores.to_csv('gene_scores.csv', index=False)


##Determine pathogenic classification of modeled SNPs
clinsig = brdf[["dbsnp","clinsig"]]
clinsig = clinsig.drop_duplicates(subset=['dbsnp'])
scores = pd.read_csv('final_scores.csv')
clinscores = scores.merge(clinsig, how="left",on="dbsnp")
clinscores.to_csv('clinscores.csv', index=False)
