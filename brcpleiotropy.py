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

#Complete disease data frame from all clinvar snps associated with the listed secondary diseases
diseases = pd.concat(map(pd.read_csv, glob.glob(os.path.join('./diseases', "*.csv"))), ignore_index=True)
diseases.columns = ['name', 'genes', 'proteins', 'conditions', 'clinsig', 'review', 'accession', 'grch37chr', 'grch37loc', 'grch38chr', 'grch38loc','varid', 'alleleid', 'dbsnp', 'spdi']
diseases = diseases.drop_duplicates()

#Dataframe from clinvar snps associated with breast cancer
brdf = pd.read_csv('./diseases/clinvar_breastcancer.csv')
brdf.columns = ['name', 'genes', 'proteins', 'conditions', 'clinsig', 'review', 'accession', 'grch37chr', 'grch37loc', 'grch38chr', 'grch38loc','varid', 'alleleid', 'dbsnp', 'spdi']
brdf = brdf[(brdf['conditions'].str.contains(r'(?=.*breast)', regex=True, case=False, na=False))]
brcgenes = brdf['genes'].to_list()
brcgenes = set(brcgenes)
len(brdf)

#Complete list of genes associated with breast cancer (split on multiple listings)
splitgs = set()
for gene in brcgenes:
    gene = str(gene)
    if "|" in gene:
        gsplit = gene.split("|")
        splitgs.update(gsplit)
    else:
        splitgs.add(gene)
len(splitgs)

#number of genes associated with breast cancer (not split on multiple listings - grouping = 1 gene)
nosplits = set()
for gene in brcgenes:
    gene = str(gene)
    if "|" in gene:
        continue
    else:
        nosplits.add(gene)
len(nosplits)

#Load disease list
disfile = open('./data/finalconditions.txt', "r")
content = disfile.read()
dis_set = content.split("\n")

#Disease set with Breast cancer and other cancer
disfileb = open('./data/finalconditionsbr.txt', "r")
contentb = disfileb.read()
dis_setb = contentb.split("\n")

#Disease set with Breast cancer and other cancer
disfile1 = open('./data/finalconditionswcs.txt', "r")
content1 = disfile1.read()
dis_set1 = content1.split("\n")

#Identify SNPs associated with both breast cancer and disease
dis = '|'.join(dis_set)
pleios = brdf[brdf['conditions'].str.contains(dis, case=False)]
len(pleios)

#Identify genes associted with pleiotropic snps (not split on multiple listings - grouping = 1 gene)
pgenes = pleios['genes'].unique()
nopgs = set()
for gene in pgenes:
    gene = str(gene)
    if "|" in gene:
        continue
    else:
        nopgs.add(gene)
len(nopgs)


##Pleiotropic SNPs split on multiple listings
pgs = set()
for gene in pgenes:
    gene = str(gene)
    if "|" in gene:
        pgss = gene.split("|")
        pgs.update(pgss)
    else:
        pgs.add(gene)
len(pgs)

###list of genes for heatmap = pgs

####GENERATE DENSITY PLOTS FOR ALL SNPS####
#Using diseases dataframe

#Create a new dataframe with only the genes in 'genes' 
genes_re = re.compile("|".join(pgs))
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


for i in pgs:
    bygene = densitydf[densitydf['gene'].str.contains(i)]
    bygene = bygene.dropna()
    countdf = pd.DataFrame(columns = ['location', 'condition'])
    for index,row in bygene.iterrows():
        for dis in dis_setb:
            if dis == 'Other Cancer':
                matches = ["cancer", "Cancer", "carcinoma", "Carcinoma", "tumor", "Tumor", "oma"]
                if any(x in row['conditions'] for x in matches) and not (brc.casefold() or hcp.casefold() in row['conditions']):
                    newdf = pd.DataFrame({"location":[row['location']], "condition":'Other Cancer'})
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
newallgenecounts = pd.DataFrame(pgs, columns=['Gene'])
for i in dis_set1:
    if i == "Breast Cancer":
        genecount = brdf['genes'].value_counts().rename_axis('Gene').reset_index(name='Breast Cancer')
        newallgenecounts = newallgenecounts.merge(genecount, how='left', on=['Gene'])
    elif i == "Other Cancer":
        modcons = list()
        for entry in brdf['conditions']:
            entry = str(entry)
            lst = entry.split("|")
            nob = [x for x in lst if "breast" not in x]
            nob = [x for x in nob if "Breast" not in x]
            nohcp = [x for x in nob if "Hereditary cancer-predisposing" not in x]
            nohcp_join = ''.join(nohcp)
            modcons.append(nohcp_join)
        brdf['modcons'] = modcons
        subset = brdf[(brdf['modcons'].str.contains(r'(?=.*cancer)', regex=True, case=False, na=False)) | (brdf['modcons'].str.contains(r'(?=.*carcinoma)', regex=True, case=False, na=False)) | (brdf['modcons'].str.contains(r'(?=.*tumor)', regex=True, case=False, na=False)) | (brdf['modcons'].str.contains(r'(?=.*oma)', regex=True, case=False, na=False))]
        genecount = subset['genes'].value_counts().rename_axis('Gene').reset_index(name='Other Cancer')
        newallgenecounts = newallgenecounts.merge(genecount, how='left', on=['Gene'])
    else:
        subset = brdf[brdf['conditions'].str.contains(i, case=False, na=False)]
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
agc = sns.heatmap(newallgenecounts, cmap="coolwarm", robust=True, annot=True, fmt="1g")
plt.title('Heatmap of All Disease SNPs on Key Genes')
#agc.figure.savefig('agcall.png', dpi=100) ##To save heatmap file

###List of SNPs and locations for MendelVar
###Mendelvar - coassociated SNPs, location interval
#pleios df = snps assocaited with brc and additional condition

pleiosnps = pleios.drop_duplicates(subset="dbsnp")

##Generate location interval
loc_list = pleiosnps['grch38loc'].tolist()
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

chr_list = pleiosnps['grch38chr'].tolist()
string = "chr"
chrs = [string + x for x in chr_list]
db_list = pleiosnps['dbsnp'].tolist()
snpdf = pd.DataFrame(list(zip(chr_list,loc1,loc2,db_list)), columns = ['chromosome', 'location1', 'location2', 'dbsnp'])
snpdf.to_csv('All SNPs', header=False, index=False, sep=" ")

###PROTEIN SEQUENCE ALIGNMENTS

### GENERATION OF DATA SET FOR SNPS WITH BREAST CANCER + DISEASE : CO-LOCALIZED SNPS
#Data set = pleiosnps

sfp = pleiosnps[['genes','grch38loc','proteins','dbsnp']]
sfp = sfp.replace(['ATM|C11orf65', 'C11orf65|ATM'], 'ATM')
sfp = sfp.replace(['FANCC|AOPEP'], 'FANCC')
sfp = sfp.replace(['FANCD2|LOC107303338'], 'FANCD2')
sfp = sfp.replace(['HRAS|LRRC56'], 'HRAS')
sfp = sfp.replace(['TH2LCRR|RAD50'], 'RAD50')
sfp = sfp.replace(['BRCA2|LOC106721785'], 'BRCA2')
sfp.dropna(subset = ['proteins'], inplace=True)
sfp = sfp[~sfp['proteins'].str.contains('[a-z]')]
sfp = sfp.reset_index(drop=True)
sfp = sfp.sort_values('genes')

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

snp_res = pd.DataFrame(list(zip(sfp_g, sfp_d, AA_pos, AA_old, AA_new)), columns = ['gene', 'dbsnp','AA_pos', 'wt_AA', 'mut_AA'])
snp_res.to_csv('snp_residues.csv', index=False)

for filename in os.listdir('./fastas'):
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
                    print('WT AA does not match isoform entry for SNP:')
                    print(row['gene'],row['dbsnp'])
                    
    saveas = "coassoc_" + filename
    SeqIO.write(snp_proteins, filename, "fasta")

### RESIDUE CONSERVATION CALCULATIONS
# Read in CO-LOCALIZED snps csv file with gene (gene name), dbsnp (ID#), and AA_pos (residue position)
snpres = pd.read_csv('snp_residues.csv')

# Generate list of all genes represented in file
genes = snpres.gene.unique()

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
    for index, row in snpres.iterrows(): #loop through rows in snpres df
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
                alldf = alldf.merge(df, how='left', on='AA') #merge the new df onto the alldf by AA position
                alldf[key] = alldf[key].replace(np.nan,0)
                
    saveas = gene + ' Co-Localized SNP Amino Acid Residue Conservation Counts.csv'
    alldf.to_csv(saveas, index=False) # save to .csv file
    
    # generate percentage counts table
    for column in alldf: #for each column in the alldf
        if column != 'AA': #not the AA column
            alldf[column] = alldf[column].transform(lambda x: (x / x.sum())*100) #transform counts into percentages of total
    alldf = alldf.round(2) #round to 2 decimal places
    saveas = gene + ' Co-Localized SNP Amino Acid Residue Conservation Percentages.csv'
    alldf.to_csv(saveas, index=False) # save to .csv file
    
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
    #genedf = genedf.drop(columns=['gene'])
    genedf = genedf.rename(columns={"dbsnp": "dbSNP", "AA_pos":"Residue", "wt_AA": "Isoform 1 AA", "mut_AA": "SNP AA"})
    
    summarydf = genedf.merge(topaadf, how='left', on="dbSNP")
    summarydf.loc[(summarydf['3rd Highest AA:% as Res'] == 'G : 0.0'), '3rd Highest AA:% as Res'] = 'NaN'
    summarydf.loc[(summarydf['3rd Highest AA:% as Res'] == 'P : 0.0'), '3rd Highest AA:% as Res'] = 'NaN'
            
    
    saveas = gene + ' Co-Localized SNP Residue Conservation Summary.csv'
    summarydf.to_csv(saveas, index=False) # save to .csv file
    
    finalsum = pd.concat([finalsum,summarydf], ignore_index=True)

finalsum.to_csv('All Genes Summary Table.csv', index=False) # complete summary table save to .csv file


###PyRosetta Score Generation
from rosetta.protocols.relax import FastRelax
sfxn = get_score_function(True)
fr = FastRelax()
fr.set_scorefxn(sfxn)
fr.max_iter(100)

#Call in File with conserved snps for mutagenesis
prsnps = pd.read_csv('./data/pyrosettasnps.csv')

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
                    else:
                        t1 = three_to_one(tres)
                if row['i1aa'] == t1:                
                    mutate_residue(original_pose, fres, row['snpaa'])
                    mscore = sfxn(original_pose)
                    netscore = mscore - ascore
                    scoredf = scoredf.append({'dbsnp': row['dbsnp'], 'score': netscore}, ignore_index=True)
                else:
                    print(row['dbsnp'] + ' no match')

finalscores = prsnps.merge(scoredf, on="dbsnp")
finalscores.to_csv('final_scores.csv', index=False)
genescores.to_csv('gene_scores.csv', index=False)
