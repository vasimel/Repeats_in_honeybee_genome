import pandas as pd

path_to_gff = r"GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
gff = pd.read_csv(
    path_to_gff,
    sep='\t',
    comment='#',
    names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrib'],
    usecols=[2,8]
)
genes = gff.query('type=="gene"')
print(len([j.split(" ")[-1].strip('"') for i in genes['attrib'].str.split(';').to_list() for j in i if "gene_id" in j]))
genes['gene_name'] = [j.split(" ")[-1].strip('"') for i in genes['attrib'].str.split(';').to_list() for j in i if "gene_id" in j]


print(len([j.split(" ")[-1].strip('"') for i in genes['attrib'].str.split(';').to_list() for j in i if 'db_xref "GeneID' in j]))
genes['gene_id'] = [j.split(" ")[-1].strip('"').replace(':','_').upper() for i in genes['attrib'].str.split(';').to_list() for j in i if 'db_xref "GeneID' in j]
print(genes)
genes[['gene_name', 'gene_id']].to_csv("name_id.csv",sep='\t',header=None,index=None)