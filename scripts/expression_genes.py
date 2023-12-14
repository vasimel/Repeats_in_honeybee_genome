import os
import pandas as pd
caste = "DRONE"
df = pd.read_csv("GCF_003254395.2_Amel_HAv3.1_genomic.gtf", header=None, comment='#', sep='\t',
                 usecols=[0, 2, 3, 4, 6, 8])
df.rename({0: "chr", 2: "type", 3: "start", 4: "stop", 6: "dir", 8: "attrib"}, axis=1, inplace=True)
df = df.query("type=='gene' or type=='exon' or type=='transcript'")
df['gene_id'] = [j.split(" ")[-1].strip('"') for i in df['attrib'].str.split(';').to_list() for j in i if
                 "gene_id" in j]
df['transcript_id'] = [j.split(" ")[-1].strip('"') for i in df['attrib'].str.split(';').to_list() for j in i if
                       "transcript_id" in j]
df.drop(columns=['attrib'], inplace=True)
dx = pd.read_csv(f"../{caste}/expressed_genes.calc.bed", sep='\t', names=["chr", "start", "stop", "median"], header=None)
merged = pd.merge(left=df, right=dx, how='left', on=["chr", "start", "stop"])
x = merged.groupby(['gene_id'])
print(dx)
print(merged[:20])

if os.path.exists(f"../{caste}/features/"):
    if os.path.exists(f"../{caste}/features/expressed_genes.bee.bed"):
        os.remove(f"../{caste}/features/expressed_genes.bee.bed")

    if os.path.exists(f"../{caste}/features/expressed_tss.bee.bed"):
        os.remove(f"../{caste}/features/expressed_tss.bee.bed")

    if os.path.exists(f"../{caste}/features/expressed_tes.bee.bed"):
        os.remove(f"../{caste}/features/expressed_tes.bee.bed")

    if os.path.exists(f"../{caste}/features/expressed_introns.bee.bed"):
        os.remove(f"../{caste}/features/expressed_introns.bee.bed")

    if os.path.exists(f"../{caste}/features/expressed_exons.bee.bed"):
        os.remove(f"../{caste}/features/expressed_exons.bee.bed")


    if os.path.exists(f"../{caste}/features/repressed_genes.bee.bed"):
        os.remove(f"../{caste}/features/repressed_genes.bee.bed")

    if os.path.exists(f"../{caste}/features/repressed_tss.bee.bed"):
        os.remove(f"../{caste}/features/repressed_tss.bee.bed")

    if os.path.exists(f"../{caste}/features/repressed_tes.bee.bed"):
        os.remove(f"../{caste}/features/repressed_tes.bee.bed")

    if os.path.exists(f"../{caste}/features/repressed_introns.bee.bed"):
        os.remove(f"../{caste}/features/repressed_introns.bee.bed")

    if os.path.exists(f"../{caste}/features/repressed_exons.bee.bed"):
        os.remove(f"../{caste}/features/repressed_exons.bee.bed")
else:
    os.mkdir(f"../{caste}/features/")

for gene_name, gene_data in x:
    dir = gene_data['dir'].unique().item()
    is_expressed = gene_data[gene_data['type'] == 'gene']['median'].item() >= 1.0
    gene = gene_data[gene_data['type'] == 'gene'].values[0]
    transcript = gene_data.groupby(['transcript_id'])
    for transcript_id, transcript_data in transcript:
        if len(transcript_id) > 0:
            exons_ranges = []
            intron_ranges = []
            for k in transcript_data.to_numpy()[1:]:
                exons_ranges.append([k[0], k[2], k[3]])
            if len(exons_ranges) > 1:
                intron_ranges = [[exons_ranges[0][0], exons_ranges[j][2], exons_ranges[j+1][1]] for j in range(len(exons_ranges)-1)]

            with open(f"../{caste}/features/expressed_exons.bee.bed" if is_expressed else f"../{caste}/features/repressed_exons.bee.bed", 'a') as file:
                file.write('\n'.join(['\t'.join(list(map(str, i))) for i in exons_ranges])+'\n')

            if len(intron_ranges) > 0:
                with open(f"../{caste}/features/expressed_introns.bee.bed" if is_expressed else f"../{caste}/features/repressed_introns.bee.bed", 'a') as file:
                    file.write('\n'.join(['\t'.join(list(map(str,i))) for i in intron_ranges])+'\n')

    with open(f"../{caste}/features/expressed_genes.bee.bed" if is_expressed else f"../{caste}/features/repressed_genes.bee.bed", 'a') as file:
        file.write('\t'.join([gene[0], str(gene[2]), str(gene[3])])+'\n')
    if dir == '+':
        tss = [gene[0], str(gene[2]-1000), str(gene[2]+1000)]
        tes = [gene[0], str(gene[3]-1000), str(gene[3]+1000)]
    else:
        tss = [gene[0], str(gene[3]-1000), str(gene[3]+1000)]
        tes = [gene[0], str(gene[2]-1000), str(gene[2]+1000)]

    with open(f"../{caste}/features/expressed_tss.bee.bed" if is_expressed else f"../{caste}/features/repressed_tss.bee.bed", 'a') as file:
        file.write('\t'.join(tss)+'\n')
    with open(f"../{caste}/features/expressed_tes.bee.bed" if is_expressed else f"../{caste}/features/repressed_tes.bee.bed", 'a') as file:
        file.write('\t'.join(tes)+'\n')