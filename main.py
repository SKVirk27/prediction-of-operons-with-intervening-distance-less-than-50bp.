import os
import pandas as pd
import gffpandas.gffpandas as gffpd
os.chdir(os.path.dirname(__file__))

# function declaration for calculating intergenic distance, distance< 50, strand same and on same contig.
def Operon(s):
    if (s['intergenic_dist'] < 50) and (s['Strand1'] == s['Strand2']):
          return 'OP'

    else:
          return 'NOP'

def MetaOperon(M):
    if (M['intergenic_dist'] < 50) and (M['Strand1'] == M['Strand2']) and (M['contig1'] == M['contig2']):
          return 'OP'

    else:
          return 'NOP'

#******************************************************************************************************
#Creating database for B_subtilis from ptt file
Data_B_subtilis = pd.read_csv('B_subtilis_168.ptt.gz',skiprows=2, delimiter="\t",lineterminator='\n')
Loc_split= Data_B_subtilis['Location']
Data_B_subtilis[['loc_start','loc_end']] =Loc_split.str.split("[..]", n=1, expand=True) # spiltting the location column into start and end
Loc_end=Data_B_subtilis['loc_end']
Data_B_subtilis['loc_end']=Loc_end.str.replace("[.]","",regex=True)
Data_B_subtilis = Data_B_subtilis.iloc[: , 1:]
Data_B_subtilis=Data_B_subtilis.iloc[:, [8,9,0,1,2,3,4,5,6,7]] # rearranging the columns.
Data_B_subtilis[['loc_start','loc_end']]= Data_B_subtilis[['loc_start','loc_end']].apply(pd.to_numeric)# String to numeric
#Data_B_subtilis.to_csv("B_subtilis_168.csv") # file saved as csv

# Finding Intergenic distance between Gene pairs for operon prediction.
a=Data_B_subtilis['loc_start'].values # assigning the variable for start location
b=Data_B_subtilis['loc_end'].values # assigning the variable for end location
start= a[1:].copy() # created a copy for start and end.
end=b[:-1].copy()
intergenic_dist=start-(end+1)  # Calculating intergenic distance between two genes.
Strand=Data_B_subtilis['Strand'].iloc[:-1]
B_subtilis_Operon = pd.DataFrame(intergenic_dist,columns = ['intergenic_dist'])
# created the database for operon prediction.
B_subtilis_Operon['Gene1']=Data_B_subtilis['Gene'].iloc[:-1] # copied the genes from original dataframe and droping end value.
B_subtilis_Operon['Strand1']=Data_B_subtilis['Strand'].iloc[:-1] # same for column strand
B_subtilis_Operon['Gene2'] = Data_B_subtilis['Gene'].shift(-1)# copied the genes and shifted them by one.
B_subtilis_Operon['Strand2'] = Data_B_subtilis['Strand'].shift(-1)# same for column strand
B_subtilis_Operon['OP/NOP'] = B_subtilis_Operon.apply(Operon, axis=1) # applied the operon defined function.
B_subtilis_Operon.to_csv("B_subtilis_OperonPrediction.csv")# Operon prediction database created.
#***********************************************************************************************************************

#Creating database for E_coli from ptt file
Data_E_coli = pd.read_csv('E_coli_K12_MG1655.ptt.gz',skiprows=2, delimiter="\t",lineterminator='\n')
Loc_split= Data_E_coli['Location']
Data_E_coli[['loc_start','loc_end']] =Loc_split.str.split("[..]", n=1, expand=True)# spiltting the location column into start and end
Loc_end=Data_E_coli['loc_end']
Data_E_coli['loc_end']=Loc_end.str.replace("[.]","",regex=True)
Data_E_coli = Data_E_coli.iloc[: , 1:]
Data_E_coli=Data_E_coli.iloc[:, [8,9,0,1,2,3,4,5,6,7]]# rearranging the columns.
Data_E_coli[['loc_start','loc_end']]= Data_E_coli[['loc_start','loc_end']].apply(pd.to_numeric)# String to numeric
#Data_E_coli.to_csv("E_coli_K12_MG1655.csv") # file saved as csv

# Finding Intergenic distance between Gene pairs of E_coli for operon prediction.
c=Data_E_coli['loc_start'].values
d=Data_E_coli['loc_end'].values
start1= c[1:].copy()
end1=d[:-1].copy()
intergenic_dist1=start1-(end1+1)  # Calculating intergenic distance between two genes.
E_coli_Operon = pd.DataFrame(intergenic_dist1,columns = ['intergenic_dist'])
E_coli_Operon['Gene1']=Data_B_subtilis['Gene'].iloc[:-1]# copied the genes from original dataframe and droping end value.
E_coli_Operon['Strand1']=Data_E_coli['Strand'].iloc[:-1]# same for column strand
E_coli_Operon['Gene2'] = Data_E_coli['Gene'].shift(-1) #copied the genes and shifted them by one.
E_coli_Operon['Strand2'] = Data_E_coli['Strand'].shift(-1)
E_coli_Operon['OP/NOP'] = E_coli_Operon.apply(Operon, axis=1)
E_coli_Operon.to_csv("E_coli_OperonPrediction.csv")

#***********************************************************************************************************************
#Creating database for Halobacterium from ptt file
Data_Halob_NRC1 = pd.read_csv('Halobacterium_NRC1.ptt.gz',skiprows=2, delimiter="\t",lineterminator='\n')
Loc_split= Data_Halob_NRC1['Location']
Data_Halob_NRC1[['loc_start','loc_end']] =Loc_split.str.split("[..]", n=1, expand=True)
Loc_end=Data_Halob_NRC1['loc_end']
Data_Halob_NRC1['loc_end']=Loc_end.str.replace("[.]","",regex=True)
Data_Halob_NRC1 = Data_Halob_NRC1.iloc[: , 1:]
Data_Halob_NRC1=Data_Halob_NRC1.iloc[:, [8,9,0,1,2,3,4,5,6,7]]
Data_Halob_NRC1[['loc_start','loc_end']]= Data_Halob_NRC1[['loc_start','loc_end']].apply(pd.to_numeric)# String to numeric
#Data_Halob_NRC1.to_csv("Halobacterium_NRC1.csv") # file saved as csv

# Finding Intergenic distance between Gene pairs of Halobacterium_NRC1 for operon prediction.
e=Data_Halob_NRC1['loc_start'].values
f=Data_Halob_NRC1['loc_end'].values
start2= e[1:].copy()
end2=f[:-1].copy()
intergenic_dist2=start2-(end2+1)  # Calculating intergenic distance between two genes.
Halob_NRC1_Operon= pd.DataFrame(intergenic_dist2,columns = ['intergenic_dist'])
Halob_NRC1_Operon['Gene1']=Data_Halob_NRC1['Gene'].iloc[:-1]# copied the genes from original dataframe and droping end value.
Halob_NRC1_Operon['Strand1']=Data_Halob_NRC1['Strand'].iloc[:-1]# same for column strand
Halob_NRC1_Operon['Gene2'] = Data_Halob_NRC1['Gene'].shift(-1)# copied the genes and shifted them by one.
Halob_NRC1_Operon['Strand2'] = Data_Halob_NRC1['Strand'].shift(-1)
Halob_NRC1_Operon['OP/NOP'] = Halob_NRC1_Operon.apply(Operon, axis=1)
Halob_NRC1_Operon.to_csv("Halob_NRC1_OperonPrediction.csv")
#************************************************************************************************************************

#Creating database for Synechocystis from ptt file
Data_Synec_PCC = pd.read_csv('Synechocystis_PCC6803_uid159873.ptt.gz',skiprows=2, delimiter="\t",lineterminator='\n')
Loc_split= Data_Synec_PCC['Location']
Data_Synec_PCC[['loc_start','loc_end']] =Loc_split.str.split("[..]", n=1, expand=True)
Loc_end=Data_Synec_PCC['loc_end']
Data_Synec_PCC['loc_end']=Loc_end.str.replace("[.]","",regex=True)
Data_Synec_PCC = Data_Synec_PCC.iloc[: , 1:]
Data_Synec_PCC=Data_Synec_PCC.iloc[:, [8,9,0,1,2,3,4,5,6,7]]
Data_Synec_PCC[['loc_start','loc_end']]= Data_Synec_PCC[['loc_start','loc_end']].apply(pd.to_numeric)# String to numeric
#Data_Synec_PCC.to_csv("Synechocystis_PCC6803.csv") # file saved as csv

# Finding Intergenic distance between Gene pairs of Synechocystis_PCC6803 for operon prediction.
g=Data_Synec_PCC['loc_start'].values
h=Data_Synec_PCC['loc_end'].values
start3= g[1:].copy()
end3=h[:-1].copy()
intergenic_dist3=start3-(end3+1)  # Calculating intergenic distance between two genes.
Synec_PCC_Operon= pd.DataFrame(intergenic_dist3,columns = ['intergenic_dist'])
Synec_PCC_Operon['Gene1']=Data_Synec_PCC['Gene'].iloc[:-1]# copied the genes from original dataframe and droping end value.
Synec_PCC_Operon['Strand1']=Data_Synec_PCC['Strand'].iloc[:-1]
Synec_PCC_Operon['Gene2'] = Data_Synec_PCC['Gene'].shift(-1)# copied the genes and shifted them by one.
Synec_PCC_Operon['Strand2'] = Data_Synec_PCC['Strand'].shift(-1)
Synec_PCC_Operon['OP/NOP'] = Synec_PCC_Operon.apply(Operon, axis=1)
Synec_PCC_Operon.to_csv("Synec_PCC_OperonPrediction.csv")

#***********************************************************************************************************************
annotation = gffpd.read_gff3('2088090036 (1).gff',)
attr_to_columns = annotation.attributes_to_columns()
annotation.to_csv('temp.csv')
Data_Metag=pd.read_csv('temp.csv')
Data_Metag[['GeneID', 'Locus', 'Product']] = Data_Metag['attributes'].str.split(';', expand=True)
Data_Metag['GeneID'] = Data_Metag['GeneID'].map(lambda x: str(x)[2:])
Data_Metag=Data_Metag.drop(columns=['attributes'],axis=1)
#Data_Metag.to_csv('Data_Metag.csv')

# Finding Intergenic distance between Gene pairs of Synechocystis_PCC6803 for operon prediction.
I=Data_Metag['start'].values
H=Data_Metag['end'].values
start4= I[1:].copy()
end4=H[:-1].copy()
intergenic_dist4=start4-(end4+1)  # Calculating intergenic distance between two genes.
Metagenome_Operon= pd.DataFrame(intergenic_dist4,columns = ['intergenic_dist'])
Metagenome_Operon['contig1'] = Data_Metag['seq_id'].iloc[:-1]
Metagenome_Operon['Gene1']=Data_Metag['GeneID'].iloc[:-1]# copied the genes from original dataframe and droping end value.
Metagenome_Operon['Strand1']=Data_Metag['strand'].iloc[:-1]#same as above step
Metagenome_Operon['Gene2'] = Data_Metag['GeneID'].shift(-1) # copied the genes and shifted them by one.
Metagenome_Operon['Strand2'] = Data_Metag['strand'].shift(-1)
Metagenome_Operon['contig2'] = Data_Metag['seq_id'].shift(-1)
Metagenome_Operon['OP/NOP'] = Metagenome_Operon.apply(MetaOperon, axis=1)
Metagenome_Operon.to_csv("Metagenome_OperonPrediction.csv")

