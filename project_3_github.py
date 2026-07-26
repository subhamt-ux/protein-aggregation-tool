import numpy as np
import pandas as pd

amino_acids = ['I','V','L','F','C','M','A','G','T','S','W','Y','P','H','E','Q','D','N','K','R']
hydrophobicity_values = [4.5,4.2,3.8,2.8,2.5,1.9,1.8,-0.4,-0.7,-0.8,-0.9,-1.3,-1.6,-3.2,-3.5,-3.5,-3.5,-3.5,-3.9,-4.5]
#Gives the hydrophobicity values according to the Kyle-Dolittle scale for each amino acid (each index gives the amino acid's hydrophobicity value at the 
#same index of amino_acids

def aggregation_scores(filename_tsv:str, window_length:int):
    '''
    Given an organism's proteome (the file), the top 1% of aggregate-prone proteins are given. 
    These particular proteins are found by summing two calculated values for each protein:
        1. Overall Hydrophobicity Value: Sums up every amino acid's hydrophobicity value in a protein then divides by the number of proteins, finding the average.
        2. Highest Aggregation Window Value: Gets the average hydrophobicity value for every n amino acids (variable is window_length)
    The sums, also the Composite Score for each protein, are then normalized by being converted
    into z-scores, and the top 1% sums are found, which is the top 1% of aggregate-prone proteins in an organism.
    '''
    df = pd.read_csv(filename_tsv,sep = "\t")
    df["Sequence"] = df['Sequence'].fillna("")
    df["Sequence"] = df["Sequence"].str.strip()
    df["string_lengths"] = df['Sequence'].str.len()
    max = df["string_lengths"].max()
    df["Sequence List"] = df["Sequence"].str.ljust(max," ")
    df["Sequence List"] = df["Sequence List"].apply(list)
    sequence_array = np.array(df['Sequence List'].tolist())

    letters = np.array(amino_acids)
    hydrophobic_scores = np.array(hydrophobicity_values)
    sequence_array = sequence_array[:,np.newaxis,:]
    letters = letters[np.newaxis,:,np.newaxis]
    hydrophobic_scores = hydrophobic_scores[np.newaxis,:,np.newaxis]

    h_scores_per_sequence = np.where(sequence_array == letters,hydrophobic_scores,0)
    h_score_per_sequence = np.sum(h_scores_per_sequence,-1)
    h_score_per_sequence = np.sum(h_score_per_sequence,-1)

    df['Sum'] = h_score_per_sequence
    df["Avg Hydrophobic Value"] = df['Sum'] / df["string_lengths"]

    df["Avg Hydrophobic Value"] = df['Sum'] / df["string_lengths"]
    #The Overall Hydrophobicity Value for each protein (1.)

    h_scores_per_sequence = np.sum(h_scores_per_sequence,1)

    index_0 = np.arange(df["Entry"].size)
    index_1_1 = np.arange(max-window_length+1)
    index_1_2 = np.arange(window_length)
    index_1 = index_1_1[:,np.newaxis] + index_1_2[np.newaxis,:]

    windows_per_sequence = h_scores_per_sequence[index_0[:,np.newaxis,np.newaxis],index_1[np.newaxis,:,:]]
    window_sums_per_sequence = np.sum(windows_per_sequence,-1)
    max_per_sequence = np.max(window_sums_per_sequence,-1)
    df['Max Agg.'] = max_per_sequence
    #Highest Aggregation Window Value for each protein (2.)

    df["Composite Score"] = df["Avg Hydrophobic Value"] + df['Max Agg.']
    #The composite score for each protein

    max_score = df["Composite Score"].max()
    df["Composite Score"] = (df["Composite Score"] - df["Composite Score"].mean())/df["Composite Score"].std()
    cutoff = df["Composite Score"].quantile(0.99)
    most_aggregation_prone = df.loc[df["Composite Score"]>=cutoff]
    return most_aggregation_prone

if __name__ == "__main__":
    df = aggregation_scores(filename_tsv = "small_enough.tsv", window_length = 7)
    print("Top 1% Most Aggregation Prone Proteins")
    print(df.loc[:,["Protein names","Avg Hydrophobic Value",'Max Agg.',"Composite Score"]])


