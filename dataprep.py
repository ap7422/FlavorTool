
import pandas as pd

from rdkit.Chem import AllChem
from rdkit import Chem

import numpy as np
from sklearn import decomposition
from math import sqrt


def search_cluster(smile, cluster_list, smiles_list):
    smile_index = smiles_list.index(smile)
    cluster_id = None
    for count, cluster_Def in enumerate(cluster_list, 0):
        if smile_index in cluster_Def:
            cluster_id = count
    if cluster_id == None:
        print("error")
    return cluster_id


def ClusterFps(fps, cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs

def get_similarities(flv):

    ester_data = pd.read_excel("fol//ester.xlsx")
    name_db = pd.read_csv("fol//out.csv")

    name_dic = {}
    for row in name_db.iterrows():
        row
        name_dic[row[1][1]] = row[1][2]
    name_dic
    ms = []
    smiles = []
    for row in ester_data.iterrows():
        if row[1][4] > -1:
            m = Chem.MolFromSmiles(row[1][1])
            ms.append(m)
            smiles.append(row[1][1])


    data = pd.read_excel("fol//flv_data.xlsx")
    for column in data.columns[1:]:
        if data[column][0] not in smiles and type(data[column][0]) != np.float64 :

            m = Chem.MolFromSmiles(data[column][0])
            ms.append(m)
            smiles.append(data[column][0])

    for sml in flv.keys():
        if sml not in smiles:
            m = Chem.MolFromSmiles(sml)
            ms.append(m)
            smiles.append(sml)


    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in ms]
    w = ClusterFps(fps, cutoff=0.1)
    number_flv = 25
    array = [0]*len(w)
    large_array = []
    for small in range(0, number_flv):
        large_array.append(array.copy())

    flv_comp_dict = {}

    for index, row in data.iterrows():
        flv_comp_dict[row["Unnamed: 0"]] = {}
        if index > 2:
            for column in data.columns[1:]:
                if type(data[column][0]) != np.float64 and row[column] > 0:
                    cluster = search_cluster(data[column][0], w, smiles)
                    if data[column][0] in flv_comp_dict[row["Unnamed: 0"]]:
                        flv_comp_dict[row["Unnamed: 0"]][data[column][0]] += row[column]
                    else:
                        flv_comp_dict[row["Unnamed: 0"]][data[column][0]] = row[column]
                    large_array[index][cluster] += row[column]

    flv_ar = [0]*len(w)
    for sml in flv.keys():
        cluster = search_cluster(sml, w, smiles)
        flv_ar[cluster] = flv[sml]

    large_array_2 = large_array[2:]
    large_array_2[0] = flv_ar

    df = pd.DataFrame.from_records(large_array_2)

    useless = []
    for column in df.columns:
        if df[column].sum() < .01:
            useless.append(column)


    pca = decomposition.PCA()
    pca.fit(df)

    pca.n_components = 15
    x_reduced = pca.fit_transform(df)

    name_list = data["Unnamed: 0"]
    name_list_trimmed = name_list[2:]
    score_array = []
    other = 2
    for row in df.iterrows():
        n = 1
        score = 0
        sum_of_other = 0
        for count, col in enumerate(df.columns, 0):
            score += abs(df.iloc[[0]][count][0] - row[1][col])
            sum_of_other += row[1][col]
            if df.iloc[[0]][count][0] > 0 or row[1][col] > 0:
                n += 1

        score_array.append(score / (sum_of_other))
        other += 1

    name_dict = {}
    for count, score in enumerate(score_array, 2):
        name_dict[name_list_trimmed[count]] = score

    sorter_ver = sorted(name_dict, key=name_dict.get, reverse=False)
    #sorter_ver.remove(" flv green tea cy07066")
    val_list = [[], [], []]
    name_list = [[], [], []]
    #TODO add check name
    for count, flv_name in enumerate(sorter_ver[1:4], 0):
        for smile in flv_comp_dict[flv_name].keys():

            val_list[count].append(flv_comp_dict[flv_name][smile])
            name_list[count].append(name_dic[smile])
    print(sorter_ver[0])
    print(sorter_ver[1])
    print(sorter_ver[2])
    print(name_dict[sorter_ver[0]])
    print(name_dict[sorter_ver[1]])
    print(name_dict[sorter_ver[2]])

    return val_list, name_list, sorter_ver[2], name_dict[sorter_ver[2]]


