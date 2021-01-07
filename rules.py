# TODO remove iteration
# clealy just adding the rules
import numpy as np
import pandas as pd
def percent_ester(mols, framework):
    total = 0
    for mol in mols:
        if mol.prop['CountOfEsters'] > 0:
            total += mol.prop['CountOfEsters']*mol.mprc

    return "Percent Ester", total, [.2, .5], 0


def ald_count(mols, framework):
    total = 0
    for mol in mols:
        if mol.prop['CountOfAldehydes'] > 0:
            total += mol.prop['CountOfAldehydes']*mol.mprc

    return "Mol Percent Aldehydes", total, [.5, .20], 2


def amine_stability(mols, framework):
    total = 0
    for mol in mols:
        if mol.prop["Perfume Color Stability in Liquid Formulations"] > 0:
            total += mol.prop["Perfume Color Stability in Liquid Formulations"]

    return "Amine Color Stability", total, [20, 65], "http://hebe.na.pg.com/molprop/docs/PRM_Color_Stability_EXP-19-BC4529-1.pdf"


def mol_percent(mols, framework):
    total = 0

    for mol in mols:
        total += mol.mprc*100

    return "Mol Percent", total, [101, 101], 6


def weight_percent(mols, framework):
    total = 0

    for mol in mols:
        total += mol.prc

    return "Weight Percent", total, [101, 101], 7


# just a place to store a bunch of random data

rulelist = [weight_percent, mol_percent, percent_ester, ald_count, amine_stability]
namelist = ['Normal Boiling Point (ACD, deg C) ()',
            'LDL NP Liq-Air Partition Coefficient log(1/K)',
            "CountOfEsters",
            "Perfume Color Stability in Liquid Formulations",
            "CountOfAldehydes",
            "General Molecular Physicochemical Properties (molecular weight)",
            "Fragrance Intensity Response Factor in Hair Dyes",
            ""
            ]
smilies_dict_select = {"Vanillin": "c1(C=O)cc(OC)c(O)cc1", "Methanol": "CO",
                       "Carvone": "O=C1C[C@@H](C\C=C1\C)C(C)=C", "Triacetin": "CC(=O)OC(COC(=O)C)COC(C)=O"}
form_dict_select = {"Callen's Flavor": {"c1(C=O)cc(OC)c(O)cc1": 50, "CO": 25, "CC(=O)OC(COC(=O)C)COC(C)=O": 20, "CCNC(=O)C1CC(CCC1C(C)C)C": 5}}
oil_dict_select = form_dict_select


data = pd.read_excel("fol//flv_data.xlsx")

oil_dict_select = {"Callen's Flavor": {"c1(C=O)cc(OC)c(O)cc1": 50, "CO": 25, "CC(=O)OC(COC(=O)C)COC(C)=O": 20, "CCNC(=O)C1CC(CCC1C(C)C)C": 5}}
flv_comp_dict = {}
for index, row in data.iterrows():
    flv_comp_dict[row["Unnamed: 0"]] = {}
    if index > 2:
        for column in data.columns[1:]:
            if type(data[column][0]) != np.float64 and row[column] > 0:

                if data[column][0] in flv_comp_dict[row["Unnamed: 0"]]:
                    flv_comp_dict[row["Unnamed: 0"]][data[column][0]] += row[column]
                else:
                    flv_comp_dict[row["Unnamed: 0"]][data[column][0]] = row[column]

flv_comp_dict.pop(list(flv_comp_dict)[0])
flv_comp_dict.pop(list(flv_comp_dict)[0])
form_dict_select = {}
for thing in list(flv_comp_dict):
    form_dict_select[str(thing)] = flv_comp_dict[thing]

oil_dict_select = form_dict_select
