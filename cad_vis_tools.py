
from bokeh.models import ColumnDataSource, TextInput, Button, Paragraph, Select
from bokeh.palettes import Category20c
from bokeh.plotting import figure
from bokeh.transform import cumsum
import math
import pandas as pd
from math import pi
import requests
from xml.etree import ElementTree
from rules import rulelist, namelist, smilies_dict_select, form_dict_select, oil_dict_select

from picture import saveimg
import string
import random


# generic object that holds data about a single molecule once retrieved from the cadmode site
class MolObj:
    def __init__(self, name, smile, val_dict, percent):
        self.name = name
        self.smile = smile
        self.prop = val_dict
        self.prc = percent
        self.mprc = None
        self.rule ={}
        self.prob_lvl = None
   # def __cmp__(self, other):
      #  return cmp(self.prob_lvl, other.prob_lvl)

    def __repr__(self):
        return "Class that hold a molecule's data"


# inputs the dictonary and retusn colors
# the random bit is due to the fact that the scheme is only 20 colors
def get_color(x):
    if len(x) == 1:
        return ['#3182bd']
    elif len(x) == 2:
        return ['#3182bd', '#6baed6']
    elif len(x) > 20:
        final_list = []
        cololor_list = Category20c[20]

        for thing in range(0, len(x)):
            random_number =random.randint(0, 19)
            random_color = cololor_list[random_number]
            final_list.append(random_color)
        return final_list
    else:
        return Category20c[len(x)]

# just a simple function that parses the result and color vector
# just to remove redundant code that would have to exist at the end of every rule
def get_color_v(result, vector):
    if result > vector[1]:
        color = 2
    elif vector[0] < result < vector[1]:
        color = 1
    else:
        color = 0
    return color


# for flavors that contain an oil this tool transforms it into its components and adds it to the main dictionary
# this should only be used right
def smiles_parser(x):
    for flv_comp in list(x):
        if flv_comp in oil_dict_select.keys():

            for smile in oil_dict_select[flv_comp].keys():
                if smile in list(x):
                    x[smile] += x[flv_comp]*oil_dict_select[flv_comp][smile]/100
                else:
                    x[smile] = x[flv_comp]*oil_dict_select[flv_comp][smile]/100
            x.pop(flv_comp)
    return x


# returns the pie chat seen
# input is a dictionary that has the smile/oil name as the key and percentage in the flavor
def get_plot(x):
    data = pd.Series(x).reset_index(name='value').rename(columns={'index': 'country'})
    data['angle'] = data['value'] / 100 * 2 * pi
    data['color'] = get_color(x)
    source = ColumnDataSource(data)

    p = figure(plot_height=350, title="Pie Chart", toolbar_location=None,
               tools="hover", tooltips="@country: @value", x_range=(-0.5, 1.0))

    p.wedge(x=0, y=1, radius=0.4,
            start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
            line_color="white", fill_color='color',  source=source)

    p.axis.axis_label = None
    p.axis.visible = False
    p.grid.grid_line_color = None
    return p


# inputs smiles dictionary
# list of ids you want calculated
# dict of the ids and name
# and the framework you want to calculate props from
# output is the mlc list
def get_value(smiles, ids, iddict, framework):
    percent = list(smiles.values())
    smilestrg = dict_to_str(smiles)

    files = {
        'user': (None, 'stein.bm'),
        'smiles': ('input.smi', smilestrg.encode('ascii')),
        'prop_ids': (None, ids),
    }
    response = requests.post('http://cadmol.na.pg.com/molprop/calc.cgi', files=files, auth=('testing', 'system'))
    tree = ElementTree.fromstring(response.text)
    mlclist = []
    count = 0
    for molecule in tree.findall('molecule'):
        mol_name = molecule.find('name').text
        mol_smile = molecule.find('smiles').text
        prop_dict = {}
        for property_blob in molecule.find('calc_values').findall('property'):
            mol_id = property_blob.find('property_id').text
            value = property_blob.find('value').text
            prop_dict[iddict[mol_id]] = int(float(value))

        mol = MolObj(mol_name, mol_smile, prop_dict, percent[count])
        count += 1
        mlclist.append(mol)

    # calculating mol percent
    sumofmass = 0
    for mol in mlclist:
        sumofmass += mol.prc/mol.prop["General Molecular Physicochemical Properties (molecular weight)"]
    for mol in mlclist:
        mol.mprc = (mol.prc/mol.prop["General Molecular Physicochemical Properties (molecular weight)"])/sumofmass
    # adding the rules value for the single mol
    for mol in mlclist:
        for rule in rulelist:
            rul_val = rule([mol], framework)
            mol.rule[rul_val[0]] = rul_val
    # calculating problem level
    for mol in mlclist:
        total = 0
        for value in mol.rule.values():
            total += value[1]*mol.mprc

        mol.prob_lvl = total
    for mol in mlclist:
        if mol.name == None:
            mol.name = " Unknown"


    mlclist.sort(key=lambda x: x.prob_lvl, reverse=True)

    return mlclist


# inputs the list of mlc returned
# optional is a name of a mlc a user wants to see
# outputs are the two lists used to make the two tables seen below the pictures
def get_table_values(mlclist, usr_pref=None):
    # generating sum list
    sum_val = []
    sum_color = []
    name = []
    name_link = []
    if usr_pref is None:
        table = mlclist[:3]
    else:
        table = usr_pref
    for rule_name in mlclist[0].rule.keys():
        total = 0
        for mol in mlclist:

            total += mol.rule[rule_name][1]
        sum_val.append("{:05.2f}".format(total))
        color = get_color_v(total, mlclist[0].rule[rule_name][2])
        sum_color.append(color)
        name.append(rule_name)
        name_link.append(mlclist[0].rule[rule_name][3])


    other = []
    other_color = []
    for count, rule_name in enumerate(mlclist[0].rule.keys(), 0):
        total = 0
        for mol in table:
            total += mol.rule[rule_name][1]
        other_val = float(sum_val[count]) - total
        other.append("{:05.2f}".format(other_val))
        other_color.append(get_color_v(other_val, mlclist[0].rule[rule_name][2]))
    names = []
    total_list = [[], [], []]
    color_list = [[], [], []]
    for c, mol in enumerate(table, 0):
        names.append(mol.name)
        for rule_name in mlclist[0].rule.keys():
            value = mol.rule[rule_name][1]
            total_list[c].append("{:05.2f}".format(value))
            color_list[c].append(get_color_v(value, mlclist[0].rule[rule_name][2]))

    other_half = [other, other_color, sum_val, sum_color, name, name_link]
    final = total_list+color_list + other_half
    return final, names


# input mlclist and possibly a list of names preferred by the user
def get_pictures(mlclist, usr_pref=None):
    plot_list = []
    letters = string.digits

    if usr_pref is None:
        table = mlclist[:3]
    else:
        table = usr_pref
    for count, mol in enumerate(table, 0):
        plot_list.append(saveimg(mol.smile, ''.join(random.choice(letters) for i in range(10))))

    return plot_list

# input is the list of the names of the ids you require in the models
# output is the id_list and the id_dict required by get_value
def get_ids(name_list):
    response = requests.get('http://cadmol.na.pg.com/molprop/property_map.cgi', auth=('testing', 'system'))
    tree = ElementTree.fromstring(response.text)
    id_list = []
    id_dict = {}
    for property_blob in tree.findall('property'):
        id_val = property_blob.find('id').text
        name = property_blob.find('name').text
        state = property_blob.find('cba').text
        if name in name_list and state == 'Y':
            id_list.append(id_val)
            id_dict[id_val] = name
    return id_list, id_dict



def dict_to_str(moldict):
    string = ""
    for key in moldict.keys():
        string += key+"\n"
    return string


def list_to_str(idlist):
    string = ""
    for iteam in idlist:
        string += iteam+","
    return string[:len(string)-1]


