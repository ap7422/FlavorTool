from bokeh.io import curdoc
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, TextInput, Button, Paragraph, Select

from bokeh.events import ButtonClick
from bokeh.models.widgets import DataTable, TableColumn, HTMLTemplateFormatter, Slider, Select
from rules import rulelist, namelist, smilies_dict_select, form_dict_select, oil_dict_select
from bokeh.layouts import layout

from cad_vis_tools import MolObj, get_color, smiles_parser, get_plot, get_value, get_table_values, get_pictures, get_ids, dict_to_str, list_to_str, get_color_v
from dataprep import get_similarities


# updates plot and appends to bokeh
def update_title(event):
    if text.value != "" or text.value != 'enter smile here':
        x[text.value] = conc.value
        x['Empty'] = 0
        sumval = sum(dict.values(x))
        x["Empty"] = 100-sumval
        if x["Empty"] < 5:
            x.pop('Empty')
        text.value = ""
        p = get_plot(x)
        col2.children.remove(col2.children[0])
        col2.children.append(p)


# pretty important function that parses
def run_rules(event):
    from bokeh.layouts import layout
    global mol_list, slider, data, names, x, sim_tab, para
    if "Empty" in x:
        x.pop('Empty')
    x = smiles_parser(x)
    mol_list = get_props(x)
    data, names = get_table_values(mol_list)
    table = get_table(data, names)
    plot = get_pictures(mol_list)
    slider = get_slider(mol_list)
    name_table = get_name_table(mol_list)
    value, name, aflv_name, error = get_similarities(x)
    comp_table = get_sim_table(name[0], value[0])
    sim_tab = comp_table
    para = paragraph(aflv_name, error)
    cur_layout = layout([[col1, col2, name_table, [comp_table, para]], plot, [table, slider]])
    curdoc().clear()
    curdoc().add_root(cur_layout)


def get_slider(mol_list):
    namelist = []
    for mol in mol_list:
        namelist.append(mol.name)
    slider = Select(title="Molecules", options=namelist)
    slider.on_change("value", change_table)
    return slider


def change_table(attr, old, new):
    from bokeh.layouts import layout
    global mol_list, slider
    mol_dict = {}
    for mol in mol_list:
        mol_dict[mol.name] = mol

    prefrence = mol_list[:2] + [mol_dict[slider.value]]
    name_table = get_name_table(mol_list)
    data, names = get_table_values(mol_list, prefrence)
    table = get_table(data, names)
    plot = get_pictures(mol_list, prefrence)
    slider = get_slider(mol_list)
    cur_layout = layout([[col1, col2, name_table, [sim_tab, para]], plot, [table, slider]])
    curdoc().clear()
    curdoc().add_root(cur_layout)


def get_table(data, names):
    data = dict(Mol1=data[0], Mol2=data[1], Mol3=data[2],
                 other=data[6], sum=data[8], s=data[9], name=data[10], link=data[11])

    source = ColumnDataSource(data)
    color_template = """
          <div style="background:<%= 
              (function colorfromint(){
                  if(s == 2){
                      return("rgba(255, 0, 0, .4)")}
                  else if (s  == 1 ){
                       return("rgba(255, 255, 51, .4)")}
                  else 
                  {return("rgba(0, 128, 0, .4)")}
                  }()) %>; 
              color: black"> 
          <%= value %>
          </div>
          """
    link_template = """
                    <a href= <%= link %> > <%= value %> </a>
                    
                    
                    """
    columns = []
    color_formatter = HTMLTemplateFormatter(template=color_template)
    link_formatter = HTMLTemplateFormatter(template=link_template)
    columns.append(TableColumn(field='Mol1', title=names[0],  width=215))
    columns.append(TableColumn(field='Mol2', title=names[1], width=215))
    columns.append(TableColumn(field='Mol3', title=names[2], width=215))
    columns.append(TableColumn(field='other', title="Other", width=100))
    columns.append(TableColumn(field='sum', title='Value', formatter=color_formatter, width=100))
    columns.append(TableColumn(field="name", title="Model", width=150, formatter=link_formatter))
    data_table = DataTable(source=source,
                           columns=columns,
                           index_position=None,
                           fit_columns=True,
                           selectable=True,
                           sortable=False,
                           width=(3*215 + 250), height=450)
    return data_table


def change_text(attr, old, new):
    text.value = smilies_dict_select[select.value]


def get_props(smiles_dict):
    idstrg = list_to_str(idlist)
    mol_list = get_value(smiles_dict, idstrg, iddict, 1)
    return mol_list


def change_form(attr, old, new):
    if "Empty" in x.keys():
        x.pop("Empty")
    tot = 0
    for key in form_dict_select[flavor_select.value]:
        x[key] = form_dict_select[flavor_select.value][key]
        tot += form_dict_select[flavor_select.value][key]
    if tot < 100:
        empt = 100-tot
    x["Empty"] = empt
    p = get_plot(x)
    col2.children.remove(col2.children[0])
    col2.children.append(p)


def add_oil(attr, old, new):
    text.value = oil_add.value


def paragraph(name, error):
    if error < .5:
        add_text = "Which means they are pretty similar"
    elif .5 < error < 1:
        add_text = "Which means they are sort of similar"

    else:
        add_text = ". Which means there are not many similar flavors"
    front = "Your flavor is closest to: "
    second = ". This has and error of "
    third = ("{:05.2f}".format(error))

    text = front + str(name) + second + third +add_text
    para = Paragraph(text=text, width=250)
    return para


def get_name_table(mol_list):
    name_list = []
    comp_list = []

    mol_list.sort(key=lambda x: x.prc, reverse=True)
    for mol in mol_list:
        name_list.append(mol.name)
        comp_list.append("{:05.2f}".format(mol.prc))

    data = dict(Names=name_list, percent=comp_list)
    source = ColumnDataSource(data)
    columns = []
    columns.append(TableColumn(field='Names', title="Constituent", width=175))
    columns.append(TableColumn(field='percent', title="Value", width=50))
    data_table = DataTable(source=source,
                           columns=columns,
                           index_position=None,
                           fit_columns=True,
                           selectable=True,
                           sortable=False,
                           width=225, height=450)
    return data_table


def get_sim_table(name_list, val_list):
    new_val = []
    for val in val_list:
        new_val.append("{:05.2f}".format(val))
    data = dict(Names=name_list, percent=new_val)
    source = ColumnDataSource(data)
    columns = []
    columns.append(TableColumn(field='Names', title="Consitituent", width=200))
    columns.append(TableColumn(field='percent', title="Value", width=50))
    data_table = DataTable(source=source,
                           columns=columns,
                           index_position=None,
                           fit_columns=True,
                           selectable=True,
                           sortable=False,
                           width=250, height=450)
    return data_table

idlist, iddict = get_ids(namelist)
x = {"Empty": 100}
p = get_plot(x)
# Set up widgets
text = TextInput(title="Enter SMILE", value='')
conc = Slider(title="Concentration", value=0.0, start=0, end=100, step=1)
button = Button(label="Add", button_type="success")
button2 = Button(label="Run", button_type="success")
select = Select(title="Pick a Constituent:", value="foo", options=list(smilies_dict_select.keys()))
select.on_change('value', change_text)

flavor_select = Select(title="Pick a Flavor:", value="foo", options=list(form_dict_select.keys()))
flavor_select.on_change('value', change_form)

oil_add = Select(title="Add Oil:", value="foo", options=list(oil_dict_select.keys()))
oil_add.on_change('value', add_oil)


button.on_event(ButtonClick, update_title)
button2.on_event(ButtonClick, run_rules)
# Set up layouts and add to document
col1 = column(text, conc,  select,oil_add, flavor_select, button, button2)
col2 = column(p)
col3 = column()
layout = layout([[col1, col2, col3]])
slider = None
mol_list = None
data = None
para = None
names = None
sim_tab = None
curdoc().add_root(layout)
curdoc().title = "Conc"
