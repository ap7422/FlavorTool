from rdkit import Chem
from rdkit.Chem import rdEHTTools
from rdkit.Chem import rdDistGeom
from bokeh.plotting import figure
from bokeh.models.tools import WheelZoomTool
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
import io
from PIL import Image


def show_png(data):
    bio = io.BytesIO(data)
    img = Image.open(bio)
    return img


def saveimg(smile, pltnm):
    atorvastatin = Chem.MolFromSmiles(smile)
    mh = Chem.AddHs(atorvastatin)
    rdDistGeom.EmbedMolecule(mh)
    _,res = rdEHTTools.RunMol(mh)
    static_chgs = res.GetAtomicCharges()[:atorvastatin.GetNumAtoms()]
    d = Draw.MolDraw2DCairo(400, 400)
    SimilarityMaps.GetSimilarityMapFromWeights(atorvastatin, list(static_chgs),draw2d=d)
    d.FinishDrawing()
    thing = show_png(d.GetDrawingText())
    name = "http://localhost:5006/fol/static/" + pltnm + ".png"
    thing.save("C:\\Users\\patil.py\\Documents\\11F-Drive\\PFastWebLocalApp\\FlavorTool\\fol\\static\\" +pltnm+".png")

    p = figure(x_range=(0, 1), y_range=(0,1), toolbar_location=None, plot_width=200, plot_height=200)
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.axis.visible = False
    thing = WheelZoomTool()
    p.add_tools(thing)
    p.toolbar.active_scroll = thing

    p.image_url(url=[name], x=0, y=1, w=1, h=1)
    return p
