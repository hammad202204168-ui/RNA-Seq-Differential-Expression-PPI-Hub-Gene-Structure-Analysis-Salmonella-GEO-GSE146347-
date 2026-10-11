# PyMOL: cartoon renders coloured by AlphaFold confidence (pLDDT is stored in the B-factor column).
# NOT tested in the environment that wrote it (no PyMOL there). Run from the project folder:  pymol -cq scripts/render_alphafold_cartoons.pml
# One PNG per protein into results/figures/cartoons/ ; edit the glob if your files live elsewhere.
python
import glob, os
from pymol import cmd
os.makedirs("results/figures/cartoons", exist_ok=True)
for f in sorted(glob.glob("data/alphafold/AF-*-F1-model_v*.pdb")):
    name = os.path.basename(f).split("-")[1]
    cmd.reinitialize(); cmd.load(f, name)
    cmd.hide("everything"); cmd.show("cartoon"); cmd.bg_color("white")
    cmd.color("0x0053D6", "b > 90")
    cmd.color("0x65CBF3", "b < 90 and b > 70")
    cmd.color("0xFFDB13", "b < 70 and b > 50")
    cmd.color("0xFF7D45", "b < 50")
    cmd.orient(); cmd.set("ray_opaque_background", 1); cmd.set("cartoon_fancy_helices", 1)
    cmd.png(f"results/figures/cartoons/{name}.png", width=1600, height=1600, dpi=300, ray=1)
python end
