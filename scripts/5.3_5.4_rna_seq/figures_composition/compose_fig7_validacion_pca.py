import fitz
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

base = Path("/home/noe/work/RNA-seq_test/results/basal_aquaporins")
out_pdf = base / "figura7_validacion.pdf"
out_png = base / "figura7_validacion.png"

paneles = [
    ("pca_aquaporins.pdf", "A) PCA de las 22 muestras paired-end"),
    ("correlation_samples.pdf", "B) Correlacion Pearson muestra-a-muestra"),
]

def pdf_to_array(path, zoom=3.0):
    doc = fitz.open(str(path))
    page = doc[0]
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, alpha=False)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    doc.close()
    return img

fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=130)
for ax, (pdf_name, label) in zip(axes, paneles):
    img = pdf_to_array(base / pdf_name)
    ax.imshow(img)
    ax.set_title(label, fontsize=13, loc="left", fontweight="bold", pad=8)
    ax.axis("off")
fig.tight_layout(pad=1.2)
fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
print("Fig7 PDF:", out_pdf.stat().st_size)
print("Fig7 PNG:", out_png.stat().st_size)
