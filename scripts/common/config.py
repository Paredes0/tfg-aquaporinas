"""
Configuración centralizada de rutas.

Reemplaza las rutas hardcoded que aparecían en los scripts originales
(c:\\Users\\Lab.Micaela VI\\Desktop\\Noe Paredes\\...).

Uso desde un script:
    from scripts.common.config import paths

    df = pd.read_csv(paths.tabla_traduccion(), sep='\\t')

Override mediante variable de entorno:
    export TFG_DATA_ROOT="C:/Users/Usuario/Desktop/resultados finales"
    export TFG_RNA_SEQ_ROOT="Z:/work/RNA-seq_test"
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(var_name: str, default: str) -> Path:
    """Devuelve Path desde variable de entorno o default."""
    return Path(os.environ.get(var_name, default))


# ── Raíces principales ──────────────────────────────────────────────────────
# Donde viven los datos primarios del TFG (genoma, FASTA, tablas, etc.).
DATA_ROOT: Path = _env_path(
    'TFG_DATA_ROOT',
    r'C:\Users\Usuario\Desktop\resultados finales'
)

# Donde vive el pipeline de RNA-seq (puede estar en otra unidad)
RNA_SEQ_ROOT: Path = _env_path(
    'TFG_RNA_SEQ_ROOT',
    r'Z:\work\RNA-seq_test'
)

# Carpeta GDR (Genome Database for Rosaceae) auxiliar
GDR_ROOT: Path = _env_path(
    'TFG_GDR_ROOT',
    str(DATA_ROOT / 'GDR_fxa')
)


# ── Subdirectorios de DATA_ROOT ─────────────────────────────────────────────
class Paths:
    """Accesores a archivos individuales. Cada método devuelve una Path concreta."""

    # 5.1 — Curaduría
    @staticmethod
    def proteinas_dir() -> Path:
        return DATA_ROOT / 'analisis proteinas aquaporina'

    @staticmethod
    def tabla_traduccion() -> Path:
        return Paths.proteinas_dir() / 'tabla_aquaporinas_traduccion.tabular'

    @staticmethod
    def clasificacion_simple() -> Path:
        return Paths.proteinas_dir() / 'clasificacion_filogenetica_simple.csv'

    @staticmethod
    def fasta_gff3_peptidos() -> Path:
        return Paths.proteinas_dir() / 'aquaporin_peptides.fasta'

    @staticmethod
    def fasta_exonerate() -> Path:
        return Paths.proteinas_dir() / 'exonerate_genes_aqp.fasta'

    @staticmethod
    def meme_combined() -> Path:
        return Paths.proteinas_dir() / 'MEME_exonerate_gff3_aqp.txt'

    @staticmethod
    def topologias_gff3() -> Path:
        return Paths.proteinas_dir() / 'predicted_topologies_gff3.3line'

    @staticmethod
    def topologias_exonerate() -> Path:
        return Paths.proteinas_dir() / 'predicted_topologies_exonerate.3line'

    @staticmethod
    def pepstats_gff3() -> Path:
        return Paths.proteinas_dir() / 'pepstats_gff3.txt'

    @staticmethod
    def pepstats_exonerate() -> Path:
        return Paths.proteinas_dir() / 'pepstats_exonerate.txt'

    @staticmethod
    def deeploc_gff3() -> Path:
        return Paths.proteinas_dir() / 'deeploc_gff3.csv'

    @staticmethod
    def deeploc_exonerate() -> Path:
        return Paths.proteinas_dir() / 'deeploc_exonerate.csv'

    # 5.2 — Filogenia
    @staticmethod
    def treefile_gff3() -> Path:
        return Paths.proteinas_dir() / 'fxa_aqp_gff3_129_clipkit.fasta.treefile'

    @staticmethod
    def iqtree_gff3() -> Path:
        return Paths.proteinas_dir() / 'fxa_aqp_gff3_129_clipkit.fasta.iqtree'

    @staticmethod
    def treefile_exonerate() -> Path:
        return Paths.proteinas_dir() / 'exonerate_aqp.treefile'

    @staticmethod
    def iqtree_exonerate() -> Path:
        return Paths.proteinas_dir() / 'exonerate_aqp.iqtree'

    @staticmethod
    def treefile_final() -> Path:
        """Árbol final BUENO de 281 secuencias funcionales."""
        return RNA_SEQ_ROOT / 'arbol_acuaporinas_2_bueno_sin_parciales.treefile'

    # 5.3–5.4 — RNA-seq
    @staticmethod
    def rna_seq_dir() -> Path:
        return DATA_ROOT / 'RNA-seq'

    @staticmethod
    def basal_aquaporins() -> Path:
        return Paths.rna_seq_dir() / 'basal_aquaporins'

    @staticmethod
    def basal_tpm() -> Path:
        return Paths.basal_aquaporins() / 'basal_aquaporins_tpm.csv'

    @staticmethod
    def de_leaf_dir() -> Path:
        return Paths.rna_seq_dir() / 'de_leaf'

    @staticmethod
    def de_roots_dir() -> Path:
        return Paths.rna_seq_dir() / 'de_roots'

    # 5.5 — Reanotación
    @staticmethod
    def reanotacion_candidates() -> Path:
        return Paths.basal_aquaporins() / 'reannotation_candidates.tsv'

    # 5.6 — Homeólogos
    @staticmethod
    def homeolog_groups_summary() -> Path:
        return DATA_ROOT / 'aqp_finales' / 'homeolog_groups_summary.tsv'

    @staticmethod
    def homeolog_analysis_dir() -> Path:
        return Paths.rna_seq_dir() / 'homeolog_analysis'

    @staticmethod
    def dominant_subgenome() -> Path:
        return Paths.homeolog_analysis_dir() / 'dominant_subgenome.csv'


paths = Paths()


def verify_paths(strict: bool = False) -> dict[str, bool]:
    """
    Comprueba qué archivos existen en la configuración actual.
    Útil para debugging y para los tests de smoke.

    Args:
        strict: si True, lanza FileNotFoundError ante el primer archivo faltante.

    Returns:
        Dict {nombre_metodo: existe}
    """
    result = {}
    for name in dir(Paths):
        if name.startswith('_'):
            continue
        method = getattr(Paths, name)
        if not callable(method):
            continue
        try:
            p = method()
            exists = p.exists()
            result[name] = exists
            if strict and not exists:
                raise FileNotFoundError(f"No existe: {p} (método paths.{name}())")
        except Exception as e:
            result[name] = False
            if strict:
                raise
    return result


if __name__ == '__main__':
    print(f"DATA_ROOT     = {DATA_ROOT}")
    print(f"RNA_SEQ_ROOT  = {RNA_SEQ_ROOT}")
    print(f"GDR_ROOT      = {GDR_ROOT}")
    print()
    print("Verificación de paths:")
    results = verify_paths(strict=False)
    for name, exists in sorted(results.items()):
        mark = '✓' if exists else '✗'
        print(f"  {mark} paths.{name}()")
    n_ok = sum(results.values())
    print(f"\n  Total: {n_ok}/{len(results)} archivos encontrados.")
