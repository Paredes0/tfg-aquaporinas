# data/ — Datos derivados para reproducir el análisis

Datos derivados (ligeros) que permiten reproducir la filogenia, el Anexo I y las figuras de RNA-seq sin necesidad de re-ejecutar el pipeline completo ni descargar los datos primarios pesados.

> Los **datos primarios** (genoma 'Benihoppe' FASTA/GFF3, 22 muestras RNA-seq paired-end) **no** se incluyen por tamaño: se referencian con sus identificadores públicos en `annexes/A_repo_overview/INPUTS_CONJUNTO_DATOS.md`.

## Estructura

```
data/
├── filogenia/                          Reconstrucción filogenética definitiva (§5.4 / 6.2)
│   ├── arbol_acuaporinas.treefile      Árbol final con soportes SH-aLRT/aBayes/UFboot (282 hojas, 277 nodos)
│   ├── arbol_acuaporinas.iqtree        Informe IQ-TREE: modelo Q.PLANT+R6, log L = −45.149,26
│   ├── arbol_acuaporinas.contree       Árbol consenso de ultrafast bootstrap
│   └── fxa_without_partials.clipkit.fasta   Alineamiento final (MAFFT E-INS-i + ClipKIT, 430 sitios)
└── rna_seq/                            Matrices derivadas del RNA-seq (§5.5 / 6.3)
    ├── basal/                          Expresión basal (6.3.1)
    │   ├── basal_aquaporins_tpm.csv          TPM por gen y muestra (121 funcionales + 8 parciales = 129)
    │   ├── basal_aquaporins_summary.csv      Long-format gen × tejido
    │   ├── basal_aquaporins_detection.csv    Booleano TPM > 1 por tejido
    │   └── basal_aquaporins_normalized.csv   Cuentas normalizadas DESeq2
    ├── de/                             Expresión diferencial control vs estrés (6.3.2)
    │   ├── de_aquaporins_leaf.csv            DESeq2 hoja, restringido a acuaporinas
    │   └── de_aquaporins_roots.csv           DESeq2 raíz, restringido a acuaporinas
    └── homeologos/                     Homeólogos y dominancia subgenómica (6.3.3)
        ├── homeolog_groups.tsv               121 genes → 32 grupos homeólogos
        ├── homeolog_groups_summary.tsv       32 grupos × métricas
        ├── collapsed_tpm.csv / collapsed_counts.csv   Expresión agregada por grupo
        ├── dominance_overall.csv / dominance_by_tissue.csv   Proporción por subgenoma
        ├── dominant_subgenome.csv            Subgenoma dominante por grupo
        └── summary_statistics.csv            Métricas clave
```

## Notas

- El DESeq2 completo del transcriptoma (`results_leaf.csv` / `results_roots.csv`, ~17 MB cada uno, ~109.000 genes) **no se incluye**: las figuras y el análisis usan solo el subconjunto de acuaporinas (`de_aquaporins_*.csv`). El resultado completo se regenera con `scripts/5.5_rna_seq/5.5.2_de_abundancia/07_de_analysis.R`.
- Los datos de curado por herramienta (DeepTMHMM, MEME, Pepstats, DeepLoc, FASTAs, consenso, GFF3 corregido) están en `annexes/B_curado_gff3_vs_exonerate/`.
- El árbol incluido es el de 282 hojas / 277 nodos (con la rama parcial `Fxa6Dg03790`), que es el que respalda las cifras citadas en §6.2.2 del TFG. Detalle en `annexes/I_soportes_filo/README.md`.
