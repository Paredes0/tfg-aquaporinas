# Anexo H — Parámetros fisicoquímicos por acuaporina

Tabla por gen de las **121 acuaporinas funcionales** de *F.* x *ananassa* 'Benihoppe'
con los parámetros citados en el cuerpo del TFG (apartados 6.1 y 6.2.3).

## Contenido

| Archivo | Contenido |
|---|---|
| `Anexo_B_parametros_fisicoquimicos.csv` | Una fila por acuaporina (121), ordenada por subfamilia y sub-subfamilia. |
| `Anexo_H_script_reproducible.py` | Genera la tabla desde `data/curado/tabla_aqp_ordenada.csv` vía `scripts/common/config.py`. |

### Columnas

- `gene_id` — identificador del locus (p. ej. Fxa5Ag00335).
- `subfamilia` — asignación filogenética (PIP, TIP, NIP, SIP, XIP).
- `isoforma` — sub-subfamilia con nomenclatura completa (p. ej. SIP1;3, PIP2;1, SIP2;2), tomada del agrupamiento homeólogo (`homeolog_groups.tsv`). Distingue, p. ej., las SIP1;1 (NPT) de las SIP1;3 (NPS).
- `modelo` — fuente de la secuencia usada: **«GFF3 oficial»** o **«Exonerate (sustituido)»**. 20 de las 121 emplean el modelo de Exonerate por su mayor integridad estructural; las otras 101, el del GFF3 oficial.
- `mRNA_modelo` — identificador del transcrito realmente empleado (mRNA del GFF3 o de Exonerate, según `modelo`).
- `longitud_aa`, `pI`, `Mw_kDa`, `TMHs` — longitud, punto isoeléctrico, peso molecular (kDa) y nº de hélices transmembrana (DeepTMHMM).
- `NPA_bucle_B` / `NPA_bucle_E` — variante del motivo NPA en cada bucle (NPA canónico o variantes NPL/NPT/NPS/NPI).
- `localizacion` — localización subcelular predicha (DeepLoc 2.0).

## Cifras de control (uniformes con el TFG)

- 121 acuaporinas (PIP 32, TIP 34, NIP 37, SIP 12, XIP 6).
- Modelo de secuencia: **101 GFF3 oficial + 20 Exonerate** (sustituidos por mayor integridad; mismos 20 del GFF3 corregido del Anexo B).
- Motivo NPA del **bucle B** en SIP: **NPL ×3** (las tres SIP2;2), **NPT ×2** (SIP1;1, ambas de Exonerate), **NPS ×7** (las siete SIP1;3) — ningún NPA canónico.
- Motivo NPA del **bucle E** en SIP: **NPA ×8 + NPI ×4**.

Ambos recuentos coinciden exactamente con el apartado 6.2.3 del TFG. Este anexo es el
referenciado como soporte de los parámetros fisicoquímicos por gen y de las variantes NPA.
