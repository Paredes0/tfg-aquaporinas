# Anexo J — Parámetros fisicoquímicos por acuaporina

Tabla por gen de las **121 acuaporinas funcionales** de *F.* x *ananassa* 'Benihoppe'
con los parámetros citados en el cuerpo del TFG (apartados 6.1 y 6.2.3).

## Contenido

| Archivo | Contenido |
|---|---|
| `Anexo_J_parametros_fisicoquimicos.csv` | Una fila por acuaporina (121), ordenada por subfamilia y sub-subfamilia. |
| `Anexo_J_script_reproducible.py` | Genera la tabla desde `data/curado/tabla_aqp_ordenada.csv` vía `scripts/common/config.py`. |

### Columnas

- `gene_id` — identificador del locus (p. ej. Fxa5Ag00335).
- `subfamilia` / `sub_subfamilia` — asignación filogenética (PIP, TIP, NIP, SIP, XIP) y sub-subfamilia (FaPIP1, FaSIP1;3, …).
- `longitud_aa`, `pI`, `Mw_kDa`, `TMHs` — longitud, punto isoeléctrico, peso molecular (kDa) y nº de hélices transmembrana (DeepTMHMM).
- `NPA_bucle_B` / `NPA_bucle_E` — variante del motivo NPA en cada bucle (NPA canónico o variantes NPL/NPT/NPS, relevantes en SIP).
- `localizacion` — localización subcelular predicha (DeepLoc 2.0).

## Cifras de control

- 121 acuaporinas (PIP 32, TIP 34, NIP 37, SIP 12, XIP 6).
- Motivo NPA del bucle B en SIP: **NPT ×2, NPS ×7, NPL ×3** — coincide con el texto del TFG (NPL en las tres SIP2;2, NPT en dos SIP1, NPS en las siete SIP1;3).

Este anexo es el referenciado en el TFG como soporte de los parámetros fisicoquímicos
por gen y de las variantes NPA por SIP.
