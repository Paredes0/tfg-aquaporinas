# Cómo obtener un DOI Zenodo para este repositorio

Pasos a ejecutar **una sola vez**, tras el primer push del repo a GitHub.

## 1. Crear cuenta Zenodo (si no la tienes)

Entrar en [zenodo.org](https://zenodo.org) y registrarse usando la opción "Sign in with GitHub" (es lo más sencillo).

## 2. Conectar GitHub ↔ Zenodo

1. Entrar en [zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/).
2. Autorizar Zenodo (te pedirá permisos para listar tus repos).
3. Buscar `tfg-aquaporinas` en la lista y **activar el switch** (ON).

A partir de ese momento, cada **release tag** que se publique en GitHub generará automáticamente un depósito en Zenodo con DOI propio.

## 3. Hacer el primer release

Desde la línea de comandos local:

```powershell
cd 'C:\Users\Usuario\Desktop\tfg-aquaporinas'
git checkout main
git tag -a v1.0.0-tfg-entregado -m "TFG entregado para defensa"
git push origin main --tags
```

GitHub detectará el tag y mostrará un botón "Create release" en la página del repo. Pulsarlo y publicar el release (puedes dejarlo con título "TFG entregado").

## 4. Recibir el DOI

A los pocos minutos, Zenodo:
1. Crea el depósito.
2. Genera el DOI (formato `10.5281/zenodo.XXXXXXX`).
3. Te envía un email con el enlace.

## 5. Insertar el DOI en el TFG y el README

- En `README.md` de este repo, sustituir `XXXXXXX` por el DOI real (commit post-release).
- En la sección §10 del TFG.docx, sustituir `10.5281/zenodo.XXXXXXX` por el DOI real.

## Política de versiones

- **v1.0.0-tfg-entregado**: versión entregada en mayo de 2026.
- Futuras revisiones (post-defensa, correcciones) → `v1.1.0`, `v1.2.0`, etc. Cada tag genera un DOI propio + un "concept DOI" que apunta a la última versión.

## Inputs externos a subir manualmente a Zenodo

Ver `annexes/A_repo_overview/INPUTS_CONJUNTO_DATOS.md` para la lista de datos derivados que conviene subir como dataset separado en Zenodo (~250 KB). Pueden ir en el mismo depósito que el código o en uno aparte; depende de tu preferencia.
