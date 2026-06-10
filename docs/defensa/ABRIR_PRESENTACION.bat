@echo off
cd /d "%~dp0"
echo ============================================================
echo   Presentacion de defensa  -  servidor local
echo.
echo   Se abrira sola en el navegador:  http://localhost:8000
echo   DEJA ESTA VENTANA ABIERTA mientras la usas.
echo   Cierra la ventana (o Ctrl+C) para detener el servidor.
echo   (Si la pagina sale en blanco al abrir, pulsa F5 para refrescar.)
echo ============================================================
echo.
start "" http://localhost:8000
python -m http.server 8000 2>nul || py -m http.server 8000 2>nul || (
  echo.
  echo  [ERROR] No se encontro Python en este equipo.
  echo  - Opcion A: instala Python desde python.org y vuelve a ejecutar este archivo.
  echo  - Opcion B: usa la version online:
  echo              https://paredes0.github.io/tfg-aquaporinas/defensa/
  echo.
  pause
)
