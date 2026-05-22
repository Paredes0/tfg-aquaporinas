#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generar visualizaciones profesionales para TFG
Distribución de longitudes de aquaporinas de Rosaceae
"""

import matplotlib.pyplot as plt
import numpy as np
import statistics
from matplotlib.patches import Rectangle
import seaborn as sns

# Configurar estilo para publicación científica
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def parse_fasta(filename):
    """Parse FASTA file and return list of sequence lengths"""
    lengths = []
    current_seq = []
    
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    lengths.append(len(''.join(current_seq)))
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq:
            lengths.append(len(''.join(current_seq)))
    
    return lengths

def calculate_statistics(lengths):
    """Calcular estadísticas descriptivas"""
    return {
        'n': len(lengths),
        'min': min(lengths),
        'max': max(lengths),
        'mean': statistics.mean(lengths),
        'median': statistics.median(lengths),
        'stdev': statistics.stdev(lengths) if len(lengths) > 1 else 0,
        'p25': np.percentile(lengths, 25),
        'p75': np.percentile(lengths, 75),
    }

def create_histogram(lengths, output_file='histograma_longitudes_aqp.png'):
    """Crear histograma profesional para TFG"""
    
    stats = calculate_statistics(lengths)
    
    # Crear figura con tamaño apropiado para publicación
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Crear histograma
    n, bins, patches = ax.hist(lengths, bins=30, color='#3498db', 
                               edgecolor='black', alpha=0.7, linewidth=1.2)
    
    # Añadir líneas verticales en 140 y 380 aa
    ax.axvline(x=140, color='red', linestyle='--', linewidth=2.5, 
               label='Umbral inferior (140 aa)', alpha=0.8)
    ax.axvline(x=380, color='red', linestyle='--', linewidth=2.5, 
               label='Umbral superior (380 aa)', alpha=0.8)
    
    # Sombrear la región entre 140 y 380
    ax.axvspan(140, 380, alpha=0.1, color='green', 
               label=f'Rango seleccionado (140-380 aa)')
    
    # Añadir línea vertical para la media
    ax.axvline(x=stats['mean'], color='orange', linestyle=':', 
               linewidth=2, label=f'Media ({stats["mean"]:.1f} aa)', alpha=0.8)
    
    # Etiquetas y títulos
    ax.set_xlabel('Longitud de la secuencia (aminoácidos)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frecuencia (número de secuencias)', fontsize=14, fontweight='bold')
    # Título en el pie de figura (norma APA/UCAM), no embebido en la imagen
    
    # Leyenda
    ax.legend(loc='upper right', fontsize=11, framealpha=0.95)
    
    # Grid más suave
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Mejorar apariencia de los ejes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Añadir texto con estadísticas clave
    textstr = f'Media ± DE: {stats["mean"]:.1f} ± {stats["stdev"]:.1f} aa\n'
    textstr += f'Rango: {stats["min"]} - {stats["max"]} aa\n'
    textstr += f'Mediana: {stats["median"]:.1f} aa'
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    # Ajustar layout
    plt.tight_layout()
    
    # Guardar figura
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] Histograma guardado: {output_file}")
    
    # También guardar en formato vectorial (mejor calidad para publicación)
    output_pdf = output_file.replace('.png', '.pdf')
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] Version vectorial guardada: {output_pdf}")
    
    plt.close()

def create_boxplot(lengths, output_file='boxplot_longitudes_aqp.png'):
    """Crear boxplot adicional"""
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    # Crear boxplot horizontal
    bp = ax.boxplot([lengths], vert=False, widths=0.6, patch_artist=True,
                     boxprops=dict(facecolor='#3498db', alpha=0.7),
                     medianprops=dict(color='red', linewidth=2),
                     whiskerprops=dict(linewidth=1.5),
                     capprops=dict(linewidth=1.5))
    
    # Añadir líneas verticales en 140 y 380
    ax.axvline(x=140, color='red', linestyle='--', linewidth=2, 
               label='Umbral 140 aa', alpha=0.8)
    ax.axvline(x=380, color='red', linestyle='--', linewidth=2, 
               label='Umbral 380 aa', alpha=0.8)
    
    # Sombrear región
    ax.axvspan(140, 380, alpha=0.1, color='green')
    
    ax.set_xlabel('Longitud de la secuencia (aminoácidos)', fontsize=14, fontweight='bold')
    # Título en el pie de figura (norma APA/UCAM), no embebido en la imagen
    ax.set_yticks([])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] Boxplot guardado: {output_file}")
    plt.close()

def create_statistics_table(lengths, output_file='tabla_estadisticas.txt'):
    """Crear tabla de estadísticas para el TFG"""
    
    stats = calculate_statistics(lengths)
    
    # Calcular porcentajes en rangos
    in_range_140_380 = sum(1 for l in lengths if 140 <= l <= 380)
    below_140 = sum(1 for l in lengths if l < 140)
    above_380 = sum(1 for l in lengths if l > 380)
    
    percent_in_range = (in_range_140_380 / stats['n']) * 100
    percent_below = (below_140 / stats['n']) * 100
    percent_above = (above_380 / stats['n']) * 100
    
    # Crear tabla en formato texto
    table_text = "=" * 80 + "\n"
    table_text += "TABLA RESUMEN: ESTADISTICAS DE LONGITUDES DE AQUAPORINAS DE ROSACEAE\n"
    table_text += "=" * 80 + "\n\n"
    
    table_text += "+----------------------------------------+-------------------+\n"
    table_text += "| Parametro Estadistico                  | Valor             |\n"
    table_text += "+----------------------------------------+-------------------+\n"
    table_text += f"| N (Total secuencias)                   | {stats['n']:>17} |\n"
    table_text += f"| Longitud Minima (aa)                   | {stats['min']:>17} |\n"
    table_text += f"| Longitud Maxima (aa)                   | {stats['max']:>17} |\n"
    table_text += f"| Media (aa)                             | {stats['mean']:>17.2f} |\n"
    table_text += f"| Desviacion Estandar (aa)               | {stats['stdev']:>17.2f} |\n"
    table_text += f"| Mediana (aa)                           | {stats['median']:>17.2f} |\n"
    table_text += f"| Percentil 25 (aa)                      | {stats['p25']:>17.2f} |\n"
    table_text += f"| Percentil 75 (aa)                      | {stats['p75']:>17.2f} |\n"
    table_text += "+----------------------------------------+-------------------+\n"
    table_text += f"| % dentro de 140-380 aa                 | {percent_in_range:>16.2f}% |\n"
    table_text += f"| Secuencias en rango 140-380            | {in_range_140_380:>17} |\n"
    table_text += f"| Secuencias < 140 aa                    | {below_140:>17} |\n"
    table_text += f"| Secuencias > 380 aa                    | {above_380:>17} |\n"
    table_text += "+----------------------------------------+-------------------+\n"
    
    # Guardar tabla
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(table_text)
    
    print(f"[OK] Tabla guardada: {output_file}")
    print("\n" + table_text)
    
    # También crear versión CSV para Excel
    csv_file = output_file.replace('.txt', '.csv')
    with open(csv_file, 'w', encoding='utf-8') as f:
        f.write("Parámetro Estadístico,Valor\n")
        f.write(f"N (Total secuencias),{stats['n']}\n")
        f.write(f"Longitud Mínima (aa),{stats['min']}\n")
        f.write(f"Longitud Máxima (aa),{stats['max']}\n")
        f.write(f"Media (aa),{stats['mean']:.2f}\n")
        f.write(f"Desviación Estándar (aa),{stats['stdev']:.2f}\n")
        f.write(f"Mediana (aa),{stats['median']:.2f}\n")
        f.write(f"Percentil 25 (aa),{stats['p25']:.2f}\n")
        f.write(f"Percentil 75 (aa),{stats['p75']:.2f}\n")
        f.write(f"% dentro de 140-380 aa,{percent_in_range:.2f}\n")
        f.write(f"Secuencias en rango 140-380,{in_range_140_380}\n")
        f.write(f"Secuencias < 140 aa,{below_140}\n")
        f.write(f"Secuencias > 380 aa,{above_380}\n")
    
    print(f"[OK] Version CSV guardada: {csv_file}")
    
    return table_text

def create_pie_chart(lengths, output_file='grafico_distribucion_rangos.png'):
    """Crear gráfico de pastel con distribución de rangos"""
    
    # Calcular distribuciones
    below_140 = sum(1 for l in lengths if l < 140)
    in_range = sum(1 for l in lengths if 140 <= l <= 380)
    above_380 = sum(1 for l in lengths if l > 380)
    
    sizes = [below_140, in_range, above_380]
    labels = [f'< 140 aa\n({below_140} seq.)', 
              f'140-380 aa\n({in_range} seq.)', 
              f'> 380 aa\n({above_380} seq.)']
    colors = ['#ff9999', '#90EE90', '#ffcc99']
    explode = (0.05, 0.1, 0.05)  # Resaltar el rango principal
    
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    
    wedges, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, 
                                       colors=colors, autopct='%1.1f%%',
                                       shadow=True, startangle=90,
                                       textprops={'fontsize': 12, 'fontweight': 'bold'})
    
    # Mejorar formato de porcentajes
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontsize(14)
        autotext.set_fontweight('bold')
    
    # Título en el pie de figura (norma APA/UCAM), no embebido en la imagen
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] Grafico de pastel guardado: {output_file}")
    plt.close()

def main():
    print("=" * 80)
    print("GENERANDO VISUALIZACIONES PROFESIONALES PARA TFG")
    print("=" * 80)
    print()
    
    # Archivo de entrada (ruta portable: override con $TFG_DATA_ROOT)
    import os
    BASE_DIR = os.environ.get('TFG_DATA_ROOT', r'C:\Users\Usuario\Desktop\resultados finales')
    filename = os.path.join(BASE_DIR, 'dataset fragaria ananassa', 'ncbi_data',
                            'TODAS_ROSACEAE_aqp_prot_FINAL_V13.fasta')

    print(f"Leyendo archivo: {filename}")
    lengths = parse_fasta(filename)
    print(f"Total de secuencias leídas: {len(lengths)}\n")

    # Crear directorio de salida
    output_dir = os.path.join(BASE_DIR, 'visualizaciones_tfg')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Directorio de salida: {output_dir}\n")
    
    # Generar visualizaciones
    print("Generando visualizaciones...\n")
    
    # 1. Histograma principal
    create_histogram(lengths, 
                     os.path.join(output_dir, 'histograma_longitudes_aqp.png'))
    
    # 2. Boxplot
    create_boxplot(lengths, 
                   os.path.join(output_dir, 'boxplot_longitudes_aqp.png'))
    
    # 3. Gráfico de pastel
    create_pie_chart(lengths, 
                     os.path.join(output_dir, 'grafico_distribucion_rangos.png'))
    
    # 4. Tabla de estadísticas
    create_statistics_table(lengths, 
                           os.path.join(output_dir, 'tabla_estadisticas.txt'))
    
    # 5. Exportar datos crudos para Excel
    excel_file = os.path.join(output_dir, 'datos_longitudes.csv')
    with open(excel_file, 'w', encoding='utf-8') as f:
        f.write("Secuencia,Longitud_aa\n")
        for i, length in enumerate(sorted(lengths), 1):
            f.write(f"{i},{length}\n")
    print(f"[OK] Datos exportados a CSV: {excel_file}")
    
    print("\n" + "=" * 80)
    print("[OK] PROCESO COMPLETADO")
    print("=" * 80)
    print(f"\nTodos los archivos están en: {output_dir}")
    print("\nArchivos generados:")
    print("  - histograma_longitudes_aqp.png (y .pdf)")
    print("  - boxplot_longitudes_aqp.png")
    print("  - grafico_distribucion_rangos.png")
    print("  - tabla_estadisticas.txt (y .csv)")
    print("  - datos_longitudes.csv")

if __name__ == "__main__":
    main()
