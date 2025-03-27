import os
import gzip
import shutil

# Iterar sobre cada archivo .fastq
for filename in os.listdir('.'):
    if filename.startswith("Galaxy") and filename.endswith(".fastqsanger.gz"):
        
        # Renombrar archivo
        renamed = filename.split('[')[-1].replace(']', '')
        os.rename(filename, renamed)
        print(f"Renombrado: {filename} -> {renamed}")

        # Cambiar extensión final a .fastq
        decompressed_name = renamed.replace('.fastqsanger.gz', '.fastq')

        # Descomprimir
        with gzip.open(renamed, 'rb') as f_in, open(decompressed_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Descomprimido: {renamed} -> {decompressed_name}")

        # Eliminar el archivo .gz original
        os.remove(renamed)
