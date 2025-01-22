# organism-assembly
Este repositorio contiene los scripts necesarios para ensamblar genomas de organismos a partir de lecturas generadas por oxford nanopore.

## Requisitos
- [Nextflow](https://www.nextflow.io/)
- [Dockers](https://www.docker.com/)

## Uso
Para correr el pipeline se debe ejecutar el siguiente comando:
```bash
nextflow run run_assembly.nf --input_file <input_file> --db_name <db_name> --min_length <min_length> --max_length <max_length> --min_quality <min_quality> --coverage <coverage> --organism_size <organism_size>
```

Donde:
- `input_file` es el archivo de lecturas en formato fastq.
- `db_name` es el nombre de la base de datos contaminados a utilizar.
- `min_length` es la longitud mínima de las lecturas a utilizar.
- `max_length` es la longitud máxima de las lecturas a utilizar.
- `min_quality` es la calidad mínima de las lecturas a utilizar.
- `coverage` es la cobertura mínima a utilizar.
- `organism_size` es el tamaño del organismo a ensamblar.

## Ejemplos
Preferentemente, se recomienda tener los archivos organizados de la siguiente manera:
```
organism-assembly/
├── data/
│   ├── SE.gz
│   └── pseudomonas-aeruginosa.fasta.gz
├── db/
│   ├── klebsiella_db.fasta.gz
│   ├── echerichia_db.00.nhr
│   ├── echerichia_db.00.nin
│   ├── echerichia_db.00.nsq
│   ├── echerichia_db.01.nhr
│   ├── echerichia_db.01.nin
│   ├── echerichia_db.01.nsq
│   ├── ..
│   └── echerichia_db.nto
├── results/
├── run_assembly.nf
└── Dockerfile
```
La ejecución puede realizarse con una base de datos ya obtenida de un archivo fasta utilizando
```bash
nextflow run run_assembly.nf -resume --input_file data/pseudomonas-aeruginosa.fasta.gz --db_name echerichia_db
```

También puede realizarse con una base de datos obtenida de un archivo fasta utilizando
```bash
nextflow run run_assembly.nf -resume --input_file data/SE.gz --db_name klebsiella_db.fasta.gz
```
