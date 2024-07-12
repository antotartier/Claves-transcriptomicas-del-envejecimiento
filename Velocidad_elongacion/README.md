# Cálculo de la velocidad de la elongación de la RNA polimerasa II
Cálculo del ritmo de la elongación a partir de los archivos BAM ordenados e indexados del experimento de RNAseq. 
## Anotiación de los intrones
El primero de los pasos es generar un archivo de anotación de los intrones en los que se calculará la velocidad de la elongación. En este sentido, los scripts llamados "recursive_annot" generan dicha anotación partiendo de los archivos SJ de STAR, del GTF y el archivo del index de STAR "chrNameLength.txt". 
## Cálculo de la pendiente de los intrones
Una vez generada la anotación de los intrones, se calcula la pendiente de su cobertura en las diferentes muestras utilizando el script "txspeed.py". Antes de ejecutar este script es necesario crear un entorno python con HTSeq, numpy y scikit-learn instalados.
```
conda create -n deepseq -c conda-forge HTSeq numpy scikit-learn
conda activate deepseq
```
Una vez activado el entorno, el script se corre desde la terminal (una vez se le han otorgado los permisos de ejecución) y sus argumentos están disponibles ejecutando el comando `./txspeed.py --help`.
## Análisis de las diferencias en la velocidad de la elongación
Tras calcular la pendiente de los intrones, está se transformó en la velocidad mediante la fórmula `v=-1/s`, donde "v" es la velocidad y "s" la pendiente. Una vez calculado el ritmo de la elongación se analizaron los cambios en el mismo utilizando los scritps llamados "speed_analysis".
## Análisis de la expresión de los reguladores de la velocidad de la elongación
En el caso de los ratones LAKI, se llevó a cabo un análisis de la expresión de los reguladores positivos y negativos del ritmo de la elongación de la Pol II. Para ello, se utilizó el script "reg_elong_LAKI" y los sets de genes de dichos reguladores descargados de la base de datos MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp).