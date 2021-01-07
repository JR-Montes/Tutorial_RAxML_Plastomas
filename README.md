## Tutoria de Maxima Verosimilitud con Plastomas en RAxML ver. 8

Este es un tutorial para correr un análisis de máxima verosimilitud con genomas completos de cloropasto de pinos
 

by José Rubén Montes

___


En RaxML se pueden implementar una gran cantidad de análisis, además de buscar árboles y calcular los valores de soporte, este programa puede hacer consensos convencionales o específicos para identificar taxones tramposos (rogue taxa) o aplica una prueba de Shimodaira-Hasegawa entre otras cosas.

En el presente tutorial se describen brevemente los análisis más comunes de máxima verosimilitud que se pueden hacer con RAxML v. 8.2.10 (2017), que corresponde a una versión de línea de comandos. Esta versión también se utiliza para correr datos multilocus. 

Recuerda que el tipo de RAxML específico que uses depende de la computadora que tengas, en este caso estoy empleando el `raxmlHPC-PTHREADS` pero el tuyo puede ser diferente. 

El comando `-T 2` hace referencia al número de núcleos del procesador que se van a emplear, si tu computadora solo tiene un núcleo no debes usar ese comando

**Los archivos que necesitas son:**

| Archivos | inputs |  
| ------   | ------ |
|**i.** Matriz molecular   |(secuencias.phy)| 
|**ii.** Archivo de particiones molecular  | (particion_adn.txt)|
|**iii.** Matriz mixta       | (total.phy)   |
|**iv.** Archivo de particiones mixto |(morfo_sec.txt) |


**ii**. Archivo de particiones

```
DNA, cod = 
DNA, nocod = 
```

___

**1.** Búsqueda del mejor árbol de verosimilitud y cálculo de los valores de soporte, con particiones moleculares. Este análisis se realiza en tres partes, dos implican algo de tiempo y la tercera es solo para resumir los resultados. Para correr el análsis se utilizara un archivo en PHYLIP con 63 taxa y 117,333 sitos del genoma de cloroplasto:

`Cembroides_Plastome_63t.phy`


___

El programa de `raxmlHPC-PTHREADS` debe o tiene que estar instalado en la carpeta `bin`dentro de la carpeta `local`en la ruta del usuario `usr` para poder correr el análisis desde cualquier lugar del servidor. `raxmlHPC-PTHREADS` corre con lenguaje de programación `perl`. En este sentido la ruta para llamar a `raxmlHPC-PTHREADS` sería la siguiente:

```
lenguaje ../executable ./ "path[...] -f d -m [modelo] -s [archivo de entrada.phy] -q [archivo de particiones.txt] -# [número de búsquedas] -n [archivo de salida.phy] -T [número de procesadores] -p [número de raíz]"
```

**NOTA**: Es muy importante que escribas las comillas para poder correr el análisis


Ejemplo:

**a)** Búsqueda heurística del mejor árbol con 1000 replicas

```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS -f d -m GTRGAMMA -s Cembroides_Plastome_t63.phy -q Particion_genes.txt -# 1000 -n Heuristica -T 22 -p 12345"
``` 
 

Uno de los archivos de salida será:
`RaxML_bestTree.Analysis.Cembroides_Plastome_63t.phy`

Este es el mejor árbol de las 1000 replicas, puedes verlo en FigTree o Dendroscope. En ese árbol vas a apuntar los resultados del bootstrap, debes verificar que el árbol está en formato NEWICK

**b)** Análisis de bootstrap con 1000 replicas.

```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS -f d -m GTRGAMMA -s Cembroides_Plastome_t63.phy -q Particion_genes.txt -# 1000 -b 12345 -n bootstrap -T 2 -p 12345"
```

**NOTA**: En algunos casos `RAxML` impre un error en la pantalla que indica que el archivo.phy debe ser convertido a `FASTA` para poder seguir con el análisis. Si esto sucede transforma el `archivo.phy` a un `archivo.fasta`

Los árboles de bootsrap (1000 árboles) se guardan en el archivo: `RaxML_bootstrap.Analysis.Cembroides_Plastome_63t.fasta`

Estos árboles puedes conservarlos y en caso de ser necesario juntarlos con otro archivo diferente que contenga —digamos— otras 1000 replicas y así juntar 2000. Esta forma de hacer el bootstrap te permite adicionar replicas que vayas haciendo en diferentes máquinas y no tener que emplear una sola computadora para hacer miles de replicas.


**c)** Resumen final, anotar el árbol con los valores de soporte.

```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS -f b -m GTRGAMMA -s Cembroides_Plastome_t63.fasta -q Particion_genes.txt -z RAxML_bootstrap.*.fasta -t RAxML_bestTree.*.phy -n BS_TREE -T 22"
```

Se genera un árbol llamado:
`RaxML_bipartitions.Analysis.RAxML_info.Analysis.RAxML_bootstrap*.fasta`

Ese árbol lo puedes abrir con FigTree o Dendroscope y ya contiene los valores de soporte anotados sobre el mejor árbol encontrado en la búsqueda heurística  de 1000 replicas.


___

## 2. Análisis de evidencia total con particiones

**a)** Búsqueda del mejor árbol de verosimilitud con particiones moleculares y morfológica (binaria y multiestado).

___

```
raxmlHPC-pthreads-sse3 -f d -m ASC_MULTICAT --asc-corr=lewis -K MK -s total.phy -q morfo_sec.txt -# 1000 -n combinado -T 2 -p 12345
```

Este comando declara dos particiones una para secuencias que por omisión será analizada con `GTR` y otra de morfología que será analizada con el Mkv, calculando verosimilitud condicional a NO tener datos invariantes y considerando CAT para la heterogeneidad de tasas (ASC_MULTICAT). La corrección se hará con el método de Lewis (2001) `(--asc-corr=lewis)` y para la matriz de transición de caracteres multiestado se usara la matriz tipo MK `(-K MK)`. Es importante señalar que en el caso de analizar caracteres binarios el comando `-K` ya no es necesario.

**b)** Análisis de bootstrap con 1000 replicas.

```
raxmlHPC-pthreads-sse3 -f d -m ASC_MULTICAT --asc-corr=lewis -K MK -s total.phy -q morfo_sec.txt -# 1000 -b 12345 -n bootstrap -T 2 -p 12345
```

**c)** Resumen final, anotar el árbol con los valores de soporte.

```
raxmlHPC-pthreads-sse3 -f b -m ASC_MULTICAT --asc-corr=lewis -K MK -s total.phy -q morfo_sec.txt -z RAxML_bootstrap.* -t RAxML_bestTree.* -n BS_TREE -T 2
```

___

## 3. Búsqueda rápida exploratoria con valores de soporte, caracteres moleculares sin particiones.

___

`raxmlHPC-pthreads-sse3 ­-f a ­-p 12345 ­-s secuencias.phy ­-x 12345 ­-# 500 ­-m GTRCAT ­-n exploratoria`

El árbol resultante será
`RaxML_bipartitions.exploratoria`

Ese árbol resultó de un bootstrap de 1000 replicas, con algoritmo rápido, al final de cual se hizo una búsqueda heurística con 1 replica. Sobre el árbol encontrado en la replica heurística se anotaron los valores de soporte.
