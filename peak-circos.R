## contexto: https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/07_handling_peaks_bedtools.html

# cargar paquetes
pacman::p_load( "dplyr", "vroom", "circlize" )

# cargamos los datos .bed de la primera muestra
peaks1 <- vroom( file = "https://data.biofreelancer.com/sample1-cromatin", col_names = FALSE )

# cargamos los datos .bed de la segunda muestra
peaks2 <- vroom( file = "https://data.biofreelancer.com/sample2-cromatin", col_names = FALSE )

# renombrar columnas
peaks1 <- peaks1 %>%           # a partir de la data original 
  rename( chromosome = 1,      # nombre_nuevo = numero de columna
          start = 2,
          end = 3,
          name = 4,
          score = 5,
          strand = 6,
          signal = 7,
          pvalue = 8,
          qvalue = 9,
          peak = 10 )

# rename the cols
peaks2 <- peaks2 %>%           # a partir de la data original 
  rename( chromosome = 1,      # nombre_nuevo = numero de columna
          start = 2,
          end = 3,
          name = 4,
          score = 5,
          strand = 6,
          signal = 7,
          pvalue = 8,
          qvalue = 9,
          peak = 10 )

# Reordenamos los picos porque el circos se dibuja en el orden de la tabla
peaks1 <- peaks1 %>% 
  arrange(chromosome, start)    # Se ordena primero por chromosoma, luego por posicion de inicio

peaks2 <- peaks2 %>%
  arrange(chromosome, start)

# Find the highest signal
max_peak <- c( peaks1$signal, peaks2$signal ) %>%     # juntamos todas las seniales
  max( na.rm = TRUE ) %>%                             # encontramos el valor maximo, descartando NA's
  ceiling( )                                          # redondeamos al sigiuente entero

# Iniciamos el cromosoma con la especie humana 
# en version hg38 (GRCh38) del genoma humano
circos.initializeWithIdeogram( species = "hg38",
                               track.height = 0.2,        # este valor es un % (0.2 = 20%) del radio del circos
                               axis.labels.cex = 0.1 )    # reducimos el tamanio de los numeros en el circos

# Cada nuevo track se debe inciar
circos.track(
  ylim = c( 0, max_peak ),          # cual es el "alto" de este track
  bg.border = NA,                   # de que color es el borde del track
  track.height = 0.1                # este valor es un % (0.2 = 20%) del radio del circos
)

# Vamos a dibujar cada segmento (cada cromosoma), uno por uno
# Se hace con un loop for
# se dibuja para cada valor unico en el primer dataframe
for ( chrom in unique( peaks1$chromosome ) ) {
  
  chrom_data <- peaks1 %>%
    filter( chromosome == chrom )    # sacamos los valores para el cromosoma que toca
  
  circos.lines(                      # dibujamos lineas en el circos
    x = chrom_data$start,            # x lo define el inicio del peak
    y = chrom_data$signal,           # y lo define la senial de open cromatin
    sector.index = chrom,            # aqui indicamos que cromosoma estamos dibujando
    type = "l",                      # dibujamos lineas continuas
    col = "green4"                   # la linea de color verde
  )
  
}

# ...Cada nuevo track se debe inciar
circos.track(
  ylim = c( 0, max_peak ),          # cual es el "alto" de este track
  bg.border = NA,                   # de que color es el borde del track
  track.height = 0.1                # este valor es un % (0.2 = 20%) del radio del circos
)

# Repetimos el loop, pero para el dataset 2
for ( chrom in unique( peaks2$chromosome ) ) {     # esta vez es con peaks2
  
  chrom_data <- peaks2 %>%
    filter(chromosome == chrom)     # sacamos los valores para el cromosoma que toca
  
  # Esto es lo mismo que en el loop anterior, solo cambia el color
  circos.lines(
    x = chrom_data$start,
    y = chrom_data$signal,
    sector.index = chrom,   
    type = "l",             
    col = "red4"
  )
  
}

# Agregamos un titulo al plot
title( main = "Open Chromatin in 2 samples",
       cex.main = 1.5 )                          # este valor es el tamanio de la letra

# Cerramos el circos
circos.clear()

# No olvides guardar el plot en la ventana "Plots"
# Como un 800 x 800 PNG
# FIN DE LA CLASE!
