# Al llarg d'aquest codi s'utilitzaran diferents paquets, per la qual cosa, 
# aquells que no estiguen instal·lats, s'hauran d'instal·lar prèviament.

#Codi exemple d'instal·lació:

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("SummarizedExperiment")


# Creació d'un objecte .Rda:

# Carreguem els paquets necessaris:
library(SummarizedExperiment)
# Desem les dades com a data frame, indicant que les dades estan separades 
#per tabuladors:
data <- read.csv("ST000002_AN000002_clean.csv", sep="\t", stringsAsFactors=FALSE)
data <- as.data.frame(data)

# Fem algunes modificacions.
# Volem afegir els noms dels metabòlits (en la primera columna) 
# com a noms de les files
rownames(data) <- data[,1]
#Eliminem la primera columna:
data <- data[,-1]

# Veiem les primeres entrades del data frame:
head(data)


# A continuació, procedim a determinar cada un dels elements per 
# a crear l'objecte:

# 1. Dades experimental: matriu amb les dades de l'experiment
# En primer lloc eliminem la fila amb els grups
# i ho desem a expr_data
expr_data <- data[2:nrow(data), 1:ncol(data)]
# La matriu sols contindrà els valors numèrics
expr_data[] <- lapply(expr_data, as.numeric)
str(expr_data)
# Convertim el dataframe en matriu:
counts <- as.matrix(expr_data)


# 2. colData: DataFrame amb les metadades de cada mostra, en aquest cas,
# cada una de les columnes:

# Establim diferent informació: l'ID de la mostra (sample_id) i un factor
# amb el grup al qual pertanyen les mesures (transplantation).

sample_id <- colnames(data[1:ncol(data)]) # Nom de les columnes
transplantation <- factor(data[1, 1:ncol(data)], levels = c("After", "Before"), labels = c("After", "Before"))

# Addicionalment afegim l'ID de l'estudi, el tipus de mostra i l'origen de la 
# mostra. En aquest cas, totes les mostres tenen les mateixes:
study_id <- 'ST000002'
sample_type <- 'Tissue'
sample_source <- 'Intestine'

#Creem colData amb la informació que hem creat:
colData <- DataFrame(study_id, sample_id, sample_type, sample_source, transplantation)

# rowData: DataFrame amb informació sobre cada una de les entrades 
# (en aquest cas, metabòlits)

metabolite_name <- rownames(data[-1,]) # Nom de les files
# Creem rowData
rowData <- DataFrame(metabolite_name)

# Afegim metadades addicionals en referència al projecte i l'estudi:
data_source <- "Metabolomics Workbench"
project_title <- "Intestinal Samples II pre/post transplantation"
study_type <- "GC-MS Analysis"
units <- 'Peak height'
# Aquestes s'afegiran després de crear l'objecte.

# Creem l'objecte SummarizedExperiment:
se <- SummarizedExperiment(assays = SimpleList(counts), 
                           colData = colData, rowData = rowData)

# Afegim les metadades:
metadata(se)$data_source <- data_source
metadata(se)$project_title <- project_title
metadata(se)$study_type <- study_type
metadata(se)$units <- units


# Desem l'objecte en format .Rda.
save(se, file = "se-PAC1-MicoFrauEstella.Rda")


# Anàlisi exploratòria de les dades:
# Carreguem els paquets necessaris:
library(dplyr)
library(stats)
library(ggplot2)
library("devtools")
# Per instal·lar factoextra: install_github("kassambara/factoextra")
library(factoextra)
library(corrplot)
library(countdata)
library(dplyr)

# Abans de començar amb l'anàlisi, fem alguns canvis en les dades:

# Trasposem les files i columnes, així serà més fàcil operar amb els grups.
# Per fer-ho trasposem la matrriu de les dades i la reconvertim a dataframe.
data_t <- as.data.frame(t(as.matrix(data)))
# Fem dels grups, factors.
data_t$Groups <- factor(data_t$Groups, levels = c("After", "Before"), labels = c("After", "Before"))
# Ara ja tenim les variables (cada un dels metabòlits) com a columnes i els grups com a factors


# De nou,  hem de convertir les dades a valors numèrics.
data_t[, names(data_t) != "Groups"] <- lapply(
  data_t[, names(data_t) != "Groups"], 
  as.numeric
)

#En primer lloc veiem un resum de les dades. Aquí veiem les mitjanes 
#dels valors per a cada metabòlit.
summary(data_t)
# A partir de les dades obtingudes
# l'apartat anterior també podem comparar entre mostres:
summary(counts)
# Necessitarem els grups com a factors:
grp <- as.factor(data_t$Groups)
# Utilitzem prcomp() per fer l'anàlisi de components
data_t.pca <- prcomp(data_t[,-1], center = TRUE, scale = TRUE)
# Representa el percentatge de les variàncies explicat per cada component:
# En cas de tindre error: dev.off()
fviz_eig(data_t.pca, addlabels = TRUE)
# En aquest cas, a partir del component 6 s'explica quasi un 80% de la variància.

# Utilitzem aquesta funció per representar el resultat de la PCA. S'agrupen 
#els individus (mostres) que tenen un perfil similar.
fviz_pca_ind(data_t.pca, col.ind = grp, addEllipses = FALSE)

# També es poden representar les variables. Les variables positivament 
#correlacionades apunten cap al mateix quadrant del gràfic i estan agrupades, les que no, apunten 
#a quadrants oposats.
# col.var = "contrib" ens permet representar amb colors la contribució
# d'aquestes variables. Com major és aquest valor, més contribueix la variable al component

fviz_pca_var(data_t.pca, col.var = "contrib", gradient.cols = c("blue", "yellow", "red"))

# Aquest gràfic representa els metabòlits que més contribueixen als dos components.
fviz_contrib(data_t.pca, choice = "var", axes = 1:2, top = 20)

#A més a més, també podem fer una anàlisi de fold change.
# per a cada metabòlit, calculem la mitjana per a after i before i el fold change
names <- colnames(data_t[,-1])
fca <- data.frame()
for (var in names){
  mean_a <- mean(data_t[data_t$Groups == "After", var], na.rm = TRUE)
  mean_b <- mean(data_t[data_t$Groups == "Before", var], na.rm = TRUE)
  fold_change <- mean_b/mean_a
  fca <- rbind(fca,data.frame(variable = var, fold_change = fold_change))
}
# s'ha generat un data frame, el qual modificarem perquè 
#tinga dues variables, meatbòlit i fold change:

rownames(fca) <- fca[,1]
colnames(fca) <- c("metabolite","fold_change")
fca <- subset(fca, select = "fold_change")

# podem representar el foldchange
plot(fca$fold_change)

# Ordenem els metabòlits per valor de fold change
fca_top <- arrange(fca, desc(fold_change))
# Seleccionem els 5 primers i els 5 últims:
head(fca_top, 5)
tail(fca_top, 5)
