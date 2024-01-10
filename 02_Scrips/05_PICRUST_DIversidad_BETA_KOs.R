#Vamos a hace run codigo para analiisis senceillo de los datos que obtuvimos de PICRUTS2 


library(tidyverse)
library(nortest)
library(vegan)
library(gridExtra)
library(compositions)


#### Importamos datos de PICRUST pred_metagenome_unstrat_descrip ######
Datos<- read.delim("./01_Datos/picrust2/pred_metagenome_unstrat_descrip.tsv",
                   header=TRUE, sep="\t")

#### Importamos los metadatos#####

metadatos <- read.csv("./01_Datos/sample-metadata.tsv", sep = '\t' )
metadatos<- metadatos[2:9,]
#LOs debemos ordenar para que se lean bien 
metadatos <- metadatos[order(metadatos$sample_name), ]
metadatos<-metadatos[,c(1:4)]


#Cambio nombre de columna de Sample_name a name 
metadatos<-dplyr::rename(metadatos,name=sample_name)


#vamos a cambiar los metadastos de la siguientes forma 
#Este <- Dzilam-Bocas
#Oeste <- El Palmar 

#Este cambio se hace para evitar confusiones en el anañisis y descripción de datos ya que tambien hay sitios de muestreo con esos nombres 

metadatos$Region <- ifelse(metadatos$Region == 'DZILAM-BOCAS', 'Este', 'Oeste')
metadatos$Condicion <- ifelse(metadatos$Condicion == 'alterado', 'Alterado', 'Conservado')


#### Damos formato a los datos de PIICRUT2  ####

#Vamos a eliminar la columna de descripción 
Datos<-Datos[,-2]
#damos rownames 
rownames(Datos)<-Datos[,1]
Datos<-Datos[,-1]
#Transponemos lo datos 
Datos<-as.data.frame(t(Datos))


#### Hacemos transformación CLR #####
DatosCLR<-as.data.frame(clr(Datos))

#### Agregamos metadatos #######

#Generamos columna de sitios 
Datos$name<-row.names(Datos)

#Vamos a unir metadatos 
Datos<-left_join(Datos,metadatos)


#### Sacamos matrices de Aitchison #####



#Se calcula la matriz de bray 
bray<-vegdist(DatosCLR, method = "euclidean")

#Se realiza un NMDS 
set.seed(1234)
nm <- metaMDS(bray)
nm

#Se extraen los puntos del NMDS par graficar 
point<- as.data.frame(nm$points)
point
point <-rownames_to_column(point) 

#Cambiamos nombres de las columnas de df de puntos
colnames(point)<- c("name", "MDS1","MDS2")

point<- merge(point,metadatos)

#Vamos a realizar un PERMANOVA con Adonis para encontrar diferencias entre los sitios de muestreo, condiciones de conservación y region de muestreo para nuesta matriz 

#permanova para Region
permaR<-adonis2(bray~metadatos$Region, permutations = 9999,by = NULL)
p_valueR <- permaR$`Pr(>F)`[1]

#Permanova para estatus de conservacion 
permaC<-adonis2(bray~metadatos$Condicion, permutations = 9999,by = NULL)
p_valueC <- permaC$`Pr(>F)`[1]

#permanova para Sitio de muestreo 
permaS<-adonis2(bray~metadatos$Sitio, permutations = 9999, by = NULL)
p_valueS <- permaC$`Pr(>F)`[1]


#Se pueden ver si la dispersión de los puntos es significantes, lo cual nos puede ayudar a saber si los resultados del permanova son influenciados por la dispersión de los datos 

disperR<-betadisper(bray,metadatos$Region)
anova(disperR)
disperR<-betadisper(bray,metadatos$Condicion)
anova(disperR)

#### Graficamos nuestro NMDS y agregamos los valores obtenidos de la PERMANOVA#######

R<-ggplot(point, aes(x=MDS1,y=MDS2, color=Region))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR KOs-Región", subtitle = paste("stress=", nm$stress), caption = paste("p-value=",p_valueR))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()


C<-ggplot(point, aes(x=MDS1,y=MDS2, color=Condicion))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR KOs-Condición", subtitle = paste("stress=", nm$stress), caption = paste("p-value=",p_valueC))+
  scale_color_manual(values=c("#fdb462", "#b3de69"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()

#Grafica
grid.arrange(C,R, ncol=2)







