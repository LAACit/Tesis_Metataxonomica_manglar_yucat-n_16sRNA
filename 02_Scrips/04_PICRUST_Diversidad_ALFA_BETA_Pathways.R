#Vamos a hace run codigo para analiisis senceillo de los datos que obtuvimos de PICRUTS2 

library(gridExtra)
library(tidyverse)
library(readr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(compositions)
library(nortest)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)

####### Importamos datos###################################

#Primero importamos nuestros datos. Vamos a utilizar dos datos. Primero la salida de PICRUTs2 de KO a nivel KO y segundo a nivel L1-L2 

#Importamos los nivesles L1-L2
datosNiveles<- read.delim("./01_Datos/picrust2/pred_metagenome_unstrat_descrip_Niveles.csv", header = T)

#Importamos los metadatos 

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





####### Selección Niveles KEGG y normalizacion clr ###################


#agrupamos por niveles 

Datos_sumaL1<-aggregate(datosNiveles[c("DZ14","DZ15","DZ24","DZ25","PM14","PM15","SIS24","SIS25")], by =list(datosNiveles$Nivel1), FUN=sum)

Datos_sumaL2<-aggregate(datosNiveles[c("DZ14","DZ15","DZ24","DZ25","PM14","PM15","SIS24","SIS25")], by =list(datosNiveles$Nivel2), FUN=sum)

Datos_sumaL3<-datosNiveles[,c(3,5:12)]


#Sabemos que nuestros datos están muy sesgados a la izquierda por los ceros que contienen. Vamos a normalizarlo utilizando una clr 

Datos_sumaL1_CLR<-as.data.frame(clr(Datos_sumaL1[,c(2:9)]))
Datos_sumaL2_CLR<-as.data.frame(clr(Datos_sumaL2[,c(2:9)]))
Datos_sumaL3_CLR<-as.data.frame(clr(Datos_sumaL3[,c(2:9)]))

#Ahora unimos nuestro resultado de CLR con las columnas que tienen los datos 
Datos_sumaL1_CLR<-cbind(Datos_sumaL1[,1],Datos_sumaL1_CLR)
Datos_sumaL2_CLR<-cbind(Datos_sumaL2[,1],Datos_sumaL2_CLR)
Datos_sumaL3_CLR<-cbind(Datos_sumaL3[,1],Datos_sumaL3_CLR)

#El nombre de la primera columna es dificil de sustituir (Tinen encabezado con nombre no adecuado), hacemos un vector que contenga los nombres orignlaes de las columnas y lo sustituimos en los  df Datos_sumaL1_CLR que hicimos 
colna<-colnames(datosNiveles)

#Asignamos los nombres de columnas a los df que creamos 
colnames(Datos_sumaL1_CLR)<-colna[c(1,5:12)]
colnames(Datos_sumaL2_CLR)<-colna[c(2,5:12)]
colnames(Datos_sumaL3_CLR)<-colna[c(3,5:12)]

#Los hacemos una tabla larga 
L1_longer<-pivot_longer(Datos_sumaL1_CLR,cols = !Nivel1)
L2_longer<-pivot_longer(Datos_sumaL2_CLR,cols = !Nivel2)
L3_longer<-pivot_longer(Datos_sumaL3_CLR,cols= !pathway)
#Les agregamos los metadatos 
L1_longer<-left_join(L1_longer,metadatos)
L2_longer<-left_join(L2_longer,metadatos)
L3_longer<-left_join(L3_longer,metadatos)


#Vamos a generar un vector que tenga ordenado de este a oeste los nombres de los sitios de muestreo, esto nos va a ayudar al momento de graficar en orden geografico 
level_order <- c("DZ14","DZ15","DZ24","DZ25","SIS24","SIS25","PM14","PM15") 

####### Podemos Graficar ######
ggplot()+
  geom_tile(data=L1_longer,aes(x=factor(name,level_order),y=Nivel1,fill=value))+
  scale_fill_distiller(palette="RdYlBu")+
  ylab("Pathway")+
  xlab("Muestra")+
  labs(colour = "PP")+
  facet_wrap(~Region,scales = "free_x")+
  theme_get()+
  labs(fill='Valor CLR') 

ggplot()+
  geom_tile(data=L2_longer,aes(x=factor(name,level_order),y=Nivel2,fill=value))+
  scale_fill_distiller(palette="RdYlBu")+
  ylab("Pathway")+
  xlab("Muestra")+
  labs(colour = "PP")+
  theme_gray()+
  facet_wrap(~Region,scales = "free_x")+
  theme_get()+
  labs(fill='Valor CLR') 


ggplot()+
  geom_tile(data=L3_longer,aes(x=factor(name,level_order),y=pathway,fill=value))+
  scale_fill_distiller(palette="RdYlBu")+
  ylab("Pathway")+
  xlab("Muestra")+
  labs(colour = "PP")+
  theme_gray()+
  facet_wrap(~Region,scales = "free_x")+
  theme_get()+
  labs(fill='Valor CLR')+
  theme(axis.text.y = element_text(size = 4))  

##### Beta diversidad #######
####### Pruebas de aitchison ###########

#Vamos a calcular las matrices de aitchison para nuestros datos. 
#Como ya calculamos la CLR para todo el df original ahora sólo 
#es necesario utilizar una verdist euclideana para obtener resultados 

#Primero pasamos las columna que son descripciones a rows

rownames(Datos_sumaL1_CLR)<-Datos_sumaL1_CLR$Nivel1
Datos_sumaL1_CLR<-Datos_sumaL1_CLR[,c(2:9)]
Datos_sumaL1_CLR<-as.data.frame(t(Datos_sumaL1_CLR))

rownames(Datos_sumaL2_CLR)<-Datos_sumaL2_CLR$Nivel2
Datos_sumaL2_CLR<-Datos_sumaL2_CLR[,c(2:9)]
Datos_sumaL2_CLR<-as.data.frame(t(Datos_sumaL2_CLR))

rownames(Datos_sumaL3_CLR)<-Datos_sumaL3_CLR$pathway
Datos_sumaL3_CLR<-Datos_sumaL3_CLR[,c(2:9)]
Datos_sumaL3_CLR<-as.data.frame(t(Datos_sumaL3_CLR))


####Se calcula la matriz de bray#### 
bray<-vegdist(Datos_sumaL3_CLR, method = "euclidean")

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
permaR<-adonis2(bray~metadatos$Region, permutations = 9999)
p_valueR <- permaR$`Pr(>F)`[1]

#Permanova para estatus de conservacion 
permaC<-adonis2(bray~metadatos$Condicion, permutations = 9999)
p_valueC <- permaC$`Pr(>F)`[1]

#permanova para Sitio de muestreo 
permaS<-adonis2(bray~metadatos$Sitio, permutations = 9999)
p_valueS <- permaC$`Pr(>F)`[1]


#Se pueden ver si la dispersión de los puntos es significantes, lo cual nos puede ayudar a saber si los resultados del permanova son influenciados por la dispersión de los datos 

disperR<-betadisper(bray,metadatos$Region)
anova(disperR)
disperR<-betadisper(bray,metadatos$Condicion)
anova(disperR)

#summary(ASV_1) %>% view()


####### Graficamos nuestro NMDS y agregamos los valores obtenidos de la PERMANOVA#######

R<-ggplot(point, aes(x=MDS1,y=MDS2, color=Region))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR Euclidean Region Nivel 2", subtitle = paste("stress=", nm$stress), caption = paste("p-value",p_valueR))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()


C<-ggplot(point, aes(x=MDS1,y=MDS2, color=Condicion))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR Euclidean Condicion Nivel 2  ", subtitle = paste("stress=", nm$stress), caption = paste("p-value",p_valueC))+
  scale_color_manual(values=c("#fdb462", "#b3de69"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()

#Grafica
grid.arrange(C,R, ncol=2)






#### Ahora vamos a probar a todos los Nivel 2 con Wilcox #####

#Ahora hacemos las pruebas para todos 

L2Wc<-L2_longer%>%
  nest(data = -Nivel2) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Condicion, data=.x) 
                    %>% tidy))%>%unnest(test)
L2Wr<-L2_longer%>%
  nest(data = -Nivel2) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Region, data=.x) 
                    %>% tidy))%>%unnest(test)












#
#
##
#
#
####### Seleccion de Grupos de inteŕes"#######################

#Ahora vamos a tomar sólo aquellas rutas que pertenecen a dos categoria de interés la primera los Xenobiotics 

Xenos<-Datos_sumaL3_CLR[Datos_sumaL3_CLR$Nivel2=="Xenobiotics_biodegradation_and_metabolism",]%>%.[,c(3:11)]

Xenos_Longer<-pivot_longer(Xenos,cols = !description)%>%left_join(.,metadatos)

ggplot()+
  geom_tile(data=Xenos_Longer,aes(x=factor(name,level_order) ,y=description,fill=value))+
  scale_fill_distiller(palette="RdBu")+
  ylab("Pathway")+
  xlab("Sitio de muestreo")+
  labs(colour = "PP")+
  facet_wrap(~Region,scales = "free_x")

#Y ahora la ruta de energías ya que incluye CBQ 

CiclosBGQ<-Datos_sumaL3_CLR[Datos_sumaL3_CLR$Nivel2=="Energy_metabolism",]%>%.[,c(3:11)]
CiclosBGQ_longer<-pivot_longer(CiclosBGQ,cols = !description)%>%left_join(.,metadatos)

ggplot()+
  geom_tile(data=CiclosBGQ_longer,aes(x=factor(name,level_order) ,y=description,fill=value))+
  scale_fill_distiller(palette="RdBu")+
  ylab("Pathway")+
  xlab("Sitio de muestreo")+
  labs(colour = "PP")+
  facet_wrap(~Region,scales = "free_x")


#### Prueba para Xenobioticos con Wilcox #####

# se toma sólo la fila de ruta de degradación de xenobiotico del df 

XenoGral<-Datos_sumaL2_CLR[,31]

#vamos su normalidad
lillie.test(XenoGral)
hist(XenoGral)
#como p-value = 0.06098 entonces ES normal, pero el histograma se ve medio rato, entonces consideramos que no es normal 

XenoGral<-as.data.frame(XenoGral)
 
#Agregamos nombres de columnas 
rownames(XenoGral)<-rownames(Datos_sumaL2_CLR)


#Ahora vamos probar si hay diferencias en el valor de Xenobioticos entre regiones y condiciones

XenoGralRegion<-wilcox.test(XenoGral$XenoGral~metadatos$Region)
XenoGralRegion$p.valueX

XenoGralCond<-wilcox.test(XenoGral$XenoGral~metadatos$Condicion)
XenoGralCond$p.value

XenoGralSit<-pairwise.wilcox.test(XenoGral$XenoGral,metadatos$Condicion)
XenoGralSit$p.value

#No hay diferencias entre sitios  regiones y condiciones 

#### Ahora vamos a probar a todos los grupos#####

#tomamos la matriz de Xenos CLR y ponenmos nombre de rows

X<-Xenos
rownames(X)<-X%>%.[,1]
X<-X[,-1]
X<-as.data.frame(t(X))

#Ahora hacemos las pruebas para todos 

XXc<-Xenos_Longer%>%
  nest(data = -description) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Condicion, data=.x) 
                    %>% tidy))%>%unnest(test)

XXr<-Xenos_Longer%>%
  nest(data = -description) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Region, data=.x) 
                    %>% tidy))%>%unnest(test)

