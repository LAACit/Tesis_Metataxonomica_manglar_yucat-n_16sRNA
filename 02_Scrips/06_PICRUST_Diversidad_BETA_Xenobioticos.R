#Analisis de diversidad para predicción de rutas de biodegradaciónde Xenobioticos 

library(tidyverse)
library(vegan)
library(compositions)
library(nortest)
library(tibble)
library(broom)
library(gridExtra)
library(viridis)

#### Importamos datos y metadatos #####
#Abrimos los datos de PICRUTs2

datos<-read.table("./01_Datos/picrust2/pred_metagenome_unstrat_descrip_Niveles.csv", header = T)

#Selecionamos los que son pathways de Xenobiotics
xenos<-datos[datos$Nivel2=="Xenobiotics_biodegradation_and_metabolism",]
xenos<-xenos[,-c(1,2,3)]

datos<-xenos
rm(xenos)

#Hace transformación CLR 
datos[,c(2:9)]<-clr(datos[,c(2:9)])


#Abrimos metadatos
metadatos <- read.csv("./01_Datos/sample-metadata.tsv", sep = '\t' )
metadatos<- metadatos[2:9,]
#Losordenar para que se lean bien 
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
#### Tabla larga datos Xenos  #######################
Xenos_Longer<-pivot_longer(datos,cols = !description)%>%left_join(.,metadatos)

#### Grafica de heatmap para Xenos  ############


level_order <- c("DZ14","DZ15","DZ24","DZ25","SIS24","SIS25","PM14","PM15") 


ggplot()+
  geom_tile(data=Xenos_Longer,aes(x=factor(name,level_order),y=description,fill=value))+
  scale_fill_distiller(palette="RdBu",name="Valor CLR")+
  ylab("Pathway")+
  xlab("Muestra")+
  labs(colour = "PP")+
  facet_wrap(~Region,scales = "free_x")+
  theme_get()



#### Prueba para Xenobioticos con Wilcox #####

totalXeno<-as.data.frame(colSums(datos[,c(2:9)]))
colnames(totalXeno)<-c("total_xeno")


XenoGralRegion<-wilcox.test(totalXeno$total_xeno~metadatos$Region)
XenoGralRegion$p.valueX

XenoGralCond<-wilcox.test(totalXeno$total_xeno~metadatos$Condicion)
XenoGralCond$p.value

XenoGralSit<-pairwise.wilcox.test(totalXeno$total_xeno,metadatos$Sitio)
XenoGralSit$p.value

#No hay diferencias entre sitios  regiones y condiciones 


#### Ahora vamos a probar a todos los grupos de Xenos #####

#Ahora hacemos las pruebas para todos 

XXc<-Xenos_Longer%>%
  nest(data = -description) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Condicion, data=.x) 
                    %>% tidy))%>%unnest(test)
XXr<-Xenos_Longer%>%
  nest(data = -description) %>%
  mutate(test = map(.x=data, ~wilcox.test(value~Region, data=.x) 
                    %>% tidy))%>%unnest(test)


#### Hacemos un NMDS para los datos ######

rownames(datos)<-datos[,1]                   
datos<-datos[,-1]
datos<-as.data.frame(t(datos))

#Se calcula la matriz de bray 
bray<-vegdist(datos,method = "euclidean")

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


#### Graficamos nuestro NMDS y agregamos los valores obtenidos de la PERMANOVA#######

R<-ggplot(point, aes(x=MDS1,y=MDS2, color=Region))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR Xenobiotics", subtitle = paste("stress=", nm$stress), caption = paste("p-value=",p_valueR))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()


C<-ggplot(point, aes(x=MDS1,y=MDS2, color=Condicion))+
  geom_point(size=6)+
  theme_minimal()+
  labs(title="NMDS CLR Xenobiotics  ", subtitle = paste("stress=", nm$stress), caption = paste("p-value=",p_valueC))+
  scale_color_manual(values=c("#fdb462", "#b3de69"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse()

#Grafica
grid.arrange(C,R, ncol=2)


