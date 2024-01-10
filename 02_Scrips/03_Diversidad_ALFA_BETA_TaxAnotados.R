#importar los dato de taxonomia 

#Librerias necearias 
library(tidyverse)
library(nortest)
library(vegan)
library(gridExtra)
  

#####  Importamos datos ####
taxonomy<-as.data.frame(read.csv("./01_Datos/qiime/level-2.csv"))
  
#preparamos los datos. La tabla eliminando filas que no son de nuestro interes 
row.names(taxonomy)<-taxonomy[,1]
taxonomy<-taxonomy[,c(2:length(taxonomy))]
str(taxonomy)
  
#importamos los metadatos. Estoa deben ser los utilizados en Qiime2
metadata<-read.csv("./01_Datos/sample-metadata.tsv", sep = '\t')
#Aqui corregimos los metadatos porque decartamos unos valores de fila 
metadata <- metadata[c(2,3,4,5,6,7,8,9),]
#Los debemos ordenar en el mismo orden en el que están los sitios de la tabla de taxonomia para evitar errores
metadata <- metadata[order(metadata$sample_name), ]


#vamos a cambiar los metadastos de la siguientes forma 
#Este <- Dzilam-Bocas
#Oeste <- El Palmar 

#Este cambio se hace para evitar confusiones en el anañisis y descripción de datos ya que tambien hay sitios de muestreo con esos nombres 

metadata$Region <- ifelse(metadata$Region == 'DZILAM-BOCAS', 'Este', 'Oeste')


#####  Calcular los indices de diversidad####

#Shannon;How difﬁcult would it be to predict correctly the species of the next individual collected?
#Chao1;Index is based upon the number of rare classes (i.e., OTUs) found in a sample
#Simpson;Probability that two individuals picked at random belong to the same species
#Evenness


#Calculamos indice de shannon 
shannon<-diversity(taxonomy, index="shannon")
  
#Calculamos Chao1
chao1<-estimateR(taxonomy)[2,]
chao1<-as.data.frame(chao1)%>%rownames_to_column()
  
#Calculamos indice de GINI-simpson 
#En este caso se calcula 1-D 
simpson<-diversity(taxonomy, index="simpson")


simpson<-as.data.frame(simpson)%>%rownames_to_column()
rm(simpeson)
  
#Ahora sacamos evenness #
spec <-as.data.frame(log(specnumber(taxonomy)))
str(spec)
#Usamos la formula  para hacer el calculo de evenness
evenness<- shannon/spec$`log(specnumber(taxonomy))`
evenness<-as.data.frame(evenness)%>%rownames_to_column()


  
#Ahora ponemos los cuatro indices juntos en un sólo dataframe  
  
#primero pasamos shanon a df 
shannon<-as.data.frame(shannon)%>%rownames_to_column()
  
#Unimos en un df
Alfa<-left_join(shannon,evenness)%>%left_join(.,chao1)%>%left_join(.,simpson)
Alfa<-dplyr::rename(Alfa,sample_name=rowname)%>%left_join(.,metadata)
  
#Podemos guardar valores de alfa diversidad. 
#write_csv(Alfa,"/home/laac/Desktop/AlfaLORden.csv")
  
  
#




#####  Prueba  resultados de alfa diversidad. #####
  
#Shapiro test p>0.05 es normal. Si p<0.05 no normal.
#lillie test p>0.05 es normal. Si p<0.05 no es normal. 
  
lillie.test(Alfa$shannon)
#Es normal 

lillie.test(Alfa$simpson)
#No normal

lillie.test(Alfa$evenness)
#Es normal 

lillie.test(Alfa$chao1)
#Es normal  


  
#Los resultados fueron normales entonces podemos hacer t.test y para probar las diferencias en los indices de diversidad y un wilcox test para los datos no normales. 
  


#####  Hacemos los test t test y wilcox test para diversidad alfa ######

#primero probamos shannon por region 
ShaTTestRegion<-pairwise.t.test(Alfa$shannon,Alfa$Region)
ShaTTestRegion<- ShaTTestRegion$p.value
ShaTTestRegion<-round(ShaTTestRegion[1,1],digits = 3)
  
#por condicion shannon 
ShaTTestCondicion<-pairwise.t.test(Alfa$shannon,Alfa$Condicion)
ShaTTestCondion<- ShaTTestCondicion$p.value
ShaTTestCondion<-round(ShaTTestCondion[1,1],digits = 3)
  
ShaTTestSitio<-pairwise.t.test(Alfa$shannon,Alfa$Sitio)
ShaTTestSitio
  
#Ahora probamos indice de Simpson 
SimpWilcoxCondicion<-wilcox.test(Alfa$simpson~Alfa$Condicion)
SimpWilcoxCondion<- SimpWilcoxCondicion$p.value
SimpWilcoxCondion
  
SimpWilcoxRegion<-wilcox.test(Alfa$simpson~Alfa$Region)
SimpWilcoxRegion<- SimpWilcoxRegion$p.value
SimpWilcoxRegion
  
SimpWilcoxSitio<-pairwise.wilcox.test(Alfa$simpson,Alfa$Sitio)
SimpWilcoxSitio
  

#evenness
EvenTTesCondicion<-pairwise.t.test(Alfa$evenness,Alfa$Condicion)
EvenTTesCondion<- EvenTTesCondicion$p.value
EvenTTesCondion<-round(EvenTTesCondion[1,1],digits = 3)
  
EvenTTesRegion<-pairwise.t.test(Alfa$evenness,Alfa$Region)
EvenTTesRegion<- EvenTTesRegion$p.value
EvenTTesRegion<-round(EvenTTesRegion[1,1],digits = 3)
  
EvenTTesSitio<-pairwise.t.test(Alfa$evenness,Alfa$Sitio)
EvenTTesSitio

#Ahora Chao1   

ChaoTTestCondicion<-pairwise.t.test(Alfa$chao1,Alfa$Condicion)
ChaoTTestCondicion<-ChaoTTestCondicion$p.value
ChaoTTestCondicion<-round(ChaoTTestCondicion,digits = 3)
  
ChaoTTestRegion<-pairwise.t.test(Alfa$chao1,Alfa$Region)
ChaoTTestRegion<-ChaoTTestRegion$p.value
ChaoTTestRegion<-round(ChaoTTestRegion, digits = 3)
  




#####  Grafica alfa ################

s2<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Condicion,y=shannon, fill=Condicion))+
  ggtitle("Indice de Shannon")+# adjust y-axis
  xlab("Estatus de conservación") +
  ylab("Indice de Shanon") +
  theme_light() +
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#fdb462","#b3de69"))+
  labs(caption = paste("p-value Condicion ",ShaTTestCondion))+ 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

t2<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Condicion,y=evenness, fill=Condicion))+
  scale_fill_manual(values=c("#fdb462","#b3de69")) + 
  ggtitle("Indice de Evenness")+# adjust y-axis
  xlab("Estatus de conservación") +
  ylab("Indice de Evennes") +
  theme_light() +
  labs(caption = paste("p-value Condicion ",EvenTTesCondion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))


r2<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Condicion,y=chao1, fill=Condicion))+
  scale_fill_manual(values=c("#fdb462","#b3de69")) + 
  ggtitle("Indice de Chao1")+# adjust y-axis
  xlab("Estatus de conservación") +
  ylab("Indice de Chao1") +
  theme_light() +
  labs(caption = paste("p-value Condicion ",ChaoTTestCondicion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

c2<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Condicion,y=simpson, fill=Condicion))+
  scale_fill_manual(values=c("#fdb462","#b3de69")) + 
  ggtitle("Indice de Simpson")+# adjust y-axis
  xlab("Estatus de conservación") +
  ylab("Indice de Simpson") +
  theme_light() +
  labs(caption = paste("p-value Condicion ",SimpWilcoxCondion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))


s3<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Region,y=shannon, fill=Region))+
  xlab("Region") +
  ylab("Indice de Shanon") +
  theme_light() +
  theme(legend.position = "none")+
  labs(caption = paste("p-value Region ",ShaTTestRegion))+ 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

t3<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Region,y=evenness, fill=Region))+
  xlab("Region") +
  ylab("Indice de Evennes") +
  theme_light() +
  labs(caption = paste("p-value Region ",EvenTTesRegion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))


r3<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Region,y=chao1, fill=Region))+
  xlab("Region") +
  ylab("Indice de Chao1") +
  theme_light() +
  labs(caption = paste("p-value Region ",ChaoTTestRegion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

c3<-ggplot()+
  geom_boxplot(data=Alfa,aes(x=Region,y=simpson, fill=Region))+
  xlab("Region") +
  ylab("Indice de Simpson") +
  theme_light() +
  labs(caption = paste("p-value Region ",SimpWilcoxRegion))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
        axis.text = element_text(angle = 0),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

grid.arrange(s2,r2,t2,c2,s3,r3,t3,c3,ncol=4)








#####  Calculamos  Beta diversidad ########################


#Agregamos un PsudoCount 
taxonomy<-taxonomy+0.01
  
#Se calcula la matriz de bray 
bray<-vegdist(taxonomy, method = "robust.aitchison")
  
#Se realiza un NMDS 
set.seed(1234)
nm <- metaMDS(bray)
nm
  
#Se extraen los puntos del NMDS par graficar 
point<- as.data.frame(nm$points)
point
point <-rownames_to_column(point) 
  
#Cambiamos nombres de las columnas de df de puntos
colnames(point)<- c("sample_name", "MDS1","MDS2")
  
point<- merge(point,metadata)
  
#Vamos a realizar un PERMANOVA con Adonis para encontrar diferencias entre los sitios de muestreo, condiciones de conservación y region de muestreo para nuesta matriz 
  
#permanova para Region
permaR<-adonis2(bray~metadata$Region, permutations = 9999)
p_valueR <- permaR$`Pr(>F)`[1]
  
#Permanova para estatus de conservacion 
permaC<-adonis2(bray~metadata$Condicion, permutations = 9999)
p_valueC <- permaC$`Pr(>F)`[1]
  
#permanova para Sitio de muestreo 
permaS<-adonis2(bray~metadata$Sitio, permutations = 9999)
p_valueS <- permaC$`Pr(>F)`[1]
  
  
#Se pueden ver si la dispersión de los puntos es significantes, lo cual nos puede ayudar a saber si los resultados del permanova son influenciados por la dispersión de los datos 
  
disperR<-betadisper(bray,metadata$Region)
anova(disperR)
disperR<-betadisper(bray,metadata$Condicion)
anova(disperR)
  
#summary(ASV_1) %>% view()
  
#Graficamos nuestro NMDS y agregamos los valores obtenidos de la PERMANOVA
  
R<-ggplot(point, aes(x=MDS1,y=MDS2, color=Region))+
    geom_point(size=6)+
    theme_minimal()+
    labs(title="NMDS Robust Aitchison Region", subtitle = paste("stress=", nm$stress), caption = paste("p-value",p_valueR))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_light()+
    stat_ellipse()
  
  
C<-ggplot(point, aes(x=MDS1,y=MDS2, color=Condicion))+
    geom_point(size=6)+
    theme_minimal()+
    labs(title="NMDS Robust Aitchison Condicion  ", subtitle = paste("stress=", nm$stress), caption = paste("p-value",p_valueC))+
    scale_color_manual(values=c("#fdb462", "#b3de69"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_light()+
    stat_ellipse()
  
#Grafica
grid.arrange(C,R, ncol=2)





