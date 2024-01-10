#Vamos a calcular al Diversidad alfa de ASV utilizando el vegan y tomando como referencia los indices calculados en qiime 

library(qiime2R)
library(tidyverse)
library(gridExtra)
library(nortest)

######## Importamos las ASVs y los metadatos de QIIME2######

SVs<-read_qza("./01_Datos/qiime/table.qza")
metadata<-read.csv("./01_Datos/sample-metadata.tsv", sep = '\t')

#Aqui corregimos los metadatos porque decartamos unos valor al hacer la rarefacción 
metadata <- metadata[c(2,4,5,6,7,8,9),]

#LOs debemos ordenar para que se lean bien 

metadata <- metadata[order(metadata$sample_name), ]

#vamos a cambiar los metadastos de la siguientes forma 
#Este <- Dzilam-Bocas
#Oeste <- El Palmar 

#Este cambio se hace para evitar confusiones en el anañisis y descripción de datos ya que tambien hay sitios de muestreo con esos nombres 

metadata$Region <- ifelse(metadata$Region == 'DZILAM-BOCAS', 'Este', 'Oeste')





taxonomy<-read_qza("./01_Datos/qiime/taxonomyDM6.qza")
head(taxonomy$data)

taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)


######## Importamos los  indices de diversidad de QIIME2#####

shannon<-read_qza("./01_Datos/qiime/shannon_vector.qza")
shannon_data<-shannon$data  %>% rownames_to_column("sample_name")
eveness<- read_qza("./01_Datos/qiime/evenness_vector.qza")
eveness_data <- eveness$data  %>% rownames_to_column("sample_name")
Chao1<- read_qza("./01_Datos/qiime/Chao1.qza")
Chao1_data <- Chao1$data %>% rownames_to_column("sample_name")

simpson<-read_qza("./01_Datos/qiime/simpson.qza")
simpson<-simpson$data%>% rownames_to_column("sample_name")

#Cuando se hizo la rarefaccion se perdio un muestra, por esto se elimina una muestra SIS24, en este caso se debe remover 
Chao1_data <- Chao1_data[c(1,2,3,4,5,6,8),]
simpson <- simpson[c(1,2,3,4,5,6,8),]

#Agrupamos los resultados en un df
DiversidadAlfaASV<-
  metadata %>% 
  left_join(shannon_data) %>% left_join(eveness_data) %>% left_join(Chao1_data)%>%left_join(simpson)

 

######## Primero vamos a revisar la normalidad de diversidad alfa #####


#Shapiro test p>0.05 es normal. Si p<0.05 no normal.
#lillie test p>0.05 es normal. Si p<0.05 no es normal.

shapiroShannonASV<-lillie.test(DiversidadAlfaASV$shannon_entropy)
shapiroShannonASV<-shapiroShannonASV$p.value
#hist(DiversidadAlfaASV$shannon_entropy)

shapiroEvennessASV<-lillie.test(DiversidadAlfaASV$pielou_evenness)
shapiroEvennessASV<-shapiroEvennessASV$p.value
#hist(DiversidadAlfaASV$pielou_evenness)

shapiroChao1ASV<-lillie.test(DiversidadAlfaASV$chao1)
shapiroChao1ASV<-shapiroChao1ASV$p.value
#hist(DiversidadAlfaASV$chao1)

shapiroSimpASV<-lillie.test(DiversidadAlfaASV$simpson)
shapiroSimpASV<-shapiroSimpASV$p.value

#Todos los resultados fueron para datos normales 



######## Pruebas alfa diversidad ######

#En niniguno de los casos anteriores se obtuvieron p>0.05 por lo tanto se pueden manejar como datos normales. Se procede a realizar un ttest para todos los casos 


#Primero probamos las condiciones por Shannon por estatus de conservacion
ttestShanonCondicionesASV <-t.test(DiversidadAlfaASV$shannon_entropy~DiversidadAlfaASV$Condicion)
PttestShanonCondicionesASV<-ttestShanonCondicionesASV$p.value
PttestShanonCondicionesASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice se Shannon entre estatus de conservacion

ttestShanonRegionASV <-t.test(DiversidadAlfaASV$shannon_entropy~DiversidadAlfaASV$Region)
PttestShanonRegionASV<-ttestShanonRegionASV$p.value
PttestShanonRegionASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice se Shannon entre regiones


#Ahora Evenness
ttestEvennessCondicionesASV <-t.test(DiversidadAlfaASV$pielou_evenness~DiversidadAlfaASV$Condicion)
PttestEvennessCondicionesASV<-ttestEvennessCondicionesASV$p.value
PttestEvennessCondicionesASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice se Evenness entre estatus de conservacion

ttestEvennessRegionASV <-t.test(DiversidadAlfaASV$pielou_evenness~DiversidadAlfaASV$Region)
PttestEvennessRegionASV<-ttestEvennessRegionASV$p.value
PttestEvennessRegionASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice se evenness entre regiones

#Ahora Chao1
ttestChao1CondicionesASV <-t.test(DiversidadAlfaASV$chao1~DiversidadAlfaASV$Condicion)
PttestChao1CondicionesASV<-ttestChao1CondicionesASV$p.value
PttestChao1CondicionesASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice se Chao1 entre estatus de conservacion

ttestChao1RegionASV <-t.test(DiversidadAlfaASV$chao1~DiversidadAlfaASV$Region)
PttestChao1RegionASV<-ttestChao1RegionASV$p.value
PttestChao1RegionASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice de Chao1 entre regiones


#Ahora Simpson
ttestSimpsonCondicionesASV <-t.test(DiversidadAlfaASV$simpson~DiversidadAlfaASV$Condicion)
PttestSimpsonCondicionesASV<-ttestSimpsonCondicionesASV$p.value
PttestSimpsonCondicionesASV
 
#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice de Simpson entre condiciones

#Ahora Simpson
ttestSimpsonRegionASV <-t.test(DiversidadAlfaASV$simpson~DiversidadAlfaASV$Region)
PttestSimpsonRegionASV<-ttestSimpsonRegionASV$p.value
PttestSimpsonRegionASV

#el valor parece ser mayor a 0.05, no es significativo. Inidica que no hay diferencias en el indice de Simpson entre Regiones




######## Media indices alfa por Condicion Region y Sitio####

AlfaCondicion<-DiversidadAlfaASV%>%group_by(Condicion) %>% 
  summarise_at(vars("shannon_entropy","pielou_evenness","chao1","simpson"), median)

AlfaRegion<-DiversidadAlfaASV%>%group_by(Region) %>% 
  summarise_at(vars("shannon_entropy","pielou_evenness","chao1","simpson"), median)

AlfaSitio<-DiversidadAlfaASV%>%group_by(Sitio) %>% 
  summarise_at(vars("shannon_entropy","pielou_evenness","chao1","simpson"), median)

AlfaSitio<-DiversidadAlfaASV%>%group_by() %>% 
  summarise_at(vars("shannon_entropy","pielou_evenness","chao1","simpson"), median)


######## Graficamos los datos alfa diversidad por Sitio #################

s<-ggplot(DiversidadAlfaASV, aes(x=`Condicion`, y=shannon_entropy, fill=`Condicion`)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0,9)) + # adjust y-axis
  facet_grid(~`Region`) + # create a panel for each body site
  xlab("Estatus de conservación") +
  ylab("Shannon Diversity") +
  theme_minimal() +
  ggtitle("Shannon")+
  scale_fill_manual(values=c("#fdb462", "#b3de69")) + #specify custom colors
  theme(legend.position="none")+ #remove the legend as it isn't needed
  theme(plot.title = element_text(hjust = 0.5))


e<-ggplot(DiversidadAlfaASV, aes(x=`Condicion`, y=pielou_evenness, fill=`Condicion`)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0,1)) + # adjust y-axis
  facet_grid(~`Region`) + # create a panel for each body site
  xlab("Estatus de conservación") +
  ylab("pielou_evenness") +
  theme_minimal() +
  ggtitle("Evenness")+
  scale_fill_manual(values=c("#fdb462", "#b3de69")) + #specify custom colors
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))

c <- ggplot(DiversidadAlfaASV, aes(x=`Condicion`, y=chao1, fill=`Condicion`)) +
  geom_boxplot() + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(0,800)) + # adjust y-axis
  facet_grid(~`Region`) + # create a panel for each body site
  xlab("Estatus de conservación") +
  ylab("Chao1") +
  theme_minimal() +
  ggtitle("Chao1")+
  scale_fill_manual(values=c("#fdb462", "#b3de69")) + #specify custom colors
  theme(legend.position="none")+ #remove the legend as it isn't needed
  theme(plot.title = element_text(hjust = 0.5))


grid.arrange(s,e,c, ncol=3)




######## Grafica datos alfa por Condición y Region ##########

sha<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Region,y=shannon_entropy, fill=Region))+
  ggtitle("Indice de Shannon")+
  ylab("Shannon")+
  labs( caption = paste("p-value region ", round(PttestShanonRegionASV,2)))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",plot.caption = element_text(hjust = 0.5))

pie<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Region,y=pielou_evenness,fill=Region))+
  ggtitle("Indice de Pielou Evenness")+
  ylab("Pielou Evenness")+
  labs(caption=paste("p-value region ", round(PttestEvennessRegionASV,4)))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",plot.caption = element_text(hjust = 0.5))


ch<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Region,y=chao1,fill=Region))+
  ggtitle("Indice Chao1")+
  ylab("Chao1")+
  labs(caption=paste("p-value region ", round(PttestChao1RegionASV,4)))+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",plot.caption = element_text(hjust = 0.5))

  
simp<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Region,y=simpson,fill=Region))+
  ggtitle("Indice de Simpson")+
  ylab("Simpson")+
  labs(caption=paste("p-value region ", round(PttestSimpsonRegionASV,4)))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",plot.caption = element_text(hjust = 0.5))



sha2<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Condicion,y=shannon_entropy,fill=Condicion))+
  labs( caption = paste("p-value condicion ", round(PttestShanonCondicionesASV,2)))+
  ylab("Shannon")+
  theme_minimal()+
  scale_fill_manual(values=c("#fdb462","#b3de69"))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),plot.caption = element_text(hjust = 0.5))

pie2<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Condicion,y=pielou_evenness, fill=Condicion))+
  labs(caption=paste("p-value condicion", round(PttestEvennessCondicionesASV,4)))+
  ylab("Pielou Evenness")+
  theme_minimal()+
  scale_fill_manual(values=c("#fdb462","#b3de69"))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),plot.caption = element_text(hjust = 0.5))

ch2<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Condicion,y=chao1,fill=Condicion))+
  labs(caption=paste("p-value condicion ", round(PttestChao1CondicionesASV,4)))+
  ylab("Choa1")+
  theme_minimal()+
  scale_fill_manual(values=c("#fdb462","#b3de69"))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),plot.caption = element_text(hjust = 0.5))

simp2<-ggplot()+
  geom_boxplot(data=DiversidadAlfaASV,aes(x=Condicion,y=simpson,fill=Condicion))+
  labs(caption=paste("p-value condicion", round(PttestSimpsonCondicionesASV,4)))+
  ylab("Shannon")+
  theme_minimal()+
  scale_fill_manual(values=c("#fdb462","#b3de69"))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),plot.caption = element_text(hjust = 0.5))

grid.arrange(sha,simp,pie,ch,sha2,simp2,pie2,ch2,ncol=4)








