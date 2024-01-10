#vamos a hacer un NMDS con la matriz de Unifraq Weight y Unweigh esto quiere decir que es un NMDs que toma en consideración la distancia phylogenitca que hay entre los ASV y lo sitios.Esto ayudara a tener la beta diversity 

#importamos las librrias importantes
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(gridExtra)

######## Importar datos #############

#Abrimos los archivos de interes. Estos los obtuvimos de qiim2 
WeightUniFrac<- read_qza("./01_Datos/qiime/weighted_unifrac_distance_matrix.qza")
UnWeightUniFrac<- read_qza("./01_Datos/qiime/unweighted_unifrac_distance_matrix.qza")

#importamos unos metadatos 

metadata <- read.csv("./01_Datos/sample-metadata.tsv", sep='\t')
metadata<-metadata[c(2,4,5,6,7,8,9),c(1:4)]

#LOs debemos ordenar para que se lean bien 

metadata <- metadata[order(metadata$sample_name), ]

#vamos a cambiar los metadastos de la siguientes forma 
#Este <- Dzilam-Bocas
#Oeste <- El Palmar 

#Este cambio se hace para evitar confusiones en el anañisis y descripción de datos ya que tambien hay sitios de muestreo con esos nombres 

#metadata$Region <- ifelse(metadata$Region == 'DZILAM-BOCAS', 'Este', 'Oeste')


######## Extracción matrices ##########

#De estas debemos extraer aquella llamada datos, que es la matriz de disimilitud

WeightUniFracM<- WeightUniFrac$data
UnWeightUniFracM<- UnWeightUniFrac$data

#Podemos remover los qza 
rm(UnWeightUniFrac)
rm(WeightUniFrac)

#Ahora que tenemos la matrices podemos hacerlos NMDS 
set.seed(12345)
NMDSWeightUniFrac<-metaMDS(WeightUniFracM)
NMDSUnWeightUniFrac<-metaMDS(UnWeightUniFracM)

#OJO, Aqui sale un mensaje que probablemente tenemos pocos datos. 

#ahora podemos hacer la extrancción de lo ejes de cada NMDS 
PointWeightUniFrac<-as.data.frame(NMDSWeightUniFrac$points)%>%rownames_to_column()
colnames(PointWeightUniFrac)<-c("sample_name","MDS1","MDS2")
PointUnWeightUniFrac<-as.data.frame(NMDSUnWeightUniFrac$points)%>%rownames_to_column()
colnames(PointUnWeightUniFrac)<-c("sample_name","MDS1","MDS2")

#Podemos unir los metadatos y los puntos 
PointWeightUniFrac<-left_join(PointWeightUniFrac,metadata)
PointUnWeightUniFrac<-left_join(PointUnWeightUniFrac,metadata)



######## PERMANOVA UNIFRAC #################

#Ahora podemos hacer las pruebas PERMANOVAS 


PERMAUnWeightUniFracMCondicion<-adonis2(UnWeightUniFracM~metadata$Condicion, permutations = 999,by = NULL )
p_vUnWeightCondicion <- PERMAUnWeightUniFracMCondicion$`Pr(>F)`[1]
p_vUnWeightCondicion

PERMAUnWeightUniFracMRegion<-adonis2(UnWeightUniFracM~metadata$Region, permutations = 999,by = NULL )
p_vUnWeightRegion <- PERMAUnWeightUniFracMRegion$`Pr(>F)`[1]
p_vUnWeightRegion


PERMAWeightUniFracMCondicion<-adonis2(WeightUniFracM~metadata$Condicion, permutations = 999,by = NULL )
p_vWeightCondicion <- PERMAWeightUniFracMCondicion$`Pr(>F)`[1]
p_vWeightCondicion


PERMAWeightUniFracMRegion<-adonis2(WeightUniFracM~metadata$Region, permutations = 999,by = NULL )
p_vWeightRegion <- PERMAWeightUniFracMRegion$`Pr(>F)`[1]
p_vWeightRegion


PERMAWeightUniFracMSitio<-adonis2(WeightUniFracM~metadata$Sitio, permutations = 999,by = NULL )
p_vWeightSitio <- PERMAWeightUniFracMSitio$`Pr(>F)`[1]
p_vWeightSitio

PERMAUnWeightUniFracMSitio<-adonis2(UnWeightUniFracM~metadata$Sitio, permutations = 999,by = NULL )
p_vUnWeightSitio <- PERMAUnWeightUniFracMSitio$`Pr(>F)`[1]
p_vUnWeightSitio


######## Grafica ##################


#Ahora podemos graficar, con los valores de stress y el p values por NDMS 

rW<-ggplot(PointWeightUniFrac,aes(x=MDS1,y=MDS2, color=Region, shape=Sitio))+
  geom_point(size=5)+
  theme_minimal()+
  ggtitle("NMDS WeightUniFrac Regiones de muestreo ", subtitle = paste0("stress",NMDSWeightUniFrac$stress))+
  labs(caption = paste("p-value",p_vWeightRegion))
  

rUw<-ggplot(PointUnWeightUniFrac,aes(x=MDS1,y=MDS2, color=Region, shape=Sitio))+
  geom_point(size=5)+
  theme_minimal()+
  ggtitle("NMDS UnWeightUniFrac Regiones de muestreo", subtitle = paste0("stress",NMDSWeightUniFrac$stress))+
  labs(caption = paste("p-value",p_vUnWeightRegion))


cW<-ggplot(PointWeightUniFrac,aes(x=MDS1,y=MDS2, color=Condicion,shape=Sitio))+
  geom_point(size=5)+
  theme_minimal()+
  ggtitle("NMDS WeightUniFrac Estatus de conservacion", subtitle = paste0("stress",NMDSUnWeightUniFrac$stress))+
  labs(caption = paste("p-value",p_vWeightCondicion))+
  scale_color_manual(values=c("#fdb462", "#b3de69"))

cUw<-ggplot(PointUnWeightUniFrac,aes(x=MDS1,y=MDS2, color=Condicion,shape=Sitio))+
  geom_point(size=5)+
  theme_minimal()+
  ggtitle("NMDS UnWeightUniFrac Estatus de conservación", subtitle = paste0("stress",NMDSUnWeightUniFrac$stress))+
  labs(caption = paste("p-value",p_vUnWeightCondicion))+
  scale_color_manual(values=c("#fdb462", "#b3de69"))
  

grid.arrange(rW,rUw,cW,cUw,nrow=2)



