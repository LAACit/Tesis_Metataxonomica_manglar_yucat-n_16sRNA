#L1<-Datos_sumaL2[,-c(1)]
#T1<-Datos_sumaL2[,c(1)]


#HEATMAP de KOS que resultan del ALDEX2

L1<-Datos
T1<-row.names(Datos)

#sumanos los datos
SumTax<-rowSums(L1)
#sumamos el total de valor de todos los datos
Total<-sum(L1)
Total <- sum(L1, na.rm = TRUE)
#sumamos para obtener total por sitio de muestreo 
SumSitio<-colSums(L1,na.rm = TRUE)


###### Calcula el Porcetaje por taxon de forma  general  #############
#Ahora calculamos el porcentaje por taxon 
SumLPorcentaje<-c()
for(i in 1:length(SumTax)){
  print(SumTax[i])
  SumLPorcentaje[i]<-(SumTax[i]*100)/Total
}


#PAra los xenos 
umLPorcentaje<-data.frame(valor=c(SumLPorcentaje), descripción=c(T1))
u<-cbind(umLPorcentaje,decrip)

#HAcemos un left

u2<-left_join(SiginificantesRBH,u)


#Heat map diferenciales 


#Importamos los metadatos 

metadatos <- read.csv("./01_Datos/sample-metadata.tsv", sep = '\t' )
metadatos<- metadatos[2:9,]
#LOs debemos ordenar para que se lean bien 
metadatos <- metadatos[order(metadatos$sample_name), ]
metadatos<-metadatos[,c(1:4)]
colnames(metadatos)<-c("Sitio","Condicion","Localidad","Region")

#Selecionamos los datos para obtener los nombres 
Datos$function.<-rownames(Datos)
#Unimos los datos de ALDEX2 y Datos 
u3<-left_join(SiginificantesRBH,Datos)
u3<-left_join(u3,decrip)
u3<-u3[,c(12:20)]

u4<-pivot_longer(data = u3,
                 cols = !function.,
                 values_to ="Valor_Counts",
                 names_to = "Sitio")

u5<-left_join(u4,metadatos)

level_order <- c("DZ14","DZ15","DZ24","DZ25","SIS24","SIS25","PM14","PM15") 

ggplot()+
  geom_tile(data=u5,aes(x=factor(Sitio,level_order),y=function.,fill=Valor_Counts))+
  scale_fill_distiller(palette="RdYlBu")+
  ylab("KO")+
  xlab("Muestra")+
  labs(colour = "PP")+
  facet_wrap(~Region,scales = "free_x")+
  theme_get()+ 
  guides(count = guide_legend(title = "Valor counts"))

write.csv(SiginificantesRBH,"./03_Resultados/ResultadosSignificantesRegionalKOSALDEX2.csv")

