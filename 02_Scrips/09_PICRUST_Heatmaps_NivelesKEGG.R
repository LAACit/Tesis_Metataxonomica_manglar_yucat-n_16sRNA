library(tidyverse)
library(ggplot2)
library(viridis)
library(compositions)
library(reshape)
library(pheatmap)

datos<-read.table("./01_Datos/picrust2/pred_metagenome_unstrat_descrip_Niveles.csv", header = T)

#Ahora vamos a eliminar los que son NoID
#datos=datos[datos$Nivel2!="NoID",]

####procedemos a hacer un agregate para nivel #####
nivel1<- aggregate(datos[c("DZ14","DZ15","DZ24","DZ25","SIS24","SIS25","PM14","PM15")], by =list(datos$Nivel1), FUN=sum)

row.names(nivel1)<-nivel1$Group.1
nivel1<-nivel1[,-1]
nivel1clr<-as.data.frame(clr(nivel1))
nivel1clr<-arrange(nivel1clr, -nivel1clr$SIS25)
nivel1clr$description<-row.names(nivel1clr)

nivel1clrL<-melt(nivel1clr, id.var="description",
             variable.names= "variable",
             value.name = "abundancia",
             variable.factor = FALSE)

ggplot(nivel1clrL,aes(x=variable,fill=value,y=description))+
  geom_tile()+
  scale_fill_distiller(palette="RdBu")

#### nivel 2 #############


nivel2<- aggregate(datos[c("DZ14","DZ15","DZ24","DZ25","SIS24","SIS25","PM14","PM15")], by =list(datos$Nivel2), FUN=sum)


row.names(nivel2)<-nivel2$Group.1
nivel2<-nivel2[,-1]
nivel2clr<-as.data.frame(clr(nivel2))
nivel2clr<-arrange(nivel2clr, -nivel2clr$SIS25)
nivel2clr$description<-row.names(nivel2clr)

nivel2clrL<-melt(nivel2clr, id.var="description",
                 variable.names= "variable",
                 value.name = "abundancia",
                 variable.factor = FALSE)

ggplot(nivel2clrL,aes(x=variable,fill=value,y=description))+
  geom_tile(color = "gray",
            lwd = 0.1,
            linetype = 1)+
  scale_fill_viridis_c(option = "plasma")+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("muestra")+
  scale_fill_continuous(name = "Valor CLR")



#####Ahora sólo Xenos###### 

#Hacemos el CLR de todos los datos 
datosCLR<-as.data.frame(clr(datos[,5:12]))
columna<-datos[,c(1,2,3,4)]
datosCLR<-cbind(columna,datosCLR)

xenos<-datosCLR[datosCLR$Nivel2=="Xenobiotics_biodegradation_and_metabolism",]
xenos<-xenos[,-c(1,2,3)]

xenosL<-melt(xenos, id.var="description",
             variable.names= "variable",
             value.name = "abundancia",
             variable.factor = FALSE)

ggplot(xenosL,aes(x=variable,fill=value,y=description))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle = 90))

