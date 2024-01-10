#Importamos las librerias necesarias 

library(tidyverse)
library(ALDEx2)

#OJO SE DEBE ESCOJER QUE IMPORTAR AL COMIENZO 

##### Importamos para todas las rutas ####
datos<-read.table("./01_Datos/picrust2/pred_metagenome_unstrat_descrip_Niveles.csv", header = T)

datos<-datos[,-c(1,2,4)]
rows<-datos[,1]
datos<-datos[,-1]
datos<- as.data.frame(lapply(datos, function(col) as.integer(col)))
rownames(datos)<-rows

### Guardamos las descripciones y los KO para identificar posteriormete los resultados significantes 

selex.sub<-datos
#Verificamos que los datos sean int. Esto es necesario par ausar ALDEx2
str(selex.sub)


##### Importamos datos de Rutas de Xenobioticos #### 
datos<-read.table("./01_Datos/picrust2/pred_metagenome_unstrat_descrip_Niveles.csv", header = T)

#Ahora vamos a eliminar los que son NoID
datos=datos[datos$Nivel2!="NoID",]

#Selecionamos los que son pathways de Xenobiotics
datos<-datos[datos$Nivel2=="Xenobiotics_biodegradation_and_metabolism",]
datos<-datos[,-c(1,2,3)]

rows<-datos[,1]
datos<-datos[,-1]

datos<- as.data.frame(lapply(datos, function(col) as.integer(col)))

rownames(datos)<-rows

### Guardamos las descripciones y los KO para identificar posteriormete los resultados significantes 

decrip<-datos[,c(1,2)]

selex.sub<-datos
#Verificamos que los datos sean int. Esto es necesario par ausar ALDEx2
str(selex.sub)





##### Importamos datos csv de QIIME 2  ###############
Datos<-as.data.frame(read.csv("./01_Datos/qiime/level-3.csv"))
l<-length(Datos)

rownames(Datos)<-Datos[,1]
Datos<-Datos[,c(2:l)]
Datos<-as.data.frame(t(Datos))

selex.sub<-Datos
#Verificamos que los datos sean int. Esto es necesario par ausar ALDEx2
str(selex.sub)

decrip<-rownames(Datos)


##### Importamos datos de PICRUSt2 el pred_metagenome_unstrat_descrip ######

Datos<- read.delim("./01_Datos/picrust2/pred_metagenome_unstrat_descrip.tsv",
                   header=TRUE, sep="\t")
rows<-Datos[,1]

### Guardamos las descripciones y los KO para identificar posteriormete los resultados significantes 

decrip<-Datos[,c(1,2)]

Datos<-Datos[,-c(1,2)]
Datos<- as.data.frame(lapply(Datos, function(col) as.integer(col)))
rownames(Datos)<-rows

selex.sub<-Datos
#Verificamos que los datos sean int. Esto es necesario par ausar ALDEx2
str(selex.sub)


##### Generamos vector con orden de los sitios por variables #####
condsReg<-c("Este","Este","Este","Este","Oeste","Oeste","Oeste","Oeste")

condsEst<-c("Conservado","Conservado","Alterado","Alterado","Conservado","Conservado","Alterado","Alterado")


##### Aplicamos la función ALDEX ######

x.allC <- aldex(selex.sub, condsEst, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)

x.allR <- aldex(selex.sub, condsReg, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)



##### Graficamos ######

#Regiones 
##Diagrama MW
par(mfrow=c(1,3))
aldex.plot(x.allR, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference",all.cex =1, called.cex =2,
           rare.cex=2)
#DIagrama Ma
aldex.plot(x.allR, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference",all.cex = 2,called.cex =2,
           rare.cex=2)
#VOlcano plot
aldex.plot(x.allR, type="volcano", test="welch", main='volcano plot')



#Estatus de conservación 
## Diagramas MA

par(mfrow=c(1,3))
#Diagrama Mw
aldex.plot(x.allC, type="MW", test="welch", xlab="Log-ratio abundance",
           ylab="Difference",all.cex = 2,called.cex =2
           ,rare.cex=2)
#Diagrama Ma
aldex.plot(x.allC, type="MA", test="welch", xlab="Dispersion",
           ylab="Difference",all.cex = 1,called.cex =2
           ,rare.cex=2)
## Volcano plot
aldex.plot(x.allC, type="volcano", test="welch", main='volcano plot')


##### Extraemos los items que fueron significativos ####
SiginificantesR <- x.allR[x.allR$we.ep<=0.05,]
SiginificantesR$function.<-row.names(SiginificantesR)
SiginificantesC <- x.allC[x.allC$we.ep<=0.05,]
SiginificantesC$function.<-row.names(SiginificantesC)

##### Significantes bonferroni p 0.06 #####
SiginificantesRBH <- x.allR[x.allR$we.eBH<0.05,]
SiginificantesRBH$function.<-row.names(SiginificantesRBH)

SiginificantesCBH <- x.allC[x.allC$we.eBH<0.05,]
SiginificantesCBH$function.<-row.names(SiginificantesCBH)




