#le as imagens geradas pelo momocs.R
#salva os resultados para uso depois
rm(list=ls())
library(EBImage)
library(Momocs)
library(shape)
mytrim <- function (x) {gsub("^\\s+|\\s+$", "", x)}
hasblack <- function(ln) { length(which(ln==0)) }

##PEGA A ESCALA DE CADA FOLHA , OU SEJA QUANTOS PIXELS CORRESPONDEM A 1 CM
lf <- list.files('reticulata_com_base_nir_shape/nir-1folha_Alves2014/escala', full.names=TRUE)
lf2 <- list.files('reticulata_com_base_nir_shape/nir-1folha_Alves2014/escala', full.names=FALSE)
sizes = NULL
for(l in 1:length(lf)) {
  img = readImage(files=lf[l])
  xt = img@.Data
  rn = apply(xt,1,hasblack)
  ln = apply(xt,2,hasblack)
  rn = rn[rn>0]
  ln = ln[ln>0]
  ll = c(rn,ln)
  pix = median(ll)
  unit = 1
  #which are the left most column having a black value
  sizes = rbind(sizes,data.frame(FILE=lf2[l],pix=pix,unit=unit,stringsAsFactors = F))
}

##PEGA AS IMAGENS E IMPORTA AS FORMAS NA FORMA DE COORDENADAS X E Y
lf <- list.files('reticulata_com_base_nir_shape/nir-1folha_Alves2014/folha', full.names=TRUE)
coo = import_jpg(lf)
coob = coo
coo = Out(coo)
coof = coo
class(coof)
##PEGA O APICE E A BASE DE CADA FOLHA E ADICIONA ISSO COMO LANDMARK
mlkds = NULL
for(i in 1:length(coof[[1]])) {
  print(i)
  x = coof[[1]][[i]]
  g1 = which(x[,1]==min(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx1 = g1[d2==min(d2)][1]

  g1 = which(x[,1]==max(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx2 = g1[d2==min(d2)][1]
  g1 = c(idx1,idx2)
  plot(x)
  points(x[g1,1],x[g1,2],pch=21,col='red',cex=3)
  #return(g1)
  mlkds[[i]] = g1
}
coof[['ldk']] = mlkds

#RESCALA OS OBJETOS PARA CM SEGUNDO AS ESCALAS ESTRAIDAS
#PEGA O IDENTFICADOR DAS FOLHAS
sizes$FILE

#ADICIONA UM INDICE para mag (magnification) - fake
ff = as.factor(sizes$pix)
lv = levels(ff)
levels(ff) = 1:length(lv)
ss = data.frame(mag=ff,sizes[,-1])
colnames(ss)[3] ='cm'
rownames(ss) = sizes$FILE
#coof = coofb
names(coof)
coof[['fac']] = ss

##FAZ A CONVERSAO DE UNIDADE DE PIXEL PARA UNIDADE EM CM
ss = unique(ss)
rownames(ss) = NULL
names(coof)
cooescaled = rescale(coof,scale_mapping =ss,magnification_col ='mag')

#mostra o que aconteceu os eixos agora sao cm
#cooescaled %T>% stack()

#FILTRA APENAS PDBFF PARA ESTA ETAPA
#coofpdb = filter(cooescaled, epdb==1)
#coofpdb %T>% stack()

#NAO FAZ DIFERENCA SELECIONANDO POUCAS AMOSTRAS POR ESPECIE O PADRAO E O MESMO
#alinha usando os landmarks (apice e base) adicionados
#ou seja assume homologia entre a base e o apice
###IMPORTANTE###
cooescaled.aligned = fgProcrustes(cooescaled)
#IMPORTANTE RESPONDEU: cannot apply fgProcrustes on 
#less than three landmarks. coo_bookstein is 
#returned (usa entao coo_bookstein).. IGNORE a mensagem e prossiga


##incluindo o tamanho junto com os pcas
x = cooescaled[[1]][[1]]
pegalonglarg <- function(x) {
  maxx = diff(range(x[,1]))
  mayy = diff(range(x[,2]))
  area = sum(x)
    c(maxx,mayy,area)
}
tamanho = sapply(cooescaled[[1]],pegalonglarg)
tamanho = t(tamanho)
rownames(tamanho) = NULL
colnames(tamanho) = c("comprimento.cm","largura.cm","area.foliar")
tamanho = data.frame(tamanho)
tamanho$RatioLargura.Comprimento = tamanho[,2]/tamanho[,1]
rownames(tamanho) = names(cooescaled)
head(tamanho)


#pega os fourier (o numero de harmonicas e decidido automaticamente)
cooescaled.fourier = efourier(cooescaled.aligned, norm=F)
#respondeu: 'nb.h' set to NUM (99% harmonic power)

#lets look the amplitude of fitted coefficients
#hist(cooescaled.fourier,drop=0)

#calcula o PCA que descreve a forma
pcaforma = PCA(cooescaled.fourier)
summary(pcaforma)
dir.create("Resultados")
#imprime os eixos para entender quais os importantes (precisa fazer isso visualmente)
#eixo 1 vs. os demais. Todas as figuras tem o eixo PCA1
fn = "Resultados/shapereticulata-nir-1folha_Alves2014.pdf"
#vamos gerar um pdf com 18 x 18 cm, que é um tamanho que cabe num artigo
#da para melhorar isso
pdf(fn,width=(18/2.54),height=(18/2.54))
for(i in c(2:6)) {
    par(mar=c(1,1,1,1))
    xl = pcaforma$x[,1]
    xa = abs(diff(range(xl)))*0.2
    xl = range(xl)+c(xa*(-1),xa)
    yl = pcaforma$x[,i]
    ya = abs(diff(range(yl)))*0.2
    yl = range(yl)+c(ya*(-1),ya)
    plot(pcaforma,xax = 1, yax = i,pts.shp =500,chullfilled=T,eigen=T,zoom=0.9,col.shp='lightgreen',amp.shp = 0.6, size.shp = 1.4,pos.shp = c("range", "full", "circle", "xy", "range_axes", "full_axes")[1],pch=21,cex=0.3,col='red',nb.grids=5,title="",unit=TRUE,xlim=xl,ylim=yl)
}
dev.off()

#SALVA A IMAGE POR GARANTIA OU PARA USO DEPOIS
save.image(file="Resultados/extracaoShapereticulata-nir-1folha_Alves2014.Rdata")


#salvando os dados para análises LDA etc.
#salvar todos os eixos (filtrar os informativos depois)
final = pcaforma$x
rn = rownames(final)
final = data.frame(final)
rn = rownames(final)
#junta as metricas de tamanho da folha
g = match(rn,rownames(tamanho))
final = data.frame(final,tamanho[g,])
write.table(final,file=paste('Resultados/dadosShapereticulata-nir-1folha_Alves2014_',
                             Sys.Date(),'.csv',sep=""),sep="\t",na='',quote=T,
            row.names = F)
#se colocar T em row.names aparecerão os nomes dos arquivos que estão em ImagensSimple
# notifica que a analise acabou porque nao somos obrigados a ficar conferindo toda hora :)
library(beepr)
beepr::beep(8)

