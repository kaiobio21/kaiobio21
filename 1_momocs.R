#install.packages("Momocs",dependencies = T)
#remove.packages("Momocs")
#library(devtools)
#devtools::install_github("MomX/Momocs",dependencies = T)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("EBImage", force = T)


library(EBImage)
library(Momocs)
rm(list=ls())
dir()

arqs = paste("paraMomocs/",dir("paraMomocs"),sep="")
#dir.create("ImagensSimple")
#dir.create("Escalas")
arqs.fn = dir("paraMomocs")


for(a in 1:length(arqs.fn)) {
  print(paste(a,length(arqs),sep="/"))
  img = readImage(files=arqs[a])
  colorMode(img) = Grayscale
  z = getFrame(img,1)
  z2 = getFrame(img,2)
  z3 = getFrame(img,3)
  zz = z*z2*z3
  zz[zz< 0.9] = 0
  zz[zz>= 0.9] = 1
  z = zz
  #remove noise
  z = medianFilter(z,2)

  plot(z)
  #separa escala da imagem
  #quais celulas da matrix formam um quadrado
  x = z@.Data
  hasblack <- function(ln) {
    length(which(ln==0))
  }
  #quais colunas tem tinta preta
  lbl = apply(x,2,hasblack)

  #quais sao valores sequenciais
  lbd = diff(lbl)
  #hist(lbd)
  #abline(v=20)
  #qual é a primeira coluna que tem
  x1 = which(lbd>20)[1]+1
  g2 = which(lbd<(-20))
  x2 = g2[g2>x1][1]-1
  yy = x[,x1:x2]

  #isso basicamente pegou as colunas onde esta a escala
  plot(z)
  abline(v=c(x1,x2),col='red')

  #dentro dessas colunas quais linhas tem tinta preta
  cbl = apply(yy,1,hasblack)
  lbd = diff(cbl)
  #pega o bloco superior
  y1 = which(lbd>20)[1]+1
  g2 = which(lbd<(-20))
  y2 = g2[g2>y1][1]-1
  abline(h=c(y1,y2),col='red')


  #linhas e colunas correspondentes a escala + 2 pixels
  xs = (x1-2):(x2+2)
  ys = (y1-2):(y2+2)
  #points(xs,ys,pch=21,cex=0.1,bg='red',col='red')

  #pega a escala
  escala = z
  escala@.Data = x[ys,xs]
  escala = medianFilter(escala,1)
  #plot(escala)

  #elimina a escala da matriz sendo trabalhada
  #colocar branco no lugar da escala
  x[ys,xs] = 1

  #quantas folhas tem na imagem?
  lb1 = apply(x,1,sum)
  lb2 = apply(x,2,sum)

  #quais linhas são inteiramente brancas?
  brlinhas = which(lb2==nrow(x))
  plot(z)
  abline(h=brlinhas,col='red')
  #separa as folhas entre as linhas em branco
  #pega a diferenca dos valores
  lbd = diff(brlinhas)
  #inicio e fim das folhas
  inifols = brlinhas[which(lbd>1)]
  fimfols = brlinhas[which(lbd>1)+1]
  abline(h=inifols,lwd=2,col='blue')
  abline(h=fimfols,lwd=2,col='blue')

  #para cada folha salva o arquivo separadamente com um buffer de 2 pixels
  for(f in 1:length(inifols)) {
      #dado da folha
      linhainicial = inifols[f]
      linhafinal = fimfols[f]
      lnf = seq(linhainicial,linhafinal)
      #salva uma imagem para cada folha
      afolha = x[,lnf]
      #adiciona 2 linhas em branco de cada lado
      branco = rep(1,nrow(x))
      afolha = cbind(branco,branco,afolha,branco,branco)

      #pega a imagem original e filtra
      zf = z
      zf@.Data = afolha
      zf = medianFilter(zf,1)

      #cria um nome novo para esse arquivo
      norg = arqs.fn[a]
      add = paste("_folha-",LETTERS[f],".jpg",sep="")
      norg = gsub("\\.jpg|\\.jpeg",add,norg,ignore.case = T)
      fn = paste("ImagensSimple/",norg,sep="")

      #salva a imagem da folha
      writeImage(zf,fn, quality = 100)

      #salva uma escala para cada folha
      fn = paste("Escalas/",norg,sep="")
      writeImage(escala,fn, quality = 100)
    
  }
}

# notifica que a analise acabou porque nao somos obrigados a ficar conferindo toda hora :)
library(beepr)
beepr::beep(8)
