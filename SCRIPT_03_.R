# Instalação e Carregamento dos Pacotes Necessários para a Aula -----------

pacotes <- c("rgdal","raster","tmap","maptools","tidyverse","broom","knitr",
             "kableExtra","RColorBrewer")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}


################################## SHAPEFILES ##################################

# 1. PARTE INTRODUTÓRIA

# Carregando um shapefile -------------------------------------------------
shp_sp <- read_sf(dsn = "shapefile_sp", layer = "estado_sp")

# Características básicas do objeto shp_sp
summary(shp_sp)

# Classe e tipo do objeto carregado
class(shp_sp)

# Acessando a base de dados e outros componentes do objeto shp_sp ---------
shp_sp

# Para acessar as variáveis da base de dados atrelada ao shapefile, utilizaremos
# o operador $:
shp_sp$NM_MUNICIP
shp_sp$CD_GEOCMU

# Plotagem básica de um shapefile -----------------------------------------
plot(shp_sp)


# Introdução à manipulação de dados em shapefiles -------------------------

#  caso haja a necessidade da inserção de dados externos? ---------------

# Carregando uma base de dados real a respeito dos municípios de SP:
load("dados_sp.RData")

# Observando a base de dados carregada
dados_sp %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = TRUE, 
                font_size = 12)

# Para combinar os dados do objeto dados_sp com a base de dados de nosso 
# shapefile, podemos utilizar a função merge():
shp_dados_sp <- merge(x = shp_sp,
                      y = dados_sp,
                      by.x = "CD_GEOCMU",
                      by.y = "codigo")



# Salvando nosso shapefile:
write_sf(obj = shp_dados_sp, 
         layer = "nosso_novo_shapefile", 
         driver = "ESRI Shapefile", 
         dsn = "shp_novo")

# 2. VISUALIZAÇÃO DE DADOS ESPACIAIS

# Utilizando a tmap: ------------------------------------------------------
tm_shape(shp = shp_dados_sp) +
  tm_fill(col = "idh", palette = "Blues")


# Como saber quais paletas de cores podem ser utilizadas? -----------------
display.brewer.all()

# Vamos reconstruir o último mapa, utilizando uma nova paleta de cor e propondo
# 4 variações de cores:
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "Greens")

# Adicionando um histograma ao mapa anterior
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "Greens", 
          legend.hist = TRUE)

# Reposicionando o histograma
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "Greens", 
          legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)

# Posicionando manualmente o histograma, e adicionando um título
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "BuPu", 
          legend.hist = TRUE) +
  tm_layout(legend.text.size = 0.7,
            legend.title.size = 0.9,
            legend.hist.size = 0.5,
            legend.hist.height = 0.2,
            legend.hist.width = 0.3,
            frame = FALSE,
            main.title = "A Distribuição do IDH nos Municípios de SP")

# Adicionando uma bússola e bordas aos polígonos
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "Reds", 
          legend.hist = TRUE) +
  tm_layout(legend.text.size = 0.7,
            legend.title.size = 0.9,
            legend.hist.size = 0.5,
            legend.hist.height = 0.2,
            legend.hist.width = 0.3,
            frame = F,
            main.title = "A Distribuição do IDH nos Municípios de SP") +
  tm_borders(alpha = 0.8) +
  tm_compass(type = "8star", 
             show.labels = 3)

# Adicionando escala
tm_shape(shp = shp_dados_sp) + 
  tm_fill(col = "idh", 
          style = "quantile", 
          n = 4, 
          palette = "Reds", 
          legend.hist = TRUE) +
  tm_layout(legend.text.size = 0.7,
            legend.title.size = 0.9,
            legend.hist.size = 0.5,
            legend.hist.height = 0.2,
            legend.hist.width = 0.3,
            frame = F,
            main.title = "A Distribuição do IDH nos Municípios de SP") +
  tm_borders(alpha = 0.8) +
  tm_compass(type = "8star", 
             show.labels = 3)+
  tm_scale_bar()


# Instalação e Carregamento dos Pacotes Necessários para a Aula -----------

pacotes <- c("tidyverse","sf","tmap","rgeos","adehabitatHR","knitr",
             "kableExtra")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}


# 1. CRIANDO UM OBJETO SF A PARTIR DE UMA BASE DE DADOS

# Carregando a base de dados
load("shoppings.RData")

# Observando a classe do objeto shoppings:
class(shoppings)

# Observando as variáveis da base de dados shoppings:
shoppings %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = TRUE, 
                font_size = 12)

# Criando um objeto do tipo sf a partir de um data frame:
sf_shoppings <- st_as_sf(x = shoppings, 
                         coords = c("longitude", "latitude"), 
                         crs = 4326)

# Observando a classe do objeto sf_shoppings:
class(sf_shoppings)

# Um componente interessante do objeto sf_shoppings é chamado geometry:
sf_shoppings$geometry

# Note que um objeto sf é, de fato, um data frame georreferenciado. Não há
# polígonos atrelados a ele, nem a necessidade de se utilizar o operador @!

# Plotando o objeto sf_shoppings de forma espacial:
tm_shape(shp = sf_shoppings) + 
  tm_dots(size = 1)

# Adicionando uma camada de um mapa do Leafleet que considere a bounding box do 
# objeto sf_shoppings:
tmap_mode("view")
tm_shape(shp = sf_shoppings) + 
  tm_dots(col = "deepskyblue4", 
          border.col = "black", 
          size = 0.2, 
          alpha = 0.8)

tmap_mode("plot")

#tmap_mode("plot") #Para desativar as camadas de mapas on-line


# 2. COMBINANDO UM OBJETO SIMPLE FEATURE COM UM SHAPEFILE

# Carregando um shapefile do município de São Paulo
shp_saopaulo <- read_sf("shapefile_municipio", "municipio_sp")

# Visualização gráfica do objeto shp_saopaulo:
tm_shape(shp = shp_saopaulo) + 
  tm_borders()

# Combinando o objeto shp_saopaulo com o objeto sf_shoppings:
tm_shape(shp = shp_saopaulo) + 
  tm_borders(alpha = 0.5) +
  tm_shape(shp = sf_shoppings) + 
  tm_dots(col = "regiao", 
          size = 0.02)

# 3. KERNEL DENSITIES

# A técnica de kernel densities calcula a densidade da presença de pontos de
# interesse em determinada área geográfica.

# O primeiro passo será criar um objeto sp com a base de dados atrelada a ele:
shoppings_sp_df <- SpatialPointsDataFrame(data = shoppings,
                                          coords = coordenadas_shoppings,
                                          proj4string = CRS("+proj=longlat"))


# Note como a função SpatialPointsDataFrame() permite a existência de um data
# frame junto a nosso objeto de classe sp:
shoppings_sp_df@data

# Para o cálculo das kernel densities, podemos utilizar a função kernelUD():
shoppings_dens <- kernelUD(xy = shoppings_sp_df,
                           h = "href",
                           grid = 1000,
                           boundary = NULL)

plot(shoppings_dens)

# Para estabelecer as zonas com maior densidade, propomos o seguinte:
zona1 <- getverticeshr(x = shoppings_dens, percent = 20) 
zona2 <- getverticeshr(x = shoppings_dens, percent = 40) 
zona3 <- getverticeshr(x = shoppings_dens, percent = 60) 
zona4 <- getverticeshr(x = shoppings_dens, percent = 80)

tmap_options(check.and.fix = TRUE) 

tm_shape(shp = shp_saopaulo) + 
  tm_fill(col = "gray90") + 
  tm_borders(col = "white", alpha = 0.5) + 
  tm_shape(shp = shoppings_sp_df) + 
  tm_dots(col = "regiao", size = 0.25) + 
  tm_shape(zona1) + 
  tm_borders(col = "firebrick4", lwd = 2.5) +
  tm_fill(alpha = 0.4, col = "firebrick4") + 
  tm_shape(zona2) + 
  tm_borders(col = "firebrick3", lwd = 2.5) + 
  tm_fill(alpha = 0.3, col = "firebrick3") + 
  tm_shape(zona3) + 
  tm_borders(col = "firebrick2", lwd = 2.5) + 
  tm_fill(alpha = 0.2, col = "firebrick2") +
  tm_shape(zona4) + 
  tm_borders(col = "firebrick1", lwd = 2.5) + 
  tm_fill(alpha = 0.1, col = "firebrick1")


# FIM ---------------------------------------------------------------------
#PARTE 2. BANCOS DE DADOS ON-LINE-----------------------------
#Instalar Pacote
install.packages("ipeadatar")
install.packages("tmap")
#Acionar pacote
library(ipeadatar)
library(tmap)
devtools::install_github("ipeaGIT/geobr", subdir = "r-package")
library(geobr)

## Vamos explorar os dados do IPEA
#analisar temas disponiveis
temas=available_subjects(language =  "br")
#Analisar series de dados disponiveis
series=available_series(language = "br")
#Filtrar s?ries de dados por temas
series_social=subset(series,theme == "Social")
#Buscar dados -> Analfabetismo
EDUC=ipeadata("ADH_T_ANALF15M", language = c("en", "br"), quiet = FALSE)
#Filtrar Dados por estado
EDU_RN=subset(EDUC,tcode >=2400000 & tcode<=2500000)
#Filtrar Dados por Pesquisa
EDURN_2010=subset(EDU_RN,date=="2010-01-01")

#Associar a Mapa usando o GEOBR

##Extrair malha do RN
RN=read_municipality(code_muni = 24,year = 2010,simplified = TRUE,showProgress = TRUE)
plot(RN)
##Juntando SHPS + informações baixadas
shp_educ_rn <- merge(x = RN,
                      y = EDURN_2010,
                      by.x = "code_muni",
                      by.y = "tcode")
#Plotagem Basica
tm_shape(shp = shp_educ_rn) +
  tm_fill(col = "value", palette = "Greens")
#------------------------------------FIM PARTE 2 ------------------------------#
#PARTE 3 voronoi e IDW--------------------------------#
# instalar pacotes e verficar diretorio de trabalho
install.packages("raster")
install.packages("sf")
install.packages("sp")

library(raster)
library(sf)
library(sp)
library(tmap)
## outros pacotes que vamos usar ao longo da aula
install.packages("dismo")
install.packages("stars")
install.packages("gstat")
install.packages("fields")
install.packages("rcompanion")
install.packages("bestNormalize")
install.packages("geoR")
install.packages("automap")
install.packages("deldir")
library(dismo)
library(stars)
library(gstat)
library(fields)
library(rcompanion)
library(bestNormalize)
library(geoR)
library(automap)
library(ggplot2)
library(dplyr)
library(deldir)
## Importando dados
st_layers("aula3.gpkg")
##Vamos criar um objeto contendo a camanda com os limites municipais
mun <- st_read("aula3.gpkg", layer="mun_abc")
#Vamos associa-la a crs SIRGAS 2000 23S
st_crs(mun)<-31983
#Agora vamos extrair o layer com as estações pluviométricas e definir o mesmo CRS
estacoes <- st_read("aula3.gpkg", layer="estacoes")
st_crs(estacoes)<-31983
View(estacoes)
#Plotagem Básica dos pontos com contornos
tm_shape(estacoes) + tm_dots(col="chuva", pal="Blues") +
  tm_shape(mun) + tm_borders(col="black")

# Calculando o Vizinho mais proximo
estacoes_sp <- as(estacoes,"Spatial")

#Calcular poligono do vizinho mais perto
voronoi <- voronoi(estacoes_sp)
voronoi_sf <- st_as_sf(voronoi)
View(voronoi_sf)
plot(st_geometry(voronoi_sf))
plot(st_geometry(estacoes), pch=20, cex=0.4, col="blue", add=TRUE)
tm_shape(voronoi_sf) + 
  tm_fill(col="chuva", style="quantile", pal="Blues") + tm_borders() +
  tm_shape(mun) + tm_borders(col="black")
tm_shape(voronoi_sf) + 
  tm_fill(col="temperat", style="quantile", pal="Reds") + tm_borders() +
  tm_shape(mun) + tm_borders(col="black")

# Inverso da Distancia - IDW
##Usaremos um raster como base, transformaremos sua resolução em 10x
modelo_raster <- raster("srtm_abc.tif", values=FALSE)
modelo_raster
modelo_raster <- aggregate(modelo_raster, fact = 10)
st_crs(estacoes) <- crs(modelo_raster)

library(stars)
library(gstat)

##IDP = quociente de distancia #Alterar pesos para fazer testes
gs_idw_2 <- gstat(formula=chuva~1, locations=estacoes, set=list(idp=4))

##Validação Cruzada -> Erros do modelo
cv_idw_2 <- gstat.cv(gs_idw_2)
View(as.data.frame(cv_idw_2))

##Raiz quadrada do erro m?dio
sqrt(mean(cv_idw_2$residual^2))                        # raiz do erro medio quadrado

##Percentual da variabilidade espacial explicada
1-(var(cv_idw_2$residual)/var(cv_idw_2$observed))      # estimativa da porcentagem da variacao explicada (R2)
chuva_idw_2 <- interpolate(object=modelo_raster, model=gs_idw_2) # Funcao do pacote raster - demora um pouco 
plot(chuva_idw_2)
plot(st_geometry(mun), add=TRUE)


#PARTE 4 SUPERFICIE DE TENDENCIA E SPLINE
## Modelo de Interpolação Superfície de Tendencia -> Alterar Graus do polinomio
gs_superficie_1 <- gstat(formula=chuva~1, location=estacoes, degree=1)

## Cross validation -> Erro e Explicação
cv_superficie_1 <- gstat.cv(gs_superficie_1)
sqrt(mean(cv_superficie_1$residual^2))
1-(var(cv_superficie_1$residual)/var(estacoes$chuva))

##Criar Raster
superficie_1 <- interpolate(modelo_raster, gs_superficie_1)
plot(superficie_1)
plot(st_geometry(mun), add=TRUE)

##Grafico tridimensional para representar a dependencia espacial
persp(superficie_1,border=NA, col="blue", theta = 320, phi=40)


## SPLINE
##Instalar Pacote fields
install.packages("fields")
library(fields)

##Criar modelo SPLINE
modelo_tps <- Tps(x=st_coordinates(estacoes), Y=estacoes$chuva)

##Calculando o erro e o R2
sqrt(mean(modelo_tps$residual^2))
1-(var(modelo_tps$residual)/var(estacoes$chuva))
chuva_tps <- interpolate(modelo_raster, modelo_tps)
plot(chuva_tps)
plot(st_geometry(mun), add=TRUE)
##Proje??o 3D
persp(chuva_tps, border=NA, col="blue", shade=0.3, theta = 320, phi=40)


##KRIGAGEM
estacoes_sp %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(size=chuva), color="blue", alpha=3/4) + 
  ggtitle("Chuva mm)") + coord_equal() + theme_bw()


##Criar Grid
grid=rasterToPoints(modelo_raster, fun=NULL, spatial=T)
crs(estacoes_sp) <- crs(grid)

##Gerar Variograma e o modelo de interpola??o
# calculates sample variogram values 
variograma <- variogram(chuva~1, estacoes_sp, cutoff=100000)
plot(variograma,plot.numbers=F)

##Informações do Variograma
estacoes_geodata <- as.geodata(estacoes_sp[,"chuva"])
variograma_geor <- variog(estacoes_geodata)
dev.new()
eyefit_geor <- eyefit(variograma_geor)

##Variograma Original
dev.off()
eyefit_geor
plot(variograma_geor)
lines(eyefit_geor)

##Ajuste do Variograma (variofit)
variofit_geor <- variofit(variograma_geor, ini.cov.pars = eyefit_geor)
variofit_geor
plot(variograma_geor)
lines(variofit_geor)

#Parâmetros
variofit_geor$nugget          # efeito pepita
variofit_geor$cov.pars[1]    # patamar
variofit_geor$cov.pars[2]    # alcance

#Ajustar modelo, #patamar; #tipo de modelo;#alcance; #pepita 
lzn.fit <- fit.variogram(variograma, model=vgm(708771.9, "Gau", 52413.24, 0)) # fit model
plot(lzn.fit,cutoff=100000)

##gerar krigagem
lzn.kriged <- krige(chuva ~ 1, estacoes_sp, grid, model=lzn.fit)
lzn.kriged %>% as.data.frame %>%
  ggplot(aes(x=x, y=y)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous() + scale_y_continuous() +
  theme_bw()

##Plotagem
tm_shape(shp=lzn.kriged)+
  tm_dots("var1.pred",pal = "Blues",alpha=0.2)+tm_shape(mun)+tm_borders(col="black") + tm_shape(estacoes_sp) + tm_dots(col="chuva", pal="Blues")







