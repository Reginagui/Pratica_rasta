# 2.1) Shapes via geobr (dados vetoriais)
# Estado da Bahia
geobr::read_state(code_state = "BA", year = 2020, simplified = TRUE)
bahia <- geobr::read_state(code_state = "BA", year = 2020, simplified = TRUE)
bahia
library(geobr)
library(terra)
vect(bahia)
bahia_v <- vect(bahia)  # converter para SpatVector

# Município de Ilhéus (BA) como área de estudo
ilheus <- geobr::read_municipality(code_muni = 2913606, year = 2020, simplified = TRUE)
ilheus_v <- vect(ilheus)

crs(bahia_v)
crs(ilheus_v)
plot(bahia_v, col = "gray90")
plot(ilheus_v, add = TRUE, col = "tomato")
library(geodata)
# 2.2) Raster WorldClim (12 bandas de temperatura média mensal)
# Requer o pacote 'geodata' para download (gera SpatRaster já no formato do terra)
if (requireNamespace("geodata", quietly = TRUE)) {
  wc_tavg <- geodata::worldclim_global(var = "tavg", res = 10, path = "data/rasters")
  # 'wc_tavg' é multibanda (12 camadas: meses)
  wc_tavg
} else {
  # fallback: gerar um raster sintético para fins didáticos (se geodata não estiver disponível)
  set.seed(42)
  wc_tavg <- rast(ext( -41, -37, -16, -13), resolution = 0.01, crs = "EPSG:4674") # WGS84 / SIRGAS 2000
  wc_tavg <- c(wc_tavg, wc_tavg, wc_tavg) # 3 bandas de exemplo
  names(wc_tavg) <- paste0("tavg_", 1:nlyr(wc_tavg))
  values(wc_tavg) <- runif(ncell(wc_tavg) * nlyr(wc_tavg), min = 18, max = 28)
  wc_tavg
}
# Definir um CRS métrico apropriado (UTM zona 23S para Ilhéus/BA)
crs_utm23s <- "EPSG:31983"

bahia_u   <- project(bahia_v, crs_utm23s)
ilheus_u  <- project(ilheus_v, crs_utm23s)

# Reprojetar raster: alinhar com projeção UTM
project(wc_tavg, crs_utm23s, method = "bilinear")
wc_tavg_u <- project(wc_tavg, crs_utm23s, method = "bilinear")
crs(bahia_u)
crs(wc_tavg_u)

# Tabela de atributos
head(as.data.frame(ilheus_u))

# Exemplo: selecionar polígonos por atributo (no caso, ilhéus_u é 1 feição)
sel <- ilheus_u[ ilheus_u$code_muni == 2913606, ]
plot(bahia_u, col="gray95")
plot(sel, add=TRUE, col="steelblue")
# Amostra de pontos em Ilhéus (grade simples)
set.seed(1)
pts <- spatSample(ilheus_u, size = 100, method = "random")
plot(ilheus_u, col="grey90")
points(pts, pch=20, cex=.6)
# Buffer de 1 km nos pontos
buf1k <- buffer(pts, width = 1000)
plot(ilheus_u, col="grey90")
plot(buf1k, add=TRUE, border="red")
# Criar um retângulo artificial cortando a área de Ilhéus
bb <- as.polygons(ext(ilheus_u), crs(ilheus_u))
bb_cut <- buffer(centroids(bb), 15000) # círculo ~15 km
# Interseção: parte comum
ilh_inter <- intersect(ilheus_u, bb_cut)
# Diferença: parte de Ilhéus fora do círculo
ilh_diff  <- erase(ilheus_u, bb_cut)

par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(ilheus_u, col="grey90", main="Original")
plot(bb_cut, add=TRUE, border="purple")
plot(ilh_inter, col="lightgreen", main="intersect(ilhéus, círculo)")
plot(ilh_diff,  col="salmon", main="erase(ilhéus, círculo)")
par(mfrow=c(1,1))
set.seed(42)
br <- geobr::read_state(code_state = "all", year = 2020, simplified = TRUE)
br_v <- project(vect(br), crs_utm23s)
pts <- spatSample(br_v, size = 500, method = "random")
plot(br_v, col="grey90")
points(pts, pch=20, cex=.6)
# Quais pontos caem dentro da Bahia (TRUE/FALSE)?
inside <- relate(pts, bahia_u, relation = "intersects")
table(inside)
# Selecionar apenas os pontos dentro
pts_in <- pts[inside, ]
pts_out <- pts[!inside, ]
plot(bahia_u, col="grey90")
points(pts_in, pch=20, col="blue")
points(pts_out, pch=20, col="red")
# Recortar pela extensão da Bahia
wc_crop <- crop(wc_tavg_u, bahia_u)
# Mask: manter somente os pixels dentro do polígono
wc_mask <- mask(wc_crop, bahia_u)
wc_mask
plot(wc_mask[[1]], main="WorldClim (tavg mês 1) — Bahia")
plot(bahia_u, add=TRUE, border="black")
# Exemplo simples: anomalia = mês 1 - média anual
mean_annual <- app(wc_mask, fun = mean, na.rm = TRUE)
anom_m1 <- wc_mask[[1]] - mean_annual
plot(anom_m1, main="Anomalia (mês 1 - média anual)")
plot(ilheus_u, add=TRUE, border = "red")
plot(bahia_u, add=TRUE)
# Reclassificar média anual em faixas (valores ilustrativos)

rcl <- matrix(c(-Inf, 20, 1,
                20, 24, 2,
                24, 28, 3,
                28,  Inf, 4), ncol = 3, byrow = TRUE)
tzone <- classify(mean_annual, rcl = rcl, include.lowest = TRUE, right=FALSE)
tzone <- as.factor(tzone)
levels(tzone) <- data.frame(ID=1:4, zone=c("≤20", "20–24", "24–28", ">28"))
plot(tzone, main="Zonas térmicas (faixas)")
plot(ilheus_u, add=TRUE, border = "red")
plot(bahia_u, add=TRUE)
# Upscale (coarsen): reduzir resolução espacial (agrupar células) com 'aggregate'
coarse <- aggregate(wc_mask[[1]], fact = 4, fun = "mean", na.rm = TRUE)

# Downscale (refinar): aumentar resolução (dividir células) com 'disagg'
fine   <- disagg(coarse, fact = 4, method = "bilinear")

par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(wc_mask[[1]], main="Original")
plot(coarse, main="Upscale: aggregate(fact=4)")
plot(fine,   main="Downscale: disagg(fact=4)")
# Estatística zonal: média por zona (ex.: buffers de 1 km como zonas)

# Usando extract() com `fun`
buff_stats <- extract(wc_mask[[1]], buf1k, fun = mean, na.rm = TRUE)
head(buff_stats)
# Pelas zonas
zone_stats <- zonal(wc_mask[[1]], tzone, fun = mean)
zone_stats
# Extração pontual: valor do raster em cada ponto
vals_pts <- extract(wc_mask[[1]], pts_in)
head(vals_pts)
# Filtro focal (janela 3x3) - suavização
ker <- matrix(1, 3, 3)
smooth <- focal(wc_mask[[1]], w = ker, fun = mean, na.rm = TRUE)
plot(smooth, main="Focal mean (3x3)")
plot(ilheus_u, add=TRUE)
# Focal (janela 3x3) - Preenchimento
ker <- matrix(1, 5, 5)
wc_temp_na <- wc_mask[[1]]
wc_temp_na[sample(cells(wc_temp_na), 200)] <- NA
wc_temp_sem_na <- focal(wc_temp_na, w = ker, fun = mean, na.rm = TRUE, na.policy = "only")

par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(wc_temp_na, main="Com NA")
plot(bahia_u, add=TRUE)
plot(wc_temp_sem_na, main="Sem NA")
plot(bahia_u, add=TRUE)
plot(mask(wc_temp_sem_na, bahia_u), main="Sem NA recortado")
plot(bahia_u, add=TRUE)
par(mfrow=c(1,1))

# Superfície de distância até os pontos amostrados (centroides como "alvo")
# Criar um raster binário com 1 nas células dos pontos, NA no restante
rtmp <- rasterize(pts_in, wc_mask[[1]], field = 1)
rdist <- distance(rtmp)  # distância euclidiana até os "alvos"
plot(rdist, main="Distância até pontos (m)")
plot(bahia_u, add=TRUE, border="black")
plot(ilheus_u, add=TRUE, border="black")
points(pts_in, pch=20, col="red")
# Mapas rápidos com plot(), sobrepondo camadas
plot(wc_mask[[1]], main="Mapa — Bahia e tavg (mês 1)")
plot(bahia_u, add=TRUE, border="black")
plot(ilheus_u, add=TRUE, border="black")
points(pts_in, pch=20, col="blue")
# Exportar resultados
dir.create("outputs/dados", recursive = TRUE, showWarnings = FALSE)
writeVector(bahia_u,  "outputs/dados/bahia_utm23s.gpkg", overwrite=TRUE)
writeVector(ilheus_u,  "outputs/dados/ilheus_utm23s.gpkg", overwrite=TRUE)
writeRaster(wc_mask,   "outputs/dados/wc_tavg_ilheus.tif", overwrite=TRUE, gdal="COMPRESS=DEFLATE")
writeRaster(anom_m1,   "outputs/dados/wc_tavg_anom_m1.tif", overwrite=TRUE, gdal="COMPRESS=DEFLATE")
