root <- normalizePath(commandArgs(trailingOnly=TRUE)[1], winslash='/')
dir.create(file.path(root,'figures'), showWarnings=FALSE)
dir.create(file.path(root,'qa'), showWarnings=FALSE)
stopifnot(requireNamespace('svglite', quietly=TRUE), requireNamespace('ragg', quietly=TRUE), requireNamespace('jsonlite', quietly=TRUE))
library(grid)
read_source <- function(name) read.delim(file.path(root, 'source_data', name), check.names=FALSE, stringsAsFactors=FALSE)
a <- read_source('Figure3a_response_effects.tsv')
b <- read_source('Figure3b_module_means.tsv')
c <- read_source('Figure3c_genus_attribution.tsv')
cl <- read_source('Figure3c_lineage_labels.tsv')
d <- read_source('Figure3d_salinity_attenuation.tsv')
mods <- b$feature
stopifnot(nrow(a)==7, nrow(b)==12, nrow(c)==11, identical(names(c)[-1], mods))
vals <- as.matrix(c[,-1])*100
totals <- b$delta_NZ_minus_SZ*100
stopifnot(max(abs(colSums(vals)-totals))<1e-10)
focus <- mods %in% c('Salt_Trk_Ktr_core','N_urea_utilization','Salt_betaine_synthesis','N_nitrate_red_periplasmic')
labels_b <- c('Trk/Ktr (K03498)','Malate DH','Formate DH','Urea hydrolysis','Betaine synthesis','NapAB','Succinyl-CoA synthetase','NarGHI-like','Aconitase','Hydroxyectoine','Citrate synthase','Choline/betaine transport')
labels_c <- c('K03498','Malate\nDH','Formate\nDH','ureABC','betAB','napAB','Succinyl-\nCoA syn.','NarGHI-\nlike','Aconitase','Hydroxy-\nectoine','Citrate\nsynthase','Choline/\nbetaine\ntransport')
lineage <- cl$display_label[match(c$genus_key, cl$lineage_key)]
stopifnot(!anyNA(lineage))
COL <- list(ink='#252A30', muted='#626E7A', grid='#DDE1E5', nz='#416C99', sz='#E5A13D', accent='#D96846', salt='#39958B', line='#A9B1B9')
pal <- grDevices::colorRampPalette(c('#E5A13D','#FFFFFF','#416C99'))(601)
heat <- function(v) pal[pmax(1,pmin(601,round((v+30)/60*600)+1))]
tx <- function(text,x,y,size=8,just='left',face='plain',col=COL$ink,rot=0) grid.text(text,x=unit(x,'native'),y=unit(y,'native'),just=just,rot=rot,gp=gpar(fontsize=size,fontfamily='Arial',fontface=face,col=col,lineheight=.97))
ln <- function(x0,y0,x1,y1,col=COL$ink,lwd=.65) grid.lines(x=unit(c(x0,x1),'native'),y=unit(c(y0,y1),'native'),gp=gpar(col=col,lwd=lwd))
rect <- function(x,y,w,h,fill,border=NA,lwd=.4) grid.rect(x=unit(x,'native'),y=unit(y,'native'),width=unit(w,'native'),height=unit(h,'native'),gp=gpar(fill=fill,col=border,lwd=lwd))
pt <- function(x,y,pch=16,col=COL$ink,size=1.8,fill=col,lwd=.9) grid.points(x=unit(x,'native'),y=unit(y,'native'),pch=pch,size=unit(size,'mm'),gp=gpar(col=col,fill=fill,lwd=lwd))
heading <- function(letter,title,subtitle=NULL) {
  tx(letter,1,3.5,10,face='bold');tx(title,6.5,3.5,9,face='bold')
  if(!is.null(subtitle)) tx(subtitle,1,9,7.2,col=COL$muted)
}
panel_a <- function() {
  heading('a','Overall basin terms','Seven station pairs; exact P at right')
  lab <- c('Read species','MAG composition','GMP','SCA','Carrier prevalence','KO + MAG detection','KO + mapping rate')
  xp <- function(x) 34+x/.42*43
  ys <- seq(23,62,length.out=7)
  tx('P',88,16,7.2,just='right',face='italic',col=COL$muted)
  for(i in seq_len(nrow(a))) {
    tx(lab[i],32,ys[i],7.4,just='right')
    ln(xp(a$low[i]),ys[i],xp(a$high[i]),ys[i],COL$line,1.6)
    pt(xp(a$R2[i]),ys[i],if(a$P[i]<.05)16 else 21,COL$accent,2.0,if(a$P[i]<.05) COL$accent else 'white',1.1)
    tx(sprintf('%.4f',a$P[i]),88,ys[i],7.3,just='right')
  }
  ln(xp(0),66,xp(.42),66)
  for(v in seq(0,.4,.1)) {ln(xp(v),66,xp(v),67);tx(sprintf('%.1f',v),xp(v),69,7.2,just='center')}
  tx('Basin R\u00b2 (95% interval)',55.5,73.3,7.6,just='center')
}
panel_b <- function() {
  heading('b','Exploratory module contrasts','Modules ranked by absolute station-mean basin contrast')
  xp <- function(x) 43+x*43
  pt(70,15,16,COL$nz,1.65);tx('NZ',72.5,15,7.1)
  pt(80,15,15,COL$sz,1.65);tx('SZ',82.5,15,7.1)
  ys <- seq(21,63,length.out=12)
  for(i in seq_len(nrow(b))) {
    tx(labels_b[i],41,ys[i],7.1,just='right',face=if(focus[i])'bold' else 'plain')
    ln(xp(b$NZ_mean[i]),ys[i],xp(b$SZ_mean[i]),ys[i],COL$line,1.4)
    pt(xp(b$NZ_mean[i]),ys[i],16,COL$nz,1.65)
    pt(xp(b$SZ_mean[i]),ys[i],15,COL$sz,1.65)
  }
  ln(xp(0),66,xp(1),66)
  for(v in c(0,.5,1)) {ln(xp(v),66,xp(v),67);tx(sprintf('%.1f',v),xp(v),69,7.2,just='center')}
  tx('Mean GMP',64.5,73.3,7.6,just='center')
}
panel_c <- function() {
  heading('c','Lineage contributions across all 12 modules','NZ \u2212 SZ; same module order as b')
  x0 <- 40.5; width <- 140.5; cw <- width/12; xs <- x0+(seq_len(12)-.5)*cw
  for(j in seq_len(12)) {
    tx(labels_c[j],xs[j],18.3,6.9,just='center',face=if(focus[j])'bold' else 'plain')
    if(focus[j]) rect(xs[j],25.2,cw-2,.6,COL$ink)
    rect(xs[j],29.1,cw,4.8,heat(totals[j]),COL$grid)
    tx(sprintf('%+.1f',totals[j]),xs[j],29.1,7.1,just='center',face='bold')
  }
  tx('Total NZ \u2212 SZ',38.6,29.1,7.5,just='right',face='bold')
  ys <- seq(35.1,74.1,length.out=11)
  for(i in seq_len(11)) {
    tx(lineage[i],38.6,ys[i],7.5,just='right')
    for(j in seq_len(12)) {
      rect(xs[j],ys[i],cw,3.9,heat(vals[i,j]),COL$grid)
      if(abs(vals[i,j])>=.5) tx(sprintf('%+.1f',vals[i,j]),xs[j],ys[i],7.1,just='center')
    }
  }
  ln(x0,32.6,x0+width,32.6,COL$muted,.5)
  lx <- 79; lw <- 55; n <- 300
  for(k in seq_len(n)) rect(lx+(k-.5)*lw/n,81,lw/n+.01,2.6,heat(-30+(k-.5)*60/n))
  for(v in c(-30,-15,0,15,30)) {xx<-lx+(v+30)/60*lw;ln(xx,82.4,xx,83,COL$muted,.4);tx(as.character(v),xx,85.2,6.8,just='center')}
  tx('SZ higher',lx-2,81,7.1,just='right',col=COL$muted)
  tx('NZ higher',lx+lw+2,81,7.1,col=COL$muted)
  tx('Contribution (GMP percentage points)',106.5,90.2,7.3,just='center')
}
panel_d <- function() {
  heading('d','Salinity conditioning')
  pt(127,3.5,16,COL$accent,1.9);tx('Baseline',130,3.5,7.3)
  pt(153.5,3.5,17,COL$salt,2.3);tx('+ salinity',157,3.5,7.3)
  for(g in seq_len(2)) {
    xleft <- if(g==1) 29 else 123
    xp <- function(x) xleft+x/.22*54
    dd <- d[d$catalogue==c(269,414)[g],]
    tx(c('269-MAG primary','414-MAG sensitivity')[g],xleft+27,10.7,8,just='center',face='bold')
    ys <- c(20.4,30)
    for(i in seq_len(2)) {
      tx(dd$response[i],xleft-3,ys[i],7.8,just='right')
      ln(xp(dd$conditioned_R2[i]),ys[i],xp(dd$baseline_R2[i]),ys[i],COL$line,1.7)
      pt(xp(dd$baseline_R2[i]),ys[i],16,COL$accent,2.1)
      pt(xp(dd$conditioned_R2[i]),ys[i],17,COL$salt,2.5)
      tx(sprintf('\u2212%.1f%%',dd$attenuation_pct[i]),mean(c(xp(dd$baseline_R2[i]),xp(dd$conditioned_R2[i]))),ys[i]-3.8,7.4,just='center',col=COL$muted)
    }
    ln(xp(0),33.6,xp(.22),33.6)
    for(v in c(0,.1,.2)) {ln(xp(v),33.6,xp(v),34.7);tx(sprintf('%.1f',v),xp(v),36.7,7.1,just='center')}
    tx('Basin term-specific R\u00b2',xleft+27,41.2,7.5,just='center')
  }
}
panel_fns <- list(a=panel_a,b=panel_b,c=panel_c,d=panel_d)
boxes <- list(a=c(0,0,89,75), b=c(94,0,89,75), c=c(0,78,183,93), d=c(0,177,183,44))
draw_panel <- function(letter,w,h) {
  pushViewport(viewport(xscale=c(0,w),yscale=c(h,0)))
  panel_fns[[letter]]();popViewport()
}
draw_all <- function() {
  grid.newpage()
  for(letter in names(boxes)) {
    bb <- boxes[[letter]]
    pushViewport(viewport(x=unit(bb[1],'mm'),y=unit(222-bb[2],'mm'),width=unit(bb[3],'mm'),height=unit(bb[4],'mm'),just=c('left','top')))
    draw_panel(letter,bb[3],bb[4]);popViewport()
  }
  pushViewport(viewport(xscale=c(0,183),yscale=c(222,0)))
  tx('Bold labels / black bars: four modules discussed in the text; not statistical significance.',1,173.8,7.2,col=COL$muted)
  popViewport()
}
figdir <- file.path(root,'figures'); dir.create(file.path(figdir,'panels'),showWarnings=FALSE)
for(letter in names(boxes)) {
  bb<-boxes[[letter]]
  cairo_pdf(file.path(figdir,'panels',paste0('Figure3',letter,'.pdf')),width=bb[3]/25.4,height=bb[4]/25.4,family='Arial')
  grid.newpage();draw_panel(letter,bb[3],bb[4]);dev.off()
}
stem<-file.path(figdir,'Figure_3_full12')
svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=222/25.4);draw_all();dev.off()
cairo_pdf(paste0(stem,'.pdf'),width=183/25.4,height=222/25.4,family='Arial');draw_all();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=183,height=222,units='mm',res=300,background='white');draw_all();dev.off()
ragg::agg_tiff(paste0(stem,'.tif'),width=183,height=222,units='mm',res=600,compression='lzw',background='white');draw_all();dev.off()
ragg::agg_png(file.path(root,'qa','Figure_3_preview.png'),width=183,height=222,units='mm',res=300,background='white');draw_all();dev.off()
used <- do.call(rbind,lapply(seq_len(11), function(i) data.frame(lineage=lineage[i],feature=mods,contribution_pp=vals[i,],highlighted_in_text=focus)))
write.table(used,file.path(root,'source_data','Figure3c_plotted_values.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
jsonlite::write_json(list(backend=R.version.string,dimensions_mm=c(183,222),source_rows=list(a=7,b=12,c=11,d=4),total_C_cells=132,omitted_C_cells=0,max_abs_sum_difference=max(abs(colSums(vals)-totals)),highlighted_features=mods[focus],models_rerun=FALSE),file.path(root,'qa','figure_numeric_checks.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
ptmm<-72/25.4
manifest<-list(schema_version=1,figure=list(width_pt=183*ptmm,height_pt=222*ptmm),panels=list(list(id='a',bbox_pt=c(0,147,89,222)*ptmm),list(id='b',bbox_pt=c(94,147,183,222)*ptmm),list(id='c',bbox_pt=c(0,51,183,144)*ptmm),list(id='d',bbox_pt=c(0,1,183,45)*ptmm)),row_groups=list(c('a','b')),column_groups=list(),exemptions=list(list(panels=c('c','d'),checks='all',reason='Full-width matrix and conditioning panels intentionally differ from the two-column top row.')))
jsonlite::write_json(manifest,file.path(root,'qa','panel_layout.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Exported full 12-column Figure 3. All 132 lineage cells retained.\n')
