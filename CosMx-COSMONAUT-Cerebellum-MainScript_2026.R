#COSMONAUT CEREBELLUM PIPELINE
#Adapted from Chase code starting 9/23/2024, 10/25/2024
#Author: James Tyr Redden

#nCount = number of genes measured by the assay, total RNAs measured
#nFeature = UNIQUE genes measured by the assay between cells/clusters, etc
#nFeatures parameter should be adjusted, keeping default may ignore important data
# test <- readRDS("C:/Users/m322332/Documents/PROTEIN-CEREBELLUM-SeuratREFINED.RDS")
# Working Directory Setup ####
project_WD <- '//mfad/researchmn/PHYSIO/SCHAFER/SHARED/JamesRedden/R Code/CEREBELLUM/'
local_WD <- "C:/Users/m322332/Documents/CosMx Cerebellum/"
HDD_WD <- "D:/"

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Install Packages & Initialize ####
LibStart_Seurat <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  # Added 'conflicted' to manage namespace collisions explicitly
  pkgs <- c("Seurat", "sctransform", "ggplot2", "plotly", "dplyr", "sf",
            "RColorBrewer", "gridExtra", "grid", "presto", "patchwork",
            "svglite", "tcltk", "glmGamPoi", "qs", "FNN", "parallel", "RcppAnnoy",
            "data.table", "cluster", "RANN", "future", "crayon",
            "celldex", "SingleR", "circlize", "limma", "igraph",
            "gplots", "grDevices", "edgeR", "SIMLR", 
            "pheatmap", "stats", "Rtsne", "graphics", "stringr", "foreach",
            "multtest", "scran", "uwot", "Pando", "conflicted")
  
  pkgs <- unique(trimws(pkgs))
  options(repos = BiocManager::repositories())
  
  # Vectorized installation check
  installed <- rownames(installed.packages())
  missing_pkgs <- pkgs[!(pkgs %in% installed)]
  
  if(length(missing_pkgs) > 0) {
    message("Installing missing packages...")
    BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
  }
  
  # Load packages silently
  suppressWarnings(suppressPackageStartupMessages({
    for (pkg in pkgs) library(pkg, character.only = TRUE)
  }))
  
  # Explicitly resolve known namespace conflicts to prevent silent pipeline failures
  conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("rename", "dplyr", quiet = TRUE)
  
  message("--- All Seurat libraries loaded successfully ---")
}
LibStart_Seurat()
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Parallelization Options ####
Parall_Space <- 20 # Gb
options(future.globals.maxSize = Parall_Space * 1024^3) # Exact GiB calculation

# OS-aware future planning to prevent silent single-core fallback on Windows
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = 12) 
} else {
  plan("multicore", workers = 12) # for WSL/Ubuntu environments
}
options(future.seed = TRUE)

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Gene Marker Lists ####
#**Begin Gene Lists*
#formerly named senBrainList == SenBrain_Mouse_2025 in this script
SenBrain_Mouse_2025 = c(
  "Akt1","Apoe","App","Axl","Bcl2","Bcl2l1","Ccl2","Ccl3","Ccl4","Ccl5","Ccne1",
  "Cd274","Cd36","Cd74","Cd9","Cdc27","Cdc6","Cdk4","Cdk8","Cdk9","Cdkn1a","Cdkn1b",
  "Cdkn1c","Cdkn2a","Cdkn2b","Cdkn2c","Cdkn2d","Cenpb","Csf1","Csf1r","Cxcl14","Cxcl2",
  "Cyba","Dpp4","Elmo1","Fxyd5","Gadd45a","Gas2","Gdf15","Glb1","Gpnmb","H2-Aa",
  "H2-Ab1","Hmgb1","Ifitm3","Igfbp7","Il10","Il18","Il1a","Il1b","Il23a","Il23r",
  "Il33","Il6","Irf5","Irf7","Lamp2","Lgals3","Lgals3bp","Lmna","Lpl","Malat1","Mapk1",
  "Mapk10","Mapk14","Mapk8","Meox1","Mmp12","Mmp14","Mmp2","Mmp3","Mmp9","Nfkb1",
  "Nol3","Notch3","Pak1","Pak2","Pak3","Pdcd1","Pik3ca","Pik3cb","Pik3r1","Plaur",
  "Pten","Rbl2","Serpina3n","Serpinb1a","Serpine1","Serpine2","Smad2","Smad3","Sod1","Sparc","Spp1",
  "Tgfb1","Tgfb2","Tgfbr1","Sting1","Tnf","Tp53bp1","Trp53",
  "Vim"
)
#1
cellCycleDomain = c("Cdk4", "Cdk8","Cdk9", "Cdc6", "Cdc27",
                    "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                    "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19", 
                    "Tp53bp1", "Trp53", "Jun","Rb1","Rbl2")
#2&3
nuclearStressDomain = c("Cenpb","Gadd45a","Hmgb1","Irf7","Mapk14","Lmna","Nfkb1",
                        "Pak1","Pak2","Pak3","Phc2","Phc3","Smad2","Smad3","Tp53bp1","Trp53","Xcl1")
#4
scapDomain = c("Akt1","Akt3","Bcl2", "Bcl2l1","Clu","Casp3","Pik3ca", "Pik3cb", "Pik3r1", "Serpina3n", "Serpinb1a","Serpine1","Serpine2")
#5
saspDomain = c("Ccl2","Ccl3","Ccl4","Ccl5","Csf1","Cd74","Cd81","Cxcl2","Cxcl14",
               "Il1a","Il1b","Il6","Il10","Il23a","Il33","Gdf15","Mmp12", "Mmp14", "Mmp2", "Mmp3", "Mmp9",
               "Spp1","Tgfb1","Tgfb2","Tnf","Ifng","Irf5","Irf7")
#6&7 Lysomal stress, metabolic adaptations, mitochondrial ROS
organelleStressDomain = c("Cyba","Fth1","Glb1","Hmgb1","Hspa1a","Mapk1", "Mapk8", "Mapk10", "Mapk14","Lamp2","Lgals3","Lgals3bp","Sirt1","Sod1")
#8 Cell surface markers
cellSurfaceDomain = c("Cd9","CD26","Cd274","Csf1r","Gpnmb","H2-Aa","H2-Ab1","Ifitm3","Igfbp7","Il23r","Notch3","Plaur","Tgfbr1","Tmem173","Vim")
#0 Aged brain upregulated
ageSelectiveDomain = c("Apoe","App","Clu","Malat1","Lpl")

nanostringSenQuiet = c("Ago3","Akt1","Akt3","Anapc16","Atm","Atr","Btrc","Cacna1d",
                       "Calm1","Calm2","Calm3","Cdc27","Cdk6","Fbxw11","Foxo1","Foxo3",
                       "H3f3b","Hipk2","Igfbp7","Il6","Itpr1","Itpr2","Jun","Map2k1","Map2k4","Map2k5",
                       "Map2k6","Map3k5","Map4k4","Mapk1","Mapk10","Mapk14","Mapk8",
                       "Mcu","Nfkb1","Phc2","Phc3","Pik3ca","Pik3cb","Pik3r1","Ppp3ca","Ppp3cb","Ppp3cc",
                       "Pten","Raf1","Rb1","Rps6ka2","Rps6ka3","Slc25a4","Smad2","Smad3","Stat3",
                       "Terf2ip","Tfdp2","Tgfb1","Tgfb2","Tgfbr1","Tnik","Tnrc6a","Tnrc6b","Tnrc6c",
                       "Trp53","Trpm7","Ube2e1")

damList = c('Apoe', 'B2m', 'C1qa','Ctsb', 'Ctsd', 'Ctss', 'Csf1','Csf1r', 'Cst3', 'Cst7', 'Cx3cr1', "Cxcl14",
            'Lgals3', 'Ftl1', 'Fth1', 'Clec7a', 'Cd74', 'H2-Aa','H2-Ab1', 'Itgax', 'Lyz1/2','Tgfbr1','Tmem173','Trem2','Tyrobp','Spp1')
damList <- sort(damList)
mgList = c('C1qa','C1qb','C1qc','Ctss','Csf1r','Cst3','Cx3cr1','Hexb','Itgam','P2ry12','Selplg','Tmem119')

cdknList = c("Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
             "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19")

#Cd9 low, Ctsd low
#macrophage associated list
macList = c('C2','C3','C4b','Crip1','Cd74','Lyz1/2','H2-Ab1','H2-Aa','Hsp90ab1',
            'Ldha','Pltp','Bsg','Hmgb1','H3f3b',
            'Ifitm3','Fxyd5','Vim','Ctsd','Cd9','Cd63','Cd81')

#Cd81 low, Ptprc and Tgfb1 too prevalent to include
immuneinfiltrateList = c('Ccl5', 'Cd3d', 'Cd3e', 'Cd3g', 'Cd44','Cd46',
                         'Cd74', 'Cd274', 'Crip1', 'Fau','H2-Aa', 'H2-Ab1',
                         'Ifitm3', 'Il6', 'Il23a', 'Il33', 'Klre1', 'Ltbp1',
                         'Ms4a4a', 'Nfkb1', 'Nkg7', 'Pdcd1', 'Tmsb10')

astroList = c('Gja1','Gpr37l1','Plpp3','Mfge8','Gfap','Slc1a3','Bcan')
gabaergicList = c('Gad1','Gad2','Npy','Slc6a1','Slc32a1','Sst')
neuronList = c('Adcy1','Atp1b1','Camk2a','Dnm1','Epha6','Gria1', 'Grin2a','Grin2b','Nell2','Nrgn','Snap25','Syp','Syt1','Ywhag')
oligoList = c('Plp1','Mog','Mag','Mobp','Mbp','Cryab','Fa2h')
opcList = c('Cntn1','Cspg5','Gpr17','Pdgfra', 'Ptprz1','Serpine2')
vascList = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg')
vlmcList = c('Apod','Col1a1','Dcn','Il33','Ptgds','Ptn')
epenList = c('Dbi','Foxj1','Nnat','Tuba1a/b/c','Pltp','Vim') #Pdgfrb also in choroid
periList = c('Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn')
nstemnblastList = c('Ctnna2','H3f3b','Nnat','Ptprs','Tcf4','Tuba1a/b/c','Tubb5')
chorplexList = c('Ttr','Htr2c','Lgals3bp','Cab39l','Igfbp7','Ppp1r1b')

vascList2 = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg','Apod','Col1a1','Dcn','Il33','Ptgds','Ptn',
             'Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn')

GAL3enrichedModule = c('Apoe','Cd9','Clec7a','Ctsb','Ctsd','Fth1','Gpnmb','Itgax',
                       'Lyz1/2','Lilrb4a/b','Mmp12','Spp1','Tmsb4x','Tyrobp')

SenGeneSets <- list(
  #Mouse Gene Set: SenBrain
  #Schafer Lab Mayo Clinic 2025
  SenBrain_Mouse_2025 = c(
    "Akt1","Apoe","App","Axl","Bcl2","Bcl2l1","Ccl2","Ccl3","Ccl4","Ccl5","Ccne1",
    "Cd274","Cd36","Cd74","Cd9","Cdc27","Cdc6","Cdk4","Cdk8","Cdk9","Cdkn1a","Cdkn1b",
    "Cdkn1c","Cdkn2a","Cdkn2b","Cdkn2c","Cdkn2d","Cenpb","Csf1","Csf1r","Cxcl14","Cxcl2",
    "Cyba","Dpp4","Elmo1","Fxyd5","Gadd45a","Gas2","Gdf15","Glb1","Gpnmb","H2-Aa",
    "H2-Ab1","Hmgb1","Ifitm3","Igfbp7","Il10","Il18","Il1a","Il1b","Il23a","Il23r",
    "Il33","Il6","Irf5","Irf7","Lamp2","Lgals3","Lgals3bp","Lmna","Lpl","Malat1","Mapk1",
    "Mapk10","Mapk14","Mapk8","Meox1","Mmp12","Mmp14","Mmp2","Mmp3","Mmp9","Nfkb1",
    "Nol3","Notch3","Pak1","Pak2","Pak3","Pdcd1","Pik3ca","Pik3cb","Pik3r1","Plaur",
    "Pten","Rbl2","Serpina3n","Serpinb1a","Serpine1","Serpine2","Smad2","Smad3","Sod1","Sparc","Spp1",
    "Tgfb1","Tgfb2","Tgfbr1","Sting1","Tnf","Tp53bp1","Trp53",
    "Vim"
  ),
  
  #Mouse Gene Set: SenMayo
  #https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/SAUL_SEN_MAYO.html
  SenMayo_Mouse = c("Acvr1b","Ang","Angpt1","Angptl4","Areg","Axl","Bex3","Bmp2","Bmp6","C3","Ccl1","Ccl2","Ccl20","Ccl24",
                    "Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Cd55","Cd9","Csf1","Csf2","Csf2rb","Ctnnb1","Ctsb","Cxcl1",
                    "Cxcl10","Cxcl12","Cxcl16","Cxcl2","Cxcl3","Cxcr2","Dkk1","Edn1","Egf","Egfr","Ereg","Esm1","Ets2","Fas",
                    "Fgf1","Fgf2","Fgf7","Gdf15","Gem","Gmfg","Hgf","Hmgb1","Icam1","Icam5","Igf1","Igfbp1","Igfbp2",'Igfbp3',
                    "Igfbp4","Igfbp5","Igfbp6","Igfbp7","Il10","Il13","Il15","Il18","Il1a","Il1b","Il2","Il6","Il6st","Il7",
                    "Inha","Iqgap2","Itga2","Itpka","Jun","Kitl","Lcp1","Mif","Mmp13","Mmp10","Mmp12","Mmp13","Mmp14","Mmp2",
                    "Mmp3","Mmp9","Nap1l4","Nrg1","Pappa","Pecam1","Pgf","Pigf","Plat","Plau","Plaur","Ptbp1","Ptger2","Ptges",
                    "Rps6ka5","Scamp4","Selplg","Sema3f","Serpinb3c","Serpine1","Serpine2","Spp1","Spx","Timp2","Tnf",
                    "Tnfrsf11b","Tnfrsf1a","Tnfrsf1b","Tubgcp2","Vegfa","Vegfc","Vgf","Wnt16","Wnt2"),
  
  #Mouse Gene Set: REACTOME_CELLULAR_SENESCENCE
  #https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_CELLULAR_SENESCENCE.html
  Reactome_Cellular_Senescence = c('Mapk7', 'Ube2c', 'Ccne1', 'Map2k7', 'Cdkn1b', 'Rps6ka1', 
                                   'Cdk4', 'Ubc', 'Ehmt2', 'Anapc1', 'Mapkapk2', 'Suz12', 
                                   'H2bc4', 'Kdm6b', 'Map2k3', 'Ubb', 'Asf1a', 'Ube2d1', 
                                   'Anapc16', 'Mdm2', 'Cabin1', 'Fzr1', 'Mapk9', 'Rad50', 
                                   'Rps27a', 'Map2k6', 'Cdc27', 'Mink1', 'Fos', 'Ube2e1', 
                                   'Mapk8', 'Rb1', 'Hira', 'Ets2', 'Cdkn1a', 'Rps6ka2', 
                                   'Cdc23', 'Lmnb1', 'Kat5', 'Anapc11', 'Cdk2', 'Terf1', 
                                   'Map4k4', 'Anapc2', 'Tnik', 'Ccna2', 'Ccna1', 'Ccne2', 
                                   'Nbn', 'Txn1', 'Anapc4', 'Mapkapk5', 'Anapc7', 'Anapc5', 
                                   'Ep400', 'Pot1a', 'Ezh2', 'Eed', 'Anapc15', 'Rps6ka3', 
                                   'Rbbp7', 'Terf2', 'Mre11a', 'Ets1', 'Mapkapk3', 'Map2k4', 
                                   'Terf2ip', 'Atm', 'H1f2', 'Ehmt1', 'Anapc10', 'Cdkn1c', 
                                   'Acd', 'Cdc16', 'Ubn1', 'Cdk6', 'Erf', 'H2az2', 
                                   'Id1', 'Mapk10', 'Hmga1', 'H2bc6', 'H1f1', 'H2ax', 
                                   'H2bc1', 'H1f4', 'Jun', 'Mapk11', 'Mapk14', 'Mdm4', 
                                   'Cebpb', 'Hmga2', 'Rbbp4', 'H2bc8', 'H1f5', 'Trp53', 
                                   'H2aj', 'H4c1', 'H4c9', 'H4c3', 'Ube2s', 'H4c8', 
                                   'H4c4', 'H2ac4', 'H2ac10', 'H2bc12', 'H2ac15', 'Mapk3', 
                                   'Mapk1', 'H2ac19', 'H2ac18', 'H4c12', 'Cdc26', 'H2ab1', 
                                   'H4c11', 'Uba52rt', 'H2bc21', 'H2ac20', 'H4c2', 'H2bc7', 
                                   'H2ac6', 'H2ac8', 'H4c6', 'H2bc11', 'H2ac11', 'H2ac12', 
                                   'H2bc24', 'H4c18', 'H4c17', 'H2bc23', 'H2bc22', 'H2ac22', 
                                   'Map3k5', 'H2ac7', 'H2ac13', 'Cdkn2b', 'H2bc3', 'Hmga1b', 
                                   'H3f4', 'H2bc26', 'H2ab2', 'H2ab3', 'Uba52', 'H4c14', 
                                   'H2ac23', 'H2bc13', 'H2ac24', 'H2bc15', 'H4c16', 'H1f0', 
                                   'H2bc14', 'H2bc9'),
  
  #Mouse Gene Set: HALLMARK_P53_PATHWAY
  #https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_P53_PATHWAY.html
  Hallmark_P53_Pathway = c('Cdk5r1', 'F2r', 'Rack1', 'Ptpn14', 'Ccng1', 'Dgka', 
                           'Vdr', 'Zbtb16', 'Xpc', 'Ercc5', 'Cdkn1a', 'Procr', 
                           'Cd82', 'Tm4sf1', 'Fdxr', 'Cdkn2b', 'Cdkn2a', 'Ier3', 
                           'Mxd4', 'Tnni1', 'Tspyl2', 'Gpx2', 'Irak1', 'Hspa4l', 
                           'Ralgds', 'Gadd45a', 'Zfp36l1', 'Ei24', 'Btg2', 'Rgs16', 
                           'Irag2', 'Tsc22d1', 'Fgf13', 'Ddit3', 'Atf3', 'Serpinb5', 
                           'Plk3', 'Cd81', 'Bak1', 'Upp1', 'Alox8', 'Kif13b', 
                           'Jag2', 'Plk2', 'Tnfsf9', 'Zmat3', 'Ninj1', 'Traf4', 
                           'Ccnk', 'Tpd52l1', 'Apaf1', 'Sphk1', 'Prmt2', 'Hint1', 
                           'Rad9a', 'S100a4', 'Trp63', 'Prkab1', 'Vamp8', 'Ier5', 
                           'St14', 'Acvr1b', 'S100a10', 'Ndrg1', 'Klf4', 'Hdac3', 
                           'Klk8', 'Trib3', 'Def6', 'Slc7a11', 'Sdc1', 'Tob1', 
                           'Phlda3', 'Abcc5', 'Pmm1', 'Fbxw7', 'Socs1', 'Ddb2', 
                           'Ppm1d', 'Sec61a1', 'Blcap', 'Rpl36', 'Ctsf', 'Wwp1', 
                           'Pidd1', 'Txnip', 'Foxo3', 'Polh', 'Wrap73', 'Sfn', 
                           'Nupr1', 'Iscu', 'Rrp8', 'Elp1', 'Coq8a', 'Retsat', 
                           'Tm7sf3', 'Vwa5a', 'Dcxr', 'Rps27l', 'Aen', 'Rchy1', 
                           'Steap3', 'Ankra2', 'Abhd4', 'Triap1', 'Cgrrf1', 'Tprkb', 
                           'Fam162a', 'Cdkn2aip', 'Nol8', 'Dram1', 'Pitpnc1', 'Osgin1', 
                           'Rap2b', 'Ddit4', 'Rnf19b', 'Dnttip2', 'Tax1bp3', 'Trafd1', 
                           'Ip6k2', 'Cyfip2', 'Ppp1r15a', 'Slc19a2', 'Perp', 'Rrad', 
                           'H1f2', 'Baiap2', 'Pom121', 'Ifi30', 'Eps8l2', 'Clca2', 
                           'Slc35d1', 'Ccp110', 'Mapkapk3', 'Gls2', 'Zfp365', 'Rad51c', 
                           'Plxnb2', 'Sesn1', 'Tchh', 'Sertad3', 'Rab40c', 'Hexim1', 
                           'Csrnp2', 'Rhbdf2', 'Nudt15', 'Abat', 'H2ac25', 'H2aj', 
                           'Ada', 'Ak1', 'App', 'Bmp2', 'Btg1', 'Ccnd2', 
                           'Ccnd3', 'Ctsd', 'Mknk2', 'Epha2', 'Stom', 'Ephx1', 
                           'Fas', 'Fos', 'Fuca1', 'Gm2a', 'Hbegf', 'Hmox1', 
                           'Hras', 'Il1a', 'Casp1', 'Inhbb', 'Itgb4', 'Jun', 
                           'Krt17', 'Ldhb', 'Lif', 'Mxd1', 'Mdm2', 'Slc3a2', 
                           'Nhlh2', 'Notch1', 'Pcna', 'Pdgfa', 'Ptpre', 'Pvt1', 
                           'Rb1', 'Rpl18', 'Rps12', 'Rxra', 'Sat1', 'Sp1', 
                           'Tap1', 'Tcn2', 'Tgfa', 'Tgfb1', 'Trp53', 'Cebpa', 
                           'Cdh13', 'Bax'),
  
  #Mouse Gene Set: REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP
  #https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP.html
  Reactome_SASP = c('Mapk7', 'Ube2c', 'Cdkn1b', 'Rps6ka1', 'Cdk4', 'Ubc', 
                    'Ehmt2', 'Anapc1', 'Ubb', 'Ube2d1', 'Anapc16', 'Fzr1', 
                    'Rps27a', 'Cdc27', 'Ube2e1', 'Cdkn1a', 'Rps6ka2', 'Cdc23', 
                    'Anapc11', 'Cdk2', 'Anapc2', 'Ccna2', 'Ccna1', 'Anapc4', 
                    'Anapc7', 'Anapc5', 'Anapc15', 'Rps6ka3', 'Ehmt1', 'Anapc10', 
                    'Cdkn1c', 'Cdc16', 'Cdk6', 'Cebpb', 'Ube2s', 'Mapk3', 
                    'Mapk1', 'Cdc26', 'Uba52rt', 'Cdkn2b', 'Uba52'),
  
  #Mouse Gene Set: GOBP_CELLULAR_SENESCENCE
  #https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/GOBP_CELLULAR_SENESCENCE.html
  GO_Cellular_Senescence = c('Fbxo4', 'Hmga1b', 'Abl1', 'Pawr', 'Prelp', 'Ang', 
                             'Arg2', 'Bmal1', 'B2m', 'Bcl6', 'Bmpr1a', 'Brca2', 
                             'Calr', 'Cav1', 'Cdk6', 'Cdkn1a', 'Cdkn2a', 'Cdkn2b', 
                             'Comp', 'Gch1', 'H2-M3', 'Htt', 'H2ax', 'Hmga1', 
                             'Hmga2', 'Hras', 'Id2', 'Ilk', 'Kras', 'Lmna', 
                             'Mapkapk5', 'Mif', 'Mnt', 'Cited2', 'Npm1', 'Nup62', 
                             'Prkcd', 'Pla2r1', 'Pml', 'Prkdc', 'Pten', 'Rbl1', 
                             'Spi1', 'Plk2', 'Srf', 'Suv39h1', 'Tbx2', 'Tbx3', 
                             'Cgas', 'Terf2', 'Tert', 'Top2b', 'Trp53', 'Trp63', 
                             'Twist1', 'Wnt1', 'Wrn', 'Smc5', 'Ybx1', 'Zfp217', 
                             'Zmpste24', 'H2aj', 'Akt3', 'Vash1', 'Nek4', 'Kat6a', 
                             'Zfp277', 'Map2k1', 'Map3k3', 'Mapk14', 'Zmiz1', 'Morc3', 
                             'Mirlet7d', 'Mirlet7g', 'Mirlet7f-1', 'Mirlet7f-2', 'Ndufs6', 'Sirt6', 
                             'Fzr1', 'Mtor', 'Icmt', 'Nampt', 'Ypel3', 'Eef1e1', 
                             'Rsl1d1', 'Abi3', 'Fbxo5', 'Smc6', 'Nsmce2', 'Pnpt1', 
                             'Zkscan3', 'Pot1b', 'Opa1', 'Bcl2l12', 'Nuak1', 'Gins3', 
                             'Ecrg4', 'Dnaja3', 'Wnt16', 'Sirt1', 'Prmt6'),
  
  #"SenSig" Genes
  #PMID: 37079217, padj<0.05, |LFC| > 1
  SenSig_padj0.05 = c("Apoe","App","Bcl2l1","Ccl4","Cd274","Cdkn1a","Cdkn2b","Cenpb","Csf1r","Cxcl2",
                      "Dpp4","Fxyd5","Gadd45a","Gas2","Igfbp7","Il10","Il18","Il6","Irf7","Lgals3","Lmna","Lpl","Mmp3","Mmp9",
                      "Notch3","Pik3r1", "Plaur","Pten","Serpina3n", "Serpinb1a","Serpine1","Serpine2",
                      "Sparc","Tgfb1","Tgfb2")
)


# Load Raw Seurat & Pre-filter ####
#**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
#'[Load Raw Seurat & Pre-filter]
#**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
# --- 1. Load DAta
thisObj <- readRDS(file.path(HDD_WD, "ead2aca2-438d-4d75-bc0e-a077ddd96c27_seuratObject.RDS"))

thisObj@meta.data$x_global_px <- thisObj@meta.data$x_slide_mm * 8333.3333
thisObj@meta.data$y_global_px <- thisObj@meta.data$y_slide_mm * 8333.3333

#CEREBELLUM FOV-id_sets
CB.YF01 = c(28:43)
CB.YF02 = c(31:50)
CB.YF03 = c(26:38)
CB.YF04 = c(28:48)
CB.236 = c(70:83)
CB.237 = c(77:95)
CB.239 = c(65:83)
CB.243 = c(75:94)

CB.YM01 = c(28:51)
CB.YM02 = c(30:37)
CB.YM03 = c(35:54)
CB.YM04 = c(35:55)
CB.301 = c(77:93)
CB.302 = c(71:88)
CB.304 = c(87:106)
CB.305 = c(96:116)

#definitions for "Run_Tissue_name" metadata for each slide
slideF1 = c("05_Schafer_F1_236_YF01")
slideF2 = c("06_Schafer_F2_237_YF02")
slideF3 = c("07_Schafer_F3_239_YF03")
slideF4 = c("08_Schafer_F4_243_YF04")
slideM1 = c("01_Schafer_M1_301_YM01")
slideM2 = c("02_Schafer_M2_302_YM02")
slideM3 = c("03_Schafer_M3_304_YM03")
slideM4 = c("04_Schafer_M4_305_YM04")

maleSlides <- c(slideM1, slideM2, slideM3, slideM4)
femSlides <- c(slideF1, slideF2, slideF3, slideF4)
yngMice <- c("YM01","YM02","YM03","YM04","YF01","YF02","YF03","YF04")
oldMice <- c("301","302","304","305","236","237","239","243")

thisObj@meta.data$mus <- dplyr::case_when(
  thisObj$Run_Tissue_name == slideM1 & thisObj$fov %in% CB.YM01 ~ "YM01",
  thisObj$Run_Tissue_name == slideM1 & thisObj$fov %in% CB.301  ~ "301",
  thisObj$Run_Tissue_name == slideM2 & thisObj$fov %in% CB.YM02 ~ "YM02",
  thisObj$Run_Tissue_name == slideM2 & thisObj$fov %in% CB.302  ~ "302",
  thisObj$Run_Tissue_name == slideM3 & thisObj$fov %in% CB.YM03 ~ "YM03",
  thisObj$Run_Tissue_name == slideM3 & thisObj$fov %in% CB.304  ~ "304",
  thisObj$Run_Tissue_name == slideM4 & thisObj$fov %in% CB.YM04 ~ "YM04",
  thisObj$Run_Tissue_name == slideM4 & thisObj$fov %in% CB.305  ~ "305",
  thisObj$Run_Tissue_name == slideF1 & thisObj$fov %in% CB.YF01 ~ "YF01",
  thisObj$Run_Tissue_name == slideF1 & thisObj$fov %in% CB.236  ~ "236",
  thisObj$Run_Tissue_name == slideF2 & thisObj$fov %in% CB.YF02 ~ "YF02",
  thisObj$Run_Tissue_name == slideF2 & thisObj$fov %in% CB.237  ~ "237",
  thisObj$Run_Tissue_name == slideF3 & thisObj$fov %in% CB.YF03 ~ "YF03",
  thisObj$Run_Tissue_name == slideF3 & thisObj$fov %in% CB.239  ~ "239",
  thisObj$Run_Tissue_name == slideF4 & thisObj$fov %in% CB.YF04 ~ "YF04",
  thisObj$Run_Tissue_name == slideF4 & thisObj$fov %in% CB.243  ~ "243",
  TRUE ~ "unk"
)

# Initial Regional Subset: Drop non-cerebellum cells
thisObj <- subset(thisObj, subset = mus != "unk")

thisObj@meta.data$sex <- dplyr::case_when(
  thisObj$Run_Tissue_name %in% maleSlides ~ "M",
  thisObj$Run_Tissue_name %in% femSlides ~ "F",
  TRUE ~ NA_character_
)
thisObj@meta.data$age <- dplyr::case_when(
  thisObj$mus %in% yngMice ~ "YNG",
  thisObj$mus %in% oldMice ~ "OLD",
  TRUE ~ NA_character_
)
thisObj@meta.data$region <- 'cerebellum'
thisObj@meta.data$agecomp <- as.integer(thisObj@meta.data$age == "OLD")
thisObj@meta.data$sexcomp <- as.integer(thisObj@meta.data$sex == "M")


# --- 2. Extract Subcellular Metrics Safely ---
library(data.table)
message("Extracting Nuclear fractions...")

coords_dt <- setDT(thisObj@misc[["transcriptCoords"]])
valid_cells <- colnames(thisObj)
coords_dt <- coords_dt[cell_id %in% valid_cells] 
colnames(coords_dt)

#'[Note: takes a long time to run, process plots overnight]
# library(data.table)
# 
# # --- 1. Reverse-Engineer FOV Global Origins ---
# # Convert metadata to data.table
# meta_dt <- as.data.table(thisObj@meta.data, keep.rownames = "cell_id")
# 
# # NOTE: Replace 'CenterX_local_px' and 'CenterY_local_px' with the exact 
# # column names in your metadata that store the cell's local FOV position.
# # We only need one cell per FOV to find the origin offset, so we take the first [1].
# fov_offsets <- meta_dt[, .(
#   fov_origin_x = x_global_px[1] - CenterX_local_px[1],
#   fov_origin_y = y_global_px[1] - CenterY_local_px[1]
# ), by = .(mus, fov)]
# 
# 
# # --- 2. Extract and Prepare Transcripts ---
# coords_dt <- setDT(thisObj@misc[["transcriptCoords"]])
# valid_cells <- colnames(thisObj)
# coords_dt <- coords_dt[cell_id %in% valid_cells & codeclass == "Endogenous"]
# 
# # The transcript table natively has 'fov', but not 'mus'. Merge 'mus' from cells.
# coords_dt <- merge(coords_dt, meta_dt[, .(cell_id, mus, age)], by = "cell_id", all.x = TRUE)
# 
# 
# # --- 3. Compute True Global Transcript Coordinates ---
# # Map the calculated FOV origins onto the transcripts
# plot_dt <- merge(coords_dt, fov_offsets, by = c("mus", "fov"), all.x = TRUE)
# 
# # Calculate true spatial position: Transcript Global = Transcript FOV + FOV Origin
# plot_dt[, tx_global_px := x_FOV_px + fov_origin_x]
# plot_dt[, ty_global_px := y_FOV_px + fov_origin_y]
# 
# 
# # --- 4. Plot True Spatial Transcript Density ---
# p_transcript_density <- ggplot(plot_dt, aes(x = tx_global_px, y = ty_global_px)) +
#   geom_bin2d(aes(fill = after_stat(log10(count))), bins = 150) + 
#   scale_fill_viridis_c(option = "inferno", name = "log10(Transcripts)") +
#   facet_wrap(~ mus, ncol = 4, scales = "free") + 
#   ggtitle("True Spatial Transcript Density by Sample (Cerebellum)") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     strip.text = element_text(face = "bold", size = 10),
#     legend.position = "bottom"
#   )
# 
# print(p_transcript_density)
# 
# ggsave(filename = paste0(local_WD, "CosMx_CB_TranscriptDensity.png"), 
#        plot = p_transcript_density, 
#        width = 16, height = 12, dpi = 600)

nuc_summary <- coords_dt[codeclass == "Endogenous", 
                         .(nCount_nuclear_raw = sum(CellComp == "Nuclear", na.rm = TRUE),
                           nCount_total_raw = .N), 
                         by = cell_id]

# Secure Metadata Injection using AddMetaData (Prevents Rowname Corruption)
nuc_df <- as.data.frame(nuc_summary)
rownames(nuc_df) <- nuc_df$cell_id
thisObj <- AddMetaData(thisObj, metadata = nuc_df)

# Fill NAs mathematically within Seurat
thisObj$nCount_nuclear_raw[is.na(thisObj$nCount_nuclear_raw)] <- 0
thisObj$nCount_total_raw[is.na(thisObj$nCount_total_raw)] <- 1 # Prevents division by zero
thisObj$nuclear_fraction <- thisObj$nCount_nuclear_raw / thisObj$nCount_total_raw

rm(coords_dt, nuc_summary, nuc_df, valid_cells)
thisObj@misc[["transcriptCoords"]] <- NULL 
gc()


# --- 3. Compute Composite QC Scores ---
message("Computing Transcript Dominance...")
counts_mat <- GetAssayData(thisObj, assay = "RNA", layer = "counts")
col_indices <- rep(1:ncol(counts_mat), diff(counts_mat@p))
dt_counts <- data.table(col = col_indices, val = counts_mat@x)
max_vals <- dt_counts[, .(max_val = max(val)), by = col]

max_gene_vector <- numeric(ncol(counts_mat))
max_gene_vector[max_vals$col] <- max_vals$max_val

thisObj$max_gene_count <- max_gene_vector
thisObj$top_gene_fraction <- thisObj$max_gene_count / thisObj$nCount_RNA
cat("Top gene fraction filter will remove:",
    sum(thisObj$top_gene_fraction >= 0.3, na.rm = TRUE), "cells\n")

rm(counts_mat, col_indices, dt_counts, max_vals, max_gene_vector)
gc()

thisObj$complexity <- log10(thisObj$nFeature_RNA) / log10(thisObj$nCount_RNA)
thisObj$falsecode_ratio <- thisObj$nCount_falsecode / thisObj$nCount_RNA

# --- 3. Compute Composite QC Scores ---
# Transcript Dominance (Top Gene Fraction)
counts_mat <- GetAssayData(thisObj, assay = "RNA", layer = "counts")

# OPTIMIZATION: Calculate max transcript directly from sparse matrix structure 
# Prevents the 1.9 GiB RAM spike and the Seurat metadata naming error.
col_indices <- rep(1:ncol(counts_mat), diff(counts_mat@p))
dt_counts <- data.table(col = col_indices, val = counts_mat@x)
max_vals <- dt_counts[, .(max_val = max(val)), by = col]

# Map back to a clean, unnamed numeric vector (forces exact positional alignment in Seurat)
max_gene_vector <- numeric(ncol(counts_mat))
max_gene_vector[max_vals$col] <- max_vals$max_val

thisObj$max_gene_count <- max_gene_vector
thisObj$top_gene_fraction <- thisObj$max_gene_count / thisObj$nCount_RNA

# Clean up sparse extraction variables
rm(counts_mat, col_indices, dt_counts, max_vals, max_gene_vector)
gc()

# Complexity & Falsecode Ratios
thisObj$complexity <- log10(thisObj$nFeature_RNA) / log10(thisObj$nCount_RNA)
thisObj$falsecode_ratio <- thisObj$nCount_falsecode / thisObj$nCount_RNA

# Collision Score
collision_genes <- c("Calb1", "Gabra6")
if (all(collision_genes %in% rownames(thisObj))) {
  expr_data <- FetchData(thisObj, vars = collision_genes, slot = "counts")
  thisObj$collision_Purk_Gran <- expr_data$Calb1 * expr_data$Gabra6
}

# --- 4. Plot Unified QC Metrics ---
qc_data <- thisObj@meta.data

# 1. Log-Space MAD Statistics for nCount_RNA
# Filter > 0 to prevent -Inf log transformations
log_counts <- log10(qc_data$nCount_RNA[qc_data$nCount_RNA > 0])
med_log_count <- median(log_counts)
mad_log_count <- mad(log_counts)

med_count <- 10^med_log_count
mad_counts_lower <- 10^(med_log_count - (1:3 * mad_log_count))
mad_counts_upper <- 10^(med_log_count + (1:3 * mad_log_count))

# 2. Log-Space MAD Statistics for nFeature_RNA
log_features <- log10(qc_data$nFeature_RNA[qc_data$nFeature_RNA > 0])
med_log_feat <- median(log_features)
mad_log_feat <- mad(log_features)

med_feat <- 10^med_log_feat
mad_feats_lower <- 10^(med_log_feat - (1:3 * mad_log_feat))
mad_feats_upper <- 10^(med_log_feat + (1:3 * mad_log_feat))

valid_counts <- thisObj$nCount_RNA[!is.na(thisObj$nCount_RNA) & thisObj$nCount_RNA > 0]
nCount_cutoff <- 10^(median(log10(valid_counts)) - (2 * mad(log10(valid_counts))))

valid_feats <- thisObj$nFeature_RNA[!is.na(thisObj$nFeature_RNA) & thisObj$nFeature_RNA > 0]
nFeature_cutoff <- 10^(median(log10(valid_feats)) - (2 * mad(log10(valid_feats))))

# 3. Augmented nCount_RNA Distribution
p_counts <- ggplot(qc_data, aes(x = nCount_RNA)) +
  geom_density(fill = "lightgreen", alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  # Add MAD lines (Lighter to darker for 3, 2, 1 MAD)
  geom_vline(xintercept = c(mad_counts_lower[3], mad_counts_upper[3]), linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = c(mad_counts_lower[2], mad_counts_upper[2]), linetype = "dotted", color = "gray30") +
  geom_vline(xintercept = c(mad_counts_lower[1], mad_counts_upper[1]), linetype = "dotted", color = "gray10") +
  geom_vline(xintercept = med_count, linetype = "solid", color = "blue", size = 1) + # Median
  geom_vline(xintercept = nCount_cutoff, linetype = "dashed", color = "red", size = 1) + # Cutoff
  ggtitle("nCount_RNA (Log10) with Median & +/- 1,2,3 MAD") +
  theme_minimal()

# 4. Augmented nFeature_RNA Distribution
p_features <- ggplot(qc_data, aes(x = nFeature_RNA)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  geom_vline(xintercept = c(mad_feats_lower[3], mad_feats_upper[3]), linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = c(mad_feats_lower[2], mad_feats_upper[2]), linetype = "dotted", color = "gray30") +
  geom_vline(xintercept = c(mad_feats_lower[1], mad_feats_upper[1]), linetype = "dotted", color = "gray10") +
  geom_vline(xintercept = med_feat, linetype = "solid", color = "blue", size = 1) +
  geom_vline(xintercept = nFeature_cutoff, linetype = "dashed", color = "darkgreen", size = 1) + # Cutoff
  ggtitle("nFeature_RNA (Log10) with Median & +/- 1,2,3 MAD") +
  theme_minimal()

# 5. Augmented Area vs nCount_RNA (2D Bin Density)
# Use a subset filtering out 0 counts to prevent axis scaling crashes
plot_data_area <- qc_data[qc_data$nCount_RNA > 0, ]

p_area_counts_heatmap <- ggplot(plot_data_area, aes(x = Area.um2, y = nCount_RNA)) +
  # Map the log10 transformation directly to the computed count statistic
  geom_bin2d(aes(fill = after_stat(log10(count))), bins = 100) + 
  scale_fill_viridis_c(option = "inferno", name = "log10(Count)") + # Removed trans="log10"
  geom_vline(xintercept = c(10, 350), linetype = "dashed", color = "blue") + 
  geom_hline(yintercept = nCount_cutoff, linetype = "dashed", color = "blue") + 
  scale_y_log10(labels = scales::comma) +
  coord_cartesian(xlim = c(0, 800)) + 
  ggtitle("Area vs nCount_RNA (2D Bin Density)") +
  theme_minimal()


p_counts_slide <- ggplot(qc_data, aes(x = nCount_RNA, fill = Run_Tissue_name)) +
  geom_density(alpha = 0.5) + scale_x_log10(labels = scales::comma) +
  facet_wrap(~Run_Tissue_name) + geom_vline(xintercept = nCount_cutoff, linetype = "dashed", color = "red") +
  ggtitle("nCount_RNA Distribution Across Slides") + theme_minimal() + theme(legend.position = "none")

# p_counts_slide <- ggplot(qc_data, aes(x = mus, y = nCount_RNA, fill = mus)) +
#   geom_violin(alpha = 0.7, scale = "width", color = "black", linewidth = 0.2) +
#   scale_y_log10(labels = scales::comma) +
#   geom_hline(yintercept = 200, linetype = "dashed", color = "red", linewidth = 1) +
#   ggtitle("nCount_RNA Distribution Across Samples") +
#   xlab("Sample (mus)") + ylab("nCount_RNA (Log10)") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
#     legend.position = "none"
#   )

# 6. Corrected Transcript Dominance
p_dominance <- ggplot(qc_data, aes(x = top_gene_fraction)) +
  geom_density(fill = "orange", alpha = 0.7) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  ggtitle("Transcript Dominance (Fraction of Top Gene)") +
  theme_minimal()

# 7. Nuclear Fraction Distribution
# Note: Adjust 'nuclear_fraction' if your exact metadata column is named differently (e.g., 'prop_nuc')
p_nuc <- ggplot(qc_data, aes(x = nuclear_fraction)) +
  geom_density(fill = "cyan4", alpha = 0.7) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "red") +
  ggtitle("Nuclear Transcript Fraction") + theme_minimal()

combined_plot <- gridExtra::arrangeGrob(
  p_counts, 
  p_features, 
  p_area_counts_heatmap, 
  p_counts_slide, 
  p_dominance, 
  p_nuc, 
  ncol = 2
)

ggsave(
  filename = paste0(local_WD, "CosMx_CB_Initial_QC_plots.png"), 
  plot = combined_plot, 
  width = 12,
  height = 15, 
  dpi = 600
)

combined_plot

# --- 4. Dynamic MAD Calculation & Subsetting ---
message("Calculating Dynamic Filter Thresholds...")

Area_min <- 10
Area_max <- 350
Top_Gene_Max <- 0.30 
Nuclear_Min <- 0.50

message(sprintf("Dynamic Cutoffs - Counts: %.1f | Features: %.1f", nCount_cutoff, nFeature_cutoff))

# Build logical vector with NA safeguarding
thisObj$pass_QC <- (
  thisObj$nCount_RNA > nCount_cutoff &
    thisObj$nFeature_RNA > nFeature_cutoff &
    thisObj$Area.um2 > Area_min & thisObj$Area.um2 < Area_max &
    thisObj$top_gene_fraction < Top_Gene_Max &
    thisObj$nuclear_fraction > Nuclear_Min
)
thisObj$pass_QC[is.na(thisObj$pass_QC)] <- FALSE # Catch any stray NAs

message(sprintf("Total cells passing all QC filters: %d", sum(thisObj$pass_QC)))

thisObj_filtered <- subset(thisObj, subset = pass_QC == TRUE)

message(sprintf("Cells before filtering: %d\nCells after filtering: %d", ncol(thisObj), ncol(thisObj_filtered)))



#'[Figure: Initial filter QC summary table]
# --- 5. Generate and Export QC Summary Table ---
library(gridExtra)
library(grid)
library(gtable)

# Calculate individual filter failures (identifies exactly how many cells violated each rule)
fail_nCount    <- sum(thisObj$nCount_RNA <= nCount_cutoff | is.na(thisObj$nCount_RNA))
fail_nFeature  <- sum(thisObj$nFeature_RNA <= nFeature_cutoff | is.na(thisObj$nFeature_RNA))
fail_Area      <- sum(thisObj$Area.um2 <= Area_min | thisObj$Area.um2 >= Area_max | is.na(thisObj$Area.um2))
fail_Dominance <- sum(thisObj$top_gene_fraction >= Top_Gene_Max | is.na(thisObj$top_gene_fraction))
fail_Nuclear   <- sum(thisObj$nuclear_fraction <= Nuclear_Min | is.na(thisObj$nuclear_fraction))

# Construct the summary data frame
qc_summary_df <- data.frame(
  Metric = c(
    "Starting Cells",
    sprintf("Failed: nCount_RNA \u2264 %.1f", nCount_cutoff),
    sprintf("Failed: nFeature_RNA \u2264 %.1f", nFeature_cutoff),
    sprintf("Failed: Area < %d or > %d", Area_min, Area_max),
    sprintf("Failed: Top Gene Fraction \u2265 %.2f", Top_Gene_Max),
    sprintf("Failed: Nuclear Fraction \u2264 %.2f", Nuclear_Min),
    "Total Cells Removed (Accounting for overlap)",
    "Final Retained Cells"
  ),
  Count = c(
    ncol(thisObj),
    fail_nCount,
    fail_nFeature,
    fail_Area,
    fail_Dominance,
    fail_Nuclear,
    ncol(thisObj) - ncol(thisObj_filtered),
    ncol(thisObj_filtered)
  )
)

# Format numbers with commas for readability
qc_summary_df$Count <- formatC(qc_summary_df$Count, format = "d", big.mark = ",")

# Define a clean, professional theme for the table graphic
t_theme <- ttheme_default(
  core = list(bg_params = list(fill = c("white", "gray95"), col = "white")),
  colhead = list(fg_params = list(col = "white", fontface = "bold"),
                 bg_params = list(fill = "black", col = "black"))
)

# Render the table grob
p_table <- tableGrob(qc_summary_df, theme = t_theme, rows = NULL)

# Add a title header to the table graphic
title <- textGrob("Spatial Transcriptomics QC Summary", gp = gpar(fontsize = 14, fontface = "bold"))
padding <- unit(5, "mm")
table_with_title <- gtable_add_rows(p_table, heights = grobHeight(title) + padding, pos = 0)
table_with_title <- gtable_add_grob(table_with_title, title, 1, 1, 1, ncol(table_with_title))

# Save as an independent, high-resolution PNG
ggsave(
  filename = paste0(local_WD, "CosMx_CB_QC_Summary_Table.png"), 
  plot = table_with_title, 
  width = 6.5, height = 4, dpi = 600
)

#'[NOTE: DO NOT CONTINUE UNTIL FILTERED OBJECT IS FINALIZED]
thisObj <- thisObj_filtered
thisObj$pass_QC <- NULL 
rm(thisObj_filtered, valid_counts, valid_feats)
gc()

saveRDS(thisObj, file.path(HDD_WD, "COSMONAUT_CB_Seurat_preALL_newcountcutoff.RDS"))
message("--- Initialization & QC Complete ---")

# Secondary Collision QC & Filtering ####
#**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
#'[Secondary Collision QC & Filtering]
#**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
# --- Phase 2A: Composite Module Scoring & Diagnostic Plotting ---

thisObj <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_preALL_newcountcutoff.RDS"))

# 1. Define Cell-Type Modules
module_genes <- list(
  Granule     = c("Slc17a7"), 
  Purkinje    = c("Calb1", "Car8", "Pcp2", "Itpr1"),
  Astrocyte   = c("Gfap", "Slc1a3"),
  Vasculature = c("Cldn5", "Vtn"),
  GABAergic   = c("Gad1", "Gad2", "Slc6a1", "Slc32a1"),
  Oligo       = c("Mbp", "Plp1", "Mog", "Mag"),
  Neuron      = c("Snap25", "Syt1", "Syp"),
  Microglia   = c("P2ry12", "C1qa", "C1qb", "C1qc"),
  OPC         = c("Pdgfra", "Cspg5")
)

# 2. Define Collision Axes
collision_axes <- list(
  Gran_Purk    = c("Granule", "Purkinje"),
  Astro_Vasc   = c("Astrocyte", "Vasculature"),
  GABA_Astro   = c("GABAergic", "Astrocyte"),
  Oligo_Neuron = c("Oligo", "Neuron"),
  Micro_Neuron = c("Microglia", "Neuron"),
  Oligo_OPC    = c("Oligo", "OPC")
)

# 3. Extract dense matrix safely for targeted genes
all_module_genes <- unique(unlist(module_genes))
valid_genes <- all_module_genes[all_module_genes %in% rownames(thisObj)]
expr_mat <- as.data.frame(t(as.matrix(GetAssayData(thisObj, assay = "RNA", layer = "counts")[valid_genes, ])))

# Compute aggregate module sums per cell
module_scores <- data.frame(row.names = rownames(expr_mat))
for (mod_name in names(module_genes)) {
  genes_present <- module_genes[[mod_name]][module_genes[[mod_name]] %in% valid_genes]
  if (length(genes_present) > 1) {
    module_scores[[mod_name]] <- rowSums(expr_mat[, genes_present])
  } else if (length(genes_present) == 1) {
    module_scores[[mod_name]] <- expr_mat[[genes_present]]
  } else {
    module_scores[[mod_name]] <- 0
  }
}

# 4. Compute Scores and Generate QC Plots
md <- thisObj@meta.data
plot_list <- list()

for (axis_name in names(collision_axes)) {
  m1 <- collision_axes[[axis_name]][1]
  m2 <- collision_axes[[axis_name]][2]
  
  # Product of the two module scores
  raw_score <- module_scores[[m1]] * module_scores[[m2]]
  md[[paste0("score_comp_", axis_name)]] <- raw_score
  
  valid_scores <- raw_score[raw_score > 0 & !is.na(raw_score)]
  
  if (length(valid_scores) > 10) {
    log_scores <- log10(valid_scores)
    med_score <- median(log_scores)
    mad_score <- mad(log_scores)
    
    mad1.5 <- med_score + (1.5 * mad_score)
    mad2.0 <- med_score + (2.0 * mad_score)
    mad3.0 <- med_score + (3.0 * mad_score)
    
    p <- ggplot(data.frame(Score = log_scores), aes(x = Score)) +
      geom_density(fill = "mediumpurple", alpha = 0.6) +
      geom_vline(xintercept = mad1.5, linetype = "dashed", color = "green", size = 0.8) +
      geom_vline(xintercept = mad2.0, linetype = "dashed", color = "orange", size = 0.8) +
      geom_vline(xintercept = mad3.0, linetype = "dashed", color = "red", size = 0.8) +
      ggtitle(paste0(axis_name, " Collision Density")) +
      xlab("Log10(Composite Score)") + ylab("Density") +
      theme_minimal() +
      annotate("text", x = mad1.5, y = Inf, label = "+1.5", color = "green", angle = 90, vjust = -0.5, hjust = 1.5) +
      annotate("text", x = mad2.0, y = Inf, label = "+2.0", color = "orange", angle = 90, vjust = -0.5, hjust = 1.5) +
      annotate("text", x = mad3.0, y = Inf, label = "+3.0", color = "red", angle = 90, vjust = -0.5, hjust = 1.5)
    
    plot_list[[axis_name]] <- p
  }
}

# Bind scores back to Seurat metadata
thisObj@meta.data <- md

# Arrange and save the QC diagnostic grid
combined_qc_plot <- gridExtra::arrangeGrob(grobs = plot_list, ncol = 2)

ggsave(filename = paste0(local_WD, "CosMx_CB_Composite_Collision_QC.png"),
       plot = combined_qc_plot, 
       width = 9, height = 5,
       dpi = 600)

gridExtra::grid.arrange(grobs = plot_list, ncol = 2)


# --- Phase 2B: Apply Customized Adaptive Thresholds ---

# 1. Manually adjust these based on the QC density plots
mad_multipliers <- list(
  Gran_Purk    = 2, 
  Astro_Vasc   = 2,
  GABA_Astro   = 2,
  Oligo_Neuron = 2,
  Micro_Neuron = 2,
  Oligo_OPC    = 2
)

md <- thisObj@meta.data
collision_flags <- data.frame(row.names = rownames(md))

# 2. Flag cells based on customized cutoffs
for (axis_name in names(collision_axes)) {
  multiplier <- mad_multipliers[[axis_name]]
  raw_score <- md[[paste0("score_comp_", axis_name)]]
  
  valid_scores <- raw_score[raw_score > 0 & !is.na(raw_score)]
  if (length(valid_scores) > 10) {
    log_scores <- log10(valid_scores)
    cutoff <- 10^(median(log_scores) + (multiplier * mad(log_scores)))
  } else {
    cutoff <- Inf 
  }
  
  collision_flags[[paste0("flag_comp_", axis_name)]] <- raw_score > cutoff
}

# Add combined flag back to Seurat metadata
md$is_chimera <- rowSums(collision_flags, na.rm = TRUE) > 0
thisObj@meta.data <- md

message(sprintf("Total severe composite chimeras flagged: %d", sum(thisObj$is_chimera)))


# --- 3. Giotto Spatial Delaunay Network Construction (Corrected) ---
library(Giotto)

sample_ids <- unique(thisObj$mus)
full_counts <- GetAssayData(thisObj, assay = "RNA", layer = "counts")
all_edges_list <- list()

options(future.globals.maxSize = 10000 * 1024^2)

options(giotto.no_python_warn = TRUE)

for (sample in sample_ids) {
  cat(paste0("\n--- Building Spatial Network: ", sample, " ---\n"))
  
  sample_meta <- md[md$mus == sample, ]
  valid_barcodes <- base::intersect(sample_meta$cell_ID, base::colnames(full_counts))
  
  if(length(valid_barcodes) < 10) next
  
  # Updated deprecated raw_exprs to expression based on Giotto warning
  g_sample <- Giotto::createGiottoObject(
    expression = full_counts[, valid_barcodes],
    spatial_locs = data.frame(
      x = sample_meta$x_slide_mm * 1000,
      y = sample_meta$y_slide_mm * 1000,
      row.names = sample_meta$cell_ID
    ),
    cell_metadata = sample_meta
  )
  
  g_sample <- Giotto::createSpatialNetwork(gobject = g_sample, 
                                           method = 'Delaunay', 
                                           maximum_distance_delaunay = 150) 
  
  d_net <- g_sample@spatial_network$cell$Delaunay_network
  if (is.null(d_net)) {
    d_net <- g_sample@spatial_network$Delaunay_network 
  }
  
  if (!is.null(d_net)) {
    if (methods::.hasSlot(d_net, "networkDT")) {
      edge_df <- as.data.frame(d_net@networkDT)
    } else {
      edge_df <- as.data.frame(d_net)
    }
    
    # Explicitly namespace base::colnames() to bypass the conflicted package
    if ("from" %in% base::colnames(edge_df)) {
      edge_clean <- edge_df %>%
        dplyr::select(source = from, target = to, distance = distance) %>%
        dplyr::mutate(Sample = sample)
    } else if ("source" %in% base::colnames(edge_df)) {
      edge_clean <- edge_df %>%
        dplyr::select(source = source, target = target, distance = distance) %>%
        dplyr::mutate(Sample = sample)
    }
    
    all_edges_list[[sample]] <- edge_clean
  }
  
  rm(g_sample)
  gc()
}

master_spatial_network <- dplyr::bind_rows(all_edges_list)
message(sprintf("\nTotal physical cellular connections mapped: %d", nrow(master_spatial_network)))

# --- Phase 3: Spatial Validation of Chimeras (Composite Module Corrected) ---
library(data.table)

# 1. Convert to bidirectional network for complete neighborhood lookups
edges_dt <- as.data.table(master_spatial_network)
edges_bidir <- unique(rbind(
  edges_dt[, .(cell_id = source, neighbor_id = target)],
  edges_dt[, .(cell_id = target, neighbor_id = source)]
))

# 2. Extract COMPOSITE MODULE data into data.table 
module_dt <- as.data.table(module_scores, keep.rownames = "cell_id")

# 3. Map neighbor module expression 
neighbor_expr <- merge(edges_bidir, module_dt, by.x = "neighbor_id", by.y = "cell_id", all.x = TRUE)

# 4. Aggregate maximum module expression per neighborhood
neighbor_max <- neighbor_expr[, lapply(.SD, max, na.rm = TRUE), by = cell_id, .SDcols = names(module_genes)]

# 5. Validate physical bleed with Asymmetric Dominance Rescue
md <- thisObj@meta.data
validated_flags <- data.frame(row.names = rownames(md))
match_idx <- match(rownames(md), neighbor_max$cell_id)


# QC: Collision Purity (Dominance) Distributions
purity_list <- list()

for (axis_name in names(collision_axes)) {
  m1 <- collision_axes[[axis_name]][1]
  m2 <- collision_axes[[axis_name]][2]
  
  # Extract pre-calculated module scores
  cell_m1 <- module_scores[[m1]]
  cell_m2 <- module_scores[[m2]]
  
  # Calculate the combined signal of the two colliding modules
  total_signal <- cell_m1 + cell_m2
  
  # Only evaluate cells that express at least one of the modules
  valid_idx <- total_signal > 0 
  
  if(sum(valid_idx) > 10) {
    ratio_m1 <- cell_m1[valid_idx] / total_signal[valid_idx]
    ratio_m2 <- cell_m2[valid_idx] / total_signal[valid_idx]
    
    # The "Purity" is whichever identity is stronger (ranges from 0.5 to 1.0)
    max_ratio <- pmax(ratio_m1, ratio_m2)
    
    purity_list[[axis_name]] <- data.frame(
      Purity = max_ratio,
      Axis = axis_name
    )
  }
}

purity_df <- dplyr::bind_rows(purity_list)

# Plot the distributions with reference lines for 60% and 80%
p_purity <- ggplot(purity_df, aes(x = Purity, fill = Axis)) +
  geom_density(alpha = 0.7, color = "black", linewidth = 0.2) +
  facet_wrap(~ Axis, scales = "free_y") +
  scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  geom_vline(xintercept = 0.60, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "orange", linewidth = 1) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "blue", linewidth = 1) +
  xlab("Chimera Ratio (0.5 = 1:1 Chimera, 1.0 = Pure Cell)") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_purity)

ggsave(filename = paste0(local_WD, "CosMx_CB_Collision_Purity_QC.png"),
       plot = p_purity, 
       width = 9, height = 5,
       dpi = 300)

#'[Purification x Delauney Filtering]
# 1. Define Variable Dominance Rescue Thresholds
rescue_thresholds <- list(
  Gran_Purk    = 0.80,
  Astro_Vasc   = 0.80,
  GABA_Astro   = 0.80,
  Oligo_Neuron = 0.70, # Relaxed for expected oligo x neuron contamination
  Micro_Neuron = 0.80,
  Oligo_OPC    = 0.80
)

# Universal cutoff: Any flagged cell below this purity is deleted regardless of neighbors
hard_kill_floor <- 0.60 

# 2. Convert to bidirectional network for complete neighborhood lookups
edges_dt <- as.data.table(master_spatial_network)
edges_bidir <- unique(rbind(
  edges_dt[, .(cell_id = source, neighbor_id = target)],
  edges_dt[, .(cell_id = target, neighbor_id = source)]
))

# 3. Extract COMPOSITE MODULE data into data.table 
module_dt <- as.data.table(module_scores, keep.rownames = "cell_id")

# 4. Map neighbor module expression 
neighbor_expr <- merge(edges_bidir, module_dt, by.x = "neighbor_id", by.y = "cell_id", all.x = TRUE)
neighbor_max <- neighbor_expr[, lapply(.SD, max, na.rm = TRUE), by = cell_id, .SDcols = names(module_genes)]

# 5. Validate physical bleed using customized rescue logic
md <- thisObj@meta.data
validated_flags <- data.frame(row.names = base::rownames(md))
match_idx <- match(base::rownames(md), neighbor_max$cell_id)

for (axis_name in names(collision_axes)) {
  m1 <- collision_axes[[axis_name]][1]
  m2 <- collision_axes[[axis_name]][2]
  dom_thresh <- rescue_thresholds[[axis_name]]
  
  # Retrieve the exact molecular flag directly from Phase 2B's collision_flags dataframe
  orig_flag <- collision_flags[[paste0("flag_comp_", axis_name)]]
  orig_flag[is.na(orig_flag)] <- FALSE 
  
  # Extract neighborhood maximums
  n_m1 <- neighbor_max[[m1]][match_idx]
  n_m2 <- neighbor_max[[m2]][match_idx]
  n_m1[is.na(n_m1)] <- 0
  n_m2[is.na(n_m2)] <- 0
  
  # Extract cell's own module expression
  cell_m1 <- module_scores[[m1]]
  cell_m2 <- module_scores[[m2]]
  
  # Calculate Max Dominance Ratio (Purity)
  total_signal <- cell_m1 + cell_m2
  total_signal[total_signal == 0] <- 1 
  max_ratio <- pmax(cell_m1 / total_signal, cell_m2 / total_signal)
  
  # Spatial Check: Does the neighbor drive the bleed?
  neighbor_source_exists <- (n_m1 > cell_m1 | n_m1 >= 5) | (n_m2 > cell_m2 | n_m2 >= 5)
  
  # FILTER LOGIC:
  # 1. Hard cutoff: Purity is below 0.60 (severe mixed artifact, kill even if neighbor is missing)
  # 2. Spatial cutoff: Purity is in the gray zone (< dom_thresh) AND a neighbor source confirms the bleed
  is_confirmed_artifact <- (max_ratio < hard_kill_floor) | ((max_ratio < dom_thresh) & neighbor_source_exists)
  
  validated_flags[[paste0("valid_", axis_name)]] <- orig_flag & is_confirmed_artifact
}

# 6. Apply validated spatial flags to Seurat object
md$true_spatial_chimera <- rowSums(validated_flags, na.rm = TRUE) > 0
thisObj@meta.data <- md

message(sprintf("Original molecular chimeras flagged: %d", sum(md$is_chimera)))
message(sprintf("True spatial collisions confirmed: %d", sum(md$true_spatial_chimera)))
message(sprintf("False positives (ambient noise) rescued: %d", sum(md$is_chimera) - sum(md$true_spatial_chimera)))

# 7. Final Spatial Subset
thisObj_cleaned <- subset(thisObj, subset = true_spatial_chimera == FALSE)
message(sprintf("Cells remaining after removing confirmed spatial collisions: %d", ncol(thisObj_cleaned)))

#'[Figure: QC Transcript Bleed Summary Table]

# --- 8. Generate and Export Spatial Collision Summary Table ---
library(gridExtra)
library(grid)
library(gtable)

# Calculate individual confirmed collisions per axis
fail_Gran_Purk    <- sum(validated_flags$valid_Gran_Purk, na.rm = TRUE)
fail_Astro_Vasc   <- sum(validated_flags$valid_Astro_Vasc, na.rm = TRUE)
fail_GABA_Astro   <- sum(validated_flags$valid_GABA_Astro, na.rm = TRUE)
fail_Oligo_Neuron <- sum(validated_flags$valid_Oligo_Neuron, na.rm = TRUE)
fail_Micro_Neuron <- sum(validated_flags$valid_Micro_Neuron, na.rm = TRUE)
fail_Oligo_OPC    <- sum(validated_flags$valid_Oligo_OPC, na.rm = TRUE)

# Construct the summary data frame
collision_summary_df <- data.frame(
  Metric = c(
    "Starting Cells (Post-Initial QC)",
    "Initial Molecular Chimeras Flagged",
    "Total Ambient Noise Rescued",
    "Confirmed: Granule x Purkinje",
    "Confirmed: Astrocyte x Vasculature",
    "Confirmed: GABAergic x Astrocyte",
    "Confirmed: Oligo x Neuron",
    "Confirmed: Microglia x Neuron",
    "Confirmed: Oligo x OPC",
    "Total Cells Removed (Accounting for overlap)",
    "Final Retained Cells"
  ),
  Count = c(
    ncol(thisObj),
    sum(md$is_chimera, na.rm = TRUE),
    sum(md$is_chimera, na.rm = TRUE) - sum(md$true_spatial_chimera, na.rm = TRUE),
    fail_Gran_Purk,
    fail_Astro_Vasc,
    fail_GABA_Astro,
    fail_Oligo_Neuron,
    fail_Micro_Neuron,
    fail_Oligo_OPC,
    sum(md$true_spatial_chimera, na.rm = TRUE),
    ncol(thisObj_cleaned)
  )
)

# Format numbers with commas for readability
collision_summary_df$Count <- formatC(collision_summary_df$Count, format = "d", big.mark = ",")

# Define a clean, professional theme for the table graphic
t_theme_col <- ttheme_default(
  core = list(bg_params = list(fill = c("white", "gray95"), col = "white")),
  colhead = list(fg_params = list(col = "white", fontface = "bold"),
                 bg_params = list(fill = "black", col = "black"))
)

# Render the table grob
p_table_col <- tableGrob(collision_summary_df, theme = t_theme_col, rows = NULL)
title_col <- textGrob("Spatial Transcript Bleed Summary", gp = gpar(fontsize = 14, fontface = "bold"))
padding <- unit(5, "mm")
table_with_title_col <- gtable_add_rows(p_table_col, heights = grobHeight(title_col) + padding, pos = 0)
table_with_title_col <- gtable_add_grob(table_with_title_col, title_col, 1, 1, 1, ncol(table_with_title_col))

ggsave(
  filename = paste0(local_WD, "CosMx_CB_TranscriptBleed_Summary_Table.png"), 
  plot = table_with_title_col, 
  width = 6.5, height = 5, dpi = 600
)


# Spatial Re-Mapping of Confirmed Collisions (Multi-Color)
plot_meta <- thisObj@meta.data
plot_meta$Plot_Label <- "Clean Cell"

# 1. Map the specific collision axes from the validated_flags dataframe
for (axis_name in names(collision_axes)) {
  col_name <- paste0("valid_", axis_name)
  idx <- which(validated_flags[[col_name]] == TRUE)
  plot_meta$Plot_Label[idx] <- axis_name
}

# 2. Flag cells that failed multiple distinct collision checks
multi_idx <- which(rowSums(validated_flags, na.rm = TRUE) > 1)
if (length(multi_idx) > 0) {
  plot_meta$Plot_Label[multi_idx] <- "Multiple Collisions"
}

# 3. Convert to factor with a locked level order
all_levels <- c("Clean Cell", names(collision_axes), "Multiple Collisions")
plot_meta$Plot_Label <- factor(plot_meta$Plot_Label, levels = all_levels)

# 4. Sort data frame so Clean Cells are plotted first (underneath the colored dots)
plot_meta <- plot_meta[order(plot_meta$Plot_Label), ]

# 5. Define a distinct, colorblind-friendly palette
custom_colors <- c(
  "Clean Cell"          = "gray85",
  "Gran_Purk"           = "#E41A1C",  
  "Astro_Vasc"          = "#377EB8",  
  "GABA_Astro"          = "#4DAF4A",  
  "Oligo_Neuron"        = "#984EA3",  
  "Micro_Neuron"        = "#FF7F00",  
  "Oligo_OPC"           = "#A65628",  
  "Multiple Collisions" = "black"
)

p_spatial_collisions <- ggplot(plot_meta, aes(x = x_global_px, y = y_global_px, color = Plot_Label)) +
  geom_point(size = 0.3, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = custom_colors) +
  facet_wrap(~ mus, scales = "free", ncol = 4) + 
  ggtitle("Spatial Mapping of Filtered Transcript Collisions") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  # Enlarge legend dots and make them opaque for readability
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 2)) 

print(p_spatial_collisions)

ggsave(filename = paste0(local_WD, "CosMx_CB_Confirmed_Collisions_Map.png"), 
       plot = p_spatial_collisions, 
       width = 9, height = 8,
       dpi = 600)

#'[NOTE: DO NOT OVERWRITE thisObj UNTIL FILTERING IS FINALIZED]
thisObj <- thisObj_cleaned

# Cleanup environment memory
rm(thisObj_cleaned, edges_dt, edges_bidir, module_dt, neighbor_expr, neighbor_max, validated_flags)
gc()

saveRDS(thisObj, file.path(HDD_WD, "COSMONAUT_CB_Seurat_post_spatial_QC.RDS"))

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# TF Based Clustering | Pre-FastReseg ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**TF Based Clustering*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Note:* After testing gene based clustering and WNN integrations, TF based clustering
# was chosen since it provides a more stable cell type estimation by mitigating
# some transcription dropout/cell segmentation issues. 
# A cell's segmentation/boundary may be poorly defined in higher density regions
# such as the granule layer causing transcript bleed between segmentations
# These artifacts, if occurring often such as MBP/oligodendrocyte proximity,
# can bias/homogenize the UMAP/clustering. i.e. if 1 gene is strongly expressing,
# it may incorrectly weight a sample towards the wrong cell type such as a neuron
# that is contaminated with a few oligo transcripts. TF based enrichment with
# decoupleR uses transcription factor consensus.
# This assigns a TF activity score to each cell based on an unweighted least squares method.
# We look at the expression of all the gene targets. By looking at genes in the regulon
# it takes into account both gene magnitude and directionality compared to the TF.
# This increases confidence and stability of estimates and punishes cells will fewer
# transcripts if the identity isn't clear based on TF activity.
# Downsides: There is a very minor loss of cell type niches such as oligodendrocytes
# and OPCs. Since they share many of the same TFs in the lineage, they will group
# more closely than in typical gene expression based clustering.
# But even in the below data, the populations can still be delineated even in the global clustering.
# I tried WNN integration of both gene expression and TF data, but it still needs work.
# The gene expression integration to the WNN does provide slightly more granularity,
# but segmentation noise homogenizes the cell typing and some populations that are very
# clear in TF based clustering, are convoluted. Plus, WNN takes significantly longer
# to run due to the dataset sizes even after conservative feature selection.
# However, WNN integration is likely the future next step; balancing cell type estimmations
# and gene level granularity of cell populations.

#'[Subset the obj for testing]
# fraction_to_keep <- 1 #set % of cells to use in subset
# total_cells <- ncol(thisObj_raw)
# cells_to_sample <- floor(fraction_to_keep * total_cells)
# set.seed(42)
# random_cells <- sample(colnames(thisObj_raw), size = cells_to_sample)
# thisObj_subset <- subset(thisObj_raw, cells = random_cells)
# print(paste("Original object size:", ncol(thisObj_raw), "cells"))
# print(paste("Subset object size:", ncol(thisObj_subset), "cells"))
# 
# thisObj <- thisObj_subset

#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**TF Based Optimization & Comparison*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
library(Seurat)
library(decoupleR)
library(matrixStats)

thisObj <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_post_spatial_QC.RDS"))

# Query for mouse TF activity
net <- get_collectri(organism = 'mouse', split_complexes = FALSE)

# 1. Evaluate Regulon Intersection with CosMx Panel
panel_genes <- base::rownames(thisObj)
net_intersect <- net %>% dplyr::filter(target %in% panel_genes)

regulon_sizes <- net_intersect %>% 
  dplyr::group_by(source) %>% 
  dplyr::summarize(n_targets = n()) %>%
  dplyr::arrange(desc(n_targets))

# Plot Regulon Depths
p_regulon_size <- ggplot(regulon_sizes, aes(x = n_targets)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "darkgreen") +
  theme_minimal() +
  labs(title = "CollecTRI Regulon Sizes within CosMx Panel",
       subtitle = "Red = minsize 5 | Green = minsize 10",
       x = "Number of Target Genes Present in Panel", y = "Frequency (TFs)")

print(p_regulon_size)

ggsave(
  filename = paste0(local_WD, "CosMx_CB_collectri_regulon_coveraged_QC.png"), 
  plot = p_regulon_size, 
  width = 5, height = 3,
  dpi = 600
)

# 2. Original Method: LogNormalize + minsize 5
message("Running Baseline Model (LogNormalize | minsize = 5)...")
thisObj <- NormalizeData(thisObj, assay = "RNA")
norm_counts <- as.matrix(GetAssayData(object = thisObj, assay = "RNA", layer = "data"))

acts_base <- run_ulm(mat = norm_counts, net = net, .source = 'source', .target = 'target',
                     .mor = 'mor', minsize = 5)

acts_mat_base <- acts_base %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix() %>%
  t() 

max_tf_base <- matrixStats::colMaxs(abs(acts_mat_base))

# 3. Adjusted Method: ScaleData + minsize 10
message("Running Adjusted Model (ScaleData | minsize = 10)...")
thisObj <- ScaleData(thisObj, assay = "RNA", features = base::rownames(thisObj))
scale_counts <- as.matrix(GetAssayData(object = thisObj, assay = "RNA", layer = "scale.data"))

acts_adj <- run_ulm(mat = scale_counts, net = net, .source = 'source', .target = 'target',
                    .mor = 'mor', minsize = 10)

acts_mat_adj <- acts_adj %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix() %>%
  t() 

max_tf_adj <- matrixStats::colMaxs(abs(acts_mat_adj))


# 4. Comparative QC
# Align cells for comparison (scaling and minsize changes may drop cells entirely if they lack valid TFs)
common_cells <- base::intersect(base::colnames(acts_mat_base), base::colnames(acts_mat_adj))

compare_df <- data.frame(
  cell_id = common_cells,
  Baseline_Score = max_tf_base[match(common_cells, base::colnames(acts_mat_base))],
  Adjusted_Score = max_tf_adj[match(common_cells, base::colnames(acts_mat_adj))]
)

med_base <- median(compare_df$Baseline_Score, na.rm = TRUE)
med_adj <- median(compare_df$Adjusted_Score, na.rm = TRUE)

# Plot A: Score Shift Scatter
p_scatter <- ggplot(compare_df, aes(x = Baseline_Score, y = Adjusted_Score)) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(option = "inferno", trans = "log10") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 2.0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Max TF Score: Baseline vs Adjusted",
       x = "Baseline (LogNorm | minsize=5)", y = "Adjusted (Scaled | minsize=10)")

# Plot B: Distribution Overlay
compare_long <- compare_df %>% 
  dplyr::select(Baseline_Score, Adjusted_Score) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = "Method", values_to = "Score")

p_density <- ggplot(compare_long, aes(x = Score, fill = Method, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("Baseline_Score" = "cadetblue3", "Adjusted_Score" = "darkorange")) +
  scale_color_manual(values = c("Baseline_Score" = "cadetblue4", "Adjusted_Score" = "darkorange3")) +
  geom_vline(xintercept = 2.0, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "Density Distribution Shift",
       subtitle = sprintf("Median Baseline: %.2f | Median Adjusted: %.2f", med_base, med_adj),
       x = "Maximum Absolute TF Score (t-statistic)")

# Render diagnostic grid
compare_grid <- p_scatter | p_density
print(compare_grid)

ggsave(filename = paste0(local_WD, "CosMx_CB_TF_Optimization_QC.png"), 
       plot = compare_grid, 
       width = 8, height = 3,
       dpi = 600)


# Test: LogNormalize vs ScaleData (minsize = 5)

# 1. Baseline: LogNormalize (Already run, but defining explicitly)
thisObj <- NormalizeData(thisObj, assay = "RNA")
norm_counts <- as.matrix(GetAssayData(object = thisObj, assay = "RNA", layer = "data"))

acts_log <- run_ulm(mat = norm_counts, net = net, .source = 'source', .target = 'target',
                    .mor = 'mor', minsize = 5)

acts_mat_log <- acts_log %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix() %>%
  t() 

max_tf_log <- matrixStats::colMaxs(abs(acts_mat_log))

# 2. Experimental: ScaleData 
thisObj <- ScaleData(thisObj, assay = "RNA", features = base::rownames(thisObj))
scale_counts <- as.matrix(GetAssayData(object = thisObj, assay = "RNA", layer = "scale.data"))

acts_scale <- run_ulm(mat = scale_counts, net = net, .source = 'source', .target = 'target',
                      .mor = 'mor', minsize = 5)

acts_mat_scale <- acts_scale %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix() %>%
  t() 

max_tf_scale <- matrixStats::colMaxs(abs(acts_mat_scale))

# 3. Plot the true impact of the scaling method
common_cells <- base::intersect(base::colnames(acts_mat_log), base::colnames(acts_mat_scale))

compare_df_strict <- data.frame(
  cell_id = common_cells,
  LogNorm_Score = max_tf_log[match(common_cells, base::colnames(acts_mat_log))],
  Scaled_Score = max_tf_scale[match(common_cells, base::colnames(acts_mat_scale))]
)

compare_long_strict <- compare_df_strict %>% 
  dplyr::select(LogNorm_Score, Scaled_Score) %>% 
  tidyr::pivot_longer(cols = everything(), names_to = "Method", values_to = "Score")

p_density_strict <- ggplot(compare_long_strict, aes(x = Score, fill = Method, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("LogNorm_Score" = "cadetblue3", "Scaled_Score" = "purple3")) +
  scale_color_manual(values = c("LogNorm_Score" = "cadetblue4", "Scaled_Score" = "purple4")) +
  geom_vline(xintercept = 2.0, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  labs(title = "True Impact of Scaling Method on ULM (minsize = 5)",
       x = "Maximum Absolute TF Score (t-statistic)")

print(p_density_strict)

ggsave(
  filename = paste0(local_WD, "CosMx_CB_ScalevsLogNorm_TF_QC.png"), 
  plot = p_density_strict, 
  width = 5, height = 3,
  dpi = 600
)
#**NOTES: Based on the QC, minsize 5 and the scaledata result in the higher, more consistent t-statistics*

#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Finalized TF Based Clustering: ScaleData + minsize 5*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#'[Note: fresh analysis to avoid accidental carry over from optimization]

# 1. Scale Data (Required to stabilize ULM variance)
message("Scaling data for OLS regression assumptions...")
thisObj <- ScaleData(thisObj, assay = "RNA", features = base::rownames(thisObj))
scale_counts <- as.matrix(GetAssayData(object = thisObj, assay = "RNA", layer = "scale.data"))

# 2. Query CollecTRI Database
net <- get_collectri(organism = 'mouse', split_complexes = FALSE)

# 3. Run Univariate Linear Model (minsize = 5 to preserve regulons)
message("Running decoupleR ULM TF Enrichment...")
acts <- run_ulm(mat = scale_counts, net = net, .source = 'source', .target = 'target',
                .mor = 'mor', minsize = 5)

# 4. Format and inject TF Assay
acts_mat <- acts %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix() %>%
  t() # Transpose to get (TFs x Cells)

thisObj[['tf_activities']] <- CreateAssayObject(counts = acts_mat)

# 5. Calculate Confidence and Filter Low Confidence Cells (t < 2)
tf_data <- as.matrix(GetAssayData(thisObj, assay = "tf_activities", layer = "counts"))
thisObj$Max_TF_Score <- matrixStats::colMaxs(abs(tf_data))

starting_cells_tf   <- ncol(thisObj)
total_panel_genes   <- length(base::rownames(thisObj))
total_tfs_collectri <- length(unique(net$source))
tfs_evaluated       <- nrow(acts_mat) # TFs successfully scored after minsize=5
cells_dropped_tf    <- sum(thisObj$Max_TF_Score < 2)

cat(sprintf("Cells before TF activity filter: %d\n", starting_cells_tf))
thisObj <- subset(thisObj, subset = Max_TF_Score >= 2)
final_cells_tf <- ncol(thisObj)
cat(sprintf("Cells after TF activity filter: %d\n", final_cells_tf))

# 6. TF Clustering Optimization QC Summary Table
library(gridExtra)
library(grid)
library(gtable)

# Construct the summary data frame
tf_summary_df <- data.frame(
  Metric = c(
    "Starting Cells (Post-Collision QC)",
    "CosMx Panel Genes Available",
    "Total Mouse TFs (CollecTRI Network)",
    "TFs Scored (Panel Intersection & minsize 5)",
    "Cells Filtered (Max Absolute t-statistic < 2.0)",
    "Final Retained Cells (High TF Confidence)"
  ),
  Count = c(
    starting_cells_tf,
    total_panel_genes,
    total_tfs_collectri,
    tfs_evaluated,
    cells_dropped_tf,
    final_cells_tf
  )
)

# Format numbers with commas for readability
tf_summary_df$Count <- formatC(tf_summary_df$Count, format = "d", big.mark = ",")

# Define a clean, professional theme for the table graphic
t_theme_tf <- ttheme_default(
  core = list(bg_params = list(fill = c("white", "gray95"), col = "white")),
  colhead = list(fg_params = list(col = "white", fontface = "bold"),
                 bg_params = list(fill = "black", col = "black"))
)

# Render the table grob
p_table_tf <- tableGrob(tf_summary_df, theme = t_theme_tf, rows = NULL)
title_tf <- textGrob("Transcription Factor Activity QC Summary", gp = gpar(fontsize = 14, fontface = "bold"))
padding <- unit(5, "mm")

table_with_title_tf <- gtable_add_rows(p_table_tf, heights = grobHeight(title_tf) + padding, pos = 0)
table_with_title_tf <- gtable_add_grob(table_with_title_tf, title_tf, 1, 1, 1, ncol(table_with_title_tf))

ggsave(
  filename = paste0(local_WD, "CosMx_CB_TF_QC_Summary_Table.png"), 
  plot = table_with_title_tf, 
  width = 7.0, height = 3.5, dpi = 600
)

# Clean up memory
rm(scale_counts, acts, acts_mat, tf_data)
gc()

#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**TF-based Clustering - PCA/UMAP*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

DefaultAssay(thisObj) <- "tf_activities"

# --- 1. Rigorous TF Feature Selection (Variance-Based) ---
cat("Extracting raw TF t-statistics for variance calculation...\n")
tf_matrix <- as.matrix(GetAssayData(thisObj, assay = "tf_activities", layer = "data"))

# Calculate the actual statistical variance of each TF across the dataset
tf_vars <- matrixStats::rowVars(tf_matrix)
names(tf_vars) <- base::rownames(tf_matrix)

# Sort TFs by variance
tf_vars_sorted <- sort(tf_vars, decreasing = TRUE)

# Plot the variance to find the biological elbow (signal vs noise)
var_df <- data.frame(Rank = 1:length(tf_vars_sorted), Variance = tf_vars_sorted)
p_tf_var <- ggplot(var_df, aes(x = Rank, y = Variance)) +
  scale_x_continuous(breaks = seq(0, 350, by = 50))+
  geom_point(color = "grey20", size = 0.5, alpha = 0.4) +
  geom_line(color = "darkblue", alpha = 0.3) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 150, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 200, linetype = "dashed", color = "gold") +
  geom_vline(xintercept = 250, linetype = "dashed", color = "green") +
  geom_vline(xintercept = 300, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "solid", color = "red") +
  theme_minimal() +
  labs(title = "TF Variance Distribution",
       x = "TF Rank (Most to Least Variable)", y = "Variance")

print(p_tf_var)

ggsave(filename = paste0(local_WD, "CosMx_CB_TF_Variance.png"), 
       plot = p_tf_var, 
       width = 3.75, height = 2,
       dpi = 300)


# --- 2. Compute Target Gene Sparsity (Optimized for Sparse Matrices) ---
tf_matrix <- as.matrix(GetAssayData(thisObj, assay = "tf_activities", layer = "data"))
rna_counts <- GetAssayData(thisObj, assay = "RNA", layer = "counts")

tf_vars <- matrixStats::rowVars(tf_matrix)
tf_means <- rowMeans(abs(tf_matrix))

# Calculate sparsity: 1 - (number of non-zero counts / total cells)
gene_detection_rate <- Matrix::rowSums(rna_counts > 0) / ncol(rna_counts)
gene_sparsity <- 1 - gene_detection_rate

# Map gene sparsity to CollecTRI regulons
regulon_sparsity <- net %>%
  dplyr::filter(source %in% names(tf_vars) & target %in% base::rownames(rna_counts)) %>%
  dplyr::mutate(target_sparsity = gene_sparsity[target]) %>%
  dplyr::group_by(source) %>%
  dplyr::summarize(mean_regulon_sparsity = mean(target_sparsity, na.rm = TRUE))

# --- 3. Build Diagnostic Dataframe ---
diag_df <- data.frame(
  TF = names(tf_vars),
  Variance = tf_vars,
  Mean_Abs_Score = tf_means[names(tf_vars)]
) %>%
  dplyr::left_join(regulon_sparsity, by = c("TF" = "source")) %>%
  dplyr::arrange(desc(Variance)) %>%
  dplyr::mutate(
    Rank = row_number(),
    Category = ifelse(Rank > 250, "Tail (Rank > 250)", "Active (Rank 1-250)")
  )

# --- 4. Generate Diagnostic Plots ---

# Plot A: Mean Abs Score vs Variance
p_math <- ggplot(diag_df, aes(x = Mean_Abs_Score, y = Variance, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Active (Rank 1-250)" = "darkblue", "Tail (Rank > 250)" = "red")) +
  theme_minimal() +
  labs(title = "Variance Diagnostic",
       x = "Mean Absolute t-statistic",
       y = "TF Variance")+
  theme(legend.position = "none")

# Plot B: Regulon Sparsity vs Variance
p_biology <- ggplot(diag_df, aes(x = mean_regulon_sparsity, y = Variance, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Active (Rank 1-250)" = "darkblue", "Tail (Rank > 250)" = "red")) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "black") +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Cell TF Dropout Diagnostic",
       x = "Average Sparsity of Target Genes (% Zeros)",
       y = "TF Variance")+
  theme(legend.position = "left")

# Plot C: Density Distribution (Top 3 vs Bottom 3 TFs)
top_3_tfs <- diag_df$TF[1:5]
bottom_3_tfs <- diag_df$TF[(nrow(diag_df)-4):nrow(diag_df)]
selected_tfs <- c(top_3_tfs, bottom_3_tfs)

dist_df <- as.data.frame(t(tf_matrix[selected_tfs, ])) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "TF", values_to = "Score") %>%
  dplyr::mutate(
    Status = ifelse(TF %in% top_3_tfs, "Top 3 Variable", "Bottom 3 (Tail)"),
    TF = factor(TF, levels = selected_tfs) # Lock order
  )

p_dist <- ggplot(dist_df, aes(x = Score, fill = Status)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~TF, scales = "fixed", ncol = 5) +
  scale_fill_manual(values = c("Top 3 Variable" = "darkblue", "Bottom 3 (Tail)" = "red")) +
  theme_minimal() +
  labs(title = "Raw t-statistic Distributions",
       subtitle = "Top/bottom 5 TFs",
       x = "ULM t-statistic", y = "Density")+
  theme(legend.position = "bottom")

# Compile and print
diagnostic_grid <- (p_math | p_biology) / p_dist
print(diagnostic_grid)

ggsave(filename = paste0(local_WD, "CosMx_CB_TF_Tail_Diagnostics.png"), 
       plot = diagnostic_grid, 
       width = 9.75, height = 5.5,
       dpi = 300)

#'[Scree plot function]
ScreePlot <- function(seurat_object, reduction_name, dims = 50, annotate_thresholds = TRUE) {
  # Calculate variance using the specified reduction
  variance <- (Stdev(seurat_object, reduction = reduction_name))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Dynamically identify thresholds
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  # Calculate the Mathematical Elbow (first PC where variance drops by < 0.1%)
  var_diffs <- abs(diff(pc_var_explained))
  n_elbow <- which(var_diffs < 0.1)[1]
  
  # Fallback if the curve never flattens below 0.1%
  if (is.na(n_elbow)) n_elbow <- dims 
  
  cat(paste("For reduction '", reduction_name, "':\n", sep=""))
  cat("  Automated Elbow (Signal Drop-off):", n_elbow, "PCs\n")
  cat("  80% variance:", n_80, "PCs\n")
  cat("  90% variance:", n_90, "PCs\n\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = 1:dims,
    Variance = pc_var_explained[1:dims]
  )
  
  # Create the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    scale_y_continuous(breaks = seq(0, max(Scree.df$Variance) + 2, by = 2)) +
    geom_line(color = "black") +
    geom_point(color = "grey60", size = 1.5) +
    labs(title = paste("Scree Plot for", reduction_name),
         x = "Principal Components",
         y = "% Variance Explained") +
    theme_minimal()
  
  if (annotate_thresholds) {
    p <- p +
      geom_vline(xintercept = n_elbow, linetype = "solid", color = "red", linewidth = 0.8) +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      annotate("text", x = n_elbow, y = max(Scree.df$Variance) * 0.95, 
               label = paste0("Elbow (PC ", n_elbow, ")"), color = "red", hjust = -0.1) +
      annotate("text", x = n_90, y = max(Scree.df$Variance) * 0.85, 
               label = paste0("90% (PC ", n_90, ")"), color = "black", hjust = 0.2)
  }
  
  print(p)
  
  return(invisible(list(scree_plot = p, n_elbow = n_elbow, n_80 = n_80, n_90 = n_90)))
}

top_tfs <- names(tf_vars_sorted)[1:250]

cat(sprintf("Retained %d highly variable TFs for dimensionality reduction.\n", length(top_tfs)))
thisObj <- ScaleData(thisObj, features = top_tfs, verbose = FALSE)

# --- 2. Dimensionality Reduction on Complete Network ---
thisObj <- RunPCA(thisObj,
                  assay = "tf_activities",
                  features = top_tfs, 
                  reduction.name = 'pca_TFactivity',
                  npcs = 50, 
                  verbose = FALSE)

ScreeResults_tf <- ScreePlot(thisObj, reduction_name = 'pca_TFactivity', dims = 50)

print(ScreeResults_tf$scree_plot)

ggsave(filename = paste0(local_WD, "CosMx_CB_TF_ScreePlot.png"), 
       plot = ScreeResults_tf$scree_plot, 
       width = 5, height = 3,
       dpi = 300)


#'[Figure: QC PCA Dimension cutoffs]
pc_embeddings <- as.data.frame(Embeddings(thisObj, reduction = "pca_TFactivity")[, c(8, 12, 31, 40)])
colnames(pc_embeddings) <- c("PC8", "PC13", "PC31", "PC40")

plot_df <- cbind(thisObj@meta.data[, c("mus", "x_global_px", "y_global_px")], pc_embeddings)

# 2. To avoid plotting 200,000 cells x 4 plots (800,000 points), select one representative sample 
# (Update "301" to whichever sample ID has the clearest cerebellar anatomy)
sample_df <- plot_df[plot_df$mus == "YM04", ]

# 3. Helper function for consistent spatial mapping
Plot_PC_Spatial <- function(df, pc_col, title) {
  ggplot(df, aes(x = x_global_px, y = y_global_px, color = .data[[pc_col]])) +
    geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
    scale_color_viridis_c(option = "inferno") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "right",
      legend.key.height = unit(0.5, "cm")
    ) +
    ggtitle(title) +
    labs(color = "Loading")
}

# 4. Generate the diagnostic plots
p_pc8  <- Plot_PC_Spatial(sample_df, "PC8", "PC 8 (Elbow)")
p_pc12 <- Plot_PC_Spatial(sample_df, "PC13", "PC 13 ")
p_pc31 <- Plot_PC_Spatial(sample_df, "PC31", "PC 31 (80% Variance)")
p_pc40 <- Plot_PC_Spatial(sample_df, "PC40", "PC 40 (90% Variance)")

pc_grid <- (p_pc8 | p_pc12) / (p_pc31 | p_pc40)
print(pc_grid)

ggsave(filename = paste0(local_WD, "CosMx_CB_PC_Spatial_Validation.png"), 
       plot = pc_grid, 
       width = 6.5, height = 3.5,
       dpi = 900)


# DimHeatmap: Visualizing Coordinated Variance 
# Sample 1000 extreme cells per PC
# 'balanced = TRUE' plots both the positive and negative driving TFs
p_heatmaps <- DimHeatmap(thisObj, 
                         dims = c(ScreeResults_tf$n_elbow, ScreeResults_tf$n_elbow+5,
                                  ScreeResults_tf$n_elbow+10,ScreeResults_tf$n_elbow+20,
                                  ScreeResults_tf$n_80, ScreeResults_tf$n_90), 
                         cells = 1000, 
                         reduction = "pca_TFactivity", 
                         balanced = TRUE, 
                         fast = FALSE) &
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

print(p_heatmaps)

ggsave(filename = paste0(local_WD, "CosMx_CB_PC_Heatmaps.png"), 
       plot = p_heatmaps, 
       width = 6, height = 4,
       dpi = 600)

dims2use <- 1:ScreeResults_tf$n_elbow

# 3. Empirical k.param Optimization Sweep
# Sweep test multiple neighborhood sizes. 
# Smaller k (20) = Pros: higher local granularity Cons: risk of fragmentation.
# Larger k (60) = Pros: broader global structures Cons: risk of blurring rare populations.

k_test_values <- c(20, 30, 40, 60)
umap_plots <- list()

for (k in k_test_values) {
  cat(sprintf("\nTesting k.param = %d...\n", k))
  
  # Compute temporary neighbor graph
  thisObj <- FindNeighbors(thisObj, dims = dims2use, k.param = k, 
                           assay = "tf_activities", reduction = 'pca_TFactivity', 
                           verbose = FALSE)
  
  # Run SLM clustering on this specific network
  cluster_name <- paste0("tf_snn_res.0.9_k", k)
  thisObj <- FindClusters(thisObj, resolution = 0.9, algorithm = 3, random.seed = 42, verbose = FALSE)
  thisObj[[cluster_name]] <- Idents(thisObj)
  
  # Run UMAP 
  thisObj <- RunUMAP(thisObj, dims = dims2use,
                     assay = "tf_activities", reduction = 'pca_TFactivity',
                     reduction.name = paste0('umap_k', k),
                     min.dist = 0.1, spread = 1, verbose = FALSE)
  
  # Generate diagnostic plot for this k value
  p <- DimPlot(thisObj, reduction = paste0('umap_k', k), group.by = cluster_name, label = TRUE) +
    ggtitle(sprintf("k.param = %d (Clusters: %d)", k, length(unique(thisObj[[cluster_name]][[1]])))) +
    theme(legend.position = "none")
  
  umap_plots[[as.character(k)]] <- p
}

# --- 4. Evaluate the Sweep ---
k_sweep_grid <- wrap_plots(umap_plots, ncol = 2)
print(k_sweep_grid)

ggsave(filename = paste0(local_WD, "CosMx_CB_kparam_optimization.png"), 
       plot = k_sweep_grid, width = 14, height = 12, dpi = 300)


























# --- 2. Dimensionality Reduction on Variable TFs ---
thisObj <- RunPCA(thisObj,
                  assay = "tf_activities",
                  features = VariableFeatures(thisObj), 
                  reduction.name = 'pca_TFactivity',
                  npcs = 50, 
                  verbose = FALSE)

ScreeResults_tf <- ScreePlot(thisObj, reduction_name = 'pca_TFactivity', dims = 50)
dims2use <- 1:ScreeResults_tf$n_90

# --- 3. Empirical k.param Optimization Sweep ---
# We will test a range of neighborhood sizes. 
# Smaller k (20) = high local granularity (risk of micro-fragmentation).
# Larger k (60) = broader global structure (risk of blurring rare populations).

k_test_values <- c(20, 30, 40, 60)
umap_plots <- list()

for (k in k_test_values) {
  cat(sprintf("\nTesting k.param = %d...\n", k))
  
  # Compute temporary neighbor graph
  thisObj <- FindNeighbors(thisObj, dims = dims2use, k.param = k, 
                           assay = "tf_activities", reduction = 'pca_TFactivity', 
                           verbose = FALSE)
  
  # Run SLM clustering on this specific network
  cluster_name <- paste0("tf_snn_res.0.9_k", k)
  thisObj <- FindClusters(thisObj, resolution = 0.9, algorithm = 3, random.seed = 42, verbose = FALSE)
  thisObj[[cluster_name]] <- Idents(thisObj)
  
  # Run UMAP (Only needs to be run once if you fix the topology, but doing it per k shows the exact spatial embedding effect)
  thisObj <- RunUMAP(thisObj, dims = dims2use,
                     assay = "tf_activities", reduction = 'pca_TFactivity',
                     reduction.name = paste0('umap_k', k),
                     min.dist = 0.1, spread = 1, verbose = FALSE)
  
  # Generate diagnostic plot for this k value
  p <- DimPlot(thisObj, reduction = paste0('umap_k', k), group.by = cluster_name, label = TRUE) +
    ggtitle(sprintf("k.param = %d (Clusters: %d)", k, length(unique(thisObj[[cluster_name]][[1]])))) +
    theme(legend.position = "none")
  
  umap_plots[[as.character(k)]] <- p
}

# --- 4. Evaluate the Sweep ---
k_sweep_grid <- wrap_plots(umap_plots, ncol = 2)
print(k_sweep_grid)

ggsave(filename = paste0(local_WD, "CosMx_CB_kparam_optimization.png"), 
       plot = k_sweep_grid, width = 14, height = 12, dpi = 300)































#**OLD OLD OLD OLD OLD OLD*

DefaultAssay(thisObj) <- "tf_activities"
thisObj <- ScaleData(thisObj, assay = "tf_activities")
thisObj <- FindVariableFeatures(thisObj,
                                assay = "tf_activities",
                                selection.method = "vst", nfeatures = 500, verbose = TRUE)
cat("# Variable Features: ", length(Seurat::VariableFeatures(thisObj)))

thisObj<- RunPCA(thisObj,
                 assay = "tf_activities",
                 reduction.name = 'pca_TFactivity',
                 npcs = 50, 
                 verbose = TRUE)


print(ElbowPlot(thisObj, reduction = 'pca_TFactivity',  ndims = 50)) +ggtitle("TF Elbow Plot")


ScreePlot <- function(seurat_object, reduction_name, dims = 50, annotate_thresholds = TRUE) {
  # Calculate variance using the specified reduction
  variance <- (Stdev(seurat_object, reduction = reduction_name))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Find the number of PCs needed
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  cat(paste("For reduction '", reduction_name, "':\n", sep=""))
  cat("  80% variance:", n_80, "PCs\n")
  cat("  90% variance:", n_90, "PCs\n")
  cat("  95% variance:", n_95, "PCs\n\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = 1:dims,
    Variance = pc_var_explained[1:dims]
  )
  
  # Create the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "black") +
    geom_point(color = "red", size = 3) +
    labs(title = paste("Scree Plot for", reduction_name),
         x = "Principal Components",
         y = "Percentage of Variance") +
    theme_minimal()
  
  if (annotate_thresholds) {
    p <- p +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      annotate("text", x = n_90, y = max(Scree.df$Variance) * 0.9, 
               label = paste0("90% (PC ", n_90, ")"), hjust = -0.1)
  }
  
  print(p)
  
  # Return results
  return(invisible(list(n_80 = n_80, n_90 = n_90, n_95 = n_95)))
}

ScreeResults_tf <- ScreePlot(thisObj, reduction_name = 'pca_TFactivity', dims = 50)

cat("  80% variance:", ScreeResults_tf$n_80, "PCs\n")
cat("  90% variance:", ScreeResults_tf$n_90, "PCs\n")
cat("  95% variance:", ScreeResults_tf$n_95, "PCs\n\n")

dims2use <- 1:ScreeResults_tf$n_90
k2use <- 20


thisObj <- Seurat::FindNeighbors(thisObj, dims = dims2use, k.param = k2use, assay = "tf_activities", reduction = 'pca_TFactivity')

thisObj <- Seurat::FindClusters(thisObj, resolution = .9, algorithm = 3, random.seed = 42)

thisObj <- Seurat::RunUMAP(thisObj, dims = dims2use,
                           assay = "tf_activities",
                           reduction = 'pca_TFactivity',
                           reduction.name = 'tf_umap',
                           min.dist = 0.0001, spread = 1)

base_umap <- DimPlot(thisObj, 
                     reduction = "tf_umap", 
                     group.by = 'tf_activities_snn_res.0.9', 
                     label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

# Add the repelled labels
plot_final <- LabelClusters(plot = base_umap, 
                            id = 'tf_activities_snn_res.0.9', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final)


#'[Merge fragmented clusters]

thisObj$seurat_clusters_merged <- as.character(thisObj$tf_activities_snn_res.0.9) # Create a new metadata column
thisObj$seurat_clusters_merged[thisObj$seurat_clusters_merged %in% c("0", "23", "16", "30", "28", "29", "27", "25",
                                                                     "1", "33", "18", "24", "9", "26", "32", "2", "5")] <- "0"

thisObj$seurat_clusters_merged[thisObj$seurat_clusters_merged %in% c("7", "4")] <- "4"

thisObj$seurat_clusters_merged[thisObj$seurat_clusters_merged %in% c("31", "11", "12")] <- "11"
thisObj$seurat_clusters_merged[thisObj$seurat_clusters_merged %in% c("6", "8", "3")] <- "3"
thisObj$seurat_clusters_merged[thisObj$seurat_clusters_merged %in% c("14", "19", "10", "21")] <- "10"



base_umap_merged <- DimPlot(thisObj, 
                     reduction = "tf_umap", 
                     group.by = 'seurat_clusters_merged', 
                     label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

plot_final_merged <- LabelClusters(plot = base_umap_merged, 
                            id = 'seurat_clusters_merged', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final_merged)



plot_cluster_diagnostics <- function(seurat_obj, 
                                     ident_name = "seurat_clusters", 
                                     sample_id_column = "orig.ident") {
  
  #' Plot Cell Proportions, Flag Missing Types, and Detect MAD Outliers
  #'
  #' @param seurat_obj A Seurat object
  #' @param ident_name Character. The metadata column containing cluster labels.
  #' @param sample_id_column Character. The metadata column containing sample IDs.
  #' @return A combined patchwork object of a stacked bar plot and a variance boxplot.
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  require(scales)
  
  meta <- seurat_obj@meta.data
  
  # 1. Calculate Exact Counts and FORCE Structural Zeros
  comp_df <- meta %>%
    dplyr::count(.data[[ident_name]], .data[[sample_id_column]], name = "n_cells") %>%
    tidyr::complete(.data[[ident_name]], .data[[sample_id_column]], fill = list(n_cells = 0))
  
  # 2. Calculate Proportions
  comp_df <- comp_df %>%
    dplyr::group_by(.data[[sample_id_column]]) %>%
    dplyr::mutate(sample_total = sum(n_cells),
                  proportion = n_cells / sample_total) %>%
    dplyr::ungroup()
  
  # A. Flag Missing Cell Types
  missing_types <- comp_df %>% dplyr::filter(n_cells == 0)
  
  if(nrow(missing_types) > 0) {
    cat("[WARNING] MISSING CELL TYPES DETECTED:\n")
    print(as.data.frame(missing_types %>% dplyr::select(-sample_total, -proportion)))
    cat("\n")
  } else {
    cat("[OK] All clusters are represented in all samples by at least 1 cell.\n\n")
  }
  
  # B. Calculate MAD and Flag Outliers
  outlier_df <- comp_df %>%
    dplyr::group_by(.data[[ident_name]]) %>%
    dplyr::mutate(
      med_prop = median(proportion),
      mad_prop = mad(proportion),
      mad_prop = ifelse(mad_prop == 0, 1e-6, mad_prop), # Prevent division by zero
      outlier_score = abs(proportion - med_prop) / mad_prop,
      is_outlier = outlier_score > 3.0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_outlier) %>%
    dplyr::arrange(desc(outlier_score))
  
  if(nrow(outlier_df) > 0) {
    cat("[ALERT] COMPOSITIONAL OUTLIERS (>3 MAD from Median):\n")
    print(as.data.frame(
      outlier_df %>% 
        dplyr::select(!!sym(ident_name), !!sym(sample_id_column), proportion, med_prop, outlier_score) %>%
        dplyr::mutate(
          proportion = paste0(round(proportion * 100, 2), "%"),
          med_prop = paste0(round(med_prop * 100, 2), "%"),
          outlier_score = paste0(round(outlier_score, 1), " MADs")
        ) %>%
        dplyr::rename(Sample_Prop = proportion, Expected_Median = med_prop, Deviation = outlier_score)
    ))
  } else {
    cat("[OK] No samples deviate by more than 3 MADs for any cell type.\n")
  }
  cat("====================================================\n")
  
  # Calculate exact number of unique clusters to define palette size
  n_clusters <- length(unique(comp_df[[ident_name]]))
  
  # Base high-contrast palette for Clusters (p1)
  base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
                "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#1B9E77", "#D95F02", 
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  
  # Dynamically expand the palette to match the exact number of clusters
  if (n_clusters > length(base_pal)) {
    expanded_pal <- colorRampPalette(base_pal)(n_clusters)
  } else {
    expanded_pal <- base_pal[1:n_clusters]
  }
  
  # Plot 1: Standard Stacked Bar
  p1 <- ggplot(comp_df, aes(x = as.factor(.data[[sample_id_column]]), 
                            y = proportion, 
                            fill = as.factor(.data[[ident_name]]))) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = expanded_pal) + # Dynamic Cluster Palette
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    ) +
    labs(title = "Overall Sample Composition",
         x = "Sample ID", y = "Proportion of Total Cells", fill = "Cluster")
  
  # Plot 2: Variance Boxplot
  p2 <- ggplot(comp_df, aes(x = as.factor(.data[[ident_name]]), 
                            y = proportion)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
    geom_jitter(aes(color = as.factor(.data[[sample_id_column]])), width = 0.2, size = 2, alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = "turbo") + # NEW: High-contrast palette for Samples
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Variance in Cluster Proportions Across Samples",
         subtitle = "Dots represent individual samples. Look for dots pulling far away from the boxplots.",
         x = "Cluster / Cell Type", y = "Proportion of Sample", color = "Sample")
  
  # Combine plots using patchwork
  final_plot <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1.5))
  
  return(final_plot)
}

plot_cluster_diagnostics(thisObj, 
                         ident_name = "seurat_clusters_merged", 
                         sample_id_column = "mus")


# Set merged cluster ids as identity
Idents(thisObj) <- "seurat_clusters_merged"
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# General Cell Typing of Global Clusters | Pre-FastReseg ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**General Cell Typing of Global Clusters*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
library("UCell")


gene_marker_lists <- list(
  # Neuronal (Cerebellar Cortex)
  Granule_Cell    = c("Cntn2", "Gls", "Rgs7"),
  UBC             = c("Calb2", "Grm1", "Plcb4"),
  Purkinje_Cell   = c("Pcp2", "Calb1", "Car8", "Itpr1", "Grid2"),
  Golgi_Cell      = c("Nos1", "Slc6a1"),
  Basket_Cell     = c("Pvalb", "Gad1"),
  Stellate_Cell   = c("Pvalb", "Penk", "Gad2"),
  Lugaro_Cell     = c("Calb2", "Gad2"),
  
  # Neuronal (Deep Cerebellar Nuclei)
  DCN_Excitatory  = c("Bdnf"),
  DCN_Inhibitory  = c("Slc32a1", "Spp1", "Gad2"),
  
  # Glial
  Astrocyte       = c("Aldh1l1", "Gfap", "Aqp4", "Slc1a2"),
  Bergmann_Glia   = c("Slc1a3"),
  OPC             = c("Pdgfra", "Olig2", "Gpr17"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag"),
  Microglia       = c("Tmem119", "P2ry12", "Cx3cr1", "C1qa"),
  
  # Vascular/Meningeal
  Endothelial     = c("Cldn5", "Pecam1"),
  Pericyte        = c("Pdgfrb", "Rgs5"),
  VSMC            = c("Acta2", "Tagln"),
  Fibroblast_VLMC = c("Col1a1", "Dcn"),
  Choroid_Plexus  = c("Ttr"),
  
  # Immune
  immune          = c("Cd3d", "Cd3e", "Cd3g", "Ptprc", "Ms4a1")
)


# Calculate Average SCT Expression
DefaultAssay(thisObj) <- "SCT"
avg_exp <- AverageExpression(thisObj, features = unique(unlist(gene_marker_lists)), 
                             group.by = "seurat_clusters_CellTyped", assay = "SCT", slot = "data")$SCT

# Score Gene Sets per Cluster
enrichment_matrix <- sapply(gene_marker_lists, function(gs) {
  present <- intersect(gs, rownames(avg_exp))
  if(length(present) > 0) colMeans(avg_exp[present, , drop=FALSE]) else rep(0, ncol(avg_exp))
}) %>% t()

# Plot Identity Heatmap
pheatmap(enrichment_matrix, main = "Cluster Identity: SCT Marker Enrichment",
                   scale = "row", cluster_cols = TRUE, border_color = "white",
                   color = colorRampPalette(c("blue", "grey", "red"))(100))

cluster_mapping <- c(
  "0" = "GLN", #granule layer neurons
  # "2" = "Granule_alt1",
  # "5" = "Granule_alt2",
  "20" = "OPC", 
  # "7" = "Oligodendrocyte_1",
  "4" = "OLG",
  "15" = "DCN",
  # "12" = "Astrocytes_DCN",
  "11" = "ASC", #Astrocytes
  # "14" = "VSMC/VLMC",
  "10" = "VASC", #Endothelial
  # "21" = "VASC",
  # "6" = "BG", #bergmann glia
  "13" = "PKJ", #purkinje neurons
  "3" = "MLI", #molecular layer interneurons
  "22" = "CHP", #Choroid Plexus
  "17" = "MGL" #Microglia
  # "8" = "MLI/PKJ/BG_MIX"
)

# Set merged cluster ids as current ids
Idents(thisObj) <- "seurat_clusters_merged" 

# Apply new cluster name ids
thisObj <- RenameIdents(thisObj, cluster_mapping)
thisObj$GlobalCellType <- Idents(thisObj)

Idents(thisObj) <- "GlobalCellType" 

base_umap_merged <- DimPlot(thisObj, 
                            reduction = "tf_umap", 
                            group.by = 'GlobalCellType', 
                            label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

plot_final_celltype <- LabelClusters(plot = base_umap_merged, 
                                   id = 'GlobalCellType', 
                                   repel = TRUE, 
                                   box = TRUE, 
                                   size = 3, 
                                   max.overlaps = Inf, 
                                   force = 20, 
                                   min.segment.length = 0)

print(plot_final_celltype)

plot_cluster_diagnostics <- function(seurat_obj, 
                                     ident_name = "seurat_clusters", 
                                     sample_id_column = "orig.ident") {
  
  #' Plot Cell Proportions, Flag Missing Types, and Detect MAD Outliers
  #'
  #' @param seurat_obj A Seurat object
  #' @param ident_name Character. The metadata column containing cluster labels.
  #' @param sample_id_column Character. The metadata column containing sample IDs.
  #' @return A combined patchwork object of a stacked bar plot and a variance boxplot.
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  require(scales)
  
  meta <- seurat_obj@meta.data
  
  # 1. Calculate Exact Counts and FORCE Structural Zeros
  comp_df <- meta %>%
    dplyr::count(.data[[ident_name]], .data[[sample_id_column]], name = "n_cells") %>%
    tidyr::complete(.data[[ident_name]], .data[[sample_id_column]], fill = list(n_cells = 0))
  
  # 2. Calculate Proportions
  comp_df <- comp_df %>%
    dplyr::group_by(.data[[sample_id_column]]) %>%
    dplyr::mutate(sample_total = sum(n_cells),
                  proportion = n_cells / sample_total) %>%
    dplyr::ungroup()
  
  # A. Flag Missing Cell Types
  missing_types <- comp_df %>% dplyr::filter(n_cells == 0)
  
  if(nrow(missing_types) > 0) {
    cat("[WARNING] MISSING CELL TYPES DETECTED:\n")
    print(as.data.frame(missing_types %>% dplyr::select(-sample_total, -proportion)))
    cat("\n")
  } else {
    cat("[OK] All clusters are represented in all samples by at least 1 cell.\n\n")
  }
  
  # B. Calculate MAD and Flag Outliers
  outlier_df <- comp_df %>%
    dplyr::group_by(.data[[ident_name]]) %>%
    dplyr::mutate(
      med_prop = median(proportion),
      mad_prop = mad(proportion),
      mad_prop = ifelse(mad_prop == 0, 1e-6, mad_prop), # Prevent division by zero
      outlier_score = abs(proportion - med_prop) / mad_prop,
      is_outlier = outlier_score > 3.0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_outlier) %>%
    dplyr::arrange(desc(outlier_score))
  
  if(nrow(outlier_df) > 0) {
    cat("[ALERT] COMPOSITIONAL OUTLIERS (>3 MAD from Median):\n")
    print(as.data.frame(
      outlier_df %>% 
        dplyr::select(!!sym(ident_name), !!sym(sample_id_column), proportion, med_prop, outlier_score) %>%
        dplyr::mutate(
          proportion = paste0(round(proportion * 100, 2), "%"),
          med_prop = paste0(round(med_prop * 100, 2), "%"),
          outlier_score = paste0(round(outlier_score, 1), " MADs")
        ) %>%
        dplyr::rename(Sample_Prop = proportion, Expected_Median = med_prop, Deviation = outlier_score)
    ))
  } else {
    cat("[OK] No samples deviate by more than 3 MADs for any cell type.\n")
  }
  cat("====================================================\n")
  
  # Calculate exact number of unique clusters to define palette size
  n_clusters <- length(unique(comp_df[[ident_name]]))
  
  # Base high-contrast palette for Clusters (p1)
  base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
                "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#1B9E77", "#D95F02", 
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  
  # Dynamically expand the palette to match the exact number of clusters
  if (n_clusters > length(base_pal)) {
    expanded_pal <- colorRampPalette(base_pal)(n_clusters)
  } else {
    expanded_pal <- base_pal[1:n_clusters]
  }
  
  # Plot 1: Standard Stacked Bar
  p1 <- ggplot(comp_df, aes(x = as.factor(.data[[sample_id_column]]), 
                            y = proportion, 
                            fill = as.factor(.data[[ident_name]]))) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = expanded_pal) + # Dynamic Cluster Palette
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    ) +
    labs(title = "Overall Sample Composition",
         x = "Sample ID", y = "Proportion of Total Cells", fill = "Cluster")
  
  # Plot 2: Variance Boxplot
  p2 <- ggplot(comp_df, aes(x = as.factor(.data[[ident_name]]), 
                            y = proportion)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
    geom_jitter(aes(color = as.factor(.data[[sample_id_column]])), width = 0.2, size = 2, alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = "turbo") + # NEW: High-contrast palette for Samples
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Variance in Cluster Proportions Across Samples",
         subtitle = "Dots represent individual samples. Look for dots pulling far away from the boxplots.",
         x = "Cluster / Cell Type", y = "Proportion of Sample", color = "Sample")
  
  # Combine plots using patchwork
  final_plot <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1.5))
  
  return(final_plot)
}

plot_cluster_diagnostics(thisObj, 
                         ident_name = "GlobalCellType", 
                         sample_id_column = "mus")

#'[SAVE: Save cell typed object]
# saveRDS(thisObj, file.path(local_WD, "COSMONAUT_CB_Seurat_GlobalCellTyped.RDS"))

# Subcluster Splitting Grouped Celltypes & SegShadow Identification | Pre-FastReseg ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Subcluster Splitting Grouped Celltypes & SegShadow Identification*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

project_WD <- '//mfad/researchmn/PHYSIO/SCHAFER/SHARED/JamesRedden/R Code/CEREBELLUM/'
local_WD <- "C:/Users/m322332/Documents/"
HDD_WD <- "E:/"

thisObj <- readRDS(file.path(local_WD, "COSMONAUT_CB_Seurat_GlobalCellTyped.RDS"))
thisObj

#'[Subclustering Fxns]
#'[Subcluster and ID clusters]
Subcluster_SCT <- function(seurat_obj, 
                           target_ident, 
                           ident_column = "GlobalCellType", 
                           nfeatures = 3000, 
                           resolution = 0.3, 
                           k_param = 20) {
  
  # 1. Subset the object
  cat(paste0("\n--- Processing: ", target_ident, " ---\n"))
  Idents(seurat_obj) <- ident_column
  sub_obj <- subset(seurat_obj, idents = target_ident)
  
  # 2. Set Assay and Find Variable Features
  DefaultAssay(sub_obj) <- "SCT"
  cat(paste0("Finding top ", nfeatures, " variable features in SCT assay...\n"))
  
  sub_obj <- FindVariableFeatures(
    sub_obj, 
    assay = "SCT", 
    selection.method = "vst", 
    nfeatures = nfeatures, 
    verbose = FALSE
  )
  
  # --- Direct Feature Extraction Bypass ---
  # Attempt to pull directly from the var.features slot
  var_feats <- sub_obj[["SCT"]]@var.features
  
  # Fallback for Seurat v5 Assay5 objects if the slot is empty
  if (is.null(var_feats) || length(var_feats) == 0) {
    sct_meta <- sub_obj[["SCT"]][[]] # Extract assay metadata
    if ("vst.variable" %in% colnames(sct_meta)) {
      var_feats <- rownames(sct_meta)[sct_meta$vst.variable == TRUE]
    } else {
      stop("Could not extract variable features from the SCT assay.")
    }
  }
  
  # 3. Scale the new variable features
  cat("Scaling data...\n")
  sub_obj <- ScaleData(
    sub_obj, 
    assay = "SCT", 
    features = var_feats, 
    verbose = FALSE
  )
  
  # 4. PCA and dynamic PC selection (80% variance)
  cat("Running PCA...\n")
  sub_obj <- RunPCA(
    sub_obj, 
    assay = "SCT", 
    features = var_feats, 
    reduction.name = "pca_sct_sub", 
    npcs = 50, 
    verbose = FALSE
  )
  
  # Print Elbow Plot for sanity check
  print(ElbowPlot(sub_obj, ndims = 50)) + ggtitle(paste0(target_ident, " Elbow Plot"))
  
  variance <- (Stdev(sub_obj, reduction = "pca_sct_sub"))^2
  cumvar <- cumsum(variance / sum(variance) * 100)
  n_80 <- which(cumvar >= 80)[1]
  
  # Fallback if 50 PCs don't reach 80%
  if (is.na(n_80)) n_80 <- 50
  cat(paste0("  -> 80% variance captured at PC ", n_80, ".\n"))
  
  dims2use <- 1:n_80
  
  # 5. Neighbors and Clustering
  cat(paste0("Clustering at resolution ", resolution, "...\n"))
  sub_obj <- FindNeighbors(
    sub_obj, 
    dims = dims2use, 
    k.param = k_param, 
    assay = "SCT", 
    reduction = "pca_sct_sub", 
    verbose = FALSE
  )
  
  sub_obj <- FindClusters(
    sub_obj, 
    resolution = resolution, 
    algorithm = 4, # Leiden algorithm (requires leiden base pkg) or 3 (SLM)
    random.seed = 42, 
    verbose = FALSE
  )
  
  # 6. UMAP
  cat("Running UMAP...\n")
  sub_obj <- RunUMAP(
    sub_obj, 
    dims = dims2use, 
    assay = "SCT", 
    reduction = "pca_sct_sub", 
    reduction.name = "umap_sct_sub", 
    min.dist = 0.0001, 
    spread = 3, 
    verbose = FALSE
  )
  
  # 7. Plotting
  p_umap <- DimPlot(sub_obj, reduction = "umap_sct_sub", label = FALSE) + 
    ggtitle(paste0(target_ident, " SCT Subclusters (Res: ", resolution, ")")) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")
    ) +
    coord_cartesian(clip = "off")
  
  p_final <- LabelClusters(
    plot = p_umap, 
    id = 'ident', 
    repel = TRUE, 
    box = TRUE, 
    size = 3, 
    max.overlaps = Inf, 
    force = 20, 
    min.segment.length = 0
  )
  
  print(p_final)
  
  return(list(
    "sub_obj" = sub_obj,
    "sub_umap" = p_final))
}

Run_Subcluster_Check <- function(global_obj, target_ident, ident_column = "GlobalCellType", resolution = 0.8) {
  
  # 1. Run your Subclustering Function
  sub_results <- Subcluster_SCT(
    seurat_obj = global_obj,
    target_ident = target_ident,
    ident_column = ident_column,
    nfeatures = 500,
    resolution = resolution
  )
  
  target_sub_obj <- sub_results$sub_obj
  
  # 2. Find Markers (Wilcoxon)
  cat("\n--- Finding Markers for", target_ident, "---\n")
  DefaultAssay(target_sub_obj) <- "SCT"
  cluster_markers_wilcox <- FindAllMarkers(
    target_sub_obj,
    only.pos = TRUE, 
    min.pct = 0.10, 
    logfc.threshold = 0.25, 
    test.use = "wilcox",
    slot = "data",
    verbose = FALSE
  )
  
  # 3. Process Top 5 Markers & Pi Score
  top5_markers <- cluster_markers_wilcox %>%
    group_by(cluster) %>%
    mutate(
      safe_p_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
      pi_score = -log10(safe_p_adj) * avg_log2FC
    ) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
    ungroup()
  
  top10_markers <- cluster_markers_wilcox %>%
    group_by(cluster) %>%
    mutate(
      safe_p_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
      pi_score = -log10(safe_p_adj) * avg_log2FC
    ) %>%
    slice_max(n = 10, order_by = avg_log2FC) %>%
    ungroup()
  
  # 4. Generate DotPlot
  top_genes_list <- unique(top5_markers$gene)
  p_dotplot <- DotPlot(
    target_sub_obj, 
    features = top_genes_list, 
    cols = c("lightgrey", "red"),
    dot.scale = 6,
    scale = FALSE 
  ) + 
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 10, face = "italic"),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    labs(
      title = paste("Top 5 Cluster Markers:", target_ident),
      subtitle = "Wilcoxon Test (Normalized Expression)",
      x = "Marker Gene", 
      y = "Cluster Identity"
    )
  
  print(p_dotplot)
  
  # 5. Print Marker Summary to Console for Review
  for (clust in levels(top10_markers$cluster)) {
    cat("\n========================================\n")
    cat("          Cluster", as.character(clust), "Top 5 Markers\n")
    cat("========================================\n")
    clust_data <- top10_markers[top10_markers$cluster == clust, ]
    print(clust_data[, c("gene", "avg_log2FC", "pi_score")])
  }
  
  # Return the finalized object and marker data for mapping
  return(list(
    "sub_obj" = target_sub_obj,
    "markers" = cluster_markers_wilcox,
    "top10" = top10_markers
  ))
}

QuickMarkers <- function(seurat_obj, target_ident) {
  # 1. Finding Markers
  cat("\n--- Finding Markers for", target_ident, "---\n")
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  cluster_markers_wilcox <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE, 
    min.pct = 0.10, 
    logfc.threshold = 0.25, 
    test.use = "wilcox",
    slot = "data",
    verbose = FALSE
  )
  
  # 2. Process Top 5 Markers & Pi Score
  top5_markers <- cluster_markers_wilcox %>%
    group_by(cluster) %>%
    mutate(
      safe_p_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
      pi_score = -log10(safe_p_adj) * avg_log2FC
    ) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
    ungroup()
  
  # 3. Generate DotPlot
  top_genes_list <- unique(top5_markers$gene)
  
  p_dotplot <- DotPlot(
    seurat_obj, 
    features = top_genes_list, 
    cols = c("lightgrey", "red"),
    dot.scale = 6,
    scale = FALSE 
  ) + 
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 10, face = "italic"),
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    labs(
      title = paste("Top 5 Cluster Markers:", target_ident),
      subtitle = "Wilcoxon Test (Normalized Expression)",
      x = "Marker Gene", 
      y = "Cluster Identity"
    )
  
  print(p_dotplot)
  
  # Return the markers dataframe for further use
  return(cluster_markers_wilcox)
}

Plot_Canonical_Lineage_Scores <- function(seurat_obj, group_by_col = "seurat_clusters") {
  
  # 1. Define the Gold Standard Biological Signatures
  canonical_signatures <- list(
    BEC = c("Cldn5", "Vwf", "Flt1", "Pecam1", "Slco1c1"),            
    VLM = c("Col1a1", "Lum", "Dcn", "Ptgds", "Slc6a15"),             
    PER = c("Pdgfrb", "Rgs5", "Abcc9", "Kcnj8", "Vtn"),              
    VSM = c("Acta2", "Tagln", "Myh11", "Myl9"),                      
    GLN = c("Slc17a7", "Gria2", "Neurod1", "Meg3", "Cacna1a"),       
    ASC = c("Aqp4", "Gfap", "Slc1a2", "Slc6a11", "Aldoc"),           
    MGL = c("C1qa", "C1qb", "Csf1r", "Tmem119", "Hexb"),             
    OLG = c("Mog", "Mag", "Plp1", "Mbp", "Myrf"),                    
    PKJ = c("Calb1", "Itpr1", "Pcp4", "Car8"),                       
    MLI = c("Gad1", "Gad2", "Slc32a1", "Penk", "Pvalb"),
    CHP = c("Ttr", "Htr2c"),
    OPC = c("Pdgfra", "Olig2", "Gpr17"),
    DCN = c("Pax6", "Tbr2", "Tbr1"),
    BG   = c("Slc1a3")
  )
  
  # 2. Filter signatures
  available_genes <- rownames(seurat_obj[["SCT"]])
  safe_signatures <- lapply(canonical_signatures, function(x) intersect(x, available_genes))
  safe_signatures <- safe_signatures[lengths(safe_signatures) > 0]
  
  # 3. Add Module Scores
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = safe_signatures,
    name = "LineageScore_",
    assay = "SCT",
    nbin = 12,    
    ctrl = 5      
  )
  
  # 4. Clean up column names
  score_cols <- paste0("LineageScore_", seq_along(safe_signatures))
  new_cols <- paste0("Score_", names(safe_signatures))
  for (i in seq_along(score_cols)) {
    colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == score_cols[i]] <- new_cols[i]
  }
  
  # 5. Aggregate scores
  Idents(seurat_obj) <- group_by_col
  score_summary <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(group_by_col), dplyr::all_of(new_cols)) %>%
    dplyr::group_by(!!sym(group_by_col)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
    tibble::column_to_rownames(var = group_by_col)
  
  colnames(score_summary) <- gsub("Score_", "", colnames(score_summary))
  
  # 6. Row-Scale (Z-score)
  score_mat_scaled <- t(scale(t(as.matrix(score_summary))))
  
  # --- THE FIX: Center colorbar at 0 ---
  # Find the symmetric limit based on the maximum divergence from 0
  limit <- max(abs(score_mat_scaled), na.rm = TRUE)
  # Generate 100 breaks from -limit to +limit
  palette_breaks <- seq(-limit, limit, length.out = 101) 
  
  # 7. Plot the Heatmap
  p_heat <- pheatmap(
    score_mat_scaled,
    main = paste("Canonical Module Scoring:", group_by_col),
    color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
    breaks = palette_breaks, # Forces the 0 point to be white
    cluster_rows = TRUE,    
    cluster_cols = FALSE,   
    border_color = "grey90",
    angle_col = 45,
    silent = FALSE
  )
  
  return(invisible(score_summary))
}

Run_Phase1_AnchorSearch <- function(seurat_obj, group_by_col = "seurat_clusters", reduction = "umap", target_lineages = NULL, custom_signatures = NULL) {
  
  # 1. Base Canonical Signatures
  base_sigs <- list(
    BEC = c("Cldn5", "Vwf", "Flt1", "Pecam1", "Slco1c1"),            
    VLM = c("Col1a1", "Lum", "Dcn", "Ptgds", "Slc6a15"),             
    PER = c("Pdgfrb", "Rgs5", "Abcc9", "Kcnj8", "Vtn"),              
    VSM = c("Acta2", "Tagln", "Myh11", "Myl9"),                      
    GLN = c("Slc17a7", "Gria2", "Neurod1", "Meg3", "Cacna1a"),       
    ASC = c("Aqp4", "Gfap", "Slc1a2", "Slc6a11", "Aldoc"),           
    MGL = c("C1qa", "C1qb", "Csf1r", "Tmem119", "Hexb"),             
    OLG = c("Mog", "Mag", "Plp1", "Mbp", "Myrf"),                    
    PKJ = c("Calb1", "Itpr1", "Pcp4", "Car8"),                       
    MLI = c("Gad1", "Gad2", "Slc32a1", "Penk", "Pvalb")              
  )
  
  if (!is.null(custom_signatures)) base_sigs <- c(base_sigs, custom_signatures)
  if (!is.null(target_lineages)) base_sigs <- base_sigs[names(base_sigs) %in% target_lineages]
  
  # 2. Safely filter for available genes
  available_genes <- rownames(seurat_obj[["SCT"]])
  safe_sigs <- lapply(base_sigs, function(x) intersect(x, available_genes))
  safe_sigs <- safe_sigs[lengths(safe_sigs) > 0]
  
  # 3. Add Module Scores
  seurat_obj <- AddModuleScore(seurat_obj, features = safe_sigs, name = "TempScore_", assay = "SCT", nbin = 12, ctrl = 5)
  score_cols <- paste0("TempScore_", seq_along(safe_sigs))
  new_cols <- names(safe_sigs)
  for (i in seq_along(score_cols)) colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == score_cols[i]] <- new_cols[i]
  
  # 4. Aggregate Raw Scores & Column Z-Scores
  Idents(seurat_obj) <- group_by_col
  score_summary <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(group_by_col), dplyr::all_of(new_cols)) %>%
    dplyr::group_by(!!sym(group_by_col)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
    tibble::column_to_rownames(var = group_by_col)
  
  score_mat_scaled <- as.data.frame(scale(score_summary))
  
  # 5. Identify the Predicted Anchors
  z_long <- score_mat_scaled %>%
    tibble::rownames_to_column("Cluster") %>%
    tidyr::pivot_longer(cols = -Cluster, names_to = "Lineage", values_to = "Z_Score")
  
  anchors <- z_long %>%
    dplyr::group_by(Lineage) %>%
    dplyr::slice_max(n = 1, order_by = Z_Score, with_ties = FALSE) %>%
    dplyr::mutate(Is_Anchor = TRUE)
  
  z_long <- z_long %>%
    dplyr::left_join(anchors %>% dplyr::select(Cluster, Lineage, Is_Anchor), by = c("Cluster", "Lineage")) %>%
    dplyr::mutate(Is_Anchor = tidyr::replace_na(Is_Anchor, FALSE))
  
  # 6. THE FIX: Half-Matrix Targeted Correlation
  avg_exp <- AverageExpression(seurat_obj, assays = "SCT", slot = "data", features = unlist(safe_sigs))$SCT
  avg_exp_mat <- as.matrix(avg_exp)
  colnames(avg_exp_mat) <- gsub("^g", "", colnames(avg_exp_mat))
  cor_mat <- cor(avg_exp_mat, method = "pearson")
  
  # Mask the upper triangle (keep the diagonal for visual grounding)
  cor_mat[upper.tri(cor_mat, diag = FALSE)] <- NA
  
  cor_long <- cor_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "Correlation", values_drop_na = TRUE)
  
  # 7. UMAP Extraction & Anchor Flagging
  umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$Cluster <- seurat_obj@meta.data[[group_by_col]]
  umap_df$Is_Anchor <- umap_df$Cluster %in% anchors$Cluster
  
  
  
  cluster_order <- rownames(score_mat_scaled)
  
  # Ensure factors for consistent ordering across all plots
  z_long$Cluster <- factor(z_long$Cluster, levels = rev(cluster_order))
  
  # Reverse Var1 (Y-axis) so the triangle anchors perfectly to the bottom-left
  cor_long$Var1 <- factor(cor_long$Var1, levels = rev(cluster_order))
  cor_long$Var2 <- factor(cor_long$Var2, levels = cluster_order)
  
  # Plot A: UMAP
  p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Cluster), size = 0.8, alpha = 0.5) +
    geom_point(data = subset(umap_df, Is_Anchor), aes(fill = Cluster), shape = 21, color = "black", stroke = 0.8, size = 1.5) +
    theme_void() +
    labs(title = "UMAP: Pure Anchor Mapping", subtitle = "Predicted anchors outlined in black") +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  # Plot B: Z-Scores (Magma)
  p_z <- ggplot(z_long, aes(x = Lineage, y = Cluster, fill = Z_Score)) +
    geom_tile(color = "white") +
    geom_tile(data = subset(z_long, Is_Anchor), color = "black", fill = NA, linewidth = 0.8) +
    viridis::scale_fill_viridis(option = "magma", name = "Z-Score") + 
    theme_minimal() +
    labs(title = "Canonical Z-Scores", subtitle = "Max scores outlined in black", x = "", y = "Cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), panel.grid = element_blank())
  
  # Plot C: Half-Matrix Bubble Correlogram
  p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1)) +
    geom_point(aes(size = abs(Correlation), color = Correlation)) +
    scale_color_gradient2(low = "#4575b4", mid = "grey90", high = "#d73027", midpoint = 0.5, limits = c(0, 1)) +
    scale_size_continuous(range = c(1, 6), guide = "none") +
    theme_minimal() +
    labs(title = "Bubble Correlogram", subtitle = "Lower triangle (Pearson R)", x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_line(color = "grey90"))
  
  # Assemble Dashboard
  layout <- "
  AACC
  BBCC
  "
  dashboard <- p_umap + p_z + p_cor + plot_layout(design = layout)
  
  print(dashboard)
  
  cat("\n--- Phase 1: Recommended Anchor Map ---\n")
  for (lin in unique(anchors$Lineage)) {
    anchor_clust <- anchors$Cluster[anchors$Lineage == lin]
    cat(sprintf("%s -> Cluster %s\n", lin, anchor_clust))
  }
  
  return(invisible(list(Raw = score_summary, Scaled = score_mat_scaled, Recommended_Anchors = anchors, Plot = dashboard)))
}

Plot_Subcluster_UMAP <- function(seurat_obj, group_by_col, plot_title, 
                                 reduction = "umap_sct_sub", 
                                 features = NULL) {
  
  # 1. Base Categorical UMAP (Cluster Labels)
  Idents(seurat_obj) <- group_by_col 
  
  base_umap <- DimPlot(
    seurat_obj, 
    group.by = group_by_col,
    reduction = reduction, 
    label = FALSE
  ) + 
    ggtitle(plot_title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.margin = margin(20, 20, 20, 20, unit = "pt")
    ) +
    coord_cartesian(clip = "off")
  
  plot_final <- LabelClusters(
    plot = base_umap, 
    id = group_by_col, 
    repel = TRUE, 
    box = TRUE, 
    size = 3.5, 
    max.overlaps = Inf, 
    force = 10, 
    min.segment.length = 0
  )
  
  # 2. Add Feature Heatmap logic if features are provided
  if (!is.null(features)) {
    
    # Handle 1D Heatmap (Single Gene)
    if (length(features) == 1) {
      p_feat <- FeaturePlot(
        seurat_obj, 
        features = features, 
        reduction = reduction, 
        order = TRUE, # Plots high-expression cells on top
        cols = viridis::magma(100)
      ) + 
        theme_minimal() +
        labs(title = paste("Hot Spot:", features)) +
        theme(plot.title = element_text(face = "italic"))
      
      # Combine side-by-side
      plot_final <- plot_final + p_feat
      
    } else if (length(features) == 2) {
      # Handle 2D Heatmap (Bivariate blend)
      # order = TRUE is critical here to see the "Double Positive" shadows
      p_feat <- FeaturePlot(
        seurat_obj, 
        features = features, 
        reduction = reduction, 
        blend = TRUE, 
        blend.threshold = 0.1,
        cols = c("lightgrey", "#ff0000", "#0000ff"), # Red for Gen1, Blue for Gen2, Purple for Both
        order = TRUE 
      )
      
      # FeaturePlot with blend = TRUE returns a list of 4 plots (Gen1, Gen2, Blend, Legend)
      # We want the 3rd plot (the blended overlay)
      plot_final <- plot_final + p_feat[[3]] + p_feat[[4]] # Add the 2D legend too
    }
  }
  
  # 3. Print and return
  print(plot_final)
  return(invisible(plot_final))
}

Check_Global_Mapping <- function(global_obj, sub_obj, target_cluster_name, global_reduction = "tf_umap") {
  library(Seurat)
  library(ggplot2)
  library(scales)
  
  # 1. Create a local copy (will automatically be destroyed when the function ends)
  temp_obj <- global_obj
  
  # 2. Extract the cell names and their new identities from the subclustered object
  sub_cells <- Cells(sub_obj)
  sub_new_idents <- as.character(Idents(sub_obj))
  names(sub_new_idents) <- sub_cells
  
  # 3. Create a new metadata column in the global object, defaulting to "Other"
  temp_obj$Global_Highlight <- "Other"
  
  # 4. Inject the new sub-identities into the local global metadata
  temp_obj$Global_Highlight[sub_cells] <- sub_new_idents
  
  # 5. Set the active identities to this new column
  Idents(temp_obj) <- "Global_Highlight"
  
  # 6. Define a color palette dynamically based on the number of subgroups
  unique_sub_groups <- unique(sub_new_idents)
  sub_colors <- scales::hue_pal()(length(unique_sub_groups))
  names(sub_colors) <- unique_sub_groups
  
  # Add the light grey color for all "Other" cells
  custom_colors <- c(sub_colors, "Other" = "grey85")
  
  # 7. Generate the Global Highlight UMAP
  p_global_check <- DimPlot(
    temp_obj, 
    reduction = global_reduction,
    group.by = "Global_Highlight",
    cols = custom_colors,
    label = TRUE,
    label.size = 3,
    repel = TRUE,
    pt.size = 0.5
  ) + 
    ggtitle(paste("Global UMAP:", target_cluster_name, "Subcluster Mapping Check")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  print(p_global_check)
  
  # Return the plot object invisibly in case you want to save it to a variable later
  return(invisible(p_global_check))
}

Check_Spatial_Mapping <- function(global_obj, sub_obj, sample_id, sample_col = "mus", target_cluster_name, custom_col_map = NULL) {
  
  # 1. Create local copy
  temp_obj <- global_obj
  
  # 2. Extract cells and new identities
  sub_cells <- Cells(sub_obj)
  sub_new_idents <- as.character(Idents(sub_obj))
  names(sub_new_idents) <- sub_cells
  
  # 3. Inject into metadata, defaulting to "Other"
  temp_obj$Spatial_Highlight <- "Other"
  temp_obj$Spatial_Highlight[sub_cells] <- sub_new_idents
  
  # 4. Extract metadata for the specific sample slide
  plot_df <- temp_obj@meta.data %>% 
    filter(!!sym(sample_col) == sample_id)
  
  unique_sub_groups <- unique(sub_new_idents)
  
  # 5. Determine Color Palette
  if (is.null(custom_col_map)) {
    # Default behavior: dynamic hue palette
    sub_colors <- scales::hue_pal()(length(unique_sub_groups))
    names(sub_colors) <- unique_sub_groups
    final_colors <- c("Other" = "grey85", sub_colors)
  } else {
    # Use user-provided palette
    final_colors <- custom_col_map
    # Failsafe: Ensure "Other" exists in the palette so the background plots correctly
    if (!"Other" %in% names(final_colors)) {
      final_colors["Other"] <- "grey85"
    }
  }
  
  # 6. Factor levels: ensure "Other" is first so it plots underneath
  plot_df$Spatial_Highlight <- factor(
    plot_df$Spatial_Highlight, 
    levels = c("Other", unique_sub_groups)
  )
  
  # Sort the dataframe so "Other" cells are rendered first
  plot_df <- plot_df %>% arrange(Spatial_Highlight)
  
  # 7. Generate the Spatial Plot
  p_spatial <- ggplot(plot_df, aes(x = x_slide_micron, y = y_slide_micron, color = Spatial_Highlight)) +
    geom_point(size = 0.75, alpha = 0.8) +
    scale_color_manual(values = final_colors) +
    coord_fixed() + # Forces 1:1 aspect ratio so the tissue doesn't stretch
    theme_void() +  # Removes axes and background
    ggtitle(paste("Spatial Map:", target_cluster_name, "Subclusters in Sample", sample_id)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = min(length(unique_sub_groups), 6))) 
  
  print(p_spatial)
  return(invisible(p_spatial))
}

plot_cluster_diagnostics <- function(seurat_obj, 
                                     ident_name = "seurat_clusters", 
                                     sample_id_column = "orig.ident") {
  
  #' Plot Cell Proportions, Flag Missing Types, and Detect MAD Outliers
  #'
  #' @param seurat_obj A Seurat object
  #' @param ident_name Character. The metadata column containing cluster labels.
  #' @param sample_id_column Character. The metadata column containing sample IDs.
  #' @return A combined patchwork object of a stacked bar plot and a variance boxplot.
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  require(scales)
  
  meta <- seurat_obj@meta.data
  
  # 1. Calculate Exact Counts and FORCE Structural Zeros
  comp_df <- meta %>%
    dplyr::count(.data[[ident_name]], .data[[sample_id_column]], name = "n_cells") %>%
    tidyr::complete(.data[[ident_name]], .data[[sample_id_column]], fill = list(n_cells = 0))
  
  # 2. Calculate Proportions
  comp_df <- comp_df %>%
    dplyr::group_by(.data[[sample_id_column]]) %>%
    dplyr::mutate(sample_total = sum(n_cells),
                  proportion = n_cells / sample_total) %>%
    dplyr::ungroup()
  
  # A. Flag Missing Cell Types
  missing_types <- comp_df %>% dplyr::filter(n_cells == 0)
  
  if(nrow(missing_types) > 0) {
    cat("[WARNING] MISSING CELL TYPES DETECTED:\n")
    print(as.data.frame(missing_types %>% dplyr::select(-sample_total, -proportion)))
    cat("\n")
  } else {
    cat("[OK] All clusters are represented in all samples by at least 1 cell.\n\n")
  }
  
  # B. Calculate MAD and Flag Outliers
  outlier_df <- comp_df %>%
    dplyr::group_by(.data[[ident_name]]) %>%
    dplyr::mutate(
      med_prop = median(proportion),
      mad_prop = mad(proportion),
      mad_prop = ifelse(mad_prop == 0, 1e-6, mad_prop), # Prevent division by zero
      outlier_score = abs(proportion - med_prop) / mad_prop,
      is_outlier = outlier_score > 3.0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_outlier) %>%
    dplyr::arrange(desc(outlier_score))
  
  if(nrow(outlier_df) > 0) {
    cat("[ALERT] COMPOSITIONAL OUTLIERS (>3 MAD from Median):\n")
    print(as.data.frame(
      outlier_df %>% 
        dplyr::select(!!sym(ident_name), !!sym(sample_id_column), proportion, med_prop, outlier_score) %>%
        dplyr::mutate(
          proportion = paste0(round(proportion * 100, 2), "%"),
          med_prop = paste0(round(med_prop * 100, 2), "%"),
          outlier_score = paste0(round(outlier_score, 1), " MADs")
        ) %>%
        dplyr::rename(Sample_Prop = proportion, Expected_Median = med_prop, Deviation = outlier_score)
    ))
  } else {
    cat("[OK] No samples deviate by more than 3 MADs for any cell type.\n")
  }
  cat("====================================================\n")
  
  # Calculate exact number of unique clusters to define palette size
  n_clusters <- length(unique(comp_df[[ident_name]]))
  
  # Base high-contrast palette for Clusters (p1)
  base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
                "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#1B9E77", "#D95F02", 
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  
  # Dynamically expand the palette to match the exact number of clusters
  if (n_clusters > length(base_pal)) {
    expanded_pal <- colorRampPalette(base_pal)(n_clusters)
  } else {
    expanded_pal <- base_pal[1:n_clusters]
  }
  
  # Plot 1: Standard Stacked Bar
  p1 <- ggplot(comp_df, aes(x = as.factor(.data[[sample_id_column]]), 
                            y = proportion, 
                            fill = as.factor(.data[[ident_name]]))) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = expanded_pal) + # Dynamic Cluster Palette
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    ) +
    labs(title = "Overall Sample Composition",
         x = "Sample ID", y = "Proportion of Total Cells", fill = "Cluster")
  
  # Plot 2: Variance Boxplot
  p2 <- ggplot(comp_df, aes(x = as.factor(.data[[ident_name]]), 
                            y = proportion)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
    geom_jitter(aes(color = as.factor(.data[[sample_id_column]])), width = 0.2, size = 2, alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = "turbo") + # NEW: High-contrast palette for Samples
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Variance in Cluster Proportions Across Samples",
         subtitle = "Dots represent individual samples. Look for dots pulling far away from the boxplots.",
         x = "Cluster / Cell Type", y = "Proportion of Sample", color = "Sample")
  
  # Combine plots using patchwork
  final_plot <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1.5))
  
  return(final_plot)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Subcluster and ID clusters, find shadow clusters/poor segmentations]

#'[VASC: Subclustering]
VASC_subcluster_res <- Run_Subcluster_Check(
  global_obj = thisObj, 
  target_ident = "VASC", 
  resolution = 0.8
)

VASC_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = VASC_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "VASC Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(VASC_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

VASC_AnchorClusters <- Run_Phase1_AnchorSearch(seurat_obj = VASC_subcluster_res$sub_obj, 
                                  reduction = "umap_sct_sub",
                                  group_by_col = "SCT_snn_res.0.8", 
                                  target_lineages = c("BEC", "VLM", "PER", "VSM"), 
                                  custom_signatures = NULL)

VASC_subcluster_res$sub_obj$VASC_SubClusterCellType <- as.character(VASC_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Vascular Cores
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType %in% c("7")] <- "BEC"
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType %in% c("8")] <- "VLM.BEC"
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "2"] <- "VLM"        
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "6"] <- "VSM"        
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "9"] <- "PER"        

# 3. Micro-environmental Vascular potential segmentation shadows
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "10"] <- "VLM.ASC"       
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "12"] <- "VLM.MGL"       
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType == "1"]  <- "BEC.OLG"       
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType %in% c("3","5")] <- "VASC.GLN" 
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType %in% c("11")] <- "VLM.MLI"
VASC_subcluster_res$sub_obj$VASC_SubClusterCellType[VASC_subcluster_res$sub_obj$VASC_SubClusterCellType %in% c("4")] <- "BEC.MLI"


Idents(VASC_subcluster_res$sub_obj) <- "VASC_SubClusterCellType"

Plot_Subcluster_UMAP(
  seurat_obj = VASC_subcluster_res$sub_obj,
  group_by_col = "VASC_SubClusterCellType",
  plot_title = "VASC Subclustering (Merged Annotations)"
)

QuickMarkers(seurat_obj = VASC_subcluster_res$sub_obj,
             target_ident = "VASC_SubClusterCellType")

Check_Global_Mapping(
  global_obj = thisObj, 
  sub_obj = VASC_subcluster_res$sub_obj, 
  target_cluster_name = "VASC", 
  global_reduction = "tf_umap"
)
unique(VASC_subcluster_res$sub_obj$mus)
vasc_highlight_colors <- c(
  "BEC" = "red",
  "VLM" = "green",
  "VSM" = "blue",
  "PER" = "black",
  "VLM.ASC" = "grey80",
  "VLM.MGL" = "grey80",
  "BEC.OLG" = "grey80",
  "VASC.GLN" = "grey80",
  "VLM.MLI" = "grey80",
  "BEC.MLI" = "grey80",
  "VASC.PKJ" = "grey80"
)

Check_Spatial_Mapping(
  global_obj = thisObj,
  sub_obj = VASC_subcluster_res$sub_obj,
  sample_id = "YM04",
  sample_col = "mus",
  target_cluster_name = "VASC",
  custom_col_map = vasc_highlight_colors
)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[GLN: Subclustering]
GLN_subcluster_res <- Run_Subcluster_Check(
  global_obj = thisObj, 
  target_ident = "GLN", 
  resolution = 0.8
)

GLN_subcluster_res$top10

Plot_Canonical_Lineage_Scores(GLN_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

Plot_Subcluster_UMAP(
  seurat_obj = GLN_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "GLN Subclustering (pre-ID)"
)

GLN_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = GLN_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("GLN"), 
  custom_signatures = NULL
)

# 1. Create temporary merge column
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType <- as.character(GLN_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Granule Cores
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("1", "3", "2", "4")] <- "GLN"

# 3. Micro-environmental potential segmentation shadows
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("13")] <- "GLN.MGL"
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("7")] <- "GLN.BEC"
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("6", "5", "14", "12")] <- "GLN.ASC"
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("8","10")] <- "GLN.OLG"
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("9", "11")] <- "GLN.MLI"
GLN_subcluster_res$sub_obj$GLN_SubClusterCellType[GLN_subcluster_res$sub_obj$GLN_SubClusterCellType %in% c("15")] <- "GLN.VLM"



Idents(GLN_subcluster_res$sub_obj) <- "GLN_SubClusterCellType"

Plot_Subcluster_UMAP(
  seurat_obj = GLN_subcluster_res$sub_obj,
  group_by_col = "GLN_SubClusterCellType",
  plot_title = "GLN Subclustering (Merged Annotations)"
)

Plot_Canonical_Lineage_Scores(GLN_subcluster_res$sub_obj, group_by_col = "GLN_SubClusterCellType")

QuickMarkers(seurat_obj = GLN_subcluster_res$sub_obj, target_ident = "GLN_SubClusterCellType")

Check_Global_Mapping(
  global_obj = thisObj, 
  sub_obj = GLN_subcluster_res$sub_obj, 
  target_cluster_name = "GLN", 
  global_reduction = "tf_umap"
)
unique(thisObj$mus)
gln_highlight_colors <- c(
  "GLN"      = "grey80",  # Pure populations
  "GLN.ASC"  = "grey80",     # Astrocyte/Bergmann Glia bleed
  "GLN.MGL"  = "grey80",  # Microglia bleed
  "GLN.OLG"  = "red",  # Oligodendrocyte bleed
  "GLN.VASC" = "grey80",    # Vascular/Endothelial bleed
  "GLN.MLI"  = "grey80",   # Molecular Layer Interneuron bleed
  "GLN.VLM"  = "grey80"  # Leptomeningeal/Fibroblast bleed
)
Check_Spatial_Mapping(
  global_obj = thisObj,
  sub_obj = GLN_subcluster_res$sub_obj,
  sample_id = "236",
  sample_col = "mus",
  target_cluster_name = "GLN",
  custom_col_map = gln_highlight_colors
)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[ASC: Subclustering]
ASC_subcluster_res <- Run_Subcluster_Check(
  global_obj = thisObj, 
  target_ident = "ASC", 
  resolution = 0.8
)

ASC_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = ASC_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "ASC Subclustering (pre-ID)",
  features = c("Gfap", "Sox9")
)

Plot_Subcluster_UMAP(
  seurat_obj = ASC_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "ASC Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(ASC_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

ASC_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = ASC_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("ASC", "GLN", "VASC", "MGL"), 
  custom_signatures = NULL
)

# 1. Create temporary merge column
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType <- as.character(ASC_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Astrocyte Cores
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("2", "5", "3", "1")] <- "ASC"

# 3. Micro-environmental potential segmentation shadows
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("8")] <- "ASC.MGL"
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("10")] <- "ASC.OPC"
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("6")] <- "ASC.OLG"
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("4")] <- "ASC.GLN"
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("9")] <- "ASC.CHP"
ASC_subcluster_res$sub_obj$ASC_SubClusterCellType[ASC_subcluster_res$sub_obj$ASC_SubClusterCellType %in% c("7")] <- "ASC.BEC"

Idents(ASC_subcluster_res$sub_obj) <- "ASC_SubClusterCellType"

Plot_Subcluster_UMAP(
  seurat_obj = ASC_subcluster_res$sub_obj,
  group_by_col = "ASC_SubClusterCellType",
  plot_title = "ASC Subclustering (Merged Annotations)"
)

Plot_Canonical_Lineage_Scores(ASC_subcluster_res$sub_obj, group_by_col = "ASC_SubClusterCellType")

QuickMarkers(seurat_obj = ASC_subcluster_res$sub_obj, target_ident = "ASC_SubClusterCellType")

Check_Global_Mapping(global_obj = thisObj, sub_obj = ASC_subcluster_res$sub_obj, target_cluster_name = "ASC", global_reduction = "tf_umap")

asc_highlight_colors <- c(
  "ASC"   = "red",  # Core population 
  "ASC.MGL" = "grey",  # Microglia bleed
  "ASC.OPC" = "grey",    # OPC bleed
  "ASC.OLG" = "grey",  # Oligodendrocyte bleed
  "ASC.GLN" = "blue",     # Granule cell bleed
  "ASC.CHP" = "grey",   # Choroid Plexus bleed (transcript dust)
  "ASC.VASC" = "grey"    # Vascular/Endothelial bleed
)

Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = ASC_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "ASC", 
                      custom_col_map = asc_highlight_colors)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[MGL: Subclustering]
MGL_subcluster_res <- Run_Subcluster_Check(
  global_obj = thisObj, 
  target_ident = "MGL", 
  resolution = 0.8
)

MGL_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = MGL_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "MGL Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(MGL_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

MGL_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = MGL_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("MGL", "VASC", "ASC"), 
  custom_signatures = NULL
)

# 1. Create temporary merge column
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType <- as.character(MGL_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Microglia Cores
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("2", "1")] <- "MGL"
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("4")] <- "IMM"

# 3. Micro-environmental potential segmentation shadows
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("3")] <- "MGL.OLG"
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("7")] <- "MGL.ASC"
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("5")] <- "MGL.MLI"
MGL_subcluster_res$sub_obj$MGL_SubClusterCellType[MGL_subcluster_res$sub_obj$MGL_SubClusterCellType %in% c("6", "8")] <- "MGL.GLN"

Idents(MGL_subcluster_res$sub_obj) <- "MGL_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = MGL_subcluster_res$sub_obj, group_by_col = "MGL_SubClusterCellType", plot_title = "MGL Subclustering (Merged Annotations)")

Plot_Canonical_Lineage_Scores(MGL_subcluster_res$sub_obj, group_by_col = "MGL_SubClusterCellType")

QuickMarkers(seurat_obj = MGL_subcluster_res$sub_obj, target_ident = "MGL_SubClusterCellType")

Check_Global_Mapping(global_obj = thisObj, sub_obj = MGL_subcluster_res$sub_obj, target_cluster_name = "MGL", global_reduction = "tf_umap")

mgl_highlight_colors <- c(
  "MGL"      = "blue",   # Core Microglia
  "IMM"      = "grey80", # Distinct Immune/Infiltrating population
  "MGL.OLG"  = "grey80",   # Oligodendrocyte interaction/bleed
  "MGL.ASC"  = "grey80",      # Astrocyte interaction/bleed
  "MGL.MLI"  = "red",    # Molecular Layer Interneuron interaction/bleed
  "MGL.GLN"  = "green"   # Granule cell interaction
)

Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = MGL_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "MGL", 
                      custom_col_map = mgl_highlight_colors)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[OLG: Subclustering]
OLG_subcluster_res <- Run_Subcluster_Check(
  global_obj = thisObj, 
  target_ident = "OLG", 
  resolution = 0.8
)

OLG_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = OLG_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "OLG Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(OLG_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

OLG_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = OLG_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("OLG", "ASC", "MGL"), 
  custom_signatures = NULL
)

# 1. Create temporary merge column
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType <- as.character(OLG_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Oligo Cores
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("1", "4", "8", "3", "7", "2", "6", "10")] <- "OLG"
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("16")] <- "OLG.OPC"

# 3. Micro-environmental potential segmentation shadows
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("15", "14")] <- "OLG.MLI"
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("5", "13", "9")] <- "OLG.GLN"
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("11")] <- "OLG.MGL"
OLG_subcluster_res$sub_obj$OLG_SubClusterCellType[OLG_subcluster_res$sub_obj$OLG_SubClusterCellType %in% c("12")] <- "OLG.BEC"


Idents(OLG_subcluster_res$sub_obj) <- "OLG_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = OLG_subcluster_res$sub_obj, group_by_col = "OLG_SubClusterCellType", plot_title = "OLG Subclustering (Merged Annotations)")

QuickMarkers(seurat_obj = OLG_subcluster_res$sub_obj, target_ident = "OLG_SubClusterCellType")

Check_Global_Mapping(global_obj = thisObj, sub_obj = OLG_subcluster_res$sub_obj, target_cluster_name = "OLG", global_reduction = "tf_umap")

olg_highlight_colors <- c(
  "OLG"      = "purple",   # Core Oligodendrocytes
  "OLG.OPC"  = "red",     # Oligodendrocyte Precursor Cells
  "OLG.MLI"  = "grey80",    # Molecular Layer Interneuron interaction
  "OLG.GLN"  = "green",   # Granule cell interaction
  "OLG.MGL"  = "grey80",      # Microglia interaction
  "OLG.VASC" = "grey80"     # Vascular/Endothelial interaction
)

Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = OLG_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "OLG", 
                      custom_col_map = olg_highlight_colors)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[MLI: Subclustering]
MLI_subcluster_res <- Run_Subcluster_Check(global_obj = thisObj, target_ident = "MLI", resolution = 0.8)

MLI_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = MLI_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "MLI Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(MLI_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

MLI_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = MLI_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("MLI", "ASC", "MGL"), 
  custom_signatures = NULL
)


# 1. Create the temporary merge column
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType <- as.character(MLI_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Oligo Cores
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("1", "2", "3", "11")] <- "MLI"
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("7", "5", "6")] <- "BG"

# 3. Micro-environmental potential segmentation shadows
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("13")] <- "MLI.MGL"
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("15", "12", "4")] <- "MLI.ASC"
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("9")] <- "MLI.BEC"
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("14")] <- "MLI.OPC"
MLI_subcluster_res$sub_obj$MLI_SubClusterCellType[MLI_subcluster_res$sub_obj$MLI_SubClusterCellType %in% c("8", "10")] <- "MLI.OLG"

# 3. Apply to Idents and plot
Idents(MLI_subcluster_res$sub_obj) <- "MLI_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = MLI_subcluster_res$sub_obj, group_by_col = "MLI_SubClusterCellType", plot_title = "MLI Subclustering")

QuickMarkers(seurat_obj = MLI_subcluster_res$sub_obj,
             target_ident = "MLI_SubClusterCellType")


plot_cluster_diagnostics(MLI_subcluster_res$sub_obj, 
                         ident_name = "MLI_SubClusterCellType", 
                         sample_id_column = "mus")


Check_Global_Mapping(global_obj = thisObj, sub_obj = MLI_subcluster_res$sub_obj, target_cluster_name = "MLI")

# 4. Spatial Check
mli_highlight_colors <- c(
  "MLI"      = "grey",   # Core Molecular Layer Interneurons
  "MLI.ASC"  = "grey",   # Interaction with Bergmann Glia / Astrocytes
  "MLI.MGL"  = "red",   # Microglia interaction
  "MLI.OPC"  = "grey",     # OPC interaction
  "MLI.OLG"  = "grey",     # Oligodendrocyte interaction
  "MLI.VASC" = "grey" # Vascular interaction
)
Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = MLI_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "MLI", 
                      custom_col_map = mli_highlight_colors)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[PKJ: Subclustering]
PKJ_subcluster_res <- Run_Subcluster_Check(global_obj = thisObj, 
                                           target_ident = "PKJ", 
                                           resolution = 0.8)

PKJ_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = PKJ_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "PKJ Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(PKJ_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

PKJ_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = PKJ_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("PKJ", "ASC", "MGL"), 
  custom_signatures = NULL
)


# 1. Create the temporary merge column
PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType <- as.character(PKJ_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Oligo Cores
PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType[PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType %in% c("5", "1", "8")] <- "PKJ"

# 3. Micro-environmental potential segmentation shadows
PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType[PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType %in% c("7", "3")] <- "PKJ.GLN"
PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType[PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType %in% c("2")] <- "PKJ.ASC"
PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType[PKJ_subcluster_res$sub_obj$PKJ_SubClusterCellType %in% c("9", "6", "4")] <- "PKJ.OLG"

# 3. Apply to Idents and plot
Idents(PKJ_subcluster_res$sub_obj) <- "PKJ_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = PKJ_subcluster_res$sub_obj, 
                     group_by_col = "PKJ_SubClusterCellType", 
                     plot_title = "PKJ Subclustering")

QuickMarkers(seurat_obj = PKJ_subcluster_res$sub_obj,
             target_ident = "PKJ_SubClusterCellType")


plot_cluster_diagnostics(PKJ_subcluster_res$sub_obj, 
                         ident_name = "PKJ_SubClusterCellType", 
                         sample_id_column = "mus")


Check_Global_Mapping(global_obj = thisObj, sub_obj = PKJ_subcluster_res$sub_obj, target_cluster_name = "PKJ")

# 4. Spatial Check
PKJ_highlight_colors <- c(
  "PKJ"      = "grey80",   # Core PKJ
  "PKJ.ASC"  = "red",   # Interaction with Bergmann Glia / Astrocytes
  "PKJ.OLG"  = "blue",     # Oligodendrocyte interaction
  "PKJ.VASC" = "green" # Vascular interaction
)
Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = PKJ_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "PKJ", 
                      custom_col_map = PKJ_highlight_colors)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[OPC: Subclustering]

OPC_subcluster_res <- Run_Subcluster_Check(global_obj = thisObj, 
                                           target_ident = "OPC", 
                                           resolution = 0.8)

OPC_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = OPC_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "OPC Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(OPC_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

OPC_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = OPC_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("OPC", "ASC", "MGL"), 
  custom_signatures = NULL
)


# 1. Create the temporary merge column
OPC_subcluster_res$sub_obj$OPC_SubClusterCellType <- as.character(OPC_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True Oligo Cores
OPC_subcluster_res$sub_obj$OPC_SubClusterCellType[OPC_subcluster_res$sub_obj$OPC_SubClusterCellType %in% c("1", "2", "3")] <- "OPC"


# 3. Apply to Idents and plot
Idents(OPC_subcluster_res$sub_obj) <- "OPC_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = OPC_subcluster_res$sub_obj, 
                     group_by_col = "OPC_SubClusterCellType", 
                     plot_title = "OPC Subclustering")

QuickMarkers(seurat_obj = OPC_subcluster_res$sub_obj,
             target_ident = "OPC_SubClusterCellType")


plot_cluster_diagnostics(OPC_subcluster_res$sub_obj, 
                         ident_name = "OPC_SubClusterCellType", 
                         sample_id_column = "mus")


Check_Global_Mapping(global_obj = thisObj, sub_obj = OPC_subcluster_res$sub_obj, target_cluster_name = "OPC")

# 4. Spatial Check
OPC_highlight_colors <- c(
  "OPC"      = "red"

)
Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = OPC_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "OPC", 
                      custom_col_map = OPC_highlight_colors)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[DCN: Subclustering]

DCN_subcluster_res <- Run_Subcluster_Check(global_obj = thisObj, 
                                           target_ident = "DCN", 
                                           resolution = 0.8)

DCN_subcluster_res$top10

Plot_Subcluster_UMAP(
  seurat_obj = DCN_subcluster_res$sub_obj,
  group_by_col = "SCT_snn_res.0.8",
  plot_title = "DCN Subclustering (pre-ID)"
)

Plot_Canonical_Lineage_Scores(DCN_subcluster_res$sub_obj, group_by_col = "SCT_snn_res.0.8")

DCN_AnchorClusters <- Run_Phase1_AnchorSearch(
  seurat_obj = DCN_subcluster_res$sub_obj, 
  reduction = "umap_sct_sub",
  group_by_col = "SCT_snn_res.0.8", 
  target_lineages = c("DCN", "ASC", "MGL"), 
  custom_signatures = NULL
)


# 1. Create the temporary merge column
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType <- as.character(DCN_subcluster_res$sub_obj$SCT_snn_res.0.8)

# 2. True DCN Cores
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("1", "2", "7")] <- "DCN"

# 3. Micro-environmental potential segmentation shadows
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("3")] <- "DCN.MLI"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("6")] <- "DCN.OLG"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("4", "5")] <- "DCN.GLN"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("9")] <- "DCN.BEC"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("8")] <- "DCN.VLM"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("10")] <- "DCN.MGL"
DCN_subcluster_res$sub_obj$DCN_SubClusterCellType[DCN_subcluster_res$sub_obj$DCN_SubClusterCellType %in% c("11")] <- "DCN.OPC"

# 3. Apply to Idents and plot
Idents(DCN_subcluster_res$sub_obj) <- "DCN_SubClusterCellType"

Plot_Subcluster_UMAP(seurat_obj = DCN_subcluster_res$sub_obj, 
                     group_by_col = "DCN_SubClusterCellType", 
                     plot_title = "DCN Subclustering")

QuickMarkers(seurat_obj = DCN_subcluster_res$sub_obj,
             target_ident = "DCN_SubClusterCellType")


plot_cluster_diagnostics(DCN_subcluster_res$sub_obj, 
                         ident_name = "DCN_SubClusterCellType", 
                         sample_id_column = "mus")


Check_Global_Mapping(global_obj = thisObj, sub_obj = DCN_subcluster_res$sub_obj, target_cluster_name = "DCN")

# 4. Spatial Check
DCN_highlight_colors <- c(
  "DCN"      = "red",      # Core DCN Neurons
  "DCN.MLI"  = "green",    # Interaction with Molecular Layer Interneurons (or similar inhibitory)
  "DCN.OLG"  = "orange",   # Interaction with Oligodendrocytes (White Matter background)
  "DCN.GLN"  = "purple",   # Interaction with Granule cells
  "DCN.VASC" = "blue",     # Vascular/Endothelial bleed
  "DCN.VLM"  = "magenta",  # Leptomeningeal/Fibroblast bleed
  "DCN.MGL"  = "darkblue", # Microglia interaction
  "DCN.OPC"  = "cyan"      # OPC interaction
)
Check_Spatial_Mapping(global_obj = thisObj, 
                      sub_obj = DCN_subcluster_res$sub_obj, 
                      sample_id = "YM04", sample_col = "mus", 
                      target_cluster_name = "DCN", 
                      custom_col_map = DCN_highlight_colors)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[CHP: Subclustering]
#skipping
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Segmentation Shadow Reconstructions | Pre-FastReseg ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Segmentation Collision and Shadow Cell Checks*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

#'[Global Metadata Integration]

# 1. Create a list of your sub-objects and their respective metadata columns
sub_obj_list <- list(
  VASC = list(obj = VASC_subcluster_res$sub_obj, col = "VASC_SubClusterCellType"),
  GLN  = list(obj = GLN_subcluster_res$sub_obj,  col = "GLN_SubClusterCellType"),
  ASC  = list(obj = ASC_subcluster_res$sub_obj,  col = "ASC_SubClusterCellType"),
  MGL  = list(obj = MGL_subcluster_res$sub_obj,  col = "MGL_SubClusterCellType"),
  OLG  = list(obj = OLG_subcluster_res$sub_obj,  col = "OLG_SubClusterCellType"),
  MLI  = list(obj = MLI_subcluster_res$sub_obj,  col = "MLI_SubClusterCellType"),
  PKJ  = list(obj = PKJ_subcluster_res$sub_obj,  col = "PKJ_SubClusterCellType"),
  OPC  = list(obj = OPC_subcluster_res$sub_obj,  col = "OPC_SubClusterCellType"),
  DCN  = list(obj = DCN_subcluster_res$sub_obj,  col = "DCN_SubClusterCellType")
)

# 2. Initialize the new column in the global object with the base global identities
# (Assuming your base global identities are stored in "seurat_clusters" or similar)
thisObj$HighRes_Annotation <- as.character(Idents(thisObj)) 

# 3. Loop through and overwrite with the high-resolution labels
for (lineage in names(sub_obj_list)) {
  
  # Extract the specific sub-object and its target column
  current_obj <- sub_obj_list[[lineage]]$obj
  target_col <- sub_obj_list[[lineage]]$col
  
  # Extract the cell IDs and their new labels
  cell_ids <- rownames(current_obj@meta.data)
  new_labels <- as.character(current_obj@meta.data[[target_col]])
  
  # Map them into the global object
  thisObj@meta.data[cell_ids, "HighRes_Annotation"] <- new_labels
}

# Set the active identity to your new, fully mapped column
Idents(thisObj) <- "HighRes_Annotation"
cat("Global merge complete. Total unique high-res populations:", length(unique(thisObj$HighRes_Annotation)), "\n")



thisObj@meta.data$x_slide_micron <- thisObj@meta.data$x_slide_mm * 1000
thisObj@meta.data$y_slide_micron <- thisObj@meta.data$y_slide_mm * 1000


library(FNN)

Verify_Spatial_Shadows_Dynamic <- function(seurat_obj, 
                                           annot_col = "HighRes_Annotation", 
                                           spatial_cols = c("x_slide_mm", "y_slide_mm"), # Update if yours are named differently
                                           area_col = "Area.um2",
                                           sample_col = "mus",
                                           tolerance_um = 1.5, # Buffer for irregular cell shapes
                                           lineage_map = NULL) { 
  
  cat("\n=== Initiating Dynamic Collision Radius Verification ===\n")
  
  meta <- seurat_obj@meta.data
  
  # Ensure the area column exists
  if (!area_col %in% colnames(meta)) {
    stop(sprintf("Error: Area column '%s' not found in metadata.", area_col))
  }
  
  # 1. Calculate the dynamic physical radius for every cell
  # NOTE: If your spatial_cols are in millimeters (like x_slide_mm), we must convert the distance math to match the Area.um2 unit, or vice versa.
  # Assuming spatial coords are in microns (or pixels converted to microns) for this math to work natively.
  # If your spatial coords are in mm, multiply them by 1000 before running this, or convert them inside the function.
  
  meta$Radius_um <- sqrt(meta[[area_col]] / pi)
  meta$CellID <- rownames(meta)
  
  # 2. Define the Global Umbrella Mapping
  if (is.null(lineage_map)) {
    lineage_map <- list(
      VASC = c("BEC", "VLM", "VSM", "PER", "VASC"), 
      GLN  = c("GLN"),
      ASC  = c("ASC"),
      MGL  = c("MGL", "IMM"),
      OLG  = c("OLG", "OPC"),
      MLI  = c("MLI", "BG"),
      PKJ  = c("PKJ"),
      DCN  = c("DCN"),
      CHP  = c("CHP")
    )
  }
  
  reverse_map <- list()
  for (lin in names(lineage_map)) {
    for (sub_id in lineage_map[[lin]]) {
      reverse_map[[sub_id]] <- lin
    }
  }
  
  # 3. Identify Shadows
  meta <- meta %>%
    mutate(
      Is_Shadow = str_detect(.data[[annot_col]], "\\."),
      Primary_ID = ifelse(Is_Shadow, str_split_i(.data[[annot_col]], "\\.", 1), .data[[annot_col]]),
      Contaminant_ID = ifelse(Is_Shadow, str_split_i(.data[[annot_col]], "\\.", 2), NA),
      Shadow_Status = ifelse(Is_Shadow, "Unverified", "Pure Core"),
      Nearest_Contaminant_Dist = NA_real_,
      Dynamic_Collision_Limit = NA_real_
    )
  
  samples <- unique(meta[[sample_col]])
  
  # 4. Iterate and Compute Collisions
  for (s in samples) {
    sample_indices <- which(meta[[sample_col]] == s)
    sample_meta <- meta[sample_indices, ]
    
    shadow_types_in_sample <- unique(sample_meta[[annot_col]][sample_meta$Is_Shadow])
    
    for (st in shadow_types_in_sample) {
      target_cells <- sample_meta[sample_meta[[annot_col]] == st, ]
      contaminant_id <- unique(target_cells$Contaminant_ID)
      
      parent_lineage <- reverse_map[[contaminant_id]]
      acceptable_anchors <- if(is.null(parent_lineage)) contaminant_id else lineage_map[[parent_lineage]]
      
      anchor_cells <- sample_meta[sample_meta[[annot_col]] %in% acceptable_anchors, ]
      
      if (nrow(anchor_cells) == 0) {
        meta[target_cells$CellID, "Shadow_Status"] <- "False Shadow (No Local Anchors)"
        next
      }
      
      # FNN Spatial Search
      query_coords <- as.matrix(target_cells[, spatial_cols])
      ref_coords <- as.matrix(anchor_cells[, spatial_cols])
      
      nn <- get.knnx(data = ref_coords, query = query_coords, k = 1)
      
      # Extract distances and the specific radius of the nearest anchor
      dists <- nn$nn.dist[, 1]
      nearest_anchor_radii <- anchor_cells$Radius_um[nn$nn.index[, 1]]
      query_radii <- target_cells$Radius_um
      
      # Calculate the dynamic collision limit for each individual cell pair
      collision_limits <- query_radii + nearest_anchor_radii + tolerance_um
      
      meta[target_cells$CellID, "Nearest_Contaminant_Dist"] <- dists
      meta[target_cells$CellID, "Dynamic_Collision_Limit"] <- collision_limits
      
      # Evaluate collision
      true_ids <- target_cells$CellID[dists <= collision_limits]
      false_ids <- target_cells$CellID[dists > collision_limits]
      
      meta[true_ids, "Shadow_Status"] <- "True Shadow (Collision Detected)"
      meta[false_ids, "Shadow_Status"] <- "False Shadow (Isolated)"
    }
  }
  
  seurat_obj$Shadow_Status <- meta$Shadow_Status
  seurat_obj$Nearest_Contaminant_Dist <- meta$Nearest_Contaminant_Dist
  seurat_obj$Dynamic_Collision_Limit <- meta$Dynamic_Collision_Limit
  
  # ---------------------------------------------------------
  # GENERATE QC PLOTS
  # ---------------------------------------------------------
  plot_df <- meta %>% filter(Is_Shadow == TRUE)
  
  p1 <- ggplot(plot_df, aes(x = .data[[annot_col]], fill = Shadow_Status)) +
    geom_bar() +
    scale_fill_manual(values = c("True Shadow (Collision Detected)" = "#3182bd", 
                                 "False Shadow (Isolated)" = "#e6550d", 
                                 "False Shadow (No Local Anchors)" = "black")) +
    theme_classic() + labs(title = "Collision Verification Summary", x = "", y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "bottom")
  
  # Scatter Plot: Distance vs. Dynamic Threshold
  p2 <- ggplot(filter(plot_df, !is.na(Nearest_Contaminant_Dist)), 
               aes(x = Nearest_Contaminant_Dist, y = Dynamic_Collision_Limit, color = Shadow_Status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
    scale_color_manual(values = c("True Shadow (Collision Detected)" = "#3182bd", "False Shadow (Isolated)" = "#e6550d")) +
    theme_classic() + 
    labs(title = "Cell-to-Cell Collision Boundary", 
         subtitle = "Points above the dashed line are physically colliding",
         x = "Distance to Nearest Anchor (um)", 
         y = "Calculated Collision Limit (r1 + r2 + tol)") +
    theme(legend.position = "none")
  
  print(p1 + p2 + plot_layout(widths = c(1, 1.2)))
  
  cat("\n--- Shadow Verification Audit ---\n")
  print(table(meta[[annot_col]][meta$Is_Shadow], meta$Shadow_Status[meta$Is_Shadow]))
  
  return(invisible(seurat_obj))
}




thisObj.T <- Verify_Spatial_Shadows_Dynamic(
  seurat_obj = thisObj,
  annot_col = "HighRes_Annotation",
  spatial_cols = c("x_slide_micron", "y_slide_micron"), 
  sample_col = "mus",
  area_col = "Area.um2",
  tolerance_um = 5
)

#'[NOTE: this analysis confirms that the fragmented cell data was likely filtered from the original seurat]
#'[So we will use the "pure" clusters as reference data for a re-segmentation algorithim reference dataset]

# Re-segmentation of fragmented cells FastReseg Step ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Re-segmentation of fragmented cells*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
# Initial "shadow"/segmentation error analysis proves that either A) the other
# "half" of a shadow cell was likely filtered out because of their area or lack
# of transcripts or B) entirely undetected due to depth/cell overlaying each other
# We can't capture cells in the latter category, but we can use Nanostring's Reseg
# to try to reconstruct/correct the fragmented cells. We'll use the "pure" clusters
# of cells such as PER, ASC, etc as the "reference set". The fastReseg will
# use the pure reference to reconstruct the cells and try to salvage the poorly
# segmented cells. Then we'll re-run the entire analysis with the refined 
# transcript localizations and new polygonal segmentations.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table) # Using fread for massive transcript files
library(FastReseg)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Functions: FastReseg Helper Fxns]
unique(thisObj$mus) #view all sample ids

Generate_FastReseg_Dictionary <- function(seurat_obj, annot_col = "HighRes_Annotation") {
  
  cat("\n=== Building FastReseg Reference Dictionary ===\n")
  
  meta <- seurat_obj@meta.data
  
  # 1. Identify only the "Pure" clusters (No dots in the name)
  meta$Is_Shadow <- grepl("\\.", meta[[annot_col]])
  pure_cell_ids <- rownames(meta[!meta$Is_Shadow, ])
  
  cat(sprintf("Extracted %d Pure Core cells. Discarding %d Shadow/Mixed cells.\n", 
              length(pure_cell_ids), nrow(meta) - length(pure_cell_ids)))
  
  # 2. Subset the Seurat Object
  pure_obj <- subset(seurat_obj, cells = pure_cell_ids)
  Idents(pure_obj) <- annot_col
  
  # 3. Calculate Average Expression Profile
  # Using the SCT data slot gives us library-size normalized, log-transformed counts
  # AverageExpression automatically exponentiates them back to non-log space
  cat("Calculating transcriptomic averages per lineage...\n")
  avg_exp <- AverageExpression(pure_obj, assays = "SCT", slot = "data")$SCT
  
  # Clean up column names (AverageExpression prepends the assay name)
  colnames(avg_exp) <- gsub("^SCT\\.", "", colnames(avg_exp))
  
  cat("\nDictionary successfully generated for the following lineages:\n")
  print(colnames(avg_exp))
  
  return(as.data.frame(avg_exp))
}

Run_FastReseg_Optimized <- function(tx_csv_path, target_sample, seurat_obj, reference_dict) {
  
  cat(sprintf("\n=== Preparing Data Architecture for Sample: %s ===\n", target_sample))
  
  valid_fovs <- seurat_obj@meta.data %>% filter(mus == target_sample) %>% pull(fov) %>% unique()
  
  cat("Loading massive transcript file...\n")
  raw_tx <- fread(tx_csv_path)
  clean_tx <- raw_tx[raw_tx$fov %in% valid_fovs, ]
  rm(raw_tx); gc() 
  
  cat("Sanitizing Background Transcript IDs...\n")
  clean_tx$Fixed_Cell_ID <- ifelse(clean_tx$cell_ID == 0, "0", clean_tx$cell)
  
  
  cat("Applying microscopic joggle to coordinates to bypass 3D geometry constraints...\n")
  
  # Adds random noise between -0.001 and +0.001. 
  # This destroys perfect planes mathematically but preserves them biologically.
  set.seed(42) # For reproducibility
  clean_tx$x_global_px <- clean_tx$x_global_px + runif(nrow(clean_tx), min = -1e-3, max = 1e-3)
  clean_tx$y_global_px <- clean_tx$y_global_px + runif(nrow(clean_tx), min = -1e-3, max = 1e-3)
  clean_tx$z <- clean_tx$z + runif(nrow(clean_tx), min = -1e-3, max = 1e-3)
  
  
  temp_dir <- paste0("Temp_FOVs_", target_sample)
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  
  cat("Chunking transcripts to temporary FOV files...\n")
  fov_list <- unique(clean_tx$fov)
  file_paths <- c()
  
  for (f in fov_list) {
    fov_data <- clean_tx[clean_tx$fov == f, ]
    f_path <- file.path(temp_dir, paste0("FOV_", f, ".csv"))
    fwrite(fov_data, f_path)
    file_paths <- c(file_paths, f_path)
  }
  rm(clean_tx); gc() 
  
  transDF_fileInfo <- data.frame(
    file_path = file_paths,
    slide = target_sample,
    fov = fov_list,
    stage_X = 0, 
    stage_Y = 0  
  )
  
  cat("Preparing Sparse Matrices...\n")
  sample_cells <- rownames(seurat_obj@meta.data[seurat_obj$mus == target_sample, ])
  counts_fixed <- as.matrix(t(seurat_obj@assays$RNA@counts[, sample_cells])) 
  ref_matrix <- as.matrix(reference_dict)
  
  cat("\nInitiating FastReseg Pipeline (Iterative FOV Mode)...\n")
  
  reseg_results <- fastReseg_full_pipeline(
    counts = counts_fixed,
    clust = NULL, 
    refProfiles = ref_matrix,
    
    transcript_df = NULL, 
    transDF_fileInfo = transDF_fileInfo,
    filepath_coln = "file_path",
    prefix_colns = NULL,  
    
    # Restored "z" to satisfy pp3!
    spatLocs_colns = c("x_global_px", "y_global_px", "z"), 
    fovOffset_colns = c("stage_X", "stage_Y"), 
    pixel_size = 0.12,      
    zstep_size = 0.8,       
    
    transID_coln = NULL, 
    transGene_coln = "target",
    cellID_coln = "Fixed_Cell_ID", 
    extracellular_cellID = c("0", 0),
    
    molecular_distance_cutoff = 2,
    cellular_distance_cutoff = 20,
    imputeFlag_missingCTs = TRUE,
    
    flagCell_lrtest_cutoff = 5,
    svmClass_score_cutoff = -2,
    groupTranscripts_method = "dbscan",
    spatialMergeCheck_method = "leidenCut",
    cutoff_spatialMerge = 0.5,
    
    path_to_output = paste0("FastReseg_Output_", target_sample),
    save_intermediates = FALSE, 
    return_perCellData = TRUE, 
    combine_extra = FALSE 
  )
  
  unlink(temp_dir, recursive = TRUE)
  
  cat(sprintf("\n=== Pipeline Complete for %s ===\n", target_sample))
  return(reseg_results)
}


fastreseg_dictionary <- Generate_FastReseg_Dictionary(seurat_obj = thisObj, 
                                                    annot_col = "HighRes_Annotation")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Analysis: Single sample test for FastReseg]
path_301_YM01 <- "C:/Users/m322332/Documents/Unified/Aging_Tx_data/01_Schafer_M1_301_YM01_tx_file.csv"

reseg_301 <- Run_FastReseg_Optimized(
  tx_csv_path = path_301_YM01,
  target_sample = "301",
  seurat_obj = thisObj,
  reference_dict = fastreseg_dictionary
)


#'[Figure: Single sample quick QC for FastReseg results]

Plot_Initial_Rescue_QC <- function(original_seurat, reseg_results, target_sample = "301") {
  
  cat("Extracting Original vs. Resegmented Counts...\n")
  
  # 1. Extract Original Counts
  original_cells <- original_seurat@meta.data %>% filter(mus == target_sample)
  old_counts <- original_cells$nCount_RNA
  
  # 2. Extract New Counts
  # The updated_perCellExprs matrix is Gene x Cell. We sum the columns to get nCount.
  new_counts <- Matrix::colSums(reseg_results$updated_perCellExprs)
  
  # 3. Create a combined data frame
  qc_df <- data.frame(
    nCount = c(old_counts, new_counts),
    State = factor(c(rep("Original", length(old_counts)), 
                     rep("Resegmented", length(new_counts))),
                   levels = c("Original", "Resegmented"))
  )
  
  # 4. Calculate Medians for the overlay
  medians <- qc_df %>% group_by(State) %>% summarize(med = median(nCount))
  
  # 5. Plot 1: Density Overlay
  p_density <- ggplot(qc_df, aes(x = nCount, fill = State, color = State)) +
    geom_density(alpha = 0.5, linewidth = 1) +
    geom_vline(data = medians, aes(xintercept = med, color = State), linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = c("Original" = "#e6550d", "Resegmented" = "#3182bd")) +
    scale_color_manual(values = c("Original" = "#e6550d", "Resegmented" = "#3182bd")) +
    scale_x_continuous(limits = c(0, 1000)) + # Cap x-axis to zoom in on the main distribution
    theme_classic() +
    labs(title = paste("Transcript Rescue Distribution:", target_sample),
         subtitle = "Dashed line indicates median nCount",
         x = "Total Transcripts per Cell (nCount_RNA)", y = "Density")
  
  # 6. Plot 2: Action Breakdown Bar Chart
  action_df <- as.data.frame(table(reseg_results$updated_perCellDT$reSeg_action))
  colnames(action_df) <- c("Action", "Count")
  
  p_action <- ggplot(action_df, aes(x = Action, y = Count, fill = Action)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = Count), vjust = -0.5) +
    theme_classic() +
    labs(title = "FastReseg Action Breakdown",
         x = "Resegmentation Action", y = "Number of Cells") +
    theme(legend.position = "none")
  
  print(p_density | p_action)
}

# Run the QC on the results you just got
Plot_Initial_Rescue_QC(thisObj, reseg_301, target_sample = "301")


#'[Analysis: Loop for Batch FastReseg]
#**Warning: This analysis takes a lot of time and storage space.*

tx_directory <- "C:/Users/m322332/Documents/Unified/Aging_Tx_data/" #transcript csv directory
all_samples <- unique(thisObj@meta.data$mus)

cat(sprintf("Found %d total samples to process.\n", length(all_samples)))

file_map <- list(
  "301" = "01_Schafer_M1_301_YM01_tx_file.csv",
  "YM01" = "01_Schafer_M1_301_YM01_tx_file.csv", 
  "302" = "02_Schafer_M2_302_YM02_tx_file.csv",
  "YM02" = "02_Schafer_M2_302_YM02_tx_file.csv",
  "304" = "03_Schafer_M3_304_YM03_tx_file.csv",
  "YM03"  = "03_Schafer_M3_304_YM03_tx_file.csv",
  "305" = "04_Schafer_M4_305_YM04_tx_file.csv",
  "YM04" = "04_Schafer_M4_305_YM04_tx_file.csv",
  "236" = "05_Schafer_F1_236_YF01_tx_file.csv",
  "YF01" = "05_Schafer_F1_236_YF01_tx_file.csv",
  "237" = "06_Schafer_F2_237_YF02_tx_file.csv",
  "YF02" = "06_Schafer_F2_237_YF02_tx_file.csv",
  "239" = "07_Schafer_F3_239_YF03_tx_file.csv",
  "YF03" = "07_Schafer_F3_239_YF03_tx_file.csv",
  "243" = "08_Schafer_F4_243_YF04_tx_file.csv",
  "YF04" = "08_Schafer_F4_243_YF04_tx_file.csv"
)


all_reseg_results <- list()


for (current_sample in all_samples) {
  
  cat(paste0("\n======================================================\n"))
  cat(sprintf("CHECKING SAMPLE: %s\n", current_sample))
  cat(paste0("======================================================\n"))
  
  # Define the exact file name we expect to see if it finished successfully
  save_name <- paste0("Reseg_Results_", current_sample, ".rds")
  

  # Resume check, if loop crashed, this will skip samples that are already processed
  if (file.exists(save_name)) {
    cat(sprintf("[SKIP] Found existing file: %s\n", save_name))
    cat("-> Loading saved data into master list and moving to next sample...\n")
    
    # If sample has been processed, add it back to the results list for the final results obj save
    all_reseg_results[[current_sample]] <- readRDS(save_name)
    next # move to next sample
  }

  # Look up the file path samples
  file_name <- file_map[[current_sample]]
  
  if (is.null(file_name)) {
    cat(sprintf("[!] ERROR: No file path mapped for sample %s. Skipping...\n", current_sample))
    next 
  }
  
  full_tx_path <- file.path(tx_directory, file_name)
  
  # Run the FastReseg full pipeline wrapper
  sample_results <- tryCatch({
    
    Run_FastReseg_Optimized(
      tx_csv_path = full_tx_path,
      target_sample = current_sample,
      seurat_obj = thisObj,
      reference_dict = fastreseg_dictionary
    )
    
  }, error = function(e) {
    cat(sprintf("\n[!] FATAL ERROR processing %s: %s\n", current_sample, e$message))
    return(NULL) 
  })
  
  # Save successful runs immediately in case a crash occurs
  if (!is.null(sample_results)) {
    saveRDS(sample_results, file = save_name)
    cat(sprintf("\n-> Successfully computed and saved %s to disk.\n", save_name))
    
    # Store in the results list
    all_reseg_results[[current_sample]] <- sample_results
  }
}

cat("\n=== ALL SAMPLES HAVE BEEN PROCESSED AND COMPILED ===\n")

# Save the final combined results list containing all samples
saveRDS(all_reseg_results, file = paste0(local_WD, "Master_FastReseg_All_Samples.rds"))






library(Seurat)
library(Matrix)
library(dplyr)

library(Seurat)
library(Matrix)
library(dplyr)

Compile_Global_Reseg_Seurat <- function(reseg_list, nCount_cutoff = 100) {
  
  cat("\n=== Compiling Global Post-Resegmentation Seurat Object ===\n")
  
  all_counts <- list()
  all_meta <- list()
  
  samples <- names(reseg_list)
  
  for (s in samples) {
    cat(sprintf("Processing Sample %s...\n", s))
    
    res <- reseg_list[[s]]
    
    # 1. Extract the Count Matrix
    counts <- res$updated_perCellExprs
    
    # 2. Extract and format the Metadata
    meta <- as.data.frame(res$updated_perCellDT)
    
  
    # THE FIX: Make Cell IDs globally unique

    # Prepend the sample name to the cell ID (e.g., "301_c_1_28_1")
    unique_ids <- paste0(s, "_", meta$updated_cellID)
    
    # Apply unique IDs to metadata
    meta$updated_cellID <- unique_ids
    rownames(meta) <- unique_ids
    
    # Apply unique IDs to the count matrix columns
    colnames(counts) <- unique_ids

    
    # Add the sample ID column
    meta$mus <- s
    
    # Store in our holding lists
    all_counts[[s]] <- counts
    all_meta[[s]] <- meta
  }
  
  cat("\nStitching 16 Count Matrices together (this may take a moment)...\n")
  global_counts <- do.call(cbind, all_counts)
  
  cat("Stitching Metadata together...\n")
  global_meta <- do.call(rbind, all_meta)
  
  # Safety check
  rownames(global_meta) <- colnames(global_counts)
  
  cat("\nInitializing pristine Seurat Object...\n")
  global_obj <- CreateSeuratObject(counts = global_counts, 
                                   meta.data = global_meta, 
                                   project = "FastReseg_Aging")
  
  # --- THE FILTERING STEP ---
  cat(sprintf("\n=== Applying Fragment Filter (nCount_RNA > %d) ===\n", nCount_cutoff))
  
  cells_before <- ncol(global_obj)
  global_obj_filtered <- subset(global_obj, subset = nCount_RNA > nCount_cutoff)
  cells_after <- ncol(global_obj_filtered)
  
  cat(sprintf("Starting Cells: %d\n", cells_before))
  cat(sprintf("Retained Cells: %d\n", cells_after))
  cat(sprintf("Removed Fragments: %d\n", cells_before - cells_after))
  
  return(global_obj_filtered)
}

# --- Execute the Compilation ---
# Pass your master list and define your new cutoff (e.g., 100 or 150 based on your QC plots)
thisObj_reseg <- Compile_Global_Reseg_Seurat(all_reseg_results, nCount_cutoff = 200)

# Save the pristine object to disk immediately!
saveRDS(thisObj_reseg, file = paste0(local_WD, "thisObj_Resegmented_Filtered.rds"))


thisObj_reseg

updated_celltype

















# Post-Reseg | TF Based Clustering ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**TF Based Clustering*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

project_WD <- '//mfad/researchmn/PHYSIO/SCHAFER/SHARED/JamesRedden/R Code/CEREBELLUM/'
local_WD <- "C:/Users/m322332/Documents/"
HDD_WD <- "E:/"
setwd(local_WD)
#'[Subset the obj for testing]
# fraction_to_keep <- 1 #set % of cells to use in subset
# total_cells <- ncol(thisObj_raw)
# cells_to_sample <- floor(fraction_to_keep * total_cells)
# set.seed(42)
# random_cells <- sample(colnames(thisObj_raw), size = cells_to_sample)
# thisObj_subset <- subset(thisObj_raw, cells = random_cells)
# print(paste("Original object size:", ncol(thisObj_raw), "cells"))
# print(paste("Subset object size:", ncol(thisObj_subset), "cells"))
# 
# thisObj <- thisObj_subset

library(Seurat)
library(decoupleR)
library(tidyverse)
library(ggrepel)
library(matrixStats)

thisObj_reseg <- readRDS(paste0(local_WD, "thisObj_Resegmented_Filtered.rds"))

thisObj_reseg

#'[Generate TF enrichment results]
# Query for mouse TF activity
net <- get_collectri(organism = 'mouse', split_complexes=FALSE)
net

NormalizeData(thisObj_reseg) #, assay = "RNA")
Norm_counts <- as.matrix(GetAssayData(object = thisObj_reseg)) #, assay = "RNA", layer = "data"))

# Run ulm of TFs
acts <- run_ulm(mat=Norm_counts, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts

acts_mat <- acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t() # Transpose to get (TFs x Cells)

thisObj_reseg[['tf_activities']] <- CreateAssayObject(counts = acts_mat)


# Filtering out cells whose max abs TF score is too low to confidently cell-type

#' # Extract the ULM matrix
#' tf_data <- as.matrix(GetAssayData(thisObj_reseg, assay = "tf_activities", layer = "data"))
#' 
#' # Calculate absolute max t-score per cell
#' thisObj_reseg$Max_TF_Score <- matrixStats::colMaxs(abs(tf_data))
#' 
#' # Extract plot data
#' plot_data <- data.frame(Max_TF_Score = thisObj_reseg$Max_TF_Score)
#' 
#' # Calculate Median and MAD
#' med_score <- median(plot_data$Max_TF_Score, na.rm = TRUE)
#' mad_score <- mad(plot_data$Max_TF_Score, na.rm = TRUE)
#' 
#' #'[Figure: TF-activity score distribution plot]
#' p_dist <- ggplot(plot_data, aes(x = Max_TF_Score)) +
#'   geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "cadetblue3", color = "black", alpha = 0.6) +
#'   geom_density(color = "cadetblue4", linewidth = 1.2) +
#'   
#'   # Add Median and MAD markers
#'   geom_vline(xintercept = med_score, color = "darkgreen", linetype = "solid", linewidth = 1) +
#'   
#'   # Add the cutoff line
#'   geom_vline(xintercept = 2.0, color = "red", linetype = "dashed", linewidth = 1.2) +
#'   
#'   theme_minimal() +
#'   labs(title = "Distribution of Maximum TF Scores per Cell",
#'        subtitle = paste0("Median: ", round(med_score, 2), " | MAD: ", round(mad_score, 2)),
#'        x = "Maximum Absolute TF Score (t-statistic)",
#'        y = "Density") +
#'   
#'   # Annotations
#'   annotate("text", x = med_score, y = max(density(plot_data$Max_TF_Score)$y) * 0.9, 
#'            label = "Median", color = "darkgreen", angle = 90, vjust = -0.5) +
#'   annotate("text", x = 2.0, y = max(density(plot_data$Max_TF_Score)$y) * 0.75, 
#'            label = "Cutoff (t=2)", color = "red", angle = 90, vjust = -0.5)
#' 
#' print(p_dist)
#' 
#' cat("Cells before low-activity filter:", ncol(thisObj_reseg), "\n")
#' thisObj_reseg <- subset(thisObj_reseg, subset = Max_TF_Score >= 2) #t=2 is a nominal p-value of ~0.045
#' cat("Cells after low-activity filter:", ncol(thisObj_reseg), "\n")
#' 
#' rm(tf_data)
#' gc()

DefaultAssay(thisObj_reseg) <- "tf_activities"
thisObj_reseg <- ScaleData(thisObj_reseg, assay = "tf_activities")
thisObj_reseg <- FindVariableFeatures(thisObj_reseg,
                                assay = "tf_activities",
                                selection.method = "vst", nfeatures = 500, verbose = TRUE)

cat("# Variable Features: ", length(Seurat::VariableFeatures(thisObj_reseg)))

thisObj_reseg<- RunPCA(thisObj_reseg,
                 assay = "tf_activities",
                 reduction.name = 'pca_TFactivity',
                 npcs = 50, 
                 verbose = TRUE)


print(ElbowPlot(thisObj_reseg, reduction = 'pca_TFactivity',  ndims = 50)) +ggtitle("TF Elbow Plot")


ScreePlot <- function(seurat_object, reduction_name, dims = 50, annotate_thresholds = TRUE) {
  # Calculate variance using the specified reduction
  variance <- (Stdev(seurat_object, reduction = reduction_name))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Find the number of PCs needed
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  cat(paste("For reduction '", reduction_name, "':\n", sep=""))
  cat("  80% variance:", n_80, "PCs\n")
  cat("  90% variance:", n_90, "PCs\n")
  cat("  95% variance:", n_95, "PCs\n\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = 1:dims,
    Variance = pc_var_explained[1:dims]
  )
  
  # Create the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "black") +
    geom_point(color = "red", size = 3) +
    labs(title = paste("Scree Plot for", reduction_name),
         x = "Principal Components",
         y = "Percentage of Variance") +
    theme_minimal()
  
  if (annotate_thresholds) {
    p <- p +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      annotate("text", x = n_90, y = max(Scree.df$Variance) * 0.9, 
               label = paste0("90% (PC ", n_90, ")"), hjust = -0.1)
  }
  
  print(p)
  
  # Return results
  return(invisible(list(n_80 = n_80, n_90 = n_90, n_95 = n_95)))
}

ScreeResults_tf <- ScreePlot(thisObj_reseg, reduction_name = 'pca_TFactivity', dims = 50)

cat("  80% variance:", ScreeResults_tf$n_80, "PCs\n")
cat("  90% variance:", ScreeResults_tf$n_90, "PCs\n")
cat("  95% variance:", ScreeResults_tf$n_95, "PCs\n\n")

dims2use <- 1:ScreeResults_tf$n_95
k2use <- 20


thisObj_reseg <- Seurat::FindNeighbors(thisObj_reseg, dims = dims2use, k.param = k2use, assay = "tf_activities", reduction = 'pca_TFactivity')

thisObj_reseg <- Seurat::FindClusters(thisObj_reseg, resolution = .5, algorithm = 3, random.seed = 42)

thisObj_reseg <- Seurat::RunUMAP(thisObj_reseg, dims = dims2use,
                           assay = "tf_activities",
                           reduction = 'pca_TFactivity',
                           reduction.name = 'tf_umap',
                           min.dist = 0.1, spread = 2)

base_umap <- DimPlot(thisObj_reseg, 
                     reduction = "tf_umap", 
                     group.by = 'tf_activities_snn_res.0.5', 
                     label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt"))

# Add the repelled labels
plot_final <- LabelClusters(plot = base_umap, 
                            id = 'tf_activities_snn_res.0.5', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final)


base_umap <- DimPlot(thisObj_reseg, 
                     reduction = "tf_umap", 
                     group.by = 'updated_celltype', 
                     label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

# Add the repelled labels
plot_final <- LabelClusters(plot = base_umap, 
                            id = 'updated_celltype', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final)















# 1. Standard RNA-based processing (No TF smoothing)
DefaultAssay(thisObj_reseg) <- "RNA"
thisObj_reseg <- NormalizeData(thisObj_reseg)
thisObj_reseg <- FindVariableFeatures(thisObj_reseg, selection.method = "vst", nfeatures = nrow(thisObj_reseg))
thisObj_reseg <- ScaleData(thisObj_reseg)
thisObj_reseg <- RunPCA(thisObj_reseg, npcs = 50)

print(ElbowPlot(thisObj_reseg, reduction = 'pca',  ndims = 50)) +ggtitle("PCA Elbow Plot")


ScreePlot <- function(seurat_object, reduction_name, dims = 50, annotate_thresholds = TRUE) {
  # Calculate variance using the specified reduction
  variance <- (Stdev(seurat_object, reduction = reduction_name))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Find the number of PCs needed
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  cat(paste("For reduction '", reduction_name, "':\n", sep=""))
  cat("  80% variance:", n_80, "PCs\n")
  cat("  90% variance:", n_90, "PCs\n")
  cat("  95% variance:", n_95, "PCs\n\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = 1:dims,
    Variance = pc_var_explained[1:dims]
  )
  
  # Create the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "black") +
    geom_point(color = "red", size = 3) +
    labs(title = paste("Scree Plot for", reduction_name),
         x = "Principal Components",
         y = "Percentage of Variance") +
    theme_minimal()
  
  if (annotate_thresholds) {
    p <- p +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      annotate("text", x = n_90, y = max(Scree.df$Variance) * 0.9, 
               label = paste0("90% (PC ", n_90, ")"), hjust = -0.1)
  }
  
  print(p)
  
  # Return results
  return(invisible(list(n_80 = n_80, n_90 = n_90, n_95 = n_95)))
}

ScreeResults_tf <- ScreePlot(thisObj_reseg, reduction_name = 'pca', dims = 50)

cat("  80% variance:", ScreeResults_tf$n_80, "PCs\n")
cat("  90% variance:", ScreeResults_tf$n_90, "PCs\n")
cat("  95% variance:", ScreeResults_tf$n_95, "PCs\n\n")

dims2use <- 1:10 #ScreeResults_tf$n_80

# 2. Relax the UMAP to allow islands to separate
# min.dist = 0.3 is the Seurat default and is much better for visualizing lineage islands
thisObj_reseg <- RunUMAP(thisObj_reseg, dims = dims2use, min.dist = 0.1, spread = 3)

# 3. THE VISUAL PROOF: Color by FastReseg's purified labels
library(ggplot2)
library(patchwork)

# Plot A: The Major Lineages
p1 <- DimPlot(thisObj_reseg, group.by = "updated_celltype", label = TRUE, repel = TRUE) +
  ggtitle("Post-FastReseg UMAP (RNA-based)") +
  theme(legend.position = "none")
p1
# Plot B: The Marker Check (Did the bleeding stop?)
# Pick two markers that used to "bleed" together in your raw data
# (e.g., Slc17a7 for Granule cells, Mbp for Oligos)
p2 <- FeaturePlot(thisObj_reseg, features = c("Slc17a7", "Mbp"), blend = TRUE)

print(p1 | p2)









































#'[Gene expression based clustering]
# We now revert to standard gene expression based clustering since we've cleaned/mitigated some segmentation issues.

thisObj_reseg <- NormalizeData(thisObj_reseg)
thisObj_reseg <- FindVariableFeatures(thisObj_reseg, selection.method = "vst", nfeatures = 500)
cat("# Variable Features: ", length(Seurat::VariableFeatures(thisObj_reseg)), "\n")

# 2. Extract the top 10 most variable genes to label them
top10_genes <- head(VariableFeatures(thisObj_reseg), 10)

# 3. Create the variance plot
var_plot <- VariableFeaturePlot(thisObj_reseg)

# 4. Add labels to the top 10 genes
labeled_plot <- LabelPoints(plot = var_plot, points = top10_genes, repel = TRUE)

# 5. Add a threshold line to see what nfeatures = 500 would actually cut off
# (Requires extracting the variance data to find the 500th gene's variance)
var_data <- HVFInfo(thisObj_reseg)
var_data <- var_data[order(var_data$variance.standardized, decreasing = TRUE), ]
threshold_variance <- var_data$variance.standardized[500]

labeled_plot <- labeled_plot + 
  geom_hline(yintercept = threshold_variance, linetype = "dashed", color = "red") +
  ggtitle("Variance vs. Mean Expression (Red Line = Top 500 Cutoff)") +
  theme(plot.title = element_text(hjust = 0.5))

print(labeled_plot)






thisObj_reseg <- ScaleData(thisObj_reseg)
thisObj_reseg <- RunPCA(thisObj_reseg, npcs = 50)

print(ElbowPlot(thisObj_reseg, reduction = 'pca',  ndims = 50)) +ggtitle("TF Elbow Plot")

ScreePlot <- function(seurat_object, reduction_name, dims = 50, annotate_thresholds = TRUE) {
  # Calculate variance using the specified reduction
  variance <- (Stdev(seurat_object, reduction = reduction_name))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Find the number of PCs needed
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  cat(paste("For reduction '", reduction_name, "':\n", sep=""))
  cat("  80% variance:", n_80, "PCs\n")
  cat("  90% variance:", n_90, "PCs\n")
  cat("  95% variance:", n_95, "PCs\n\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = 1:dims,
    Variance = pc_var_explained[1:dims]
  )
  
  # Create the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "black") +
    geom_point(color = "red", size = 3) +
    labs(title = paste("Scree Plot for", reduction_name),
         x = "Principal Components",
         y = "Percentage of Variance") +
    theme_minimal()
  
  if (annotate_thresholds) {
    p <- p +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      annotate("text", x = n_90, y = max(Scree.df$Variance) * 0.9, 
               label = paste0("90% (PC ", n_90, ")"), hjust = -0.1)
  }
  
  print(p)
  
  # Return results
  return(invisible(list(n_80 = n_80, n_90 = n_90, n_95 = n_95)))
}

ScreeResults_tf <- ScreePlot(thisObj_reseg, reduction_name = 'pca', dims = 50)

cat("  80% variance:", ScreeResults_tf$n_80, "PCs\n")
cat("  90% variance:", ScreeResults_tf$n_90, "PCs\n")
cat("  95% variance:", ScreeResults_tf$n_95, "PCs\n\n")

dims2use <- 1:ScreeResults_tf$n_80
k2use <- 20


thisObj_reseg <- Seurat::FindNeighbors(thisObj_reseg, dims = dims2use, k.param = k2use, reduction = 'pca')

# swap to default resolution (0.8) and optimized algo 4
thisObj_reseg <- Seurat::FindClusters(thisObj_reseg, resolution = .8, algorithm = 4, random.seed = 42)

thisObj_reseg <- Seurat::RunUMAP(thisObj_reseg, dims = dims2use,
                           # assay = "counts",
                           # reduction = 'pca',
                           reduction.name = 'global_umap',
                           min.dist = 0.0001, spread = 1)

base_umap <- DimPlot(thisObj_reseg, 
                     reduction = "global_umap", 
                     group.by = 'updated_celltype', 
                     label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

# Add the repelled labels
plot_final <- LabelClusters(plot = base_umap, 
                            id = 'updated_celltype', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final)






thisObj_reseg2 <- readRDS(paste0(local_WD, "thisObj_Resegmented_Filtered.rds"))

thisObj_reseg2 <- NormalizeData(thisObj_reseg2)
thisObj_reseg2 <- FindVariableFeatures(thisObj_reseg2, selection.method = "vst", nfeatures = 500)
thisObj_reseg2 <- ScaleData(thisObj_reseg2)
thisObj_reseg2 <- RunPCA(thisObj_reseg2, npcs = 50)

# You may want to run Harmony here if you have batch effects between your 16 samples!
# thisObj_reseg2 <- RunHarmony(thisObj_reseg2, group.by.vars = "mus")

thisObj_reseg2 <- RunUMAP(thisObj_reseg2, dims = 1:30,
                          min.dist = 0.1, spread = 1)
thisObj_reseg2 <- FindNeighbors(thisObj_reseg2, dims = 1:30)
thisObj_reseg2 <- FindClusters(thisObj_reseg2, resolution = 0.5)



base_umap <- DimPlot(thisObj_reseg2, 
                     reduction = "umap", 
                     group.by = 'updated_celltype',
                     label = FALSE) + 
  ggtitle("exp cluster") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

# Add the repelled labels
plot_final <- LabelClusters(plot = base_umap, 
                            id = 'updated_celltype', 
                            repel = TRUE, 
                            box = TRUE, 
                            size = 3, 
                            max.overlaps = Inf, 
                            force = 20, 
                            min.segment.length = 0)

print(plot_final)











library(pheatmap)
library("UCell")


gene_marker_lists <- list(
  # Neuronal (Cerebellar Cortex)
  Granule_Cell    = c("Cntn2", "Gls", "Rgs7"),
  UBC             = c("Calb2", "Grm1", "Plcb4"),
  Purkinje_Cell   = c("Pcp2", "Calb1", "Car8", "Itpr1", "Grid2"),
  Golgi_Cell      = c("Nos1", "Slc6a1"),
  Basket_Cell     = c("Pvalb", "Gad1"),
  Stellate_Cell   = c("Pvalb", "Penk", "Gad2"),
  Lugaro_Cell     = c("Calb2", "Gad2"),
  
  # Neuronal (Deep Cerebellar Nuclei)
  DCN_Excitatory  = c("Bdnf"),
  DCN_Inhibitory  = c("Slc32a1", "Spp1", "Gad2"),
  
  # Glial
  Astrocyte       = c("Aldh1l1", "Gfap", "Aqp4", "Slc1a2"),
  Bergmann_Glia   = c("Slc1a3"),
  OPC             = c("Pdgfra", "Olig2", "Gpr17"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag"),
  Microglia       = c("Tmem119", "P2ry12", "Cx3cr1", "C1qa"),
  
  # Vascular/Meningeal
  Endothelial     = c("Cldn5", "Pecam1"),
  Pericyte        = c("Pdgfrb", "Rgs5"),
  VSMC            = c("Acta2", "Tagln"),
  Fibroblast_VLMC = c("Col1a1", "Dcn"),
  Choroid_Plexus  = c("Ttr"),
  
  # Immune
  immune          = c("Cd3d", "Cd3e", "Cd3g", "Ptprc", "Ms4a1")
)


# Calculate Average SCT Expression
DefaultAssay(thisObj) <- "SCT"
avg_exp <- AverageExpression(thisObj, features = unique(unlist(gene_marker_lists)), 
                             group.by = "tf_activities_snn_res.0.9", assay = "SCT", slot = "data")$SCT

# Score Gene Sets per Cluster
enrichment_matrix <- sapply(gene_marker_lists, function(gs) {
  present <- intersect(gs, rownames(avg_exp))
  if(length(present) > 0) colMeans(avg_exp[present, , drop=FALSE]) else rep(0, ncol(avg_exp))
}) %>% t()

# Plot Identity Heatmap
pheatmap(enrichment_matrix, main = "Cluster Identity: SCT Marker Enrichment",
         scale = "row", cluster_cols = TRUE, border_color = "white",
         color = colorRampPalette(c("blue", "grey", "red"))(100))



















#'[Merge fragmented clusters]

thisObj_reseg$seurat_clusters_merged <- as.character(thisObj_reseg$tf_activities_snn_res.0.9) # Create a new metadata column
thisObj_reseg$seurat_clusters_merged[thisObj_reseg$seurat_clusters_merged %in% c("0", "23", "16", "30", "28", "29", "27", "25",
                                                                     "1", "33", "18", "24", "9", "26", "32", "2", "5")] <- "0"

thisObj_reseg$seurat_clusters_merged[thisObj_reseg$seurat_clusters_merged %in% c("7", "4")] <- "4"

thisObj_reseg$seurat_clusters_merged[thisObj_reseg$seurat_clusters_merged %in% c("31", "11", "12")] <- "11"
thisObj_reseg$seurat_clusters_merged[thisObj_reseg$seurat_clusters_merged %in% c("6", "8", "3")] <- "3"
thisObj_reseg$seurat_clusters_merged[thisObj_reseg$seurat_clusters_merged %in% c("14", "19", "10", "21")] <- "10"



base_umap_merged <- DimPlot(thisObj_reseg, 
                            reduction = "tf_umap", 
                            group.by = 'seurat_clusters_merged', 
                            label = FALSE) + 
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

plot_final_merged <- LabelClusters(plot = base_umap_merged, 
                                   id = 'seurat_clusters_merged', 
                                   repel = TRUE, 
                                   box = TRUE, 
                                   size = 3, 
                                   max.overlaps = Inf, 
                                   force = 20, 
                                   min.segment.length = 0)

print(plot_final_merged)



plot_cluster_diagnostics <- function(seurat_obj, 
                                     ident_name = "seurat_clusters", 
                                     sample_id_column = "orig.ident") {
  
  #' Plot Cell Proportions, Flag Missing Types, and Detect MAD Outliers
  #'
  #' @param seurat_obj A Seurat object
  #' @param ident_name Character. The metadata column containing cluster labels.
  #' @param sample_id_column Character. The metadata column containing sample IDs.
  #' @return A combined patchwork object of a stacked bar plot and a variance boxplot.
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  require(scales)
  
  meta <- seurat_obj@meta.data
  
  # 1. Calculate Exact Counts and FORCE Structural Zeros
  comp_df <- meta %>%
    dplyr::count(.data[[ident_name]], .data[[sample_id_column]], name = "n_cells") %>%
    tidyr::complete(.data[[ident_name]], .data[[sample_id_column]], fill = list(n_cells = 0))
  
  # 2. Calculate Proportions
  comp_df <- comp_df %>%
    dplyr::group_by(.data[[sample_id_column]]) %>%
    dplyr::mutate(sample_total = sum(n_cells),
                  proportion = n_cells / sample_total) %>%
    dplyr::ungroup()
  
  # A. Flag Missing Cell Types
  missing_types <- comp_df %>% dplyr::filter(n_cells == 0)
  
  if(nrow(missing_types) > 0) {
    cat("[WARNING] MISSING CELL TYPES DETECTED:\n")
    print(as.data.frame(missing_types %>% dplyr::select(-sample_total, -proportion)))
    cat("\n")
  } else {
    cat("[OK] All clusters are represented in all samples by at least 1 cell.\n\n")
  }
  
  # B. Calculate MAD and Flag Outliers
  outlier_df <- comp_df %>%
    dplyr::group_by(.data[[ident_name]]) %>%
    dplyr::mutate(
      med_prop = median(proportion),
      mad_prop = mad(proportion),
      mad_prop = ifelse(mad_prop == 0, 1e-6, mad_prop), # Prevent division by zero
      outlier_score = abs(proportion - med_prop) / mad_prop,
      is_outlier = outlier_score > 3.0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_outlier) %>%
    dplyr::arrange(desc(outlier_score))
  
  if(nrow(outlier_df) > 0) {
    cat("[ALERT] COMPOSITIONAL OUTLIERS (>3 MAD from Median):\n")
    print(as.data.frame(
      outlier_df %>% 
        dplyr::select(!!sym(ident_name), !!sym(sample_id_column), proportion, med_prop, outlier_score) %>%
        dplyr::mutate(
          proportion = paste0(round(proportion * 100, 2), "%"),
          med_prop = paste0(round(med_prop * 100, 2), "%"),
          outlier_score = paste0(round(outlier_score, 1), " MADs")
        ) %>%
        dplyr::rename(Sample_Prop = proportion, Expected_Median = med_prop, Deviation = outlier_score)
    ))
  } else {
    cat("[OK] No samples deviate by more than 3 MADs for any cell type.\n")
  }
  cat("====================================================\n")
  
  # Calculate exact number of unique clusters to define palette size
  n_clusters <- length(unique(comp_df[[ident_name]]))
  
  # Base high-contrast palette for Clusters (p1)
  base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
                "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#1B9E77", "#D95F02", 
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  
  # Dynamically expand the palette to match the exact number of clusters
  if (n_clusters > length(base_pal)) {
    expanded_pal <- colorRampPalette(base_pal)(n_clusters)
  } else {
    expanded_pal <- base_pal[1:n_clusters]
  }
  
  # Plot 1: Standard Stacked Bar
  p1 <- ggplot(comp_df, aes(x = as.factor(.data[[sample_id_column]]), 
                            y = proportion, 
                            fill = as.factor(.data[[ident_name]]))) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = expanded_pal) + # Dynamic Cluster Palette
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    ) +
    labs(title = "Overall Sample Composition",
         x = "Sample ID", y = "Proportion of Total Cells", fill = "Cluster")
  
  # Plot 2: Variance Boxplot
  p2 <- ggplot(comp_df, aes(x = as.factor(.data[[ident_name]]), 
                            y = proportion)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey40") +
    geom_jitter(aes(color = as.factor(.data[[sample_id_column]])), width = 0.2, size = 2, alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = "turbo") + # NEW: High-contrast palette for Samples
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Variance in Cluster Proportions Across Samples",
         subtitle = "Dots represent individual samples. Look for dots pulling far away from the boxplots.",
         x = "Cluster / Cell Type", y = "Proportion of Sample", color = "Sample")
  
  # Combine plots using patchwork
  final_plot <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1.5))
  
  return(final_plot)
}

plot_cluster_diagnostics(thisObj_reseg, 
                         ident_name = "seurat_clusters_merged", 
                         sample_id_column = "mus")


# Set merged cluster ids as identity
Idents(thisObj_reseg) <- "seurat_clusters_merged"








# TRASH_TRASH_TRASH_TRASH ####





#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*


cluster_mapping <- c(
  "0" = "GLN", #granule layer neurons
  # "2" = "Granule_alt1",
  # "5" = "Granule_alt2",
  "20" = "OPC", 
  # "7" = "Oligodendrocyte_1",
  "4" = "OLG",
  "15" = "DCN",
  # "12" = "Astrocytes_DCN",
  "11" = "ASC", #Astrocytes
  # "14" = "VSMC/VLMC",
  # "10" = "VASC", #Endothelial
  # "21" = "VASC",
  # "6" = "BG", #bergmann glia
  "13" = "PKJ", #purkinje neurons
  "3" = "MLI", #molecular layer interneurons
  "22" = "CHP", #Choroid Plexus
  "17" = "MGL" #Microglia
  # "8" = "MLI/PKJ/BG_MIX"
)

custom_col <- c(
  # --- Neurons ---
  "GLN"        = "navy",
  # "Granule_alt1"   = "royalblue",  
  # "Granule_alt2"   = "steelblue",  
  "PKJ"       = "red",
  "DCN"            = "red2",
  "MLI"            = "limegreen", 
  
  # --- Glia ---
  "ASC"     = "orange",     
  # "Astrocytes_DCN"  = "darkorange",  
  "BG"  = "gold",
  "OLG" = "darkorchid4",
  # "Oligodendrocyte_2" = "purple",
  "OPC"            = "mediumorchid",
  
  # --- Vasculature & Stroma ---
  "VASC"    = "cyan",
  # "Pericyte"       = "cyan3",
  # "VSMC/VLMC"      = "cadetblue",
  # "Choroid Plexus" = "pink",
  
  # --- Immune ---
  "MGL"      = "black"
  
  # --- Mixed / Artifacts ---
  # "MLI/PKJ/BG_MIX" = "grey50"
)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Cell typed UMAP]

DefaultAssay(thisObj) <- 'tf_activities'
umap_p <- DimPlot(thisObj,
                  reduction = "tf_umap",
                  group.by = "seurat_clusters_CellTyped",
                  raster = TRUE,
                  label = FALSE,
                  pt.size = 0.5,
                  repel = TRUE) +
  scale_color_manual(values = custom_col) +
  ggtitle("CellTyped TF UMAP") +
  theme(plot.title = element_text(hjust = 0.5))
umap_p

cellTyped_UMAP <- LabelClusters(plot = umap_p, 
                                   id = 'seurat_clusters_CellTyped', 
                                   repel = TRUE, 
                                   box = TRUE, 
                                   size = 3, 
                                   max.overlaps = Inf, 
                                   force = 20, 
                                   min.segment.length = 0)

print(cellTyped_UMAP)


base_umap_merged <- DimPlot(thisObj, 
                            reduction = "tf_umap", 
                            group.by = 'seurat_clusters_CellTyped', 
                            label = FALSE) + 
  scale_color_manual(values = custom_col) +
  ggtitle("TF-Based Regulatory Clusters") +
  theme(legend.position = "none",
        plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")) +
  coord_cartesian(clip = "off")

plot_final_merged <- LabelClusters(plot = base_umap_merged, 
                                   id = 'seurat_clusters_CellTyped', 
                                   repel = TRUE, 
                                   box = TRUE, 
                                   size = 3, 
                                   max.overlaps = Inf, 
                                   force = 20, 
                                   min.segment.length = 0)

print(plot_final_merged)


ggsave(filename = paste0(local_WD, "TF_UMAP_v2.png"), plot = cellTyped_UMAP, dpi = 600, width = 10, height = 8)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Cell type enrichment heatmap w/ new cluster names]
gene_marker_lists <- list(
  # Neuronal (Cerebellar Cortex)
  Granule_Cell    = c("Cntn2", "Gls", "Rgs7"),
  UBC             = c("Calb2", "Grm1", "Plcb4"),
  Purkinje_Cell   = c("Pcp2", "Calb1", "Car8", "Itpr1", "Grid2"),
  Golgi_Cell      = c("Nos1", "Slc6a1"),
  Basket_Cell     = c("Pvalb", "Gad1"),
  Stellate_Cell   = c("Pvalb", "Penk", "Gad2"),
  Lugaro_Cell     = c("Calb2", "Gad2"),
  
  # Neuronal (Deep Cerebellar Nuclei)
  DCN_Excitatory  = c("Bdnf"),
  DCN_Inhibitory  = c("Slc32a1", "Spp1", "Gad2"),
  
  # Glial
  Astrocyte       = c("Aldh1l1", "Gfap", "Aqp4", "Slc1a2"),
  Bergmann_Glia   = c("Slc1a3"),
  OPC             = c("Pdgfra", "Olig2", "Gpr17"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag"),
  Microglia       = c("Tmem119", "P2ry12", "Cx3cr1", "C1qa"),
  
  # Vascular/Meningeal
  Endothelial     = c("Cldn5", "Pecam1"),
  Pericyte        = c("Pdgfrb", "Rgs5"),
  VSMC            = c("Acta2", "Tagln"),
  Fibroblast_VLMC = c("Col1a1", "Dcn"),
  Choroid_Plexus  = c("Ttr"),
  
  # Immune
  immune          = c("Cd3d", "Cd3e", "Cd3g", "Ptprc", "Ms4a1")
)
# re-print cell type heatmap w/ cluster labels
avg_exp <- AverageExpression(thisObj, features = unique(unlist(gene_marker_lists)), 
                             group.by = "seurat_clusters_CellTyped", assay = "SCT", slot = "data")$SCT

# Score Gene Sets per Cluster
enrichment_matrix <- sapply(gene_marker_lists, function(gs) {
  present <- intersect(gs, rownames(avg_exp))
  if(length(present) > 0) colMeans(avg_exp[present, , drop=FALSE]) else rep(0, ncol(avg_exp))
}) %>% t()

pheatmap(enrichment_matrix, main = "Cluster Identity: SCT Marker Enrichment",
         scale = "row", cluster_cols = TRUE, border_color = "white",
         color = colorRampPalette(c("blue", "grey", "red"))(100))


# Load Cell Typed Seurat Obj ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Load Cell Typed Seurat Obj*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
local_WD <- "C:/Users/m322332/Documents/"
thisObj <- readRDS(file.path(local_WD, "COSMONAUT_CB_Seurat_GlobalCellTyped.RDS"))

custom_col <- c(
  # --- Neurons ---
  "Granule"        = "navy",
  "Granule_alt1"   = "royalblue",  
  "Granule_alt2"   = "steelblue",  
  "Purkinje"       = "red",
  "DCN"            = "red2",
  "MLI"            = "limegreen", 
  
  # --- Glia ---
  "Astrocytes"     = "orange",     
  "Astrocytes_DCN"  = "darkorange",  
  "Bergmann Glia"  = "gold",
  "Oligodendrocyte_1" = "darkorchid4",
  "Oligodendrocyte_2" = "purple",
  "OPC"            = "mediumorchid",
  
  # --- Vasculature & Stroma ---
  "Endothelial"    = "cyan",
  "Pericyte"       = "cyan3",
  "VSMC/VLMC"      = "cadetblue",
  "Choroid Plexus" = "pink",
  
  # --- Immune ---
  "Microglia"      = "black",
  
  # --- Mixed / Artifacts ---
  "MLI/PKJ/BG_MIX" = "grey50"
)
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Plot Spatial Cell Types ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Plot Spatial Cell Types*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

thisObj@meta.data$x_slide_micron <- thisObj@meta.data$x_slide_mm * 1000
thisObj@meta.data$y_slide_micron <- thisObj@meta.data$y_slide_mm * 1000

# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"


plotSpatialSample <- function(seurat_obj,
                              sample_id,
                              sample_col,
                              coord_x_micron = "x_slide_micron",
                              coord_y_micron = "y_slide_micron",
                              color_by = "ident",
                              custom_palette = NULL,
                              pt_size = 1,
                              alpha = 0.6,
                              plot_title = NULL,
                              # "ggplot" or "plotly" (interactive plot)
                              return_type = "ggplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (return_type == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' needed for return_type='plotly'. Please install it.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Sample column '", sample_col, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_x_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("X coordinate column '", coord_x_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_y_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Y coordinate column '", coord_y_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!(color_by == "ident" || color_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Coloring column/variable '", color_by, "' not found in Seurat object metadata or is not 'ident'."), call. = FALSE)
  }
  if (!sample_id %in% unique(seurat_obj@meta.data[[sample_col]])) {
    stop(paste("Sample ID '", sample_id, "' not found in column '", sample_col, "'."), call. = FALSE)
  }
  
  print(paste("Subsetting data for sample:", sample_id))
  subset_data <- subset(seurat_obj, subset = !!sym(sample_col) == sample_id)
  
  if (nrow(subset_data@meta.data) == 0) {
    stop(paste("No cells found after subsetting for sample ID '", sample_id, "'."), call. = FALSE)
  }
  
  plot_df <- subset_data@meta.data
  if (color_by == "ident") {
    plot_df$color_values <- as.character(Idents(subset_data))
    color_lab <- "Cluster ID"
  } else {
    plot_df$color_values <- as.character(plot_df[[color_by]])
    color_lab <- color_by
  }
  plot_df$color_values <- factor(plot_df$color_values) # Ensure factor for levels
  
  x_range <- range(plot_df[[coord_x_micron]], na.rm = TRUE)
  y_range <- range(plot_df[[coord_y_micron]], na.rm = TRUE)
  buffer <- 10
  xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
  ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Sample:", sample_id, "- Colored by:", color_lab)
  }
  
  # Determine Color Palette
  color_levels <- levels(plot_df$color_values)
  num_colors <- length(color_levels)
  
  if (is.null(custom_palette)) {
    # Generate default palette if none provided
    if (num_colors <= 12) {
      color_palette_values <- hcl.colors(num_colors, palette = "Set 3")
    } else {
      color_palette_values <- scales::hue_pal()(num_colors)
    }
    # Name palette for consistent mapping
    names(color_palette_values) <- color_levels
  } else {
    # Use the custom palette
    color_palette_values <- custom_palette
  }
  
  
  print("Generating ggplot...")
  p <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]], color = color_values)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_palette_values) +
    labs(title = plot_title,
         x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
    theme_minimal() +
    coord_fixed() +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
  
  if (return_type == "plotly") {
    print("Generating plotly object...")
    plot_df$tooltip_text <- paste("Cell:", rownames(plot_df), "<br>", color_lab, ":", plot_df$color_values)
    p_plotly <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]],
                                    color = color_values, text=tooltip_text)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_manual(values = color_palette_values) + # Use the chosen palette
      labs(title = plot_title, x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
      theme_minimal() + coord_fixed() + xlim(xlim_vals) + ylim(ylim_vals) +
      theme( panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.box = "horizontal") +
      guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
    
    suppressWarnings({
      plotly_obj <- plotly::ggplotly(p_plotly, tooltip = "text")
    })
    return(plotly_obj)
    
  } else {
    print("Returning ggplot object.")
    return(p)
  }
}

levels(thisObj$seurat_clusters_merged)
YM04_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                                         sample_id = "YF02",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "seurat_clusters_CellTyped",
                                         custom_palette = custom_col,
                                         pt_size = 0.75,
                                         alpha = 0.6,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(YM04_spatial_ggplot)


highlight_custom_col <- c(
  # --- Neurons ---
  "Granule"        = "grey",
  "Granule_alt1"   = "grey",  
  "Granule_alt2"   = "grey",  
  "Purkinje"       = "grey",
  "DCN"            = "red",
  "MLI"            = "grey", 
  
  # --- Glia ---
  "Astrocytes"     = "grey",     
  "Astrocytes_DCN"  = "grey",  
  "Bergmann Glia"  = "grey",
  "Oligodendrocyte_1" = "grey",
  "Oligodendrocyte_2" = "grey",
  "OPC"            = "grey",
  
  # --- Vasculature & Stroma ---
  "Endothelial"    = "grey",
  "Pericyte"       = "grey",
  "VSMC/VLMC"      = "grey",
  "Choroid Plexus" = "grey",
  
  # --- Immune ---
  "Microglia"      = "grey",
  
  # --- Mixed / Artifacts ---
  "MLI/PKJ/BG_MIX" = "grey"
)

highlight_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                                         sample_id = "YM04",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "seurat_clusters_CellTyped",
                                         custom_palette = highlight_custom_col,
                                         pt_size = 0.75,
                                         alpha = 0.6,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")+ 
  theme_void() + # Completely removes axes, backgrounds, and grid lines
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  # Insert the cluster name inside the plot
  annotate("text", x = -Inf, y = Inf, label = "Microglia", 
           hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold", color = "black")

print(highlight_spatial_ggplot)

ggsave(
  filename = paste0(local_WD, "239_microglia.png"), 
  plot = highlight_spatial_ggplot, 
  width = 6,  
  height = 6,  
  dpi = 600  
)


library(ggplot2)
library(patchwork)

# Define list of clusters
custom_col <- c(
  # --- Neurons ---
  "Granule",
  "Granule_alt1",   
  "Granule_alt2",   
  "Purkinje",
  "DCN",
  "MLI", 
  
  # --- Glia ---
  "Astrocytes",     
  "Astrocyte_DCN", 
  "Bergmann Glia",
  "Oligodendrocytes",
  "OPC",
  
  # --- Vasculature & Stroma ---
  "Endothelial",
  "Pericyte",
  "VSMC/VLMC",
  "Choroid Plexus",
  
  # --- Immune ---
  "Microglia",
  
  # --- Mixed / Artifacts ---
  "MLI/PKJ/BG_MIX" 
)


# Initialize an empty list to store the generated plots
plot_list <- list()

cat("Generating spatial highlight plots...\n")

# Loop through each cluster
for (clust in all_clusters) {
  
  # Dynamically build the palette: Set all to a light grey, then the target to red
  current_pal <- setNames(rep("grey85", length(all_clusters)), all_clusters)
  current_pal[clust] <- "red"
  
  # Generate plot
  p <- plotSpatialSample(
    seurat_obj = thisObj,
    sample_id = "305",
    sample_col = "mus",
    coord_x_micron = "x_slide_micron",
    coord_y_micron = "y_slide_micron",
    color_by = "seurat_clusters_CellTyped",
    custom_palette = current_pal,
    pt_size = 0.5, # Slightly smaller points prevent the grey from burying the red
    alpha = 0.8,
    plot_title = NULL,
    return_type = "ggplot"
  )
  
  # Apply formatting
  p <- p + 
    theme_void() + # Completely removes axes, backgrounds, and grid lines
    theme(
      legend.position = "none",
      plot.title = element_blank()
    ) +
    # Insert the cluster name inside the plot so you know what you are looking at
    annotate("text", x = -Inf, y = Inf, label = clust, 
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold", color = "black")
  
  safe_clust_name <- gsub("/", "_", clust)
  
  # ggsave(
  #   filename = paste0(local_WD, "305", safe_clust_name, ".png"), 
  #   plot = p, 
  #   width = 4,  
  #   height = 3,  
  #   dpi = 600  
  # )

  # Save plot to the list
  plot_list[[clust]] <- p
}

plot_list$Granule
plot_list$


#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Celltype Populations Size/Composition Analysis/Plots ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Celltype Populations Size/Composition Analysis/Plots*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

#'[Figure: Cluster population proportions]
# Extract metadata
meta <- thisObj@meta.data

# Calculate raw cell counts per sample, age, sex, and cluster
raw_counts <- meta %>%
  group_by(mus, age, sex, seurat_clusters_CellTyped) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate the total number of cells in each sample
sample_totals <- meta %>%
  group_by(mus) %>%
  summarise(TotalCells = n(), .groups = "drop")

# Create a complete grid to catch "Structural Zeros"
# This ensures every sample has an entry for every cluster, even if the count is 0
all_samples <- meta %>% select(mus, age, sex) %>% distinct()
all_clusters <- unique(meta$seurat_clusters_CellTyped)

complete_grid <- expand_grid(
  mus = all_samples$mus,
  seurat_clusters_CellTyped = all_clusters
) %>% left_join(all_samples, by = "mus")

# Merge data and calculate population proportions
prop_data <- complete_grid %>%
  left_join(raw_counts, by = c("mus", "age", "sex", "seurat_clusters_CellTyped")) %>%
  replace_na(list(Count = 0)) %>%
  left_join(sample_totals, by = "mus") %>%
  mutate(Percent = (Count / TotalCells) * 100)

# Factor and re-order by age
prop_data$age <- factor(prop_data$age, levels = c("YNG", "OLD"))


stat_results <- prop_data %>%
  group_by(seurat_clusters_CellTyped) %>%
  summarise(
    # Run the t-test (suppressing warnings if a cluster has 0 variance in a group)
    p_val = tryCatch(t.test(Percent ~ age)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_val, method = "fdr"),
    # Create significance stars for the plot
    sig_label = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print("--- Differential Abundance T-Test Results ---")
print(stat_results %>% arrange(p_adj))


# Merge the significance labels back into the plotting data for easy annotation
prop_data <- prop_data %>% left_join(stat_results, by = "seurat_clusters_CellTyped")

# Calculate the maximum Y value per cluster to position the significance stars cleanly
max_y_data <- prop_data %>%
  group_by(seurat_clusters_CellTyped, sig_label) %>%
  summarise(Max_Percent = max(Percent) * 1.05, .groups = "drop")

# Generate the plot
p_pop <- ggplot(prop_data, aes(x = seurat_clusters_CellTyped, y = Percent)) +
  # Split Violins mapping to Age
  geom_violin(aes(fill = age), position = position_dodge(width = 0.8), 
              alpha = 0.6, scale = "width", color = "gray30") +
  
  # Sample Data Points mapping Shape to Sex and Grouping to Age
  geom_point(aes(shape = sex, group = age), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
             size = 2.5, color = "black", alpha = 0.5) +
  
  # Significance Stars mapping
  geom_text(data = max_y_data, aes(x = seurat_clusters_CellTyped, y = Max_Percent, label = sig_label), 
            size = 5, vjust = 0) +
  scale_fill_manual(values = c("YNG" = "skyblue", "OLD" = "salmon")) +
  scale_shape_manual(values = c("M" = 16, "F" = 17)) + # 16 = circle, 17 = triangle
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Relative Population Frequencies by Age and Sex",
    subtitle = "Normalized to sample-level capture. Significance: * p<0.05, ** p<0.01, *** p<0.001 (FDR)",
    x = "Cell Identity", 
    y = "Proportion of Total Sample Cells (%)",
    fill = "Age", 
    shape = "Sex"
  )

print(p_pop)

ggsave(
  filename = paste0("CB_Population_Changes.png"), 
  plot = p_pop, 
  width = 10,  
  height = 5,  
  dpi = 600  
)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Top n DE genes per cluster dotplot]
Idents(thisObj) <- "seurat_clusters_CellTyped"
DefaultAssay(thisObj) <- "SCT"

cluster_markers_wilcox <- FindAllMarkers(
  thisObj,
  only.pos = TRUE, 
  min.pct = 0.10, 
  logfc.threshold = 0.25, 
  test.use = "wilcox",
  slot = "data"
)

library(DESeq2)
cluster_markers_deseq <- FindAllMarkers(
  thisObj,
  only.pos = TRUE, 
  min.pct = 0.10, 
  # logfc.threshold = 0.25, 
  test.use = "DESeq2"
  # slot = "data"
)

top10_markers <- cluster_markers_wilcox %>%
  group_by(cluster) %>%
  mutate(
    # Add a tiny constant to prevent -log10(0) returning Inf
    safe_p_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
    # Calculate the pi score; a stand-in for a t value since the fxn doesn't output a t value.
    pi_score = -log10(safe_p_adj) * avg_log2FC
  ) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

head(top10_markers)

cluster_peri <- top10_markers[top10_markers$cluster == "Pericyte", ]
head(cluster_peri)

top3_markers <- cluster_markers_wilcox %>%
  group_by(cluster) %>%
  mutate(
    # Add a tiny constant to prevent -log10(0) returning Inf
    safe_p_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
    # Calculate the pi score; a stand-in for a t value since the fxn doesn't output a t value.
    pi_score = -log10(safe_p_adj) * avg_log2FC
  ) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  ungroup()

# Extract just the unique gene names for the DotPlot
top_genes_list <- unique(top3_markers$gene)
head(top_genes_list)


p_dotplot <- DotPlot(
  thisObj, 
  features = top_genes_list, 
  cols = c("lightgrey", "red"),
  dot.scale = 6,
  scale = FALSE 
) + 
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 10, face = "italic"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90")
  ) +
  labs(
    title = "Top 3 Cluster Markers",
    subtitle = "Wilcoxon Test (Normalized Expression)",
    x = "Marker Gene", 
    y = "Cluster Identity"
  )

print(p_dotplot)

ggsave(
  filename = paste0("CB_ClusterTop3DotPlot.png"), 
  plot = p_dotplot, 
  width = 8,  
  height = 10,  
  dpi = 600  
)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Violin plot of cluster population proportions]
# Ensure active identity is correct and pull metadata
Idents(thisObj) <- "seurat_clusters_CellTyped"
meta_data <- thisObj@meta.data
meta_data$ClusterID <- Idents(thisObj)

sample_col <- "mus"
age_col <- "age"
sex_col <- "sex" 
cluster_col <- "ClusterID"

# Calculate sample totals
sample_totals <- meta_data %>%
  group_by(!!sym(sample_col)) %>%
  summarise(TotalCellsInSample = n(), .groups = 'drop')

# Create a complete grid to catch structural zeros
all_samples <- meta_data %>% select(all_of(c(sample_col, age_col, sex_col))) %>% distinct()
all_clusters <- unique(meta_data$ClusterID)

complete_grid <- expand_grid(
  !!sym(sample_col) := all_samples[[sample_col]],
  ClusterID = all_clusters
) %>% left_join(all_samples, by = sample_col)

# Calculate raw counts
raw_counts <- meta_data %>%
  group_by(!!sym(sample_col), ClusterID, !!sym(age_col), !!sym(sex_col)) %>%
  summarise(CellCount = n(), .groups = 'drop')

# Define your custom spatial color palette
custom_col <- c(
  # --- Neurons ---
  "Granule"        = "navy",
  "Granule_alt1"   = "royalblue",  
  "Granule_alt2"   = "steelblue",  
  "Purkinje"       = "red",
  "DCN"            = "red2",
  "MLI"            = "limegreen", 
  
  # --- Glia ---
  "Astrocytes"     = "orange",     
  "Astrocytes_DCN"  = "darkorange",  
  "Bergmann Glia"  = "gold",
  "Oligodendrocyte_1" = "darkorchid4",
  "Oligodendrocyte_2" = "purple",
  "OPC"            = "mediumorchid",
  
  # --- Vasculature & Stroma ---
  "Endothelial"    = "cyan",
  "Pericyte"       = "cyan3",
  "VSMC/VLMC"      = "cadetblue",
  "Choroid Plexus" = "pink",
  
  # --- Immune ---
  "Microglia"      = "black",
  
  # --- Mixed / Artifacts ---
  "MLI/PKJ/BG_MIX" = "grey50"
)

# Merge and calculate final proportions
plot_data <- complete_grid %>%
  left_join(raw_counts, by = c(sample_col, "ClusterID", age_col, sex_col)) %>%
  replace_na(list(CellCount = 0)) %>%
  left_join(sample_totals, by = sample_col) %>%
  mutate(Proportion = CellCount / TotalCellsInSample) %>%
  mutate(
    # Lock cluster order to match your color palette list logically
    ClusterID = factor(ClusterID, levels = names(custom_col)),
    Age = factor(!!sym(age_col), levels = c("YNG", "OLD")), 
    Sex = factor(!!sym(sex_col)),
    SampleID = factor(!!sym(sample_col))
  )

# Violin Plots
plot_all <- ggplot(plot_data, aes(x = ClusterID, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = Age), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 0.5, alpha = 0.6) +
  scale_fill_manual(name = "Age", values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  scale_color_manual(name = "Age", values = c("YNG" = "darkblue", "OLD" = "darkred")) +
  labs(
    title = "Normalized Cell Type Proportions by Age Group",
    x = "Cluster ID",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

print(plot_all)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Faceted by cluster Violin plot of cluster population proportions]
# Option 2: Facet by cluster (cleaner for large numbers of clusters)
plot_faceted <- ggplot(plot_data, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) + 
  facet_wrap(~ ClusterID, scales = "free_y", ncol = 6) + 
  labs(
    title = "Normalized Cell Type Proportions by Age Group (Faceted by Cluster)",
    x = "Age Group",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot_faceted)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Cluster specific violin of population proportions]
plot_data_mg <- plot_data %>% filter(ClusterID == "Microglia")

plot_mg <- ggplot(plot_data_mg, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width", color="gray30") +
  geom_jitter(aes(shape = Sex), width = 0.1, size = 3.0, alpha = 0.8) +
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon"), guide = "none") + 
  scale_shape_manual(values = c("F" = 16, "M" = 17)) + 
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    title = "Microglia Population Density",
    y = "Proportion of Cells per Sample",
    shape = "Sex" 
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    axis.title.y = element_text(size = 12, face = "bold"), 
    axis.text.x = element_text(face = "bold", size = 12), 
    axis.title.x = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.position = "bottom" 
  )

print(plot_mg)
ggsave(filename = paste0(local_WD, "plot_CB_Microglia_proportions.png"), plot = plot_mg, dpi = 600, width = 4, height = 5)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Stacked bar plot of cluster population proportions out of 100%]

stacked_bar_plot <- ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = ClusterID)) +
  geom_col(position = "stack") +
  facet_grid(~Age, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_col) + 
  labs(
    title = "Cell Type Composition per Sample by Age Group",
    x = "Sample ID",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face="bold"),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face="bold", size=12),
    panel.spacing.x = unit(1, "lines")
  )

print(stacked_bar_plot)

ggsave(filename = paste0(local_WD, "plot_CB_PopulationsStackedBar.png"), plot = stacked_bar_plot, dpi = 600, width = 9, height = 3)

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Senescent Signature Analysis/Plots ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Senescent Signature Analysis/Plots*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
library(UCell)
library(forcats)
library(patchwork)

SenGeneSets <- list(
  SenBrain_Mouse_2025 = c("Akt1","Apoe","App","Axl","Bcl2","Bcl2l1","Ccl2","Ccl3","Ccl4","Ccl5","Ccne1",
                          "Cd274","Cd36","Cd74","Cd9","Cdc27","Cdc6","Cdk4","Cdk8","Cdk9","Cdkn1a/p21","Cdkn1b/p27kip1",
                          "Cdkn1c/p57kip2","Cdkn2a/p16","Cdkn2b/p15ink4b","Cdkn2c/p18ink4c","Cdkn2d/p19","Cenpb","Csf1","Csf1r","Cxcl14","Cxcl2",
                          "Cyba","Dpp4","Elmo1","Fxyd5","Gadd45a","Gas2","Gdf15","Glb1","Gpnmb","H2-Aa",
                          "H2-Ab1","Hmgb1","Ifitm3","Igfbp7","Il10","Il18","Il1a","Il1b","Il23a","Il23r",
                          "Il33","Il6","Irf5","Irf7","Lamp2","Lgals3","Lgals3bp","Lmna","Lpl","Malat1","Mapk1",
                          "Mapk10","Mapk14","Mapk8","Meox1","Mmp12","Mmp14","Mmp2","Mmp3","Mmp9","Nfkb1",
                          "Nol3","Notch3","Pak1","Pak2","Pak3","Pdcd1","Pik3ca","Pik3cb","Pik3r1","Plaur",
                          "Pten","Rbl2","Serpina3n","Serpinb1a","Serpine1","Serpine2","Smad2","Smad3","Sod1","Sparc","Spp1",
                          "Tgfb1","Tgfb2","Tgfbr1","Sting1","Tnf","Tp53bp1","Trp53","Vim"),
  
  SenMayo_Mouse = c("Acvr1b","Ang","Angpt1","Angptl4","Areg","Axl","Bex3","Bmp2","Bmp6","C3","Ccl1","Ccl2","Ccl20","Ccl24",
                    "Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Cd55","Cd9","Csf1","Csf2","Csf2rb","Ctnnb1","Ctsb","Cxcl1",
                    "Cxcl10","Cxcl12","Cxcl16","Cxcl2","Cxcl3","Cxcr2","Dkk1","Edn1","Egf","Egfr","Ereg","Esm1","Ets2","Fas",
                    "Fgf1","Fgf2","Fgf7","Gdf15","Gem","Gmfg","Hgf","Hmgb1","Icam1","Icam5","Igf1","Igfbp1","Igfbp2","Igfbp3",
                    "Igfbp4","Igfbp5","Igfbp6","Igfbp7","Il10","Il13","Il15","Il18","Il1a","Il1b","Il2","Il6","Il6st","Il7",
                    "Inha","Iqgap2","Itga2","Itpka","Jun","Kitl","Lcp1","Mif","Mmp13","Mmp10","Mmp12","Mmp14","Mmp2",
                    "Mmp3","Mmp9","Nap1l4","Nrg1","Pappa","Pecam1","Pgf","Pigf","Plat","Plau","Plaur","Ptbp1","Ptger2","Ptges",
                    "Rps6ka5","Scamp4","Selplg","Sema3f","Serpinb3c","Serpine1","Serpine2","Spp1","Spx","Timp2","Tnf",
                    "Tnfrsf11b","Tnfrsf1a","Tnfrsf1b","Tubgcp2","Vegfa","Vegfc","Vgf","Wnt16","Wnt2"),
  
  Reactome_Cellular_Senescence = c('Mapk7', 'Ube2c', 'Ccne1', 'Map2k7', 'Cdkn1b/p27kip1', 'Rps6ka1', 'Cdk4', 'Ubc', 'Ehmt2', 'Anapc1', 
                                   'Mapkapk2', 'Suz12', 'H2bc4', 'Kdm6b', 'Map2k3', 'Ubb', 'Asf1a', 'Ube2d1', 'Anapc16', 'Mdm2', 
                                   'Cabin1', 'Fzr1', 'Mapk9', 'Rad50', 'Rps27a', 'Map2k6', 'Cdc27', 'Mink1', 'Fos', 'Ube2e1', 
                                   'Mapk8', 'Rb1', 'Hira', 'Ets2', 'Cdkn1a/p21', 'Rps6ka2', 'Cdc23', 'Lmnb1', 'Kat5', 'Anapc11', 
                                   'Cdk2', 'Terf1', 'Map4k4', 'Anapc2', 'Tnik', 'Ccna2', 'Ccna1', 'Ccne2', 'Nbn', 'Txn1', 
                                   'Anapc4', 'Mapkapk5', 'Anapc7', 'Anapc5', 'Ep400', 'Pot1a', 'Ezh2', 'Eed', 'Anapc15', 
                                   'Rps6ka3', 'Rbbp7', 'Terf2', 'Mre11a', 'Ets1', 'Mapkapk3', 'Map2k4', 'Terf2ip', 'Atm', 
                                   'H1f2', 'Ehmt1', 'Anapc10', 'Cdkn1c/p57kip2', 'Acd', 'Cdc16', 'Ubn1', 'Cdk6', 'Erf', 'H2az2', 
                                   'Id1', 'Mapk10', 'Hmga1', 'H2bc6', 'H1f1', 'H2ax', 'H2bc1', 'H1f4', 'Jun', 'Mapk11', 
                                   'Mapk14', 'Mdm4', 'Cebpb', 'Hmga2', 'Rbbp4', 'H2bc8', 'H1f5', 'Trp53', 'H2aj', 'H4c1', 
                                   'H4c9', 'H4c3', 'Ube2s', 'H4c8', 'H4c4', 'H2ac4', 'H2ac10', 'H2bc12', 'H2ac15', 'Mapk3', 
                                   'Mapk1', 'H2ac19', 'H2ac18', 'H4c12', 'Cdc26', 'H2ab1', 'H4c11', 'Uba52rt', 'H2bc21', 
                                   'H2ac20', 'H4c2', 'H2bc7', 'H2ac6', 'H2ac8', 'H4c6', 'H2bc11', 'H2ac11', 'H2ac12', 
                                   'H2bc24', 'H4c18', 'H4c17', 'H2bc23', 'H2bc22', 'H2ac22', 'Map3k5', 'H2ac7', 'H2ac13', 
                                   'Cdkn2b/p15ink4b', 'H2bc3', 'Hmga1b', 'H3f4', 'H2bc26', 'H2ab2', 'H2ab3', 'Uba52', 'H4c14', 
                                   'H2ac23', 'H2bc13', 'H2ac24', 'H2bc15', 'H4c16', 'H1f0', 'H2bc14', 'H2bc9'),
  
  Hallmark_P53_Pathway = c('Cdk5r1', 'F2r', 'Rack1', 'Ptpn14', 'Ccng1', 'Dgka', 'Vdr', 'Zbtb16', 'Xpc', 'Ercc5', 'Cdkn1a/p21', 
                           'Procr', 'Cd82', 'Tm4sf1', 'Fdxr', 'Cdkn2b/p15ink4b', 'Cdkn2a/p16', 'Ier3', 'Mxd4', 'Tnni1', 'Tspyl2', 'Gpx2', 
                           'Irak1', 'Hspa4l', 'Ralgds', 'Gadd45a', 'Zfp36l1', 'Ei24', 'Btg2', 'Rgs16', 'Irag2', 'Tsc22d1', 
                           'Fgf13', 'Ddit3', 'Atf3', 'Serpinb5', 'Plk3', 'Cd81', 'Bak1', 'Upp1', 'Alox8', 'Kif13b', 'Jag2', 
                           'Plk2', 'Tnfsf9', 'Zmat3', 'Ninj1', 'Traf4', 'Ccnk', 'Tpd52l1', 'Apaf1', 'Sphk1', 'Prmt2', 'Hint1', 
                           'Rad9a', 'S100a4', 'Trp63', 'Prkab1', 'Vamp8', 'Ier5', 'St14', 'Acvr1b', 'S100a10', 'Ndrg1', 
                           'Klf4', 'Hdac3', 'Klk8', 'Trib3', 'Def6', 'Slc7a11', 'Sdc1', 'Tob1', 'Phlda3', 'Abcc5', 'Pmm1', 
                           'Fbxw7', 'Socs1', 'Ddb2', 'Ppm1d', 'Sec61a1', 'Blcap', 'Rpl36', 'Ctsf', 'Wwp1', 'Pidd1', 'Txnip', 
                           'Foxo3', 'Polh', 'Wrap73', 'Sfn', 'Nupr1', 'Iscu', 'Rrp8', 'Elp1', 'Coq8a', 'Retsat', 'Tm7sf3', 
                           'Vwa5a', 'Dcxr', 'Rps27l', 'Aen', 'Rchy1', 'Steap3', 'Ankra2', 'Abhd4', 'Triap1', 'Cgrrf1', 
                           'Tprkb', 'Fam162a', 'Cdkn2aip', 'Nol8', 'Dram1', 'Pitpnc1', 'Osgin1', 'Rap2b', 'Ddit4', 'Rnf19b', 
                           'Dnttip2', 'Tax1bp3', 'Trafd1', 'Ip6k2', 'Cyfip2', 'Ppp1r15a', 'Slc19a2', 'Perp', 'Rrad', 'H1f2', 
                           'Baiap2', 'Pom121', 'Ifi30', 'Eps8l2', 'Clca2', 'Slc35d1', 'Ccp110', 'Mapkapk3', 'Gls2', 'Zfp365', 
                           'Rad51c', 'Plxnb2', 'Sesn1', 'Tchh', 'Sertad3', 'Rab40c', 'Hexim1', 'Csrnp2', 'Rhbdf2', 'Nudt15', 
                           'Abat', 'H2ac25', 'H2aj', 'Ada', 'Ak1', 'App', 'Bmp2', 'Btg1', 'Ccnd2', 'Ccnd3', 'Ctsd', 'Mknk2', 
                           'Epha2', 'Stom', 'Ephx1', 'Fas', 'Fos', 'Fuca1', 'Gm2a', 'Hbegf', 'Hmox1', 'Hras', 'Il1a', 
                           'Casp1', 'Inhbb', 'Itgb4', 'Jun', 'Krt17', 'Ldhb', 'Lif', 'Mxd1', 'Mdm2', 'Slc3a2', 'Nhlh2', 
                           'Notch1', 'Pcna', 'Pdgfa', 'Ptpre', 'Pvt1', 'Rb1', 'Rpl18', 'Rps12', 'Rxra', 'Sat1', 'Sp1', 
                           'Tap1', 'Tcn2', 'Tgfa', 'Tgfb1', 'Trp53', 'Cebpa', 'Cdh13', 'Bax'),
  
  # Reactome_SASP = c('Mapk7', 'Ube2c', 'Cdkn1b/p27kip1', 'Rps6ka1', 'Cdk4', 'Ubc', 'Ehmt2', 'Anapc1', 'Ubb', 'Ube2d1', 'Anapc16', 
  #                   'Fzr1', 'Rps27a', 'Cdc27', 'Ube2e1', 'Cdkn1a/p21', 'Rps6ka2', 'Cdc23', 'Anapc11', 'Cdk2', 'Anapc2', 
  #                   'Ccna2', 'Ccna1', 'Anapc4', 'Anapc7', 'Anapc5', 'Anapc15', 'Rps6ka3', 'Ehmt1', 'Anapc10', 'Cdkn1c/p57kip2', 
  #                   'Cdc16', 'Cdk6', 'Cebpb', 'Ube2s', 'Mapk3', 'Mapk1', 'Cdc26', 'Uba52rt', 'Cdkn2b/p15ink4b', 'Uba52'),
  # 
  GO_Cellular_Senescence = c('Fbxo4', 'Hmga1b', 'Abl1', 'Pawr', 'Prelp', 'Ang', 'Arg2', 'Bmal1', 'B2m', 'Bcl6', 'Bmpr1a', 
                             'Brca2', 'Calr', 'Cav1', 'Cdk6', 'Cdkn1a/p21', 'Cdkn2a/p16', 'Cdkn2b/p15ink4b', 'Comp', 'Gch1', 'H2-M3', 'Htt', 
                             'H2ax', 'Hmga1', 'Hmga2', 'Hras', 'Id2', 'Ilk', 'Kras', 'Lmna', 'Mapkapk5', 'Mif', 'Mnt', 
                             'Cited2', 'Npm1', 'Nup62', 'Prkcd', 'Pla2r1', 'Pml', 'Prkdc', 'Pten', 'Rbl1', 'Spi1', 'Plk2', 
                             'Srf', 'Suv39h1', 'Tbx2', 'Tbx3', 'Cgas', 'Terf2', 'Tert', 'Top2b', 'Trp53', 'Trp63', 'Twist1', 
                             'Wnt1', 'Wrn', 'Smc5', 'Ybx1', 'Zfp217', 'Zmpste24', 'H2aj', 'Akt3', 'Vash1', 'Nek4', 'Kat6a', 
                             'Zfp277', 'Map2k1', 'Map3k3', 'Mapk14', 'Zmiz1', 'Morc3', 'Mirlet7d', 'Mirlet7g', 'Mirlet7f-1', 
                             'Mirlet7f-2', 'Ndufs6', 'Sirt6', 'Fzr1', 'Mtor', 'Icmt', 'Nampt', 'Ypel3', 'Eef1e1', 'Rsl1d1', 
                             'Abi3', 'Fbxo5', 'Smc6', 'Nsmce2', 'Pnpt1', 'Zkscan3', 'Pot1b', 'Opa1', 'Bcl2l12', 'Nuak1', 
                             'Gins3', 'Ecrg4', 'Dnaja3', 'Wnt16', 'Sirt1', 'Prmt6'),
  
  SenSig_padj0.05 = c("Apoe","App","Bcl2l1","Ccl4","Cd274","Cdkn1a/p21","Cdkn2b/p15ink4b","Cenpb","Csf1r","Cxcl2",
                      "Dpp4","Fxyd5","Gadd45a","Gas2","Igfbp7","Il10","Il18","Il6","Irf7","Lgals3","Lmna","Lpl","Mmp3","Mmp9",
                      "Notch3","Pik3r1", "Plaur","Pten","Serpina3n", "Serpinb1a","Serpine1","Serpine2",
                      "Sparc","Tgfb1","Tgfb2")
)

# Calculate overlap of sen-gene lists and the CosMx 1000x panel
panel_genes <- rownames(thisObj)

# Iterate through the lists and compute overlap stats
overlap_stats <- lapply(names(SenGeneSets), function(set_name) {
  genes <- SenGeneSets[[set_name]]
  found_genes <- intersect(genes, panel_genes)
  
  data.frame(
    GeneSet = set_name,
    Total_Genes = length(genes),
    Genes_Found = length(found_genes),
    Percent_Found = (length(found_genes) / length(genes)) * 100
  )
})

overlap_df <- do.call(rbind, overlap_stats)
overlap_df <- overlap_df %>% arrange(Percent_Found)

# Print a quick summary to console
print("--- Gene Set Representation in 1000-Plex Panel ---")
print(overlap_df)

# Plot the overlap stats
p_coverage <- ggplot(overlap_df, aes(x = reorder(GeneSet, Percent_Found), y = Percent_Found, fill = Percent_Found)) +
  geom_col(color = "black") +
  # Add the overlap count/total fractions directly onto the bars for clarity
  geom_text(aes(label = paste0(Genes_Found, " / ", Total_Genes)), 
            hjust = -0.1, size = 4, fontface = "bold") +
  scale_fill_viridis_c(option = "magma", end = 0.9) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  coord_flip() +
  labs(
    title = "Gene Set Coverage in CosMx 1000-Plex Panel",
    subtitle = "Number of listed genes physically present in the Seurat object",
    x = "Gene Set",
    y = "Percentage Found (%)",
    fill = "Coverage %"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 11)
  )

print(p_coverage)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Analysis: UCell Analysis of SenBrain 2025 list]
DefaultAssay(thisObj) <- "SCT"
# Extract all genes from the active assay (SCT)
panel_genes <- rownames(thisObj)

# Iterate through your list and keep only the genes present in the panel
Filtered_SenGeneSets <- lapply(SenGeneSets, function(genes) {
  intersect(genes, panel_genes)
})

# Optional safety net: Remove any lists that ended up with 0 genes after filtering
Filtered_SenGeneSets <- Filter(function(x) length(x) > 0, Filtered_SenGeneSets)

# Calculate UCell Scores
# Computes scores for ALL filtered lists simultaneously. 
# Creates metadata columns appended with "_UCell"
thisObj <- AddModuleScore_UCell(
  obj = thisObj,
  features = Filtered_SenGeneSets, # Passing the newly filtered list
  assay = "SCT", 
  name = "_UCell",
  maxRank = 500 #added 2.26.2026 after Chase's note in lab meeting. See: https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html#32_The_maxRank_parameter
) # ^ checked median number of features using median(thisObj$nFeature_SCT)
# since median features was 311, changed maxRank to 500

# Extract generated column names for reference
signature.names <- paste0(names(Filtered_SenGeneSets), "_UCell")

# Set the active identity to your corrected spatial clusters
Idents(thisObj) <- "seurat_clusters_CellTyped"


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: UMAP overlay of UCell SenBrain 2025 results]

#**Note: Change name to plot a different list.*
primary_score <- "SenBrain_Mouse_2025_UCell"

UCell_umap_plot <- FeaturePlot(
  object = thisObj,
  features = primary_score,
  reduction = "tf_umap",
  split.by = "age", 
  raster = FALSE,
  keep.scale = "all",
  pt.size = 0.5,
  label = FALSE,
  repel = TRUE
) +
  labs(color = "Enrichment Score") +
  patchwork::plot_annotation(
    title = paste0(primary_score, " Module Scores"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  ) &
  theme(legend.position = "right") &
  scale_colour_gradient2(low = "white", mid = "grey", high = "red", midpoint = 0.2)

print(UCell_umap_plot)

ggsave(filename = paste0(local_WD, "SenScore_ageFacet_UCell_umap_plot.png"), plot = UCell_umap_plot, dpi = 900, height = 8, width = 10)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Boxplot of UCell SenBrain 2025 results]

custom_col <- c(
  "Granule"               = "navy",
  "Granule Mixed_1"       = "dodgerblue",
  "Granule Mixed_2"       = "deepskyblue",
  "Purkinje"              = "red",
  "Bergmann Glia"         = "gold",     
  "MLI_1"                 = "limegreen",  
  "MLI_2"                 = "darkgreen", 
  "Oligodendrocyte"       = "darkorchid4",
  "OPC"                   = "mediumpurple1",
  "Astrocyte"             = "darkorange",
  "DCN"                   = "red2",  
  "Oligo-Granule Mixed"   = "saddlebrown",
  "Oligo-Astrocyte Mixed" = "palevioletred",
  "Vasculature"           = "cyan",
  "Endothelial/Pericyte"  = "cyan3",
  "Vasculature Mixed"     = "cadetblue",
  "Microglia"             = "black",      
  "Choroid Plexus"        = "pink"
)

library(forcats)

# Extract age and score data
plot_df <- data.frame(
  Identity = Idents(thisObj),
  Age = thisObj$age,
  UCell_Score = thisObj[[primary_score, drop = TRUE]]
)

# Factor Age and order young/old
plot_df$Age <- factor(plot_df$Age, levels = c("YNG", "OLD"))

plot_df$Identity <- fct_reorder(plot_df$Identity, plot_df$UCell_Score, .fun = median, .desc = FALSE)

stat_results <- plot_df %>%
  group_by(Identity) %>%
  summarise(
    # Swapped to Wilcoxon Rank-Sum test for non-parametric, skewed data
    p_val = tryCatch(wilcox.test(UCell_Score ~ Age)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    # Adjust for multiple comparisons
    p_adj = p.adjust(p_val, method = "fdr"),
    # Generate significance stars
    sig_label = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

UCell_boxplot <- ggplot(plot_df, aes(x = UCell_Score, y = Identity, fill = Age)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) + 
  geom_text(data = stat_results, aes(x = -Inf, y = Identity, label = sig_label), 
            inherit.aes = FALSE, hjust = -0.5, vjust = 0.5, size = 5, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  
  labs(
    title = "UCell Enrichment",
    subtitle = paste0(primary_score, "\n Wilcoxon Rank-Sum: * FDR<0.05, ** <0.01, *** <0.001"),
    x = "UCell Enrichment Score",
    y = " ",
    fill = "Age"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
  )

print(UCell_boxplot)

ggsave(filename = paste0(local_WD, "SenBrain_UCell_boxplot.png"), plot = UCell_boxplot, dpi = 900, height = 10, width = 8)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Boxplot of UCell SenBrain 2025 results]

# Calculate YNG Baselines & Percent Change
# Establish the YNG baseline mean for each cluster
baseline_df <- plot_df %>%
  filter(Age == "YNG") %>%
  group_by(Identity) %>%
  summarise(YNG_Baseline = mean(UCell_Score, na.rm = TRUE), .groups = "drop")

# Calculate the % change for every OLD cell relative to its cluster's YNG baseline
change_df <- plot_df %>%
  filter(Age == "OLD") %>%
  left_join(baseline_df, by = "Identity") %>%
  mutate(
    # Safety net to prevent division by absolute zero
    Safe_Baseline = ifelse(YNG_Baseline == 0, .Machine$double.xmin, YNG_Baseline),
    Percent_Change = ((UCell_Score - Safe_Baseline) / Safe_Baseline) * 100
  )

# Reorder clusters based on the median percent change so the biggest increasers are at the top
change_df$Identity <- fct_reorder(change_df$Identity, change_df$Percent_Change, .fun = median, .desc = FALSE)

# Boxplot
pct_change_boxplot <- ggplot(change_df, aes(x = Percent_Change, y = Identity, fill = Identity)) +
  geom_boxplot(outlier.shape = NA) + 
  # Add a vertical dashed line at 0% to see positive vs negative shifts
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  # Map the significance stars from the previous Wilcoxon test
  geom_text(data = stat_results, aes(x = -Inf, y = Identity, label = sig_label), 
            inherit.aes = FALSE, hjust = -0.5, vjust = 0.5, size = 5, fontface = "bold") +
  # Padding to prevent the boxes from overlapping the stars
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  scale_fill_manual(values = custom_col) + 
  labs(
    title = "Magnitude of SenBrain Shift (OLD vs YNG)",
    subtitle = paste0(primary_score, "\n% Change relative to YNG Mean | * FDR<0.05, ** <0.01, *** <0.001"),
    x = "Percent Change (%)",
    y = " "
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "none" # Removed legend since Y-axis labels already show identity
  )

print(pct_change_boxplot)

ggsave(filename = paste0(local_WD, "SenBrain_Magnitude_Shift_boxplot.png"), plot = pct_change_boxplot, dpi = 900, height = 10, width = 8)

# Output a quck summary table
summary_table <- change_df %>%
  group_by(Identity) %>%
  summarise(Median_Pct_Change = median(Percent_Change, na.rm = TRUE)) %>%
  left_join(stat_results, by = "Identity") %>%
  arrange(desc(Median_Pct_Change))

cat("\n=== Median Percent Change per Cluster ===\n")
print(summary_table)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Faceted UCell Score Distribution plots]

p_full_dist <- ggplot(plot_df, aes(x = UCell_Score, fill = Age)) +
  # Use density instead of histograms for the overview to prevent visual clutter
  geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  facet_wrap(~ Identity, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "SenBrain UCell Score Distributions per Cluster",
    subtitle = "Comparing YNG vs OLD Overlays",
    x = "SenBrain UCell Score", 
    y = "Density"
  )

print(p_full_dist)

ggsave(filename = paste0(local_WD, "SenBrain_distributions.png"), plot = p_full_dist, dpi = 900, height = 7, width = 6)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Figure: Individual cluster UCell Score Distribution plots]
# includes more stat details
individual_dist_plots <- list()
unique_clusters <- levels(plot_df$Identity)

cat("\n=== Generating Detailed Distribution Plots with Stats ===\n")

for (clust in unique_clusters) {
  
  # Subset data for the specific cluster
  sub_df <- plot_df %>% filter(Identity == clust)
  
  # Skip if a cluster lacks cells for one of the ages (prevents stat crashes)
  if(length(unique(sub_df$Age)) < 2) next 
  
  # Estimate max density to scale the x-axis tick marks
  max_dens <- max(density(sub_df$UCell_Score, na.rm = TRUE)$y)
  tick_height <- max_dens * 0.08 # Ticks will be 8% of the plot height
  
  # Calculate mean and standard devs
  stats_df <- sub_df %>%
    group_by(Age) %>%
    summarize(
      Mean = mean(UCell_Score, na.rm = TRUE),
      SD = sd(UCell_Score, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create a dataframe for the SD tick marks
  sd_multipliers <- c(-2.5, -2, -1, 1, 2, 2.5)
  sd_marks <- expand_grid(Age = c("YNG", "OLD"), SD_Val = sd_multipliers) %>%
    left_join(stats_df, by = "Age") %>%
    mutate(Mark_Pos = Mean + (SD_Val * SD))
  
  # Extract mean values for stats
  y_mean <- stats_df %>% filter(Age == "YNG") %>% pull(Mean)
  o_mean <- stats_df %>% filter(Age == "OLD") %>% pull(Mean)
  
  # Calculate Stat
  # Wilcoxon rank-sum test for significance
  test_res <- tryCatch(wilcox.test(UCell_Score ~ Age, data = sub_df), error = function(e) list(p.value = NA))
  p_val <- test_res$p.value
  
  # Percent Change
  safe_y_mean <- ifelse(y_mean == 0, .Machine$double.xmin, y_mean)
  pct_change <- ((o_mean - safe_y_mean) / safe_y_mean) * 100
  
  # Format the subtitle
  dynamic_subtitle <- sprintf(
    "Wilcox p-val: %.2e | Mean Shift (OLD vs YNG): %+.1f%%\nDashed Line: Mean | Ticks on X-axis: \u00B11, 2, 2.5 SD", 
    p_val, pct_change
  )
  
  # Plot
  p_ind <- ggplot(sub_df, aes(x = UCell_Score, fill = Age, color = Age)) +
    # Overlaid Histograms
    geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.4, position = "identity") +
    # Density Curves (outline only)
    geom_density(alpha = 0, linewidth = 1) +
    # Mean Lines
    geom_vline(data = stats_df, aes(xintercept = Mean, color = Age), linetype = "dashed", linewidth = 1) +
    # Standard Deviation Marks (Small ticks on the X-axis)
    geom_segment(data = sd_marks, aes(x = Mark_Pos, xend = Mark_Pos, y = 0, yend = tick_height, color = Age), 
                 linewidth = 1, linetype = "solid") +
    # Text Annotations (Staggered vjust so YNG and OLD text don't overlap)
    # Means
    annotate("text", x = y_mean, y = Inf, label = "", vjust = 1.5, hjust = 1.1, color = "dodgerblue", fontface="bold") +
    annotate("text", x = o_mean, y = Inf, label = "", vjust = 3.0, hjust = 1.1, color = "red", fontface="bold") +
    scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) +
    scale_color_manual(values = c("YNG" = "dodgerblue", "OLD" = "red")) +
    labs(
      title = paste("SenBrain Signal Distribution:", clust),
      subtitle = dynamic_subtitle,
      x = "SenBrain UCell Score", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(face = "bold", color = "gray20")
    )
  
  # Save to list
  individual_dist_plots[[clust]] <- p_ind
}

cat("=== All Detailed Plots Generated ===\n")


print(individual_dist_plots[["Purkinje"]])

print(individual_dist_plots[["Microglia"]])

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*

# Giotto Proximity Analysis ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Giotto Proximity Analysis*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
# install.packages("terra")
# devtools::install_github("drieslab/GiottoUtils")
# if(!"pak" %in% installed.packages()) {
#   install.packages("pak")
# }
# pak::pkg_install("drieslab/Giotto")

library(Giotto)
library(dplyr)
library(stringr)

#'[Analysis: Compute Proximity/ICG Analyses]
sample_ids <- unique(thisObj@meta.data$mus) 
rownames(thisObj@meta.data) <- thisObj@meta.data$cell_ID

full_meta_df <- as.data.frame(thisObj@meta.data)
full_counts <- GetAssayData(thisObj, assay = "RNA", layer = "counts")

options(future.globals.maxSize = 3000 * 1024^2) # Setting to 3GB for safety, crashed previously due to memory


# Main Proximity/ICG Analysis loop
sample_proximity_results <- list()
sample_icg_results <- list()


for (sample in sample_ids) {
  cat(paste0("\n--- Processing Sample: ", sample, " ---\n"))
  
  barcodes_to_keep <- full_meta_df$cell_ID[full_meta_df$mus == sample]
  sample_meta <- full_meta_df[full_meta_df$mus == sample, ]
  valid_barcodes <- intersect(barcodes_to_keep, colnames(full_counts))
  
  if(length(valid_barcodes) == 0) next
  
  # Build Object
  g_sample <- createGiottoObject(
    expression = full_counts[, valid_barcodes],
    spatial_locs = data.frame(
      x = sample_meta$x_slide_mm * 1000,
      y = sample_meta$y_slide_mm * 1000,
      row.names = sample_meta$cell_ID
    ),
    cell_metadata = sample_meta
  )
  
  # Network & Proximity
  g_sample <- createSpatialNetwork(gobject = g_sample, method = 'Delaunay', maximum_distance_delaunay = 150)
  cpe_sample <- cellProximityEnrichment(g_sample, cluster_column = "seurat_clusters_CellTyped", number_of_simulations = 2000)
  
  # PROXIMITY EXTRACTION (Adding sex to the metadata here!)
  res_df <- cpe_sample$enrichm_res
  res_df$Sample_ID <- sample
  res_df$Age <- unique(sample_meta$age)
  res_df$sex <- unique(sample_meta$sex) # CRITICAL FIX
  sample_proximity_results[[sample]] <- res_df
  
  # ICG ANALYSIS
  cat("  Calculating ICGs...\n")
  g_sample <- filterGiotto(g_sample, expression_threshold = 1, feat_det_in_min_cells = 20) # Lowered for rare cells
  g_sample <- normalizeGiotto(g_sample)
  
  icg_results <- findInteractionChangedFeats(g_sample, cluster_column = "seurat_clusters_CellTyped", diff_test = "permutation", nr_permutations = 100)
  
  if(!is.null(icg_results)) {
    res_table <- icg_results$ICFscores
    res_table$Sample_ID <- sample
    res_table$Age <- unique(sample_meta$age)
    res_table$sex <- unique(sample_meta$sex) # CRITICAL FIX
    sample_icg_results[[sample]] <- res_table
  }
  rm(g_sample); gc()
}

all_proximity_df <- bind_rows(sample_proximity_results)
all_icg_df <- bind_rows(sample_icg_results)

library(dplyr)

#'[Stats: Run proximity analysis stats]
# Computes both absolute interaction/proximity changes (cell A is "closer" to cell B)
# and computes relative interaction (i.e. if we control for population sizes,
# does the frequency of interaction between cell types change computed against random null?)
# REVERTED TO TIER 1 DISCOVERY: ~ Age model to preserve degrees of freedom for N=8 cohort.
run_stat_proximity_base <- function(df) {
  df$Age <- factor(df$Age, levels = c("YNG", "OLD"))
  interactions <- unique(df$unified_int)
  
  results <- lapply(interactions, function(int) {
    subset_df <- df[df$unified_int == int, ]
    
    # Check for sufficient N
    if (length(unique(subset_df$Age)) == 2 && nrow(subset_df) >= 4) {
      m_rel <- tryCatch({ lm(enrichm ~ Age, data = subset_df) }, error = function(e) NULL)
      m_abs <- tryCatch({ lm(original ~ Age, data = subset_df) }, error = function(e) NULL)
      
      if (!is.null(m_rel) && !is.null(m_abs)) {
        s_rel <- coef(summary(m_rel))
        s_abs <- coef(summary(m_abs))
        
        row_rel <- grep("AgeOLD", rownames(s_rel), value = TRUE)
        row_abs <- grep("AgeOLD", rownames(s_abs), value = TRUE)
        
        if (length(row_rel) > 0) {
          return(data.frame(
            unified_int = int,
            # In a ~Age model, the Estimate is the exact Delta (OLD - YNG)
            delta_relative_interaction = s_rel[row_rel, "Estimate"],
            pval_relative_interaction = s_rel[row_rel, "Pr(>|t|)"],
            delta_abs_count = s_abs[row_abs, "Estimate"],
            pval_abs_count = s_abs[row_abs, "Pr(>|t|)"]
          ))
        }
      }
    }
    return(NULL)
  })
  return(bind_rows(results))
}

dual_stat_df <- run_stat_proximity_base(all_proximity_df) %>%
  mutate(FDR_relative_interaction = p.adjust(pval_relative_interaction, method = "BH"),
         FDR_abs_count = p.adjust(pval_abs_count, method = "BH"))


#'[Stats: Run Interaction Changed Gene (ICG) Analysis stats]
# controlling for Sex.
# Tells us which genes change in response to changing cell-cell proximity/interactions
get_model_stats <- function(d) {
  # 1. Variance Check: Need variation in Age to run the model
  if(length(unique(d$Age)) < 2) return(c(pval = NA, beta = NA))
  
  # 2. Try the model with fallback for single-sex interactions
  model <- tryCatch({
    if(length(unique(d$sex)) > 1) {
      lm(diff ~ Age + sex, data = d)
    } else {
      lm(diff ~ Age, data = d)
    }
  }, error = function(e) NULL)
  
  if(is.null(model)) return(c(pval = NA, beta = NA))
  
  # 3. Safe extraction
  s_stats <- coef(summary(model))
  age_row <- grep("AgeOLD", rownames(s_stats), value = TRUE)
  
  if(length(age_row) > 0) {
    return(c(pval = s_stats[age_row, "Pr(>|t|)"], 
             beta = s_stats[age_row, "Estimate"]))
  } else {
    return(c(pval = NA, beta = NA))
  }
}

run_stat_icg_sex_adj <- function(df) {
  df$Age <- factor(df$Age, levels = c("YNG", "OLD"))
  df$sex <- factor(df$sex)
  
  cat("Running statistical comparison on ICGs (Memory-Safe Adjusted)...\n")
  
  icg_summary <- df %>%
    group_by(feats, cell_type, int_cell_type) %>%
    filter(sum(Age == "YNG") >= 2 & sum(Age == "OLD") >= 2) %>%
    mutate(
      mean_induction_YNG = mean(diff[Age == "YNG"], na.rm = TRUE),
      mean_induction_OLD = mean(diff[Age == "OLD"], na.rm = TRUE),
      delta_proximity_induction = mean_induction_OLD - mean_induction_YNG
    ) %>%
    group_by(mean_induction_YNG, mean_induction_OLD, delta_proximity_induction, .add = TRUE) %>%
    group_nest(.key = "data") %>%
    mutate(stats = lapply(data, get_model_stats)) %>%
    mutate(
      pval_age_adj = sapply(stats, `[[`, "pval"),
      beta_age_adj = sapply(stats, `[[`, "beta")
    ) %>%
    select(-data, -stats) %>%
    filter(!is.na(pval_age_adj)) %>%
    mutate(FDR_age_adj = p.adjust(pval_age_adj, method = "BH")) %>%
    arrange(FDR_age_adj)
  
  return(icg_summary)
}

stat_icg_df <- run_stat_icg_sex_adj(all_icg_df)

#'[Flag interactions that were also significantly changing in abs interaction count or relative]
stat_icg_df <- stat_icg_df %>%
  mutate(unif_int = paste0(cell_type, "--", int_cell_type)) %>%
  mutate(is_prox_sig = unif_int %in% (dual_stat_df %>% 
                                        filter(FDR_relative_interaction < 0.05 | FDR_abs_count < 0.05) %>% 
                                        pull(unified_int)))

#'[Stratified FDR Calculation]
stat_icg_df <- stat_icg_df %>%
  # Group by the specific interaction pair
  group_by(unif_int) %>%
  # Recalculate FDR only within each cell-type's interaction
  mutate(FDR_stratified = p.adjust(pval_age_adj, method = "BH")) %>%
  ungroup() %>%
  arrange(FDR_stratified)

cat("Significant relative interactions:", sum(dual_stat_df$FDR_relative_interaction < 0.05, na.rm = TRUE), "\n")
cat("Significant absolute interaction counts:", sum(dual_stat_df$FDR_abs_count < 0.05, na.rm = TRUE), "\n")
cat("ICG significant genes:", sum(stat_icg_df$FDR_stratified < 0.05, na.rm = TRUE), "\n")

# View(stat_icg_df)











#'[Analysis: Compute Proximity/ICG Analyses]
sample_ids <- unique(thisObj@meta.data$mus) 
rownames(thisObj@meta.data) <- thisObj@meta.data$cell_ID

full_meta_df <- as.data.frame(thisObj@meta.data)
full_counts <- GetAssayData(thisObj, assay = "RNA", layer = "counts")

options(future.globals.maxSize = 3000 * 1024^2) # Setting to 3GB for safety, crashed previously due to memory

# Main Proximity/ICG Analysis loop
sample_proximity_results <- list()
sample_icg_results <- list()

for (sample in sample_ids) {
  cat(paste0("\n--- Processing Sample: ", sample, " ---\n"))
  
  barcodes_to_keep <- full_meta_df$cell_ID[full_meta_df$mus == sample]
  sample_meta <- full_meta_df[full_meta_df$mus == sample, ]
  valid_barcodes <- intersect(barcodes_to_keep, colnames(full_counts))
  
  if(length(valid_barcodes) == 0) next
  
  # Build Object
  g_sample <- createGiottoObject(
    expression = full_counts[, valid_barcodes],
    spatial_locs = data.frame(
      x = sample_meta$x_slide_mm * 1000,
      y = sample_meta$y_slide_mm * 1000,
      row.names = sample_meta$cell_ID
    ),
    cell_metadata = sample_meta
  )
  
  # Network & Proximity - Lowered to default 100um for tighter interaction logic
  g_sample <- createSpatialNetwork(gobject = g_sample, method = 'Delaunay', maximum_distance_delaunay = 150)
  cpe_sample <- cellProximityEnrichment(g_sample, cluster_column = "seurat_clusters_CellTyped", number_of_simulations = 2000)
  
  # PROXIMITY EXTRACTION (Adding sex to the metadata here!)
  res_df <- cpe_sample$enrichm_res
  res_df$Sample_ID <- sample
  res_df$Age <- unique(sample_meta$age)
  res_df$sex <- unique(sample_meta$sex) # CRITICAL FIX
  sample_proximity_results[[sample]] <- res_df
  
  # ICG ANALYSIS
  cat("  Calculating ICGs...\n")
  g_sample <- filterGiotto(g_sample, expression_threshold = 1, feat_det_in_min_cells = 20) 
  g_sample <- normalizeGiotto(g_sample)
  
  # Using permutation to get the initial ICG diff scores per sample
  icg_results <- findInteractionChangedFeats(g_sample, cluster_column = "seurat_clusters_CellTyped", diff_test = "permutation", nr_permutations = 100)
  
  if(!is.null(icg_results)) {
    # Applying canonical ICG filter to remove low-confidence features early
    icg_results <- filterInteractionChangedFeats(icg_results, min_cells = 5, min_int_cells = 5, min_fdr = 0.1)
    
    res_table <- icg_results$ICFscores
    res_table$Sample_ID <- sample
    res_table$Age <- unique(sample_meta$age)
    res_table$sex <- unique(sample_meta$sex)
    sample_icg_results[[sample]] <- res_table
  }
  rm(g_sample); gc()
}

all_proximity_df <- bind_rows(sample_proximity_results)
all_icg_df <- bind_rows(sample_icg_results)

library(dplyr)

#'[Helper: Safe Wilcoxon to prevent loop/dplyr crashes]
safe_wilcox <- function(x, y) {
  tryCatch({
    res <- wilcox.test(x, y)
    res$p.value
  }, error = function(e) NA)
}

#'[Stats: Run proximity analysis stats]
# Switched to Wilcoxon Rank Sum for direct YNG vs OLD comparison
run_stat_proximity_base <- function(df) {
  # Lock factor levels to guarantee directionality (OLD - YNG)
  df$Age <- factor(df$Age, levels = c("YNG", "OLD"))
  interactions <- unique(df$unified_int)
  
  results <- lapply(interactions, function(int) {
    subset_df <- df[df$unified_int == int, ]
    
    if (length(unique(subset_df$Age)) == 2) {
      
      # Extract vectors to avoid formula scoping issues
      enrich_old <- subset_df$enrichm[subset_df$Age == "OLD"]
      enrich_yng <- subset_df$enrichm[subset_df$Age == "YNG"]
      abs_old <- subset_df$original[subset_df$Age == "OLD"]
      abs_yng <- subset_df$original[subset_df$Age == "YNG"]
      
      pval_rel <- safe_wilcox(enrich_old, enrich_yng)
      pval_abs <- safe_wilcox(abs_old, abs_yng)
      
      if (!is.na(pval_rel) && !is.na(pval_abs)) {
        return(data.frame(
          unified_int = int,
          delta_relative_interaction = median(enrich_old, na.rm = TRUE) - median(enrich_yng, na.rm = TRUE),
          pval_relative_interaction = pval_rel,
          delta_abs_count = median(abs_old, na.rm = TRUE) - median(abs_yng, na.rm = TRUE),
          pval_abs_count = pval_abs
        ))
      }
    }
    return(NULL)
  })
  return(bind_rows(results))
}

dual_stat_df <- run_stat_proximity_base(all_proximity_df) %>%
  mutate(FDR_relative_interaction = p.adjust(pval_relative_interaction, method = "BH"),
         FDR_abs_count = p.adjust(pval_abs_count, method = "BH"))

#'[Stats: Run Interaction Changed Gene (ICG) Analysis stats]
# Directly comparing the 'diff' (induction) across samples via Wilcoxon
run_stat_icg_wilcox <- function(df) {
  cat("Running Wilcoxon comparison on ICGs...\n")
  
  icg_summary <- df %>%
    group_by(feats, cell_type, int_cell_type) %>%
    # Ensure enough samples per group to run a test
    filter(sum(Age == "YNG") >= 3 & sum(Age == "OLD") >= 3) %>%
    summarise(
      mean_diff_YNG = mean(diff[Age == "YNG"], na.rm = TRUE),
      mean_diff_OLD = mean(diff[Age == "OLD"], na.rm = TRUE),
      delta_induction = mean_diff_OLD - mean_diff_YNG,
      # Replaced formula with safe vector extraction
      pval_age = safe_wilcox(diff[Age == "OLD"], diff[Age == "YNG"]),
      .groups = 'drop'
    ) %>%
    filter(!is.na(pval_age)) %>% # Drop rows where test failed (e.g. 0 variance)
    mutate(FDR_global = p.adjust(pval_age, method = "BH"))
  
  return(icg_summary)
}

stat_icg_df <- run_stat_icg_wilcox(all_icg_df)

#'[Flag interactions that were also significantly changing in abs interaction count or relative]
stat_icg_df <- stat_icg_df %>%
  mutate(unif_int = paste0(cell_type, "--", int_cell_type)) %>%
  mutate(is_prox_sig = unif_int %in% (dual_stat_df %>% 
                                        filter(FDR_relative_interaction < 0.05 | FDR_abs_count < 0.05) %>% 
                                        pull(unified_int)))

#'[Stratified FDR Calculation]
# Keeping as requested, but likely redundant with the global FDR + ICG Filter
stat_icg_df <- stat_icg_df %>%
  group_by(unif_int) %>%
  mutate(FDR_stratified = p.adjust(pval_age, method = "BH")) %>%
  ungroup() %>%
  arrange(FDR_stratified)

cat("Significant relative interactions:", sum(dual_stat_df$FDR_relative_interaction < 0.05, na.rm = TRUE), "\n")
cat("Significant absolute interaction counts:", sum(dual_stat_df$FDR_abs_count < 0.05, na.rm = TRUE), "\n")
cat("ICG significant genes:", sum(stat_icg_df$FDR_stratified < 0.05, na.rm = TRUE), "\n")


#'[Figure: Quadrant plot of absolute interaction count vs relative interaction]
# interaction count is how frequent do the boundaries of the Delaunay network touch
# this is NOT population dependent. i.e. agnostic to population sizes, it simply
# checks does Cell A touch Cell B more/less often between conditions.
# This is what most people think of with "interaction" but this is agnostic to population sizes.
# The relative interaction is what giotto focuses on computing by comparing a 
# randomized scramble of cell positions, it compares random positions as the statistical null
# and compares that to how the cells actually changed in interaction while taking
# into account population sizes. i.e. If there is more of Cell A then naturally this
# interaction may increase, but by taking into account the random "baseline" which
# is generated using the population sizes per sample, this accounts for population sizes.
# This now tells us whether 2 cell types are interacting simply because there are 
# more of these cells vs are they interacting more DESPITE the population changes 
# or if the population is stable, maybe the cells are interacting more/less.

library(ggplot2)
library(ggrepel)

#'[Prep data with Wilcoxon Delta Coefficients]
plot_df <- dual_stat_df %>%
  filter(FDR_relative_interaction < 0.05 | FDR_abs_count < 0.05) %>%
  mutate(
    # Color by the MOST significant metric (min FDR)
    min_FDR = pmin(FDR_relative_interaction, FDR_abs_count, na.rm = TRUE),
    log_min_FDR = -log10(pmax(min_FDR, 1e-10)),
    
    # Weighted shifts using the new Delta values (Median OLD - Median YNG)
    weight_rel = delta_relative_interaction * -log10(pmax(FDR_relative_interaction, 1e-10)),
    weight_abs = delta_abs_count * -log10(pmax(FDR_abs_count, 1e-10))
  )

#'[Generate the Quadrant Plot]
ggplot(plot_df, aes(x = weight_rel, y = weight_abs)) +
  # Using shape 21: 'fill' handles the significance gradient, 'color' handles the outline
  geom_point(aes(fill = log_min_FDR, size = abs(weight_rel) + abs(weight_abs)),
             shape = 21, color = "black", stroke = 0.3) +
  
  # Solid, semi-transparent axis lines
  geom_hline(yintercept = 0, color = "black", alpha = 0.4) +
  geom_vline(xintercept = 0, color = "black", alpha = 0.4) +
  
  # Improved annotation: box.padding adds space, segment.color draws the leader line
  geom_text_repel(aes(label = unified_int),
                  size = 3,
                  box.padding = 1,
                  point.padding = 0.1,
                  segment.color = "grey30",
                  segment.alpha = 0.6,
                  max.overlaps = 15) +
  scale_fill_viridis_c(option = "plasma") + 
  scale_size_continuous(range = c(2, 9)) +
  
  # NOTE: I temporarily removed coord_cartesian(). 
  # Deltas (Medians) are on a different numerical scale than your old LM Betas. 
  # Let ggplot auto-scale first. If the labels are too zoomed out, add back:
  # coord_cartesian(xlim = c(-X, X), ylim = c(-Y, Y)) 
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 13, face = "bold")
  ) +
  labs(
    title = "Weighted Spatial Interaction Shift (Aging)",
    subtitle = "Relative Enrichment vs. Absolute Count (Wilcoxon Delta: OLD - YNG)",
    x = "Relative Interaction Shift (Weighted Delta)",
    y = "Absolute Interaction Count Shift (Weighted Delta)",
    size = "Magnitude of Shift",
    fill = "-log10(min_FDR)"
  )

#'[Figure: Volcano plot of gene interactions that change with proximity]

# Prep data for plotting
plot_icg <- stat_icg_df %>%
  mutate(
    # Create a categorical column for coloring (Switched to Global FDR)
    sig_category = case_when(
      FDR_stratified  < 0.05 & is_prox_sig == TRUE ~ "Sig + Physical Shift",
      FDR_stratified  < 0.05 & is_prox_sig == FALSE ~ "Sig (Static Neighborhood)",
      TRUE ~ "Not Significant"
    ),
    # Cap the -log10(FDR) for visualization if there are extreme outliers
    log_FDR = -log10(pmax(FDR_stratified , 1e-10))
  )

# Get top gene/cell-cell interactions
top_labels <- plot_icg %>%
  filter(FDR_stratified  < 0.05) %>%
  arrange(FDR_stratified ) %>%
  slice_head(n = 20)


ggplot(plot_icg, aes(x = delta_induction, y = log_FDR)) +
  # Plot non-significant points first so they sit in the background
  geom_point(data = filter(plot_icg, sig_category == "Not Significant"), 
             color = "grey85", alpha = 0.5, size = 1) +
  
  # Plot significant points, using shape to denote physical shifts
  geom_point(data = filter(plot_icg, sig_category != "Not Significant"), 
             aes(color = sig_category, shape = as.character(is_prox_sig)), 
             alpha = 0.8, size = 2.5) +
  
  # Threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +
  
  # Custom color legend
  scale_color_manual(values = c("Sig + Physical Shift" = "#D55E00", 
                                "Sig (Static Neighborhood)" = "#0072B2")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) + # Triangle for shift, Circle for static
  
  # Label the top genes
  geom_text_repel(data = top_labels, 
                  aes(label = paste0(feats, "\n(", unif_int, ")")), 
                  size = 3, 
                  box.padding = 0.5,
                  max.overlaps = 20,
                  segment.color = "grey50") +
  
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Interaction Changed Genes (ICG) in Aging Cerebellum",
    subtitle = "Wilcoxon Rank Sum (OLD - YNG) with Global FDR",
    x = "Induction Shift (Wilcoxon Delta: OLD - YNG)",
    y = "-log10(Global FDR)",
    color = "Interaction Status",
    shape = "Neighborhood physically shifted?"
  )



#'[Figure: Heatmap of changing ICG genes/changes with age]
heatmap_data <- stat_icg_df %>%
  filter(FDR_stratified < 0.05) %>%
  group_by(cell_type, int_cell_type) %>%
  summarize(icg_count = n(), .groups = 'drop')

all_cells <- unique(c(stat_icg_df$cell_type, stat_icg_df$int_cell_type))

# Create every possible combination of Influencer and Responder
complete_data <- expand_grid(Responder = all_cells, Influencer = all_cells) %>%
  # Join our actual counts into this perfect grid
  left_join(heatmap_data, by = c("Responder" = "cell_type", "Influencer" = "int_cell_type")) %>%
  # Replace the NAs (missing interactions) with 0
  replace_na(list(icg_count = 0)) %>%
  # Order alphabetically for clean axes
  mutate(Responder = fct_relevel(Responder, sort),
         Influencer = fct_relevel(Influencer, sort))

# 3. Generate the Heatmap with ggplot2
ggplot(complete_data, aes(x = Influencer, y = Responder, fill = icg_count)) +
  geom_tile(color = "grey90") + 
  
  # Add text labels inside tiles
  geom_text(aes(label = ifelse(icg_count > 0, icg_count, "")),
            color = ifelse(complete_data$icg_count > max(complete_data$icg_count, na.rm=T)/1.8, "white", "black"),
            size = 3) +
  
  scale_fill_gradient(low = "white", high = "#800080", na.value = "white") + 
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Directional Communication Landscape in Aging Cerebellum",
    subtitle = "Number of significant Interaction Changed Genes (ICGs) per cell-cell pair (Stratified FDR < 0.05)",
    x = "Influencer Cell Type (Neighbor)",
    y = "Responder Cell Type (Target)",
    fill = "# Sig. ICGs"
  )


#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# CellChat V2: Ligand-Receptor Analysis ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**CellChat V2: Ligand-Receptor Analysis*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]

library(CellChat)
library(CellChat)

# 1. Ensure you have the spatial factors and database ready
spatial_factors <- data.frame(ratio = 1, tol = 10)
CellChatDB <- CellChatDB.mouse

# 2. The Fixed Function
run_spatial_cellchat <- function(seurat_obj, name_prefix) {
  
  cat("Running CellChat for:", name_prefix, "\n")
  
  # WARNING FIX: Tell CellChat which slide each cell is on.
  # This prevents impossible cross-mouse signaling!
  seurat_obj$samples <- as.factor(seurat_obj$mus)
  
  coords <- as.matrix(seurat_obj@meta.data[, c("x_slide_micron", "y_slide_micron")])
  
  # Initialize
  cc <- createCellChat(
    object = seurat_obj,
    coordinates = coords,
    spatial.factors = spatial_factors, 
    group.by = "seurat_clusters_CellTyped",
    assay = "SCT",
    datatype = "spatial"
  )
  
  cc@DB <- CellChatDB
  
  # Standard Pre-processing
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  # ERROR FIX: Reverted back to your original 'interaction.range' argument
  cc <- computeCommunProb(
    cc, 
    type = "truncatedMean", 
    trim = 0.1,
    distance.use = TRUE, 
    interaction.range = 250, # Maximum spatial distance in microns
    scale.distance = 0.9,
    contact.dependent = TRUE, 
    contact.range = 10       # Max distance for direct physical contact
  )
  
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  return(cc)
}

# 3. Execute for Both Groups
# (Assuming obj_yng and obj_old are still in your environment)
cc_yng <- run_spatial_cellchat(obj_yng, "YNG")
cc_old <- run_spatial_cellchat(obj_old, "OLD")


# 1. Merge the objects for comparison
object.list <- list(YNG = cc_yng, OLD = cc_old)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# 2. Visualize what signaling pathways specifically turned ON or OFF in Microglia with age
netVisual_bubble(
  cellchat_merged, 
  targets.use = "Microglia", # Looking at signals received by Microglia
  comparison = c(1, 2),      # Compare YNG (1) vs OLD (2)
  angle.x = 45,
  title.name = "Differential Signaling Received by Microglia (Old vs Young)"
)





# 1. Merge the processed objects
object.list <- list(YNG = cc_yng, OLD = cc_old)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# 2. View the global difference in interaction number and strength
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight")
gg1 + gg2








library(dplyr)

# 1. Extract the merged data as a list
df_merged_list <- subsetCommunication(cellchat_merged, slot.name = "net")

# 2. Bind the list elements together and calculate the delta
mg_diff <- bind_rows(df_merged_list, .id = "group") %>%
  filter(target == "Microglia") %>%
  # We group by the specific interaction name and the cell types involved
  group_by(interaction_name, source, target) %>%
  summarize(
    prob_YNG = sum(prob[group == "YNG"]),
    prob_OLD = sum(prob[group == "OLD"]),
    .groups = 'drop'
  ) %>%
  mutate(diff = prob_OLD - prob_YNG) %>%
  arrange(desc(diff))

# 3. Extract the clean names for the bubble plot
top_3_up <- head(mg_diff$interaction_name, 3)
top_3_down <- tail(mg_diff %>% filter(prob_YNG > 0 | prob_OLD > 0) %>% arrange(diff) %>% pull(interaction_name), 3)
top_LRs <- c(top_3_up, top_3_down)

pair_df <- data.frame(interaction_name = top_LRs)


clean_LRs <- unique(as.character(top_LRs))

# 2. Re-create the data frame for the bubble plot
pair_df <- data.frame(interaction_name = clean_LRs)

# 3. Run bubble plot
netVisual_bubble(
  cellchat_merged, 
  targets.use = "Microglia", 
  comparison = c(1, 2), 
  pairLR.use = pair_df,
  title.name = "Top Unique Aging Signaling Shifts",
  angle.x = 45
)


top_pathway <- "ApoE" # Use the base pathway name

par(mfrow = c(1,2), xpd = TRUE)
netVisual_chord_gene(
  cc_yng, 
  targets.use = "Microglia", 
  signaling = top_pathway, 
  legend.pos.x = 10
)
netVisual_chord_gene(
  cc_old, 
  targets.use = "Microglia", 
  signaling = top_pathway, 
  legend.pos.x = 10
)






# A. Comparison Heatmap (Weight/Strength)
# Shows which cell-cell links are strengthening (Red) or weakening (Blue)
netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1, 2))

# B. Signaling Role Network (e.g., for the APP pathway)
# Shows who is the source vs the sink in the aging brain
pathway_of_interest <- "ApoE" 



netAnalysis_signalingRole_network(
  cellchat_merged, 
  signaling = pathway_of_interest, 
  width = 12, height = 4, 
  comparison = c(1, 2)
)


# C. Aggregate Net Comparison
# This will output two circle plots side-by-side
par(mfrow = c(1,2), xpd=TRUE)
netVisual_aggregate(cc_yng, signaling = pathway_of_interest, layout = "circle", signaling.name = "YNG: ApoE Pathway")
netVisual_aggregate(cc_old, signaling = pathway_of_interest, layout = "circle", signaling.name = "OLD: ApoE Pathway")
par(mfrow = c(1,1), xpd=TRUE)


levels(cc_old@meta$samples)

target_sample <- levels(cc_old@meta$samples)[4] 
target_sample
# 3. Run the spatial plot for that specific slide
netVisual_aggregate(
  cc_old, 
  signaling = pathway_of_interest, 
  layout = "spatial", 
  sources.use = "Microglia", 
  sample.use = target_sample, # THE FIX
  edge.width.max = 5, 
  vertex.size.max = 1,
  alpha.image = 0.2
)



# Visualize the expression of the gene pair across groups
plotGeneExpression(
  cellchat_merged, 
  features = c("Apoe", "Trem2", "Tyrobp"), 
  assay = "RNA", # Use normalized data
  type = "violin",
  group.by = "seurat_clusters_CellTyped", # X-axis
  split.by = "age", # Side-by-side comparison of YNG vs OLD
  color.use = c("OLD" =  "red", "YNG" = "blue")
)


p <- plotGeneExpression(
  cellchat_merged, 
  features = c("Apoe", "Trem2", "Tyrobp"), 
  assay = "RNA", 
  type = "violin",
  group.by = "seurat_clusters_CellTyped", 
  split.by = "age"
)

# 2. Add the statistical comparison
# This compares the 'split.by' groups (YNG vs OLD) within each X-axis category
p_stats <- p + stat_compare_means(
  aes(group = age), 
  label = "p.signif", 
  method = "wilcox.test",
  label.y.npc = 0.5,    # Places label at 90% of the plot height for each facet
  hide.ns = FALSE,      # Ensures 'ns' shows if not significant
  size = 5              # Makes the stars/text larger and visible
)

print(p_stats)












library(rstatix)
library(dplyr)
library(tidyr)

# 1. Extract the expression matrix for our genes of interest
genes_of_interest <- c("Apoe", "Trem2", "Tyrobp")
exp_data <- FetchData(thisObj, vars = c(genes_of_interest, "seurat_clusters_CellTyped", "age"))

# 2. Run the Wilcoxon Rank Sum Test (Age comparison within each Cell Type)
stat_results <- exp_data %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression") %>%
  group_by(seurat_clusters_CellTyped, Gene) %>%
  wilcox_test(Expression ~ age) %>%
  # Apply Benjamini-Hochberg (FDR) correction across all tests
  adjust_pvalue(method = "BH") %>%
  add_significance()

# 3. Filter and show the results for Microglia
print("--- Statistical Comparison: Young vs Old in Microglia ---")
stat_results %>% 
  filter(seurat_clusters_CellTyped == "Microglia") %>%
  select(Gene, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)





















# 1. Compute centrality for both objects
cc_yng <- netAnalysis_computeCentrality(cc_yng, slot.name = "netP")
cc_old <- netAnalysis_computeCentrality(cc_old, slot.name = "netP")

# 2. Generate the scatter plots
pathway.use <- "ApoE"

gg1 <- netAnalysis_signalingRole_scatter(cc_yng, signaling = pathway.use, title = "YNG")
gg2 <- netAnalysis_signalingRole_scatter(cc_old, signaling = pathway.use, title = "OLD")

# 3. Display side-by-side
library(patchwork)
gg1 + gg2





pair_df_apoe <- extractEnrichedLR(cellchat_merged, signaling = "ApoE", geneLR.return = TRUE)

plotGeneExpression(
  cellchat_merged, 
  features = c("ApoE", "Trem2", "Tyrobp"), 
  cols = c("red", "blue", "green"),
  assay = "RNA", 
  type = "dot", 
  split.by = "age"
)



#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# Segmentation Shadow Identification/Filtering ####
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
#**Segmentation Shadow Identification/Filtering*
#'[&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&]
library(FNN)
library(UCell)
# "305"
# "YM04"
# "304"
# "YM03"
# "YM02"
# "302"
# "YM01"
# "301"
# "YF04"
# "243"
# "YF03"
# "239"
# "YF02"
# "237"
# "YF01"
# "236"

FlagSpatialOutliers <- function(seurat_obj,
                                target_cluster,
                                sample_col = "mus",
                                cluster_col = "seurat_clusters_CellTyped",
                                spatial_cols = c("x_slide_mm", "y_slide_mm"),
                                dist_threshold_mm = 0.01,
                                quantile_cutoff = 0.05) {
  
  DefaultAssay(seurat_obj) <- "SCT"
  cat(sprintf("\n=== Starting Spatial QC for Cluster %s ===\n", target_cluster))
  
  # --- 1. Global Signature Generation ---
  cat("1. Calculating top markers for the cluster...\n")
  markers <- FindMarkers(seurat_obj, 
                         ident.1 = target_cluster, 
                         group.by = cluster_col,
                         only.pos = TRUE, 
                         min.pct = 0.25, 
                         logfc.threshold = 0.25)
  
  # Extract up to top 50 genes
  top_genes <- head(rownames(markers), 50)
  
  cat("2. Scoring cells using UCell...\n")
  dd_signature <- list(TargetCluster = top_genes)
  score_name <- paste0("DDScore_C", target_cluster)
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features = dd_signature, name = score_name)
  
  # Handle UCell naming quirks
  actual_score_col <- paste0("TargetCluster", score_name)
  if (!actual_score_col %in% colnames(seurat_obj@meta.data)) {
    actual_score_col <- paste0("XTargetCluster", score_name)
  }
  
  # Initialize QC column if it doesn't exist
  if (!"Spatial_QC" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$Spatial_QC <- "Pass"
  }
  
  # --- 2. Iterate Over All Samples ---
  unique_samples <- unique(seurat_obj[[sample_col, drop = TRUE]])
  plot_list <- list()
  
  for (current_sample in unique_samples) {
    cat(sprintf("\nProcessing Sample: %s\n", current_sample))
    
    # Isolate sample data
    sample_cells <- rownames(seurat_obj@meta.data[seurat_obj[[sample_col, drop=TRUE]] == current_sample, ])
    
    if (length(sample_cells) == 0) next
    
    # Isolate cluster data within the sample
    cluster_cells <- sample_cells[seurat_obj[[cluster_col, drop=TRUE]][sample_cells] == target_cluster]
    
    if (length(cluster_cells) == 0) {
      cat("  -> No cells found for this cluster in this sample. Skipping.\n")
      next
    }
    
    # --- 3. Calculate Stats and Distribution per Sample ---
    sample_scores <- seurat_obj@meta.data[cluster_cells, actual_score_col]
    
    cluster_mean <- mean(sample_scores, na.rm = TRUE)
    cluster_median <- median(sample_scores, na.rm = TRUE)
    q05 <- quantile(sample_scores, quantile_cutoff, na.rm = TRUE)
    q01 <- quantile(sample_scores, 0.01, na.rm = TRUE)
    
    score_data <- data.frame(Score = sample_scores)
    
    # Generate Distribution Plot
    p_dist <- ggplot(score_data, aes(x = Score)) +
      geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "gray80", color = "gray50") +
      geom_density(color = "blue", linewidth = 1) +
      geom_vline(xintercept = cluster_median, color = "forestgreen", linetype = "solid", linewidth = 1) +
      geom_vline(xintercept = cluster_mean, color = "red", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = q05, color = "purple", linetype = "dotted", linewidth = 1) +
      geom_vline(xintercept = q01, color = "purple", linetype = "dotted", linewidth = 1) +
      annotate("text", x = cluster_median, y = Inf, label = "Median", vjust = 1.5, hjust = -0.1, color = "forestgreen") +
      annotate("text", x = cluster_mean, y = Inf, label = "Mean", vjust = 1.5, hjust = 1.1, color = "red") +
      annotate("text", x = q05, y = 0, label = paste0(quantile_cutoff * 100, "%"), vjust = -0.5, hjust = -0.1, color = "purple") +
      annotate("text", x = q01, y = 0, label = "1%", vjust = -0.5, hjust = 1.1, color = "purple") +
      labs(
        title = paste("Core Signal Distribution: Cluster", target_cluster, "| Sample:", current_sample),
        subtitle = "Solid Green: Median | Dashed Red: Mean | Dotted Purple: 1st/5th Percentiles",
        x = paste("Cluster", target_cluster, "Core Signature Score (UCell)"), 
        y = "Density"
      ) +
      theme_minimal()
    
    # Split into Core and Outliers
    core_cells <- cluster_cells[sample_scores >= q05]
    outlier_cells <- cluster_cells[sample_scores < q05]
    
    # Tag Core cells
    seurat_obj@meta.data[core_cells, "Spatial_QC"] <- "Core"
    
    if (length(outlier_cells) == 0) {
      cat("  -> No outliers detected below the threshold.\n")
      # Still save the distribution plot even if no outliers
      plot_list[[current_sample]] <- p_dist
      next
    }
    
    # --- 4. Spatial k-NN Search for Outliers ---
    all_sample_meta <- seurat_obj@meta.data[sample_cells, c(spatial_cols, cluster_col)]
    outlier_coords <- all_sample_meta[outlier_cells, spatial_cols]
    
    # k=2 because the cell will find itself as neighbor #1
    nn_res <- get.knnx(data = all_sample_meta[, spatial_cols], 
                       query = outlier_coords, 
                       k = 2)
    
    nn_dist <- nn_res$nn.dist[, 2]
    nn_cluster <- all_sample_meta[nn_res$nn.index[, 2], cluster_col]
    
    # --- 5. Classification Rules ---
    classifications <- case_when(
      nn_cluster == target_cluster ~ "Homotypic_Outlier",
      nn_cluster != target_cluster & nn_dist < dist_threshold_mm ~ "Heterotypic_Seg_Noise",
      nn_cluster != target_cluster & nn_dist >= dist_threshold_mm ~ "Isolated_Outlier"
    )
    
    seurat_obj@meta.data[outlier_cells, "Spatial_QC"] <- classifications
    
    # --- 6. Generate Spatial Verification Plot ---
    spatial_plot_data <- seurat_obj@meta.data[cluster_cells, c(spatial_cols, "Spatial_QC")]
    spatial_plot_data$Spatial_QC <- factor(spatial_plot_data$Spatial_QC, 
                                           levels = c("Core", "Homotypic_Outlier", "Isolated_Outlier", "Heterotypic_Seg_Noise"))
    
    p_spatial <- ggplot(spatial_plot_data, aes(x = .data[[spatial_cols[1]]], y = .data[[spatial_cols[2]]], color = Spatial_QC)) +
      geom_point(data = filter(spatial_plot_data, Spatial_QC == "Core"), size = 0.5, alpha = 0.3) +
      geom_point(data = filter(spatial_plot_data, Spatial_QC == "Homotypic_Outlier"), size = 0.8, alpha = 0.8) +
      geom_point(data = filter(spatial_plot_data, Spatial_QC == "Isolated_Outlier"), size = 0.8, alpha = 0.9) +
      geom_point(data = filter(spatial_plot_data, Spatial_QC == "Heterotypic_Seg_Noise"), size = 0.8, alpha = 0.9) +
      scale_color_manual(values = c("Core" = "grey85", "Homotypic_Outlier" = "orange", 
                                    "Isolated_Outlier" = "dodgerblue", "Heterotypic_Seg_Noise" = "red")) +
      labs(title = sprintf("Spatial Segmentation QC: Cluster %s | Sample: %s", target_cluster, current_sample),
           x = "X (mm)", y = "Y (mm)") +
      theme_minimal() + coord_fixed()
    
    # Combine the plots
    plot_list[[current_sample]] <- p_dist / p_spatial
    
    # Print sample tally
    cat("  -> Tally:\n")
    print(table(seurat_obj@meta.data[cluster_cells, "Spatial_QC"]))
  }
  
  cat("\n=== QC Complete ===\n")
  return(list(object = seurat_obj, plots = plot_list))
}

qc_results_c0 <- FlagSpatialOutliers(seurat_obj = thisObj, target_cluster = "Granule")

thisObj <- qc_results_c0$object

print(qc_results_c0$plots[["YF01"]])


#stop
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
# options(future.globals.maxSize = 10000 * 1024^2)
# 
# thisObj <- FindMultiModalNeighbors(
#   thisObj,
#   reduction.list = list("pca_sct", "pca_TFactivity"),
#   knn.graph.name = "wknn_SCTxTF",
#   snn.graph.name = "wsnn_SCTxTF",
#   dims.list = list(1:25, 1:25),
#   modality.weight.name = c("SCT.weight", "TF.weight")
# )
# 
# 
# thisObj <- RunUMAP(
#   thisObj,
#   graph = "wsnn",     
#   reduction.name = "wnn.umap",
#   reduction.key = "wnnUMAP_"
# )
# 
# 
# thisObj <- RunUMAP(thisObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# 
# thisObj <- FindClusters(thisObj, graph.name = "wsnn", algorithm = 3, random.seed = 42) # resolution = 0.7
# DimPlot(thisObj, reduction = "wnn.umap", group.by = "seurat_clusters")
        # min.dist = 0.00001, spread = 1)



dims2use <- 1:38
k2use <- 20


thisObj <- Seurat::FindNeighbors(thisObj, dims = dims2use, k.param = k2use, assay = "tf_activities", reduction = 'pca_TFactivity')

thisObj <- Seurat::FindClusters(thisObj, resolution = .9, algorithm = 3, random.seed = 42)

thisObj <- Seurat::RunUMAP(thisObj, dims = dims2use, reduction = 'pca_TFactivity',
                   reduction.name = 'tf_umap',
                   min.dist = 0.00001, spread = 1)
DimPlot(thisObj, reduction = "tf_umap", group.by = 'tf_activities_snn_res.0.7', label = TRUE, pt.size = 0.5) + NoLegend()

plan("multicore", workers = 1)
plan()

#'[Merge fragmented clusters]


thisObj$seurat_clusters_merged <- as.character(thisObj$tf_activities_snn_res.0.7) # Create a new metadata column
thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("4", "34", "32", "23", "27", "28", "26", "7", "24", "20", "11", "33", "6", "0")] <- "0"

thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("19", "31", "21", "1")] <- "1"

thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("22", "30", "10")] <- "10"
thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("25", "15")] <- "15"

thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("2", "0")] <- "0"
thisObj$seurat_clusters_merged[thisObj$tf_activities_snn_res.0.7 %in% c("12", "29")] <- "12"

Idents(thisObj) <- "seurat_clusters_merged"


DimPlot(thisObj, reduction = "tf_umap", group.by = 'seurat_clusters_merged', label = TRUE, pt.size = 0.5) + NoLegend()

#'[Save new TF clustered obj]
# saveRDS(thisObj, file.path(local_WD, "COSMONAUT_CB_Seurat_postPCA_umapTFClustered_preCelltyping.RDS"))


thisObj <- readRDS(file.path(local_WD, "COSMONAUT_CB_Seurat_postPCA_umapTFClustered_preCelltyping.RDS"))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Identifying for Shadows/outliers*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

all_markers <- FindAllMarkers(thisObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

data_driven_signatures <- split(top_markers$gene, top_markers$cluster)

thisObj <- AddModuleScore_UCell(thisObj, features = data_driven_signatures, name = "_DDScore")

thisObj$seurat_clusters_merged <- factor(thisObj$seurat_clusters_merged)
levels(thisObj$seurat_clusters_merged)


library(FNN)
library(patchwork)
#'[May consider adding SD markers and Gaussian fit markers on the score distribution plot]
FindSegShadows <- function(seurat_obj,
                                cluster_of_interest,
                                umap_cols = c("tfumap_1", "tfumap_2"),
                                impurity_cutoff = 0.25,
                                inverted_score_cutoff = 0.25,
                                k_neighbors = 20,
                                cluster_col = "seurat_clusters_merged",
                                print_combined = FALSE) {
  
  
  
  #' @title Identify Core and Shadow Cells within a Seurat Cluster
  #' @description Analyzes a specific cluster to distinguish between "core" and "shadow" cells using neighborhood impurity and a UCell score.
  #'
  #' @param seurat_obj A Seurat object.
  #' @param cluster_of_interest The identifier for the cluster to analyze.
  #' @param umap_cols A character vector of length 2 specifying the column names for the UMAP coordinates.
  #' @param impurity_cutoff A numeric value (0-1). A cell is an outlier if the proportion of its neighbors from other clusters exceeds this value.
  #' @param inverted_score_cutoff A numeric value (0-1) for the inverted UCell score cutoff.
  #' @param k_neighbors An integer for the number of nearest neighbors to consider.
  #' @param cluster_col A string specifying the metadata column with cluster IDs.

  # Extract coordinates and cluster IDs for all cells from the Seurat object
  all_coords <- FetchData(seurat_obj, vars = umap_cols)
  all_ids <- seurat_obj[[cluster_col, drop = TRUE]]

  # Define the UCell score column name based on the cluster of interest
  score_col <- paste0(cluster_of_interest, "_DDScore")
  
  # Check if required columns exist in the object
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Score column '", score_col, "' not found in Seurat object metadata.", sep = ""))
  }
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Error: Cluster column '", cluster_col, "' not found in Seurat object metadata.", sep = ""))
  }
  
  # Fetch the relevant metadata for all cells
  plot_metadata <- FetchData(seurat_obj, vars = c(cluster_col, score_col))
  plot_data <- cbind(all_coords, plot_metadata)
  
  cluster_subset_data <- dplyr::filter(plot_data, .data[[cluster_col]] == cluster_of_interest)
  
  if (nrow(cluster_subset_data) == 0) {
    stop(paste("Error: No cells found for cluster '", cluster_of_interest, "'.", sep=""))
  }

  # For each cell in our target cluster, find its k nearest neighbors from the WHOLE dataset
  # The 'get.knnx' function from the FNN package is extremely fast
  neighbors <- get.knnx(data = all_coords, query = cluster_subset_data[, umap_cols], k = k_neighbors)
  
  # For each cell, calculate the "neighborhood impurity"
  # This is the fraction of its neighbors that do NOT belong to the target cluster
  impurity_scores <- apply(neighbors$nn.index, 1, function(n) mean(all_ids[n] != cluster_of_interest))
  cluster_subset_data$neighborhood_impurity <- impurity_scores
  
  # Calculate inverted score and filter for "Core" and "Shadow" cells based on cutoffs
  cluster_subset_data$inverted_score <- 1 - cluster_subset_data[[score_col]]
  
  core_cells <- rownames(filter(cluster_subset_data, neighborhood_impurity < impurity_cutoff & inverted_score < inverted_score_cutoff))
  shadow_cells <- rownames(filter(cluster_subset_data, !(rownames(cluster_subset_data) %in% core_cells)))
  
  # Assign cell status for plotting
  cluster_subset_data$cell_status <- "Shadow"
  cluster_subset_data$cell_status[rownames(cluster_subset_data) %in% core_cells] <- "Core"
  cluster_subset_data$cell_status <- factor(cluster_subset_data$cell_status, levels = c("Shadow", "Core"))
  
  
  # Helper function for 2D density calculation
  get_density <- function(x, y, ...) {
    # Calculate the 2D density distribution
    dens <- MASS::kde2d(x, y, ...)
    # Find the grid location for each point
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    # Combine the grid locations
    ii <- cbind(ix, iy)
    # Return the density value for each point from the grid
    return(dens$z[ii])
  }
  
  # Plot 1: Main UMAP showing core and shadow cells for the specified cluster
  main_plot <- ggplot(plot_data, aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]])) +
    geom_point(color = "lightgray", size = 0.5) +
    geom_point(data = cluster_subset_data, aes(fill = .data[[score_col]], color = cell_status),
               shape = 21, size = 0.75, stroke = 0.7) +
    scale_fill_viridis_c(name = "UCell Score") +
    scale_color_manual(name = "Cell Status", values = c("Shadow" = "red", "Core" = "black")) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(title = paste("Cluster", cluster_of_interest, "Identity Core")) + theme_minimal()
  
  # Plot 2: Distribution of UCell scores within the cluster
  mean_score <- mean(cluster_subset_data[[score_col]])
  sd_score <- sd(cluster_subset_data[[score_col]])
  
  # Create a data frame for annotations to keep the plot code clean
  # This makes it easy to add/remove lines (e.g., +/- 3 SD) if needed
  stats_data <- data.frame(
    intercepts = c(mean_score - 3 * sd_score, mean_score - 2 * sd_score, 
                   mean_score - sd_score, mean_score, mean_score + sd_score, 
                   mean_score + 2 * sd_score, mean_score + 3 * sd_score),
    labels = c("-3 SD", "-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD", "+3 SD")
  )
  stats_data <- stats_data[stats_data$intercepts > 0 & stats_data$intercepts < 1, ]
  
  
  # Plot 2: Distribution of UCell scores with statistical markers
  dist_plot <- ggplot(cluster_subset_data, aes(x = .data[[score_col]])) +
    geom_density(fill = "#440154FF", alpha = 0.5) +
    # Overlay a theoretical normal distribution (Gaussian fit)
    stat_function(
      fun = dnorm,
      args = list(mean = mean_score, sd = sd_score),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    # Add vertical lines for SD markers
    geom_vline(
      data = stats_data,
      aes(xintercept = intercepts),
      linetype = "dotted",
      color = "black"
    ) +
    # Add text labels for the SD markers
    geom_text(
      data = stats_data,
      aes(x = intercepts, y = Inf, label = paste(labels, "\n", round(intercepts, 2))),
      vjust = 1.1,
      hjust = 0.95,
      size = 4
    ) +
    labs(
      title = "UCell Score Distribution",
      subtitle = paste0("Mean = ", round(mean_score, 2), ", SD = ", round(sd_score, 2)),
      x = "UCell Score"
    ) +
    theme_minimal() +
    coord_cartesian(clip = "off") 
  
  # Plot 3: Gating plot showing the impurity and score cutoffs
  gating_plot <- ggplot(cluster_subset_data, aes(x = neighborhood_impurity, y = inverted_score)) +
    geom_jitter(alpha = 0.5, width = 0.02, height = 0) +
    geom_vline(xintercept = impurity_cutoff, linetype = "dashed", color = "red") +
    geom_hline(yintercept = inverted_score_cutoff, linetype = "dashed", color = "darkred") +
    labs(title = "Neighbor Gating", x = "Neighborhood Impurity Score", y = "Inverted UCell Score (1-UCell)") +
    theme_minimal()
  
  # Calculate distance from a weighted centroid for the final plot
  cluster_subset_data$distance_to_centroid <- sqrt(
    (cluster_subset_data[[umap_cols[1]]] - weighted.mean(cluster_subset_data[[umap_cols[1]]], w = cluster_subset_data[[score_col]]^100))^2 +
      (cluster_subset_data[[umap_cols[2]]] - weighted.mean(cluster_subset_data[[umap_cols[2]]], w = cluster_subset_data[[score_col]]^100))^2
  )
  cluster_subset_data$density <- get_density(x = cluster_subset_data$distance_to_centroid,
                                             y = cluster_subset_data$inverted_score, n = 100)
  
  # Plot 4: Cell identity vs. distance from centroid
  distance_plot <- ggplot(cluster_subset_data, aes(x = distance_to_centroid, y = inverted_score)) +
    geom_point(aes(fill = density, color = cell_status), shape = 21, alpha = 0.7, size = 2, stroke = 0.7) +
    scale_fill_gradientn(colors = c("blue", "green", "yellow", "red"), name = "Point Density") +
    scale_color_manual(name = "Cell Status", values = c("Shadow" = "red", "Core" = "black")) +
    geom_density_2d(color = "white", alpha = 0.5) +
    guides(color = "none") +
    labs(title = "Cell Identity-Distance (Contour)", x = "Distance from Centroid", y = "Inverted UCell Score (1-UCell)") +
    theme_minimal()

  if (print_combined) {
    print("Printing Combined Plot. This may take a moment.")
    combined_plot <- (main_plot | dist_plot) /
      (gating_plot | distance_plot)
    
    print(combined_plot)
  } else {
    print("Skipping Plot")
  }
  
  
  results <- list(
    parameters = list(
      cluster_of_interest = cluster_of_interest,
      impurity_cutoff = impurity_cutoff,
      inverted_score_cutoff = inverted_score_cutoff,
      k_neighbors = k_neighbors
    ),
    cells = list(
      core_cells = core_cells,
      shadow_cells = shadow_cells
    ),
    plots = list(
      combined_plot = combined_plot,
      main_plot = main_plot,
      distribution_plot = dist_plot,
      gating_plot = gating_plot,
      distance_identity_plot = distance_plot
    )
  )
  
  return(results)
}

ShadowRes_Cluster_0 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 0,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.25,
  inverted_score_cutoff = 0.42,
  k_neighbors = 30,
  print_combined = TRUE
)
plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_0_combined", ".png"), plot = ShadowRes_Cluster_0$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_0_UMAP", ".png"), plot = ShadowRes_Cluster_0$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_0_UCellDistr", ".png"), plot = ShadowRes_Cluster_0$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_0_NeighborGate", ".png"), plot = ShadowRes_Cluster_0$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_0_ContourMap", ".png"), plot = ShadowRes_Cluster_0$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_1 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 1,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.3,
  inverted_score_cutoff = 0.45,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_1_combined", ".png"), plot = ShadowRes_Cluster_1$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_1_UMAP", ".png"), plot = ShadowRes_Cluster_1$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_1_UCellDistr", ".png"), plot = ShadowRes_Cluster_1$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_1_NeighborGate", ".png"), plot = ShadowRes_Cluster_1$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_1_ContourMap", ".png"), plot = ShadowRes_Cluster_1$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_10 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 10,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.65,
  inverted_score_cutoff = 0.38,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_10_combined", ".png"), plot = ShadowRes_Cluster_10$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_10_UMAP", ".png"), plot = ShadowRes_Cluster_10$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_10_UCellDistr", ".png"), plot = ShadowRes_Cluster_10$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_10_NeighborGate", ".png"), plot = ShadowRes_Cluster_10$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_10_ContourMap", ".png"), plot = ShadowRes_Cluster_10$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)


ShadowRes_Cluster_12 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 12,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.5,
  inverted_score_cutoff = 0.4,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_12_combined", ".png"), plot = ShadowRes_Cluster_12$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_12_UMAP", ".png"), plot = ShadowRes_Cluster_12$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_12_UCellDistr", ".png"), plot = ShadowRes_Cluster_12$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_12_NeighborGate", ".png"), plot = ShadowRes_Cluster_12$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_12_ContourMap", ".png"), plot = ShadowRes_Cluster_12$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_13 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 13,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.3,
  inverted_score_cutoff = 0.33,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_13_combined", ".png"), plot = ShadowRes_Cluster_13$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_13_UMAP", ".png"), plot = ShadowRes_Cluster_13$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_13_UCellDistr", ".png"), plot = ShadowRes_Cluster_13$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_13_NeighborGate", ".png"), plot = ShadowRes_Cluster_13$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_13_ContourMap", ".png"), plot = ShadowRes_Cluster_13$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_14 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 14,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.25,
  inverted_score_cutoff = 0.26,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_14_combined", ".png"), plot = ShadowRes_Cluster_14$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_14_UMAP", ".png"), plot = ShadowRes_Cluster_14$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_14_UCellDistr", ".png"), plot = ShadowRes_Cluster_14$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_14_NeighborGate", ".png"), plot = ShadowRes_Cluster_14$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_14_ContourMap", ".png"), plot = ShadowRes_Cluster_14$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_15 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 15,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.5,
  inverted_score_cutoff = 0.42,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_15_combined", ".png"), plot = ShadowRes_Cluster_15$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_15_UMAP", ".png"), plot = ShadowRes_Cluster_15$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_15_UCellDistr", ".png"), plot = ShadowRes_Cluster_15$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_15_NeighborGate", ".png"), plot = ShadowRes_Cluster_15$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_15_ContourMap", ".png"), plot = ShadowRes_Cluster_15$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_16 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 16,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.4,
  inverted_score_cutoff = 0.32,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_16_combined", ".png"), plot = ShadowRes_Cluster_16$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_16_UMAP", ".png"), plot = ShadowRes_Cluster_16$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_16_UCellDistr", ".png"), plot = ShadowRes_Cluster_16$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_16_NeighborGate", ".png"), plot = ShadowRes_Cluster_16$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_16_ContourMap", ".png"), plot = ShadowRes_Cluster_16$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_17 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 17,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.3,
  inverted_score_cutoff = 0.45,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_17_combined", ".png"), plot = ShadowRes_Cluster_17$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_17_UMAP", ".png"), plot = ShadowRes_Cluster_17$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_17_UCellDistr", ".png"), plot = ShadowRes_Cluster_17$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_17_NeighborGate", ".png"), plot = ShadowRes_Cluster_17$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_17_ContourMap", ".png"), plot = ShadowRes_Cluster_17$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_18 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 18,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.2,
  inverted_score_cutoff = 0.18,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_18_combined", ".png"), plot = ShadowRes_Cluster_18$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_18_UMAP", ".png"), plot = ShadowRes_Cluster_18$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_18_UCellDistr", ".png"), plot = ShadowRes_Cluster_18$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_18_NeighborGate", ".png"), plot = ShadowRes_Cluster_18$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_18_ContourMap", ".png"), plot = ShadowRes_Cluster_18$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_3 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 3,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.5,
  inverted_score_cutoff = 0.3,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_3_combined", ".png"), plot = ShadowRes_Cluster_3$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_3_UMAP", ".png"), plot = ShadowRes_Cluster_3$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_3_UCellDistr", ".png"), plot = ShadowRes_Cluster_3$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_3_NeighborGate", ".png"), plot = ShadowRes_Cluster_3$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_3_ContourMap", ".png"), plot = ShadowRes_Cluster_3$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_5 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 5,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.5,
  inverted_score_cutoff = 0.38,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_5_combined", ".png"), plot = ShadowRes_Cluster_5$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_5_UMAP", ".png"), plot = ShadowRes_Cluster_5$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_5_UCellDistr", ".png"), plot = ShadowRes_Cluster_5$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_5_NeighborGate", ".png"), plot = ShadowRes_Cluster_5$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_5_ContourMap", ".png"), plot = ShadowRes_Cluster_5$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_8 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 8,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.45,
  inverted_score_cutoff = 0.39,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_8_combined", ".png"), plot = ShadowRes_Cluster_8$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_8_UMAP", ".png"), plot = ShadowRes_Cluster_8$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_8_UCellDistr", ".png"), plot = ShadowRes_Cluster_8$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_8_NeighborGate", ".png"), plot = ShadowRes_Cluster_8$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_8_ContourMap", ".png"), plot = ShadowRes_Cluster_8$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)

ShadowRes_Cluster_9 <- FindSegShadows(
  seurat_obj = thisObj,
  cluster_of_interest = 9,
  umap_cols = c("tfumap_1", "tfumap_2"),
  impurity_cutoff = 0.4,
  inverted_score_cutoff = 0.39,
  k_neighbors = 30,
  print_combined = TRUE
)

plot_local <- "C:/Users/m322332/Documents/CosMxPlots/"
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_9_combined", ".png"), plot = ShadowRes_Cluster_9$plots$combined_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_9_UMAP", ".png"), plot = ShadowRes_Cluster_9$plots$main_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_9_UCellDistr", ".png"), plot = ShadowRes_Cluster_9$plots$distribution_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_9_NeighborGate", ".png"), plot = ShadowRes_Cluster_9$plots$gating_plot, dpi = 900, units = "in", width = 20, height = 10)
ggsave(filename = paste0(plot_local, "ShadowRes_Cluster_9_ContourMap", ".png"), plot = ShadowRes_Cluster_9$plots$distance_identity_plot, dpi = 900, units = "in", width = 20, height = 10)


levels(thisObj$seurat_clusters_merged)
shadow_lists <- list(
  "0" = ShadowRes_Cluster_0$cells$shadow_cells,
  "1" = ShadowRes_Cluster_1$cells$shadow_cells,
  "3" = ShadowRes_Cluster_3$cells$shadow_cells,
  "5" = ShadowRes_Cluster_5$cells$shadow_cells,
  "8" = ShadowRes_Cluster_8$cells$shadow_cells,
  "9" = ShadowRes_Cluster_9$cells$shadow_cells,
  "10" = ShadowRes_Cluster_10$cells$shadow_cells,
  "12" = ShadowRes_Cluster_12$cells$shadow_cells,
  "13" = ShadowRes_Cluster_13$cells$shadow_cells,
  "14" = ShadowRes_Cluster_14$cells$shadow_cells,
  "15" = ShadowRes_Cluster_15$cells$shadow_cells,
  "16" = ShadowRes_Cluster_16$cells$shadow_cells,
  "17" = ShadowRes_Cluster_17$cells$shadow_cells,
  "18" = ShadowRes_Cluster_18$cells$shadow_cells
)
all_shadows <- unlist(shadow_lists)
length(thisObj@meta.data$mus)
length(all_shadows)

#'[Move shadows into their own cluster]
thisObj@meta.data$seurat_clusters_merged <- as.character(thisObj@meta.data$seurat_clusters_merged)

# Loop through the named list to assign cells to their unique shadow cluster
for (cluster_num in names(shadow_lists)) {
  # Define the new cluster name, e.g., "Shadow_C0"
  new_cluster_name <- paste0("Shadow_C", cluster_num)
  cells_to_reassign <- shadow_lists[[cluster_num]]
  thisObj@meta.data[cells_to_reassign, "seurat_clusters_merged"] <- new_cluster_name
}
thisObj@meta.data$seurat_clusters_merged <- as.factor(thisObj@meta.data$seurat_clusters_merged)


DimPlot(thisObj, reduction = "tf_umap", group.by = 'seurat_clusters_merged', label = TRUE, pt.size = 0.5) + NoLegend()


library(ggrepel)
# 1. Extract UMAP coordinates and cluster metadata
plot_data <- as.data.frame(Embeddings(thisObj, reduction = "tf_umap"))
plot_data$cluster <- thisObj@meta.data$seurat_clusters_merged

# 2. Identify which clusters are main clusters and which are shadows
all_levels <- levels(plot_data$cluster)
shadow_clusters <- grep("^Shadow_C", all_levels, value = TRUE)
main_clusters <- setdiff(all_levels, shadow_clusters)

# 3. Split the data for separate plotting layers
shadow_data <- plot_data %>% filter(cluster %in% shadow_clusters)
main_data <- plot_data %>% filter(cluster %in% main_clusters)%>%
  mutate(cluster = droplevels(cluster))

main_cluster_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
  "#98df8a", "#ff9896", "#c5b0d5", "#c49c94" 
)

# Important: Ensure the names of your colors match the *actual* main clusters
# This line correctly assigns names to the first N colors, assuming order matches
names(main_cluster_colors) <- levels(main_data$cluster) 

label_positions <- main_data %>%
  group_by(cluster) %>%
  summarize(
    x = median(tfumap_1),
    y = median(tfumap_2)
  )

# 5. Plotting (no changes needed here)
advanced_umap_plot <- ggplot() +
  # Layer 1: Plot the shadow cells first (grey, small, semi-transparent)
  geom_point(data = shadow_data, aes(x = tfumap_1, y = tfumap_2), 
             color = "grey80", size = 0.1, alpha = 0.2) +
  
  # Layer 2: Plot the MAIN cells on top (colored by cluster, larger size)
  geom_point(data = main_data, aes(x = tfumap_1, y = tfumap_2, color = cluster), 
             size = 0.5) +
  
  scale_color_manual(values = main_cluster_colors) +
  
  # Layer 3: Add cluster labels
  ggrepel::geom_text_repel(data = label_positions, aes(x = x, y = y, label = cluster), 
                           box.padding = 0, max.overlaps = Inf,
                           color = "black", fontface = "bold") +
  theme_classic() +
  guides(color = "none") + # A better way to remove the legend than NoLegend()
  labs(x = "TF_UMAP 1", y = "TF_UMAP 2", title = "UMAP of Merged Clusters")

print(advanced_umap_plot)



# ggsave(filename = paste0(local_WD, "TF_umap", ".png"), plot = advanced_umap_plot, dpi = 900, units = "in", width = 8, height = 9)


# saveRDS(thisObj, file.path(local_WD, "COSMONAUT_CB_Seurat_postPCA_umapTFClustered_ShadowFilt_preCelltyping.RDS"))


thisObj <- readRDS(file.path(local_WD, "COSMONAUT_CB_Seurat_postPCA_umapTFClustered_ShadowFilt_preCelltyping.RDS"))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Spatial Mapping - Checking before cell typing*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

thisObj@meta.data$x_slide_micron <- thisObj@meta.data$x_slide_mm * 1000
thisObj@meta.data$y_slide_micron <- thisObj@meta.data$y_slide_mm * 1000


# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"

plotSpatialSample <- function(seurat_obj,
                              sample_id,
                              sample_col,
                              coord_x_micron = "x_slide_micron",
                              coord_y_micron = "y_slide_micron",
                              color_by = "ident",
                              custom_palette = NULL,
                              pt_size = 1,
                              alpha = 0.6,
                              plot_title = NULL,
                              # "ggplot" or "plotly" (interactive plot)
                              return_type = "ggplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (return_type == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' needed for return_type='plotly'. Please install it.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Sample column '", sample_col, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_x_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("X coordinate column '", coord_x_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_y_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Y coordinate column '", coord_y_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!(color_by == "ident" || color_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Coloring column/variable '", color_by, "' not found in Seurat object metadata or is not 'ident'."), call. = FALSE)
  }
  if (!sample_id %in% unique(seurat_obj@meta.data[[sample_col]])) {
    stop(paste("Sample ID '", sample_id, "' not found in column '", sample_col, "'."), call. = FALSE)
  }
  
  print(paste("Subsetting data for sample:", sample_id))
  subset_data <- subset(seurat_obj, subset = !!sym(sample_col) == sample_id)
  
  if (nrow(subset_data@meta.data) == 0) {
    stop(paste("No cells found after subsetting for sample ID '", sample_id, "'."), call. = FALSE)
  }
  
  plot_df <- subset_data@meta.data
  if (color_by == "ident") {
    plot_df$color_values <- as.character(Idents(subset_data))
    color_lab <- "Cluster ID"
  } else {
    plot_df$color_values <- as.character(plot_df[[color_by]])
    color_lab <- color_by
  }
  plot_df$color_values <- factor(plot_df$color_values) # Ensure factor for levels
  
  x_range <- range(plot_df[[coord_x_micron]], na.rm = TRUE)
  y_range <- range(plot_df[[coord_y_micron]], na.rm = TRUE)
  buffer <- 10
  xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
  ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Sample:", sample_id, "- Colored by:", color_lab)
  }
  
  # --- Determine Color Palette ---
  color_levels <- levels(plot_df$color_values)
  num_colors <- length(color_levels)
  
  if (is.null(custom_palette)) {
    # Generate default palette if none provided
    if (num_colors <= 12) {
      color_palette_values <- hcl.colors(num_colors, palette = "Set 3")
    } else {
      color_palette_values <- scales::hue_pal()(num_colors)
    }
    # Name palette for consistent mapping
    names(color_palette_values) <- color_levels
  } else {
    # Use the provided custom palette
    color_palette_values <- custom_palette
    # Optional: Add checks here if needed (e.g., ensure names match levels)
    # but adhering to 'no extra checks' request for now.
  }
  
  
  print("Generating ggplot...")
  p <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]], color = color_values)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_palette_values) + # Use the chosen palette
    labs(title = plot_title,
         x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
    theme_minimal() +
    coord_fixed() +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
  
  if (return_type == "plotly") {
    print("Generating plotly object...")
    plot_df$tooltip_text <- paste("Cell:", rownames(plot_df), "<br>", color_lab, ":", plot_df$color_values)
    p_plotly <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]],
                                    color = color_values, text=tooltip_text)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_manual(values = color_palette_values) + # Use the chosen palette
      labs(title = plot_title, x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
      theme_minimal() + coord_fixed() + xlim(xlim_vals) + ylim(ylim_vals) +
      theme( panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.box = "horizontal") +
      guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
    
    suppressWarnings({
      plotly_obj <- plotly::ggplotly(p_plotly, tooltip = "text")
    })
    return(plotly_obj)
    
  } else {
    print("Returning ggplot object.")
    return(p)
  }
}

levels(thisObj$seurat_clusters_merged)
YM04_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                                         sample_id = "YF02",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "ident",
                                         custom_palette = c(
                                           "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                                           "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
                                           "#98df8a", "#ff9896", "#c5b0d5", "#c49c94" 
                                         ),
                                         pt_size = 1,
                                         alpha = 0.6,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(YM04_spatial_ggplot)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Cell Typing*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#'[Checking canonical cell type genes per cluster]
specific_marker_list <- list(
  # 1. Neuronal Cell Types
  # --- Excitatory Neurons ---
  Granule_Cell = c("Atp2b4", "Reln", "Gap43", "Grin2c"),      # Alternatives for targeted panels. Use combinatorially.
  UBC = c("Eomes", "Calb2", "Grm1", "Plcb4"),                  # Eomes (Tbr2) is the most specific marker.
  
  # --- Inhibitory Neurons ---
  Purkinje_Cell = c("Calb1", "Pcp2", "Car8", "Rora"),           # Pcp2 (L7) and Calb1 are canonical.
  Golgi_Cell = c("Gcgr", "Nrgn", "Gdnf", "Mapt"),               # Gcgr is a reliable Golgi marker.
  Basket_Cell = c("Reln", "Sln", "Pvalb"),                     # Enriched in basket cells vs. stellate cells.
  Stellate_Cell = c("Ntf3", "Penk", "Pvalb"),                  # Enriched in stellate cells. Pvalb is lower than in basket.
  Lugaro_Cell = c("Calb2", "Sox14", "Plcl1", "Vwc2"),           # Rare; Calb2 shared with UBCs.
  
  # 2. Glial Cell Types
  Astrocyte = c("Aldh1l1", "Gfap", "Aqp4", "S100b"),            # Aldh1l1 is a robust pan-astrocyte marker.
  Bergmann_Glia = c("Hopx", "Gdf10", "Slc1a3", "Sox9"),         # Hopx is highly specific to Bergmann glia.
  OPC = c("Pdgfra", "Cspg4", "Olig2", "Tcf7l2"),               # Pdgfra & Cspg4 (NG2) are canonical.
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Klk6", "Myrf"),    # Markers for mature, myelinating cells.
  Microglia = c("Tmem119", "P2ry12", "Cx3cr1", "Csf1r"),        # Tmem119/P2ry12 distinguish from macrophages.
  
  # 3. Vascular & Meningeal Cell Types
  Endothelial_Cell = c("Cldn5", "Pecam1", "Esam", "Flt1"),      # Cldn5 is a key BBB marker.
  Pericyte = c("Pdgfrb", "Rgs5", "Kcnj8", "Vtn"),               # Pdgfrb is more specific than Acta2.
  VSMC = c("Acta2", "Tagln", "Notch3", "Myh11"),                # Tagln is highly specific to smooth muscle.
  Fibroblast_VLMC = c("Col1a1", "Dcn", "Lum", "Slc6a13"),       # Collagen-producing meningeal/vascular cells.
  Choroid_Plexus = c("Ttr", "Folr1", "Aqp1", "Otx2")            # Ttr is the most definitive marker.
)
marker_lists <- specific_marker_list 

# 1. Get all unique TFs that are present in the Seurat object's TF activity assay
all_markers <- unique(unlist(marker_lists))
markers_present <- intersect(all_markers, rownames(thisObj[["SCT"]]))

# Check if we have markers to proceed
if (length(markers_present) == 0) {
  stop("None of the specified TFs were found in the 'SCT' assay.")
}

print(markers_present)
# 2. Calculate the average TF activity for these markers across your clusters
# The result is a matrix of TFs (rows) vs. clusters (columns)
avg_exp_matrix <- AverageExpression(
  thisObj,
  features = markers_present,
  group.by = "seurat_clusters_merged", # Or whichever column holds your cluster IDs
  assay = "SCT",
  slot = "data"
)$SCT # Directly extract the matrix from the list

# 3. Calculate the enrichment score for each TF list
# This loop calculates the mean activity score for each TF list in each cluster.
# It iterates through your list of cell types (e.g., Granule_Cell, Purkinje_Cell).
enrichment_scores <- list()
for (cell_type in names(marker_lists)) {
  
  # Find which TFs from the current list are in our matrix
  tfs_in_list <- intersect(marker_lists[[cell_type]], rownames(avg_exp_matrix))
  
  if (length(tfs_in_list) > 0) {
    # Calculate the mean score for this list of TFs across all clusters
    # The 'colMeans' function averages the scores for the selected TFs for each cluster
    enrichment_scores[[cell_type]] <- colMeans(avg_exp_matrix[tfs_in_list, , drop = FALSE])
  }
}

# 4. Combine the scores into a single matrix for plotting
# The result is a matrix of cell types/lists (rows) vs. clusters (columns)
enrichment_matrix <- do.call(rbind, enrichment_scores)

# At this point, 'enrichment_matrix' should be a valid matrix, not NULL.

pheatmap::pheatmap(
  enrichment_matrix,
  main = "TF List Enrichment per Cluster",
  scale = "row",         # Scale scores by row (cell type) to see which cluster is highest for each list
  cluster_cols = TRUE,   # Cluster columns (cell clusters) to group similar clusters
  cluster_rows = TRUE,   # Cluster rows (cell types) to group similar TF lists
  border_color = "white"
)


#'[Plotting top n markers per cluster]
seurat_obj = thisObj
group_by = "Global_CellTypes" #cluster #s w/ shadows: seurat_clusters_merged
n_genes = 5
shadow_pattern = "^Shadow"

Seurat::Idents(seurat_obj) <- seurat_obj[[group_by, drop = TRUE]]

all_clusters <- levels(Seurat::Idents(seurat_obj))
shadow_clusters <- grep(shadow_pattern, all_clusters, value = TRUE)
main_clusters <- setdiff(all_clusters, shadow_clusters)

seurat_subset <- subset(seurat_obj, idents = main_clusters)

message(paste("Finding markers for", length(main_clusters), "main clusters..."))
markers <- Seurat::FindAllMarkers(seurat_subset, assay = "SCT", only.pos = TRUE)


marks.df <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = n_genes, order_by = avg_log2FC) %>%
  dplyr::ungroup() %>% # It's good practice to ungroup after grouped operations
  dplyr::select(cluster, gene, avg_log2FC, p_val_adj)

top_genes <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = n_genes, order_by = avg_log2FC) %>%
  dplyr::pull(gene) %>%
  unique()

# --- New Code for Averaging and Plotting with pheatmap ---

message("Calculating average expression for top genes across clusters...")

# Get the average expression of the top genes for each cluster
# The output is a list, so we extract the matrix for the default assay (e.g., 'RNA')
avg_exp_matrix <- Seurat::AverageExpression(
  seurat_subset,
  features = top_genes,
  group.by = "ident", # Use the active identity (clusters)
  assay = "SCT",
  slot = "data"
)[[1]] # Extract the first assay's matrix


message(paste("Generating heatmap for", length(top_genes), "unique top genes..."))


custom_cluster_order <- c(
  "Granule (Neuron)",
  "Granule (Vegfa+)",
  "Granule (Mix)",
  "PKJ",
  "Bergmann Glia",
  "MLI (Neuron)",
  "MLI (Vegfa+)",
  "DCN",
  "ODC",
  "Endothelial",
  "Vascular",
  "Microglia",
  "Astrocytes",
  "Choroid Plexus"
)

CB_topGenesPerClusterHeat_grob <- pheatmap::pheatmap(
  avg_exp_matrix[, custom_cluster_order],
  scale = "row",
  main = paste(""),
  cluster_cols = FALSE, 
  cluster_rows = TRUE,  
  color = rev(RColorBrewer::brewer.pal(11, "Spectral")),
  show_rownames = TRUE, 
  fontsize_row = 8,
  angle = 315,
  silent = TRUE 
)
# Making genes italic
CB_topGenesPerClusterHeat_grob$gtable$grobs[[which(CB_topGenesPerClusterHeat_grob$gtable$layout$name == "row_names")]]$gp$fontface <- "italic"

grid::grid.newpage()
grid::grid.draw(CB_topGenesPerClusterHeat_grob)

CB_topGenesPerClusterHeat

file_name <- paste0(local_WD, "CB_topGenesPerClusterHeat", ".png")
png(file = file_name, width = 4, height = 10, units = "in", res = 900)
grid::grid.newpage() # Clears the plotting device
grid::grid.draw(CB_topGenesPerClusterHeat_grob)
dev.off()

#'[Trying to find differences between clusters]
cluster_markers <- FindMarkers(thisObj, 
                               ident.1 = 0, #1
                               ident.2 = 1, #2 
                               test.use = "roc", #default wilcox
                               assay = "SCT", 
                               min.pct = 0.1,
                               logfc.threshold = 0)

# Positive avg_log2FC means the gene is more highly expressed in cluster 1
# Negative avg_log2FC means the gene is more highly expressed in cluster 2
head(rownames(cluster_markers), 30)



# 0 = Granule (Neuron)
# 1 = Granule (Vegfa+)
# 10 = Granule (ODC/Astrocyte Mix)
# 12 = Endothelial
# 13 = MLI (Vegfa+)
# 14 = PKJ
# 15 = Vascular
# 16 = DCN
# 17 = Microglia
# 18 = Choroid Plexus
# 3 = ODC
# 5 = MLI (Neuron)
# 8 = Astrocytes
# 9 = Bergmann Glia

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Apply cell types*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cluster_Labels <- list(
  '0' = "Granule (Neuron)",
  '1' = "Granule (Vegfa+)",
  '10' = "Granule (Mix)",
  '12' = "Endothelial",
  '13' = "MLI (Vegfa+)",
  '14' = "PKJ",
  '15' = "Vascular",
  '16' = "DCN",
  '17' = "Microglia",
  '18' = "Choroid Plexus",
  '3' = "ODC",
  '5' = "MLI (Neuron)",
  '8' = "Astrocytes",
  '9' = "Bergmann Glia"
)

names(Cluster_Labels) <- levels(thisObj)
thisObj <- RenameIdents(thisObj, Cluster_Labels)
thisObj$Global_CellTypes <- Idents(thisObj)
levels(thisObj)
DimPlot(thisObj, reduction = "tf_umap", group.by = "age", raster = FALSE, label = TRUE, pt.size = 0.1) + NoLegend()

cluster_markers <- FindMarkers(thisObj, 
                               ident.1 = "OLD",
                               ident.2 = "YNG",
                               group.by = 'age',
                               features = "Vegfa",
                               subset.ident = "Granule (Vegfa+)", # <-- This is the key change
                               assay = "RNA",
                               min.pct = 0.1,
                               logfc.threshold = 0)

# View the top results
head(cluster_markers)


#'[Plotting top n markers per cluster]
seurat_obj = thisObj
group_by = "Global_CellTypes" #cluster #s w/ shadows: seurat_clusters_merged
n_genes = 5
shadow_pattern = "^Shadow"

Seurat::Idents(seurat_obj) <- seurat_obj[[group_by, drop = TRUE]]

all_clusters <- levels(Seurat::Idents(seurat_obj))
shadow_clusters <- grep(shadow_pattern, all_clusters, value = TRUE)
main_clusters <- setdiff(all_clusters, shadow_clusters)

seurat_subset <- subset(seurat_obj, idents = main_clusters)

message(paste("Finding markers for", length(main_clusters), "main clusters..."))
markers <- Seurat::FindAllMarkers(seurat_subset, assay = "SCT", only.pos = TRUE)


marks.df <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = n_genes, order_by = avg_log2FC) %>%
  dplyr::ungroup() %>% # It's good practice to ungroup after grouped operations
  dplyr::select(cluster, gene, avg_log2FC, p_val_adj)

top_genes <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = n_genes, order_by = avg_log2FC) %>%
  dplyr::pull(gene) %>%
  unique()

# --- New Code for Averaging and Plotting with pheatmap ---

message("Calculating average expression for top genes across clusters...")

# Get the average expression of the top genes for each cluster
# The output is a list, so we extract the matrix for the default assay (e.g., 'RNA')
avg_exp_matrix <- Seurat::AverageExpression(
  seurat_subset,
  features = top_genes,
  group.by = "ident", # Use the active identity (clusters)
  assay = "SCT",
  slot = "data"
)[[1]] # Extract the first assay's matrix


message(paste("Generating heatmap for", length(top_genes), "unique top genes..."))


custom_cluster_order <- c(
  "Granule (Neuron)",
  "Granule (Vegfa+)",
  "Granule (Mix)",
  "PKJ",
  "Bergmann Glia",
  "MLI (Neuron)",
  "MLI (Vegfa+)",
  "DCN",
  "ODC",
  "Endothelial",
  "Vascular",
  "Microglia",
  "Astrocytes",
  "Choroid Plexus"
)

CB_topGenesPerClusterHeat_grob <- pheatmap::pheatmap(
  avg_exp_matrix[, custom_cluster_order],
  scale = "row",
  main = paste(""),
  cluster_cols = FALSE, 
  cluster_rows = TRUE,  
  color = rev(RColorBrewer::brewer.pal(11, "Spectral")),
  show_rownames = TRUE, 
  fontsize_row = 8,
  angle = 315,
  silent = TRUE 
)
# Making genes italic
CB_topGenesPerClusterHeat_grob$gtable$grobs[[which(CB_topGenesPerClusterHeat_grob$gtable$layout$name == "row_names")]]$gp$fontface <- "italic"

grid::grid.newpage()
grid::grid.draw(CB_topGenesPerClusterHeat_grob)

CB_topGenesPerClusterHeat

file_name <- paste0(local_WD, "CB_topGenesPerClusterHeat", ".png")
png(file = file_name, width = 4, height = 10, units = "in", res = 900)
grid::grid.newpage() # Clears the plotting device
grid::grid.draw(CB_topGenesPerClusterHeat_grob)
dev.off()

# saveRDS(thisObj, file.path(local_WD, "COSMONAUT_CB_Seurat_umapTFClustered_GlobalCellTypes.RDS"))


#'[Save point with Cell type clusters]
thisObj <- readRDS(file.path(local_WD, "COSMONAUT_CB_Seurat_umapTFClustered_GlobalCellTypes.RDS"))

custom_colors <- c(
  "Granule (Neuron)" = "blue4",
  "Granule (Vegfa+)" = "lightblue",
  "Granule (Mix)" = "purple",
  "Endothelial" = "grey30",
  "MLI (Vegfa+)" = "green",
  "PKJ" = "red",
  "Vascular" = "salmon",
  "DCN" = "red3",
  "Microglia" = "black",
  "Choroid Plexus" = "cyan", 
  "ODC" = "yellow2",
  "MLI (Neuron)" = "darkgreen",
  "Astrocytes" = "orchid",
  "Bergmann Glia" = "orange"
)

plot_data <- as.data.frame(Embeddings(thisObj, reduction = "tf_umap"))
plot_data$seurat_clusters_merged <- thisObj@meta.data$seurat_clusters_merged
plot_data$cell_type <- thisObj@meta.data$Global_CellTypes

all_levels <- unique(plot_data$seurat_clusters_merged)
shadow_clusters <- grep("^Shadow_C", all_levels, value = TRUE)

shadow_data <- plot_data %>% filter(seurat_clusters_merged %in% shadow_clusters)
main_data <- plot_data %>% filter(!seurat_clusters_merged %in% shadow_clusters)%>%
  mutate(seurat_clusters_merged = droplevels(seurat_clusters_merged))

label_positions <- main_data %>%
  group_by(cell_type) %>%
  summarize(
    x = median(tfumap_1),
    y = median(tfumap_2)
  )

advanced_umap_plot <- ggplot() +
  geom_point(data = shadow_data, aes(x = tfumap_1, y = tfumap_2),
             color = "grey80", size = 0.1, alpha = 0.2) +
  geom_point(data = main_data, aes(x = tfumap_1, y = tfumap_2, color = cell_type),
             size = 0.5) +
  scale_color_manual(values = custom_colors) +
  ggrepel::geom_text_repel(
    data = label_positions,
    aes(x = x, y = y, label = cell_type),
    color = "black",
    fontface = "bold",
    size = 6,
    max.overlaps = Inf,
    box.padding = 0,     
    point.padding = 0.5,   
    nudge_y = 0, 
    nudge_x = 0,
    segment.color = 'grey50',
    segment.alpha = 0.7
  ) +
  theme_classic() +
  guides(color = "none") +
  labs(x = "TF_UMAP 1", y = "TF_UMAP 2", title = "UMAP with Updated Cell Identities")

print(advanced_umap_plot)

#'[Without labels]
advanced_umap_plot_noLabels <- ggplot() +
  geom_point(data = shadow_data, aes(x = tfumap_1, y = tfumap_2),
             color = "grey80", size = 0.1, alpha = 0.2) +
  geom_point(data = main_data, aes(x = tfumap_1, y = tfumap_2, color = cell_type),
             size = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  guides(color = "none") +
  labs(x = "TF_UMAP 1", y = "TF_UMAP 2", title = "UMAP with Updated Cell Identities")

print(advanced_umap_plot_noLabels)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped", ".png")
ggsave(filename = file_name, plot = advanced_umap_plot_noLabels, dpi = 900, units = "in", width = 10, height = 11)



umap_age <- DimPlot(thisObj,
                    group.by = "age",
                    reduction = "tf_umap",
                    raster = FALSE,
                    label = FALSE,
                    pt.size = 0.5,
                    repel = TRUE) +
  scale_color_manual(values = c("YNG" = "lightblue", "OLD" = "pink2")) +
  ggtitle("UMAP | Age") +
  theme(plot.title = element_text(hjust = 0.5)) # Center title if desired

print(umap_age)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped_Age", ".png")
ggsave(filename = file_name, plot = umap_age, dpi = 900, units = "in", width = 10, height = 11)


thisObj$age <- factor(thisObj$age, levels = c("YNG", "OLD"))
umap_age_facet <- DimPlot(thisObj,
                          group.by = "age",
                          split.by =  "age",
                          reduction = "tf_umap",
                          raster = TRUE,
                          label = FALSE,
                          pt.size = 0.5,
                          repel = TRUE) + 
  scale_color_manual(values = c("YNG" = "deepskyblue", "OLD" = "red")) +
  ggtitle("UMAP | Age") +
  theme(plot.title = element_text(hjust = 0.5))

print(umap_age_facet)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped_Age_facet", ".png")
ggsave(filename = file_name, plot = umap_age_facet, dpi = 900, units = "in", width = 20, height = 11)



# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"

plotSpatialSample <- function(seurat_obj,
                              sample_id,
                              sample_col,
                              coord_x_micron = "x_slide_micron",
                              coord_y_micron = "y_slide_micron",
                              color_by = "ident",
                              custom_palette = NULL,
                              pt_size = 1,
                              alpha = 0.6,
                              plot_title = NULL,
                              # "ggplot" or "plotly" (interactive plot)
                              return_type = "ggplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (return_type == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' needed for return_type='plotly'. Please install it.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Sample column '", sample_col, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_x_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("X coordinate column '", coord_x_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_y_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Y coordinate column '", coord_y_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!(color_by == "ident" || color_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Coloring column/variable '", color_by, "' not found in Seurat object metadata or is not 'ident'."), call. = FALSE)
  }
  if (!sample_id %in% unique(seurat_obj@meta.data[[sample_col]])) {
    stop(paste("Sample ID '", sample_id, "' not found in column '", sample_col, "'."), call. = FALSE)
  }
  
  print(paste("Subsetting data for sample:", sample_id))
  subset_data <- subset(seurat_obj, subset = !!sym(sample_col) == sample_id)
  
  if (nrow(subset_data@meta.data) == 0) {
    stop(paste("No cells found after subsetting for sample ID '", sample_id, "'."), call. = FALSE)
  }
  
  plot_df <- subset_data@meta.data
  if (color_by == "ident") {
    plot_df$color_values <- as.character(Idents(subset_data))
    color_lab <- "Cluster ID"
  } else {
    plot_df$color_values <- as.character(plot_df[[color_by]])
    color_lab <- color_by
  }
  plot_df$color_values <- factor(plot_df$color_values) # Ensure factor for levels
  
  x_range <- range(plot_df[[coord_x_micron]], na.rm = TRUE)
  y_range <- range(plot_df[[coord_y_micron]], na.rm = TRUE)
  buffer <- 10
  xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
  ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Sample:", sample_id, "- Colored by:", color_lab)
  }
  
  # --- Determine Color Palette ---
  color_levels <- levels(plot_df$color_values)
  num_colors <- length(color_levels)
  
  if (is.null(custom_palette)) {
    # Generate default palette if none provided
    if (num_colors <= 12) {
      color_palette_values <- hcl.colors(num_colors, palette = "Set 3")
    } else {
      color_palette_values <- scales::hue_pal()(num_colors)
    }
    # Name palette for consistent mapping
    names(color_palette_values) <- color_levels
  } else {
    # Use the provided custom palette
    color_palette_values <- custom_palette
    # Optional: Add checks here if needed (e.g., ensure names match levels)
    # but adhering to 'no extra checks' request for now.
  }
  
  
  print("Generating ggplot...")
  p <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]], color = color_values)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_palette_values) + # Use the chosen palette
    labs(title = plot_title,
         x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
    theme_minimal() +
    coord_fixed() +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
  
  if (return_type == "plotly") {
    print("Generating plotly object...")
    plot_df$tooltip_text <- paste("Cell:", rownames(plot_df), "<br>", color_lab, ":", plot_df$color_values)
    p_plotly <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]],
                                    color = color_values, text=tooltip_text)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_manual(values = color_palette_values) + # Use the chosen palette
      labs(title = plot_title, x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
      theme_minimal() + coord_fixed() + xlim(xlim_vals) + ylim(ylim_vals) +
      theme( panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.box = "horizontal") +
      guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
    
    suppressWarnings({
      plotly_obj <- plotly::ggplotly(p_plotly, tooltip = "text")
    })
    return(plotly_obj)
    
  } else {
    print("Returning ggplot object.")
    return(p)
  }
}


thisObj@meta.data$x_slide_micron <- thisObj@meta.data$x_slide_mm * 1000
thisObj@meta.data$y_slide_micron <- thisObj@meta.data$y_slide_mm * 1000

custom_col <- c(
  "Granule (Neuron)" = "blue4",
  "Granule (Vegfa+)" = "lightblue",
  "Granule (Mix)" = "purple",
  "Endothelial" = "grey30",
  "MLI (Vegfa+)" = "green",
  "PKJ" = "red",
  "Vascular" = "salmon",
  "DCN" = "red3",
  "Microglia" = "black",
  "Choroid Plexus" = "cyan", 
  "ODC" = "yellow2",
  "MLI (Neuron)" = "darkgreen",
  "Astrocytes" = "orchid",
  "Bergmann Glia" = "orange"
)

custom_col_highlight <- c(
  "Granule (Neuron)" = "grey",
  "Granule (Vegfa+)" = "grey",
  "Granule (Mix)" = "red",
  "Endothelial" = "grey",
  "MLI (Vegfa+)" = "grey",
  "PKJ" = "grey",
  "Vascular" = "grey",
  "DCN" = "grey",
  "Microglia" = "black",
  "Choroid Plexus" = "grey", 
  "ODC" = "grey",
  "MLI (Neuron)" = "grey",
  "Astrocytes" = "orchid",
  "Bergmann Glia" = "grey"
)


YM04_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                                         sample_id = "YM04",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "ident",
                                         custom_palette = custom_col_highlight,
                                         pt_size = 1,
                                         alpha = 1,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(YM04_spatial_ggplot)


file_name <- paste0(local_WD, "YM04_spatial_ggplot_HI", ".png")
ggsave(filename = file_name, plot = YM04_spatial_ggplot, dpi = 900, units = "in", width = 10, height = 11)



O243_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                                         sample_id = "243",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "ident",
                                         custom_palette = custom_col_highlight,
                                         pt_size = 1,
                                         alpha = 1,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(O243_spatial_ggplot)


file_name <- paste0(local_WD, "O243_spatial_ggplot_HI", ".png")
ggsave(filename = file_name, plot = O243_spatial_ggplot, dpi = 900, units = "in", width = 10, height = 11)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Cell Populations*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Pull cell count data]
library(ggpubr)
meta_data <- thisObj@meta.data
meta_data$ClusterID <- Idents(thisObj)

sample_col <- "mus"
age_col <- "age"
sex_col <- "sex" 
cluster_col <- "ClusterID"

sample_totals <- meta_data %>%
  group_by(!!sym(sample_col)) %>%
  summarise(TotalCellsInSample = n(), .groups = 'drop')

cluster_counts_per_sample <- meta_data %>%
  group_by(!!sym(sample_col), !!sym(cluster_col), !!sym(age_col), !!sym(sex_col)) %>% # Added sex_col
  summarise(CellCount = n(), .groups = 'drop')

plot_data <- cluster_counts_per_sample %>%
  left_join(sample_totals, by = sample_col) %>%
  mutate(Proportion = CellCount / TotalCellsInSample) %>%
  mutate(
    ClusterID = factor(!!sym(cluster_col)),
    Age = factor(!!sym(age_col), levels = c("YNG", "OLD")), 
    Sex = factor(!!sym(sex_col))
  )

custom_cluster_order <- c(
  "Granule (Neuron)",
  "Granule (Vegfa+)",
  "Granule (Mix)",
  "PKJ",
  "Bergmann Glia",
  "MLI (Neuron)",
  "MLI (Vegfa+)",
  "DCN",
  "ODC",
  "Endothelial",
  "Vascular",
  "Microglia",
  "Astrocytes",
  "Choroid Plexus"
)
# Apply this order to the data frame before plotting
plot_data$ClusterID <- factor(plot_data$ClusterID, levels = custom_cluster_order)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[All clusters - Main population plot]
plot_final <- ggplot(plot_data, aes(x = ClusterID, y = Proportion, fill = factor(Age, levels = c("YNG", "OLD")))) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = factor(Age, levels = c("YNG", "OLD"))), 
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), 
              size = 0.5, alpha = 0.6) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif", 
    hide.ns = TRUE,          # Remove non-significant comparisons
    label.y = 0.07,          # Manually set the y-position
    size = 7
  ) +
  scale_fill_manual(name = "Age", values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  scale_color_manual(name = "Age", values = c("YNG" = "darkblue", "OLD" = "darkred")) +
  labs(
    title = "Normalized Cell Type Proportions by Age Group with T-test",
    y = "Proportion of Cells per Sample",
    x = NULL 
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, size = 12, face = "bold"), # Increase size and bold x-axis text
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.position = "none" 
  )

print(plot_final)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#'[Faceted by cluster]
plot_faceted <- ggplot(plot_data, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  # geom_boxplot(width=0.2, outlier.shape = NA, alpha=0.5) + # Optional boxplot
  scale_fill_manual(values = c("OLD" = "salmon", "YNG" = "lightblue")) + # Adjust colors as needed
  facet_wrap(~ ClusterID, scales = "free_y", ncol = 6) + # Adjust ncol as needed
  labs(
    title = "Normalized Cell Type Proportions by Age Group (Faceted by Cluster)",
    x = "Age Group",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 0.5),
    legend.position = "bottom", # Increase size and bold x-axis text
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
  )

print(plot_faceted)



custom_cluster_order <- c(
  "Granule (Neuron)",
  "Granule (Vegfa+)",
  "Granule (Mix)",
  "PKJ",
  "MLI (Neuron)",
  "MLI (Vegfa+)",
  "DCN"
  # "Bergmann Glia",
  # "ODC",
  # "Endothelial",
  # "Vascular",
  # "Microglia",
  # "Astrocytes",
  # "Choroid Plexus"
)

plot_data_filtered <- plot_data %>%
  filter(ClusterID %in% custom_cluster_order) %>%
  mutate(ClusterID = factor(ClusterID, levels = custom_cluster_order))



# 2. Filter the data and set the factor levels to your custom order
plot_faceted <- ggplot(plot_data_filtered, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("OLD" = "salmon", "YNG" = "lightblue")) +
  facet_wrap(~ ClusterID, scales = "free_y", ncol = 4) +
  labs(
    title = "Normalized Cell Type Proportions by Age Group",
    x = "",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # --- Line to modify facet titles ---
    strip.text = element_text(size = 12, face = "bold") # <-- Added this line
  )

# Print the final plot
print(plot_faceted)



file_name <- paste0(local_WD, "CB_Neuron_PopViolins", ".png")
ggsave(filename = file_name, plot = plot_faceted, dpi = 900, units = "in", width = 8, height = 5)



custom_cluster_order <- c(
  # "Granule (Neuron)",
  # "Granule (Vegfa+)",
  # "Granule (Mix)",
  # "PKJ",
  # "MLI (Neuron)",
  # "MLI (Vegfa+)",
  # "DCN"
  "Bergmann Glia",
  "ODC",
  "Endothelial",
  "Vascular",
  "Microglia",
  "Astrocytes",
  "Choroid Plexus"
)

plot_data_filtered <- plot_data %>%
  filter(ClusterID %in% custom_cluster_order) %>%
  mutate(ClusterID = factor(ClusterID, levels = custom_cluster_order))



# 2. Filter the data and set the factor levels to your custom order
plot_faceted <- ggplot(plot_data_filtered, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("OLD" = "salmon", "YNG" = "lightblue")) +
  facet_wrap(~ ClusterID, scales = "free_y", ncol = 4) +
  labs(
    title = "Normalized Cell Type Proportions by Age Group",
    x = "",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # --- Line to modify facet titles ---
    strip.text = element_text(size = 12, face = "bold") # <-- Added this line
  )

# Print the final plot
print(plot_faceted)



file_name <- paste0(local_WD, "CB_Glia_PopViolins", ".png")
ggsave(filename = file_name, plot = plot_faceted, dpi = 900, units = "in", width = 8, height = 5)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Sen-GSEA*
# install.packages("UCell")
library(UCell)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# senBrainList
#"Apoe",, "Vim", "Malat1"

Mouse_SenBrain_2025 = c(
  "Akt1","Apoe","App","Axl","Bcl2","Bcl2l1","Ccl2","Ccl3","Ccl4","Ccl5","Ccne1",
  "Cd274","Cd36","Cd74","Cd9","Cdc27","Cdc6","Cdk4","Cdk8","Cdk9","Cdkn1a/p21","Cdkn1b/p27kip1",
  "Cdkn1c/p57kip2","Cdkn2a/p16","Cdkn2b/p15ink4b","Cdkn2c/p18ink4c","Cdkn2d/p19","Cenpb","Csf1","Csf1r","Cxcl14","Cxcl2",
  "Cyba","Dpp4","Elmo1","Fxyd5","Gadd45a","Gas2","Gdf15","Glb1","Gpnmb","H2-Aa",
  "H2-Ab1","Hmgb1","Ifitm3","Igfbp7","Il10","Il18","Il1a","Il1b","Il23a","Il23r",
  "Il33","Il6","Irf5","Irf7","Lamp2","Lgals3","Lgals3bp","Lmna","Lpl","Malat1","Mapk1",
  "Mapk10","Mapk14","Mapk8","Meox1","Mmp12","Mmp14","Mmp2","Mmp3","Mmp9","Nfkb1",
  "Nol3","Notch3","Pak1","Pak2","Pak3","Pdcd1","Pik3ca","Pik3cb","Pik3r1","Plaur",
  "Pten","Rbl2","Serpina3n","Serpinb1a","Serpine1","Serpine2","Smad2","Smad3","Sod1","Sparc","Spp1",
  "Tgfb1","Tgfb2","Tgfbr1","Sting1","Tnf","Tp53bp1","Trp53",
  "Vim")

markers <- list(my_geneset = Mouse_SenBrain_2025)

# Calculate UCell scores for the custom gene set using the SCT assay
thisObj <- AddModuleScore_UCell(
  obj = thisObj,
  features = markers,
  assay = "SCT", 
  name = "_SenUCell" 
)

signature.names <- paste0(names(markers), "_SenUCell")




#'[Sennet quick check]

colnames(thisObj@meta.data)

DefaultAssay(thisObj) <- 'tf_activities'

FeaturePlot(thisObj, features = c("Cdkn2a/p16", "my_geneset_SenUCell"), blend = TRUE)



# Generate the UMAP with GSEA scores overlaid, faceted by age
gsea_umap_plot <- FeaturePlot(
  object = thisObj,
  features = signature.names,
  reduction = "tf_umap",
  split.by = "age", 
  raster = FALSE,
  keep.scale = "all",
  pt.size = 0.5,
  label = TRUE,
  repel = TRUE
)  +
  labs(color = "Enrichment Score") +
  patchwork::plot_annotation(
    title = "SenBrain Scores | UCell Module Scores",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )&
  theme(legend.position = "bottom")&
  scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .5)

print(gsea_umap_plot)

file_name <- paste0(local_WD, "SenScore_ageFacet_gsea_umap_plot", ".png")
ggsave(filename = file_name, plot = gsea_umap_plot, dpi = 900)

custom_col <- c(
  "Granule (Neuron)" = "blue4",
  "Granule (Vegfa+)" = "lightblue",
  "Granule (Mix)" = "purple",
  "Endothelial" = "grey30",
  "MLI (Vegfa+)" = "green",
  "PKJ" = "red",
  "Vascular" = "salmon",
  "DCN" = "red3",
  "Microglia" = "black",
  "Choroid Plexus" = "cyan", 
  "ODC" = "yellow2",
  "MLI (Neuron)" = "darkgreen",
  "Astrocytes" = "orchid",
  "Bergmann Glia" = "orange"
)

plot_df <- data.frame(
  Identity = Idents(thisObj),
  UCell_Score = thisObj$my_geneset_SenUCell
)

plot_df$Identity <- fct_reorder(plot_df$Identity, plot_df$UCell_Score, .fun = 'median', .desc = TRUE)

# 4. Generate the horizontal boxplot using ggplot2
gsea_boxplot <- ggplot(plot_df, aes(x = UCell_Score, y = Identity, fill = Identity)) +
  geom_boxplot(outlier.shape = NA) + # Hides outlier points for a cleaner look
  # geom_jitter(height = 0.2, width = 0, size = 0.5, alpha = 0.2) + # Optional: add points
  scale_fill_manual(values = custom_col) +
  labs(
    title = "SenBrain Score Enrichment by Cell Type",
    x = "UCell Enrichment Score",
    y = "Cluster Identity"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# 5. Print the plot
print(gsea_boxplot)










library(rstatix)
data_for_stats <- data.frame(
  Identity = Idents(thisObj),
  UCell_Score = thisObj$my_geneset_SenUCell,
  age = thisObj$age 
)

# 3. Perform statistical analysis
# For each cluster (Identity), perform a t-test on UCell scores between age groups.
# The 'ref.group' is the baseline for comparison.
stat_results <- data_for_stats %>%
  group_by(Identity) %>%
  t_test(UCell_Score ~ age, ref.group = "YNG") %>%
  add_significance("p") # Adds a 'p.signif' column (e.g., ****, ***, ns)

# 4. Calculate Fold Change
# First, calculate the mean UCell score for each age group within each cluster.
# Then, calculate the fold change (OLD / YNG).
fc_results <- data_for_stats %>%
  group_by(Identity, age) %>%
  summarise(mean_score = mean(UCell_Score), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = age, values_from = mean_score) %>%
  mutate(fold_change = OLD / YNG)

# 5. Combine stats and fold change data for plotting
plot_data <- left_join(fc_results, stat_results, by = "Identity")

# 6. Generate the ordered bar plot
fold_change_plot <- ggplot(
  # Reorder 'Identity' based on fold_change value, so highest is on top
  plot_data, aes(x = fold_change, y = fct_reorder(Identity, fold_change), fill = Identity)
) +
  geom_col() + # Use geom_col for bar plot with pre-calculated values
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") + # Line for no change
  # Add the significance stars to the right of the bars
  geom_text(aes(label = p.signif), hjust = -0.2, size = 6, color = "black") +
  scale_fill_manual(values = custom_col) +
  labs(
    title = "Fold Change in SenBrain Score (OLD vs. YNG)",
    x = "Fold Change (Mean UCell Score OLD / YNG)",
    y = "Cluster Identity"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  # Add extra space on the right for the significance labels
  coord_cartesian(clip = 'off') + # Prevents labels from being cut off
  scale_x_continuous(expand = expansion(mult = c(0, .1)))


# 7. Print the plot
print(fold_change_plot)












data_for_stats <- data.frame(
  Identity = Idents(thisObj),
  UCell_Score = thisObj$my_geneset_SenUCell,
  age = thisObj$age 
)

# 2. Pre-calculate statistics and their positions for plotting
stat_results <- data_for_stats %>%
  group_by(Identity) %>%
  t_test(UCell_Score ~ age, ref.group = "YNG") %>%
  add_significance("p")

stat_positions <- data_for_stats %>%
  group_by(Identity) %>%
  summarise(max_score = max(UCell_Score, na.rm = TRUE))

stat_labels <- stat_results %>%
  left_join(stat_positions, by = "Identity") %>%
  mutate(y.position = max_score * 1.05) 

# 3. Determine Cluster Order
# For horizontal plots, arrange() ascending to put the highest score at the top.
cluster_order <- data_for_stats %>%
  filter(age == "OLD") %>%
  group_by(Identity) %>%
  summarise(median_score = median(UCell_Score, na.rm = TRUE)) %>%
  arrange(median_score) %>% 
  pull(Identity)

# 4. Apply the factor order to ALL relevant data frames
data_for_stats$Identity <- factor(data_for_stats$Identity, levels = cluster_order)
stat_labels$Identity <- factor(stat_labels$Identity, levels = cluster_order)

# 5. Generate the final horizontal boxplot
horizontal_boxplot <- ggplot(data_for_stats, aes(x = Identity, y = UCell_Score, fill = age)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8)) +
  
  geom_text(
    data = stat_labels, 
    aes(x = Identity, y = y.position, label = p.signif), 
    inherit.aes = FALSE,
    size = 5
  ) +
  
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  
  labs(
    title = "SenBrain Score Distribution by Age and Cell Type",
    y = "UCell Enrichment Score",
    x = NULL,
    fill = "Age"
  ) +
  
  # Set axis limits AND flip the coordinates
  coord_cartesian(ylim = c(0.3, 0.5)) +
  coord_flip() +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    # Style the y-axis text (cluster names) instead of the x-axis
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 14),
    # Remove all grid lines
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 6. Print the final plot
print(horizontal_boxplot)









stat_position_manual <- 0.41

data_for_stats <- data.frame(
  Identity = Idents(thisObj),
  UCell_Score = thisObj$my_geneset_SenUCell,
  age = thisObj$age 
)

# 2. Pre-calculate statistics and set their fixed position
# The dynamic calculation is replaced with our manual value.
stat_results <- data_for_stats %>%
  group_by(Identity) %>%
  t_test(UCell_Score ~ age, ref.group = "YNG") %>%
  add_significance("p")

stat_labels <- stat_results %>%
  mutate(y.position = stat_position_manual) # Use the manual position here

# 3. Determine Cluster Order
cluster_order <- data_for_stats %>%
  filter(age == "OLD") %>%
  group_by(Identity) %>%
  summarise(median_score = median(UCell_Score, na.rm = TRUE)) %>%
  arrange(median_score) %>% 
  pull(Identity)

# 4. Apply the factor order to ALL relevant data frames
data_for_stats$Identity <- factor(data_for_stats$Identity, levels = cluster_order)
stat_labels$Identity <- factor(stat_labels$Identity, levels = cluster_order)

# 5. Generate the final horizontal boxplot
horizontal_boxplot <- ggplot(data_for_stats, aes(x = Identity, y = UCell_Score, fill = age)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8)) +
  
  geom_text(
    data = stat_labels, 
    aes(x = Identity, y = y.position, label = p.signif), 
    inherit.aes = FALSE,
    size = 5
  ) +
  
  scale_fill_manual(values = c("YNG" = "lightblue", "OLD" = "salmon")) +
  
  labs(
    title = "",
    y = "UCell Enrichment Score",
    x = NULL,
    fill = "Age"
  ) +
  
  # Set axis limits AND flip the coordinates
  # The upper limit is now slightly larger to make room for the stat labels.
  coord_flip(ylim = c(0.40, 0.55)) +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 14),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# 6. Print the final plot
print(horizontal_boxplot)


file_name <- paste0(local_WD, "CB_SenBrainClusterBoxplots", ".png")
ggsave(filename = file_name, plot = horizontal_boxplot, dpi = 900, units = "in")#, width = 8, height = 5)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Sub-cluster Analyses*
levels(thisObj)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
SubClusterCellTyping_Singles <- function(obj, cluster_label, 
                                         resolution = 0.1, k.param = k.param, 
                                         spread = 1, min.dist = 0.1) {
  # Subset the object based on the provided cluster label
  subset_obj <- subset(obj, idents = cluster_label)
  subset_obj <- ScaleData(subset_obj, features = rownames(subset_obj))
  subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", verbose = TRUE)
  
  # Run PCA and display the elbow plot for user inspection
  subset_obj <- RunPCA(subset_obj, features = VariableFeatures(object = subset_obj), npcs = 50)
  E <- ElbowPlot(subset_obj) + ggtitle(paste("Elbow Plot of", cluster_label))
  print(E)
  
  # Loop until the user confirms that the UMAP clustering is sufficient
  repeat {
    # Request number of dimensions from the user (converted to numeric)
    dims_to_use <- as.numeric(readline(prompt = "Enter number of dims to consider: "))
    
    # Run neighbor finding, clustering, and UMAP embedding using the chosen dims
    obj_iter <- FindNeighbors(subset_obj, dims = 1:dims_to_use, k.param = k.param)
    obj_iter <- FindClusters(obj_iter, resolution = resolution)
    obj_iter <- RunUMAP(obj_iter, dims = 1:dims_to_use, min.dist = min.dist, spread = spread)
    
    # Identify all markers per cluster (no filtering)
    # all_genes <- FindAllMarkers(obj_iter, only.pos = TRUE)
    # print(head(all_genes))
    # Extract the top 10 marker genes per cluster for visualization
    # top_genes <- all_genes %>%
    #   group_by(cluster) %>%
    #   slice_max(order_by = avg_log2FC, n = 10)
    # 
    # # Ensure proper ordering and factor levels for plotting
    # top_genes <- top_genes[order(top_genes$cluster), ]
    # top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))
    # 
    # Generate and display the UMAP plot
    UmapPlot <- DimPlot(obj_iter, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) +
      ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))
    
    # # Create and display the dot plot for top marker genes per cluster
    # gene_plot <- ggplot(top_genes, aes(x = factor(cluster), y = gene, color = avg_log2FC, size = pct.1)) +
    #   geom_point() +
    #   scale_color_viridis_c() +
    #   labs(title = paste("Top 10 Genes per Sub-Cluster | Global Cluster:", cluster_label),
    #        x = "Cluster",
    #        y = "Gene",
    #        color = "Log2 Fold Change",
    #        size = "Percentage Expression") +
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # print(gene_plot)
    print(UmapPlot)
    # Ask the user if the UMAP clustering is sufficient with explicit validation
    valid_input <- FALSE
    while (!valid_input) {
      answer <- readline(prompt = "Is the UMAP clustering sufficient? (Y/YES or N/NO): ")
      answer <- toupper(trimws(answer))
      if (answer %in% c("Y", "YES", "N", "NO")) {
        valid_input <- TRUE
      } else {
        cat("Invalid input. Please enter Y/YES or N/NO.\n")
      }
    }
    
    if (answer %in% c("Y", "YES")) {
      break  # Exit the loop if the clustering is sufficient
    } else if (answer %in% c("N", "NO")) {
      # Reprint the elbow plot so the user doesn't have to scroll back
      print(E)
    }
  }
  
  # Return a list of key outputs for downstream analysis, including all genes
  return(list(
    subset_obj = obj_iter,
    # all_genes = all_genes,  # Full list of DE genes
    # top_genes = top_genes,  # Top 10 genes per cluster
    UMAP_plot = UmapPlot,
    # gene_plot = gene_plot,
    elbow_plot = E
  ))
}

library('presto')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Neuron subcluster]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster_label = "Granule (Neuron)"
SubC_GranuleNeuron_Data <- SubClusterCellTyping_Singles(
  obj = thisObj,
  cluster_label = cluster_label,
  resolution = 0.05,
  k.param = 15,
  spread = 1,
  min.dist = 0.001
)
# Dims used: 15
top_gene_list <- SubC_GranuleNeuron_Data$all_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  summarise(all_genes = list(gene)) %>%
  pull(all_genes)

print(top_gene_list)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'[Microglia subcluster]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster_label = "Microglia"

obj <- thisObj
# Subset the object based on the provided cluster label
subset_obj <- subset(obj, idents = cluster_label)
subset_obj <- ScaleData(subset_obj, features = rownames(subset_obj))
subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", verbose = TRUE)
  
# Run PCA and display the elbow plot for user inspection
subset_obj <- RunPCA(subset_obj, features = VariableFeatures(object = subset_obj), npcs = 50)
E <- ElbowPlot(subset_obj) + ggtitle(paste("Elbow Plot of", cluster_label))
print(E)
  
dims_to_use <- 7
    
# Run neighbor finding, clustering, and UMAP embedding using the chosen dims
obj_iter <- FindNeighbors(subset_obj, dims = 1:dims_to_use, k.param = 10)
obj_iter <- FindClusters(obj_iter, resolution = 1)
obj_iter <- RunUMAP(obj_iter, dims = 1:dims_to_use, min.dist = 0.1, spread = 0.1)

UmapPlot <- DimPlot(obj_iter, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) +
      ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

 print(UmapPlot)
 
 
 
 my_geneset = c("App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                "Ccne1", "CD26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
                "Cdk4", "Cdk8", "Cdk9", "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19",
                "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
                "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
                "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
                "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
                "Lmna", "Lpl", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
                "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
                "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
                "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vtn", "Xcl1")
 markers <- list(my_geneset = my_geneset)
 
 # Calculate UCell scores for the custom gene set using the SCT assay
 obj_iter <- AddModuleScore_UCell(
   obj = obj_iter,
   features = markers,
   assay = "SCT", 
   name = "_SenUCell" 
 )
 
 signature.names <- paste0(names(markers), "_SenUCell")
 
 obj_iter$age <- factor(x = obj_iter$age, levels = c("YNG", "OLD"))
 
 gsea_umap_plot <- FeaturePlot(
   object = obj_iter,
   features = signature.names,
   reduction = "umap",
   split.by = 'age',
   raster = FALSE,
   keep.scale = "all",
   pt.size = 2,
   label = FALSE,
   repel = TRUE
 )  +
   labs(color = "Enrichment Score") +
   patchwork::plot_annotation(
     title = "",
     theme = theme(plot.title = element_text(hjust = 0, size = 16))
   )&
   theme(legend.position = "none",
         axis.title = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank())&
   scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .5)
 
 print(gsea_umap_plot)
 
 file_name <- paste0(local_WD, "SenScore_microgliaSubC_umap", ".png")
 ggsave(filename = file_name, plot = gsea_umap_plot, dpi = 900)
    
    
    
    
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #'[Granule (Mix) subcluster]
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cluster_label = "Granule (Mix)"
 
 obj <- thisObj
 # Subset the object based on the provided cluster label
 subset_obj <- subset(obj, idents = cluster_label)
 subset_obj <- ScaleData(subset_obj, features = rownames(subset_obj))
 subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", verbose = TRUE)
 
 # Run PCA and display the elbow plot for user inspection
 subset_obj <- RunPCA(subset_obj, features = VariableFeatures(object = subset_obj), npcs = 50)
 E <- ElbowPlot(subset_obj) + ggtitle(paste("Elbow Plot of", cluster_label))
 print(E)
 
 dims_to_use <- 15
 
 # Run neighbor finding, clustering, and UMAP embedding using the chosen dims
 obj_iter <- FindNeighbors(subset_obj, dims = 1:dims_to_use, k.param = 10)
 obj_iter <- FindClusters(obj_iter, resolution = 1)
 obj_iter <- RunUMAP(obj_iter, dims = 1:dims_to_use, min.dist = 0.1, spread = 0.1)
 
 UmapPlot <- DimPlot(obj_iter, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) +
   ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))
 
 print(UmapPlot)
 
 
 
 my_geneset = c("App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                "Ccne1", "CD26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
                "Cdk4", "Cdk8", "Cdk9", "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19",
                "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
                "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
                "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
                "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
                "Lmna", "Lpl", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
                "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
                "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
                "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vtn", "Xcl1")
 markers <- list(my_geneset = my_geneset)
 
 # Calculate UCell scores for the custom gene set using the SCT assay
 obj_iter <- AddModuleScore_UCell(
   obj = obj_iter,
   features = markers,
   assay = "SCT", 
   name = "_SenUCell" 
 )
 
 signature.names <- paste0(names(markers), "_SenUCell")
 
 obj_iter$age <- factor(x = obj_iter$age, levels = c("YNG", "OLD"))
 
 gsea_umap_plot <- FeaturePlot(
   object = obj_iter,
   features = signature.names,
   reduction = "umap",
   split.by = 'age',
   raster = FALSE,
   keep.scale = "all",
   pt.size = 1,
   label = TRUE,
   repel = TRUE
 )  +
   labs(color = "Enrichment Score") +
   patchwork::plot_annotation(
     title = "",
     theme = theme(plot.title = element_text(hjust = 0, size = 16))
   )&
   theme(legend.position = "none",
         axis.title = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank())&
   scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .5)
 
 print(gsea_umap_plot)
 
 file_name <- paste0(local_WD, "SenScore_GranuleMix_umap", ".png")
 ggsave(filename = file_name, plot = gsea_umap_plot, dpi = 900, units = "in", width = 6.01, height = 4.77)
 
  

 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #'[Astrocyte subcluster]
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cluster_label = "Astrocytes"
 
 obj <- thisObj
 # Subset the object based on the provided cluster label
 subset_obj <- subset(obj, idents = cluster_label)
 subset_obj <- ScaleData(subset_obj, features = rownames(subset_obj))
 subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", verbose = TRUE)
 
 # Run PCA and display the elbow plot for user inspection
 subset_obj <- RunPCA(subset_obj, features = VariableFeatures(object = subset_obj), npcs = 50)
 E <- ElbowPlot(subset_obj) + ggtitle(paste("Elbow Plot of", cluster_label))
 print(E)
 
 dims_to_use <- 10
 
 # Run neighbor finding, clustering, and UMAP embedding using the chosen dims
 obj_iter <- FindNeighbors(subset_obj, dims = 1:dims_to_use, k.param = 10)
 obj_iter <- FindClusters(obj_iter, resolution = 1)
 obj_iter <- RunUMAP(obj_iter, dims = 1:dims_to_use, min.dist = 0.1, spread = 0.1)
 
 UmapPlot <- DimPlot(obj_iter, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) +
   ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))
 
 print(UmapPlot)
 
 
 
 my_geneset = c("App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                "Ccne1", "CD26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
                "Cdk4", "Cdk8", "Cdk9", "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19",
                "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
                "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
                "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
                "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
                "Lmna", "Lpl", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
                "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
                "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
                "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vtn", "Xcl1")
 markers <- list(my_geneset = my_geneset)
 
 # Calculate UCell scores for the custom gene set using the SCT assay
 obj_iter <- AddModuleScore_UCell(
   obj = obj_iter,
   features = markers,
   assay = "SCT", 
   name = "_SenUCell" 
 )
 
 signature.names <- paste0(names(markers), "_SenUCell")
 
 obj_iter$age <- factor(x = obj_iter$age, levels = c("YNG", "OLD"))
 
 gsea_umap_plot <- FeaturePlot(
   object = obj_iter,
   features = signature.names,
   reduction = "umap",
   split.by = 'age',
   raster = FALSE,
   keep.scale = "all",
   pt.size = 2,
   label = TRUE,
   repel = TRUE
 )  +
   labs(color = "Enrichment Score") +
   patchwork::plot_annotation(
     title = "",
     theme = theme(plot.title = element_text(hjust = 0, size = 16))
   )&
   theme(legend.position = "bottom",
         axis.title = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank())&
   scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .5)
 
 print(gsea_umap_plot)
 
 file_name <- paste0(local_WD, "legend", ".png")
 ggsave(filename = file_name, plot = gsea_umap_plot, dpi = 900, units = "in", width = 15, height = 20)
 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Cellchat attempt*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
 library(CellChat)
 thisObj$samples <- as.factor(thisObj$mus)
 cellchat <- createCellChat(object = thisObj,
                            coordinates = c(x_slide_micron, y_slide_micron),
                            group.by = "Global_CellTypes",
                            assay = "SCT")

 CellChatDB <- CellChatDB.mouse
 cellchat@DB <- CellChatDB
 head(CellChatDB$interaction)

 cellchat <- subsetData(cellchat)

 cellchat <- identifyOverExpressedGenes(cellchat)
 cellchat <- identifyOverExpressedInteractions(cellchat)

 cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                               distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                               contact.dependent = TRUE, contact.range = 10)
 
 cellchat <- filterCommunication(cellchat, min.cells = 10)
 cellchat <- computeCommunProbPathway(cellchat)
 cellchat <- aggregateNet(cellchat)

 groupSize <- as.numeric(table(cellchat@idents))
 
 # First, let's visualize the number of interactions between cell types
 Interaction_num <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                  label.edge= FALSE, title.name = "Number of Interactions", vertex.label.cex = 1.5, top = 0.2,
                  margin = 0.3)
 
 # Now, let's visualize the strength of those interactions
 Interaction_strength <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                  label.edge= FALSE, title.name = "Interaction Strength", vertex.label.cex = 1.5, top = 0.2,
                  margin = 0.3)

 
 file_name <- paste0(local_WD, "Interaction_num", ".png")
 png(file = file_name, width = 10, height = 10, units = "in", res = 900)
 Interaction_num
 dev.off() 
 
 file_name <- paste0(local_WD, "Interaction_strength", ".png")
 png(file = file_name, width = 10, height = 10, units = "in", res = 900)
 Interaction_strength
 dev.off() 
 
 
 netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues",
                   cluster.rows = TRUE,
                   cluster.cols = TRUE)
 netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues",
                   cluster.rows = TRUE,
                   cluster.cols = TRUE)
 
 
 netVisual_bubble(cellchat, targets.use = "Microglia", remove.isolate = FALSE)
 netVisual_chord_gene(cellchat, sources.use = "Microglia", signaling = "CSF")
 
 cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
 
 pathways.show <- c("APP") 
 APP_heat <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
 APP_heat
 
 pathways.show <- c("VEGF") 
 VEGF_heat <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
 VEGF_heat
 
 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Bisque Celltype Deconvolution*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 library(BisqueRNA)
 
 sc_counts <- SummarizedExperiment::assay(sce, "counts")
 sc_pheno <- as.data.frame(SummarizedExperiment::colData(sce))
 
 # Bisque needs a column specifically identifying the cell types. 
 # Let's ensure your "Global_CellTypes" are in a column named "label".
 sc_pheno$label <- sc_pheno$Global_CellTypes
 
 # Create the single-cell ExpressionSet object
 sc_eset <- Biobase::ExpressionSet(
   assayData = as.matrix(sc_counts),
   phenoData = new("AnnotatedDataFrame", data = sc_pheno)
 )
 
 
 # Assuming you have a matrix 'bulk_counts' and a data frame 'bulk_metadata'
 # Make sure the row names of 'bulk_metadata' match the column names of 'bulk_counts'

 bulk_eset <- Biobase::ExpressionSet(assayData = as.matrix(CB_deseq$counts),
                                     phenoData = new("AnnotatedDataFrame", data = as.data.frame(colData(CB_deseq$dds))))
 
 
 deconvolution_results <- BisqueRNA::ReferenceBasedDecomposition(
   bulk.eset = bulk_eset,
   sc.eset = sc_eset,
   cell.types = "Global_CellTypes",
   subject.names = "mus",      # The correct argument is 'subject.names'
   use.overlap = FALSE
 )
 
 # The rest of the process is the same
 cell_proportions <- deconvolution_results$bulk.props
 
 as.data.frame(CB_deseq$dds@metadata)
 
 
 
 
 microglia_proportion <- cell_proportions["Microglia", ]
 meta <- as.data.frame(colData(CB_deseq$dds))
 age_numeric <- as.numeric(as.character(age))
 
 # Now, run the correlation test using the new numeric vector
 cor_test <- cor.test(microglia_proportion, age_numeric, method = "pearson")
 
 # Print the results
 print(cor_test)
 
 # png("microglia_vs_age_plot.png", width = 800, height = 600, res = 100)
 # 
 plot_data <- data.frame(
   Age = age_numeric,
   Proportion = microglia_proportion
 )
 
 # For a boxplot, convert the numeric 'Age' to a categorical factor
 # Specifying 'levels' ensures they are ordered correctly on the x-axis
 plot_data$Age <- factor(plot_data$Age, levels = c(4, 12, 18, 24, 28))
 
 # 2. Build the ggplot object
 p_boxplot <- ggplot(plot_data, aes(x = Age, y = Proportion)) +
   geom_boxplot() +  # Create the boxplot
   geom_jitter(width = 0.2, alpha = 0.6) + # Optional: overlays the actual data points
   labs(
     title = "",
     x = "Age",
     y = "Estimated Microglia Enrichment"
   ) +
   theme_minimal()+ 
   theme(
     # 1. Remove gridlines and plot title
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     plot.title = element_blank(),
     
     # 2. Remove the x-axis title
     axis.title.x = element_blank(),
     
     # 3. Increase y-axis title font size and make it bold
     axis.title.y = element_text(size = 14, face = "bold"),
     
     # 4. (Optional but recommended) Increase axis text font size for readability
     axis.text = element_text(size = 14, face = "bold", color = "black")
   )
 
 # 3. Print the plot to view it
 print(p_boxplot)
 
 file_name <- paste0(local_WD, "MG_Deconvolution", ".png")
 ggsave(filename = file_name, plot = p_boxplot, dpi = 900, units = "in", width = 4, height = 6)
 
 
 
 
 
 
 

 
 proportions_long <- cell_proportions %>%
   as.data.frame() %>%
   rownames_to_column(var = "CellType") %>%
   pivot_longer(
     cols = -CellType,
     names_to = "SampleID",
     values_to = "Proportion"
   )
 
 # b. Prepare the metadata (assuming 'bulk_metadata' has sample IDs and ages)
 #    Make sure the sample IDs are in a column to join on.
 metadata_df <- meta %>%
   rownames_to_column(var = "SampleID")
 
 # c. Join proportions with metadata to get the age for each sample
 plot_data_full <- left_join(proportions_long, metadata_df, by = "SampleID")
 
 # d. Calculate the mean proportion for each cell type at each age
 summary_data <- plot_data_full %>%
   mutate(Age = as.numeric(as.character(Age))) %>% # Ensure age is numeric
   group_by(CellType, Age) %>%
   summarise(Mean_Proportion = mean(Proportion), .groups = 'drop')
 
 
 cell_types_to_plot <- c("Microglia", "Astrocytes", "PKJ", "Bergmann Glia")
 summary_data_subset <- summary_data %>%
   filter(CellType %in% cell_types_to_plot)
 
 
 # 2. Build the ggplot line graph
 p_line <- ggplot(summary_data, aes(x = Age, y = Mean_Proportion, group = CellType, color = CellType)) +
   geom_line(linewidth = 1) +  # Create the lines
   geom_point(size = 2.5) +     # Add points at each timepoint
   scale_color_manual(values = custom_col) + # Apply your custom color palette
   labs(
     title = "Mean Cell Type Proportion vs. Age",
     x = "Age",
     y = "Mean Estimated Proportion",
     color = "Cell Type" # Renames the legend title
   ) +
   theme_minimal()
 
 # 3. Print the plot to view it
 print(p_line)
 
#catch
#**OLD CODE* ####


gene_marker_lists <- list(
  # 1. Neurons (Excitatory: Granule + UBCs)
  Pan_Neuronal = c("Neurod1", "Neurod2", "Neurod6", "Mecp2", "Rest"), # Key neuronal development and identity TFs
  Granule = c("Pax6", "Zeb1", "Neurod1"), # Pax6 and Neurod1 are crucial for granule cell development
  
  # 2. Interneurons (Inhibitory)
  Neuron_Inhibitory = c("Dlx2", "Dlx5", "Ascl1", "Hes1"), # Dlx family are master regulators of GABAergic neurons
  
  # 3. Purkinje Cells
  Purkinje = c("Rora", "Pou4f2", "Lmx1b", "Tbx2"), # Rora is a classic and essential Purkinje cell TF
  
  # 4. Glia (Astrocytes, Bergmann Glia)
  Astrocyte_Glia = c("Sox9", "Nfia", "Nfib", "Stat3", "Hes1", "Hesx1"), # Sox9 and Nfi family are key glial regulators
  
  # 5. Fibroblasts / Mesenchymal
  Fibroblast_VLMC = c("Foxc2", "Twist1", "Twist2", "Tcf7l2"),
  
  # 6. Microglia
  Microglia = c("Spi1", "Irf8", "Irf5", "Irf9", "Runx1", "Stat6"), # Spi1 (PU.1) and the IRF family are master regulators of microglia
  
  # 7. Oligodendrocytes & OPCs
  Oligodendrocyte_Lineage = c("Olig1", "Sox10", "Tcf7l2", "Nkx6-1"), # Olig1/2 and Sox10 are master regulators of this lineage
  
  # 8. Endothelial & Vascular Cells (Endothelial, Pericyte, VSMC)
  Vascular = c("Erg", "Fli1", "Foxc2", "Sox17", "Etv6") # ERG and FLI1 are key endothelial TFs
)





specific_marker_list <- list(
  # 1. Neuronal Cell Types
  # --- Excitatory Neurons ---
  Granule_Cell = c("Atp2b4", "Reln", "Gap43", "Grin2c"),      # Alternatives for targeted panels. Use combinatorially.
  UBC = c("Eomes", "Calb2", "Grm1", "Plcb4"),                  # Eomes (Tbr2) is the most specific marker.
  
  # --- Inhibitory Neurons ---
  Purkinje_Cell = c("Calb1", "Pcp2", "Car8", "Rora"),           # Pcp2 (L7) and Calb1 are canonical.
  Golgi_Cell = c("Gcgr", "Nrgn", "Gdnf", "Mapt"),               # Gcgr is a reliable Golgi marker.
  Basket_Cell = c("Reln", "Sln", "Pvalb"),                     # Enriched in basket cells vs. stellate cells.
  Stellate_Cell = c("Ntf3", "Penk", "Pvalb"),                  # Enriched in stellate cells. Pvalb is lower than in basket.
  Lugaro_Cell = c("Calb2", "Sox14", "Plcl1", "Vwc2"),           # Rare; Calb2 shared with UBCs.
  
  # 2. Glial Cell Types
  Astrocyte = c("Aldh1l1", "Gfap", "Aqp4", "S100b"),            # Aldh1l1 is a robust pan-astrocyte marker.
  Bergmann_Glia = c("Hopx", "Gdf10", "Slc1a3", "Sox9"),         # Hopx is highly specific to Bergmann glia.
  OPC = c("Pdgfra", "Cspg4", "Olig2", "Tcf7l2"),               # Pdgfra & Cspg4 (NG2) are canonical.
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Klk6", "Myrf"),    # Markers for mature, myelinating cells.
  Microglia = c("Tmem119", "P2ry12", "Cx3cr1", "Csf1r"),        # Tmem119/P2ry12 distinguish from macrophages.
  
  # 3. Vascular & Meningeal Cell Types
  Endothelial_Cell = c("Cldn5", "Pecam1", "Esam", "Flt1"),      # Cldn5 is a key BBB marker.
  Pericyte = c("Pdgfrb", "Rgs5", "Kcnj8", "Vtn"),               # Pdgfrb is more specific than Acta2.
  VSMC = c("Acta2", "Tagln", "Notch3", "Myh11"),                # Tagln is highly specific to smooth muscle.
  Fibroblast_VLMC = c("Col1a1", "Dcn", "Lum", "Slc6a13"),       # Collagen-producing meningeal/vascular cells.
  Choroid_Plexus = c("Ttr", "Folr1", "Aqp1", "Otx2")            # Ttr is the most definitive marker.
)
marker_lists <- specific_marker_list 

# 1. Get all unique TFs that are present in the Seurat object's TF activity assay
all_markers <- unique(unlist(marker_lists))
markers_present <- intersect(all_markers, rownames(thisObj[["SCT"]]))

# Check if we have markers to proceed
if (length(markers_present) == 0) {
  stop("None of the specified TFs were found in the 'SCT' assay.")
}

print(markers_present)
# 2. Calculate the average TF activity for these markers across your clusters
# The result is a matrix of TFs (rows) vs. clusters (columns)
avg_exp_matrix <- AverageExpression(
  thisObj,
  features = markers_present,
  group.by = "seurat_clusters_merged", # Or whichever column holds your cluster IDs
  assay = "SCT",
  slot = "data"
)$SCT # Directly extract the matrix from the list

# 3. Calculate the enrichment score for each TF list
# This loop calculates the mean activity score for each TF list in each cluster.
# It iterates through your list of cell types (e.g., Granule_Cell, Purkinje_Cell).
enrichment_scores <- list()
for (cell_type in names(marker_lists)) {
  
  # Find which TFs from the current list are in our matrix
  tfs_in_list <- intersect(marker_lists[[cell_type]], rownames(avg_exp_matrix))
  
  if (length(tfs_in_list) > 0) {
    # Calculate the mean score for this list of TFs across all clusters
    # The 'colMeans' function averages the scores for the selected TFs for each cluster
    enrichment_scores[[cell_type]] <- colMeans(avg_exp_matrix[tfs_in_list, , drop = FALSE])
  }
}

# 4. Combine the scores into a single matrix for plotting
# The result is a matrix of cell types/lists (rows) vs. clusters (columns)
enrichment_matrix <- do.call(rbind, enrichment_scores)

# At this point, 'enrichment_matrix' should be a valid matrix, not NULL.

pheatmap::pheatmap(
  enrichment_matrix,
  main = "TF List Enrichment per Cluster",
  scale = "row",         # Scale scores by row (cell type) to see which cluster is highest for each list
  cluster_cols = TRUE,   # Cluster columns (cell clusters) to group similar clusters
  cluster_rows = TRUE,   # Cluster rows (cell types) to group similar TF lists
  border_color = "white"
)



tf_marker_list <- list(
Granule_Cell = c("Pax6", "Neurod1", "Zic1", "Zic2"),      # Pax6/Neurod1 are master regulators
UBC = c("Eomes"),                                        # Eomes (Tbr2) is the single most specific TF
Purkinje_Cell = c("Rora", "Foxp2", "Pou4f2"),              # Rora is the canonical Purkinje TF
Golgi_Cell = c("Etv1", "Neurod6", "Pax2"),                 # Etv1 is a key specific marker
MLI = c("Pax2", "Lhx1", "Lhx5"),                         # Pax2 defines basket/stellate interneurons

# 2. Glial Transcription Factors
Astrocyte = c("Sox9", "Nfia", "Nfib"),                     # Nfi family and Sox9 are master regulators
Bergmann_Glia = c("Hopx", "Sox9"),                       # Hopx is highly specific to this subtype
OPC = c("Olig2", "Sox10", "Ascl1"),                      # Olig2/Sox10 are foundational for the lineage
Oligodendrocyte = c("Myrf", "Olig1", "Sox10"),           # Myrf is the terminal TF for myelination
Microglia = c("Spi1", "Irf8", "Runx1"),                  # Spi1 (PU.1) and Irf8 are master regulators

# 3. Vascular & Mesenchymal Transcription Factors
Endothelial_Cell = c("Erg", "Fli1", "Sox17"),            # Erg/Fli1 are key ETS-family TFs
Pericyte = c("Foxf2", "Tbx18", "Meis2"),                 # Foxf2 is strongly associated with pericytes
Fibroblast_VLMC = c("Foxc1", "Twist1", "Tcf21")           # TFs for mesenchymal identity
)


gene_list <- unique(unlist(tf_marker_list))
gene_marker_lists <- gene_list

library(UCell)

# Use the filtered marker list created previously
marker_lists_to_use <- gene_marker_lists # Or use cerebellum_marker_lists if preferred

# Get all unique genes from the lists
all_markers_filtered <- unique(unlist(marker_lists_to_use))

# Ensure markers exist in the assay to be used (e.g., tf_activities)
all_markers_present <- intersect(all_markers_filtered, rownames(thisObj[["tf_activities"]]))
if (length(all_markers_present) == 0) {
  stop("None of the marker genes found in the tf_activities assay of the Seurat object.")
}
print(paste("Analyzing enrichment for", length(all_markers_present), "marker genes present in the tf_activities assay."))


# 2. Calculate Average Expression per Cluster
# Using the tf_activities assay data slot is recommended since tf_activities was used for normalization/PCA
avg_expression <- AverageExpression(
  thisObj,
  features = all_markers_present,
  group.by = "ident", # Use the active cluster identity
  assay = "tf_activities",      
  slot = "data"       # Use normalized data
)

if (is.list(avg_expression)) {
  avg_expression_matrix <- avg_expression[["tf_activities"]] # Maybe "tf_activities" wasn't found?
} else {
  avg_expression_matrix <- avg_expression
}
# Error occurs here if avg_expression_matrix wasn't successfully created above
scaled_avg_expression <- t(scale(t(avg_expression_matrix)))
# Replace NaN values that can result from scaling genes with 0 variance across clusters
scaled_avg_expression[is.nan(scaled_avg_expression)] <- 0

cluster_enrichment_scores <- list()

for (list_name in names(marker_lists_to_use)) {
  gene_list <- marker_lists_to_use[[list_name]]
  markers_in_avg_exp <- intersect(gene_list, rownames(scaled_avg_expression))
  
  if (length(markers_in_avg_exp) > 0) {
    # Get scaled average expression for these markers
    list_scaled_avg_exp <- scaled_avg_expression[markers_in_avg_exp, , drop = FALSE]
    # Calculate mean scaled score per cluster for this list
    scores <- colMeans(list_scaled_avg_exp, na.rm = TRUE)
    cluster_enrichment_scores[[list_name]] <- scores
  } else {
    # Handle case where no markers from the list were found/calculated
    num_clusters <- ncol(scaled_avg_expression)
    cluster_enrichment_scores[[list_name]] <- rep(0, num_clusters)
    names(cluster_enrichment_scores[[list_name]]) <- colnames(scaled_avg_expression)
    warning(paste("No markers from list", list_name, "found in average expression matrix."))
  }
}

# 4. Combine into a Matrix
cluster_enrichment_matrix <- do.call(rbind, cluster_enrichment_scores)
# Column names should already be cluster IDs from AverageExpression

# 5. Visualize with pheatmap (Similar logic as with PCs)

print("Generating Heatmap: Relative Gene List Enrichment per Cluster (Scaled by Column)")
UMAP_Cluster_EnrichHeat_ColScale <- pheatmap(avg_expression_matrix,
                                             scale = "row", # Which list(s) are prominent within each cluster?
                                             main = "Relative Gene List Enrichment per Cluster (Scaled by Column)",
                                             cluster_rows = TRUE, # Cluster rows (gene lists) by similarity
                                             cluster_cols = FALSE) # Cluster columns (clusters) by similarity

file_name <- paste0(local_WD, "UMAP_Cluster_EnrichHeat_ColScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
UMAP_Cluster_EnrichHeat_ColScale
dev.off() 













tf_marker_list <- list(
  Granule_Cell = c("Pax6", "Neurod1", "Zic1", "Zic2"),      # Pax6/Neurod1 are master regulators
  UBC = c("Eomes"),                                        # Eomes (Tbr2) is the single most specific TF
  Purkinje_Cell = c("Rora", "Foxp2", "Pou4f2"),              # Rora is the canonical Purkinje TF
  Golgi_Cell = c("Etv1", "Neurod6", "Pax2"),                 # Etv1 is a key specific marker
  MLI = c("Pax2", "Lhx1", "Lhx5"),                         # Pax2 defines basket/stellate interneurons
  
  # 2. Glial Transcription Factors
  Astrocyte = c("Sox9", "Nfia", "Nfib"),                     # Nfi family and Sox9 are master regulators
  Bergmann_Glia = c("Hopx", "Sox9"),                       # Hopx is highly specific to this subtype
  OPC = c("Olig2", "Sox10", "Ascl1"),                      # Olig2/Sox10 are foundational for the lineage
  Oligodendrocyte = c("Myrf", "Olig1", "Sox10"),           # Myrf is the terminal TF for myelination
  Microglia = c("Spi1", "Irf8", "Runx1"),                  # Spi1 (PU.1) and Irf8 are master regulators
  
  # 3. Vascular & Mesenchymal Transcription Factors
  Endothelial_Cell = c("Erg", "Fli1", "Sox17"),            # Erg/Fli1 are key ETS-family TFs
  Pericyte = c("Foxf2", "Tbx18", "Meis2"),                 # Foxf2 is strongly associated with pericytes
  Fibroblast_VLMC = c("Foxc1", "Twist1", "Tcf21")           # TFs for mesenchymal identity
)


marker_lists <- specific_marker_list 

# 1. Get all unique TFs that are present in the Seurat object's TF activity assay
all_markers <- unique(unlist(marker_lists))
markers_present <- intersect(all_markers, rownames(thisObj[["SCT"]]))

# Check if we have markers to proceed
if (length(markers_present) == 0) {
  stop("None of the specified TFs were found in the 'SCT' assay.")
}

print(markers_present)
# 2. Calculate the average TF activity for these markers across your clusters
# The result is a matrix of TFs (rows) vs. clusters (columns)
avg_exp_matrix <- AverageExpression(
  thisObj,
  features = markers_present,
  group.by = "seurat_clusters_merged", # Or whichever column holds your cluster IDs
  assay = "SCT",
  slot = "data"
)$SCT # Directly extract the matrix from the list

# 3. Calculate the enrichment score for each TF list
# This loop calculates the mean activity score for each TF list in each cluster.
# It iterates through your list of cell types (e.g., Granule_Cell, Purkinje_Cell).
enrichment_scores <- list()
for (cell_type in names(marker_lists)) {
  
  # Find which TFs from the current list are in our matrix
  tfs_in_list <- intersect(marker_lists[[cell_type]], rownames(avg_exp_matrix))
  
  if (length(tfs_in_list) > 0) {
    # Calculate the mean score for this list of TFs across all clusters
    # The 'colMeans' function averages the scores for the selected TFs for each cluster
    enrichment_scores[[cell_type]] <- colMeans(avg_exp_matrix[tfs_in_list, , drop = FALSE])
  }
}

# 4. Combine the scores into a single matrix for plotting
# The result is a matrix of cell types/lists (rows) vs. clusters (columns)
enrichment_matrix <- do.call(rbind, enrichment_scores)

# At this point, 'enrichment_matrix' should be a valid matrix, not NULL.

pheatmap::pheatmap(
  enrichment_matrix,
  main = "TF List Enrichment per Cluster",
  scale = "row",         # Scale scores by row (cell type) to see which cluster is highest for each list
  cluster_cols = TRUE,   # Cluster columns (cell clusters) to group similar clusters
  cluster_rows = TRUE,   # Cluster rows (cell types) to group similar TF lists
  border_color = "white"
)



















print("Generating Heatmap: Relative Cluster Enrichment per Gene List (Scaled by Row)")
UMAP_Cluster_EnrichHeat_RowScale <- pheatmap(cluster_enrichment_matrix,
                                             scale = "row", # Which cluster(s) are prominent for each list?
                                             main = "Relative Cluster Enrichment per Gene List (Scaled by Row)",
                                             cluster_rows = TRUE,
                                             cluster_cols = FALSE)

file_name <- paste0(local_WD, "UMAP_Cluster_EnrichHeat_RowScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
UMAP_Cluster_EnrichHeat_RowScale
dev.off() 



# Gene wise heatmap
marker_lists_to_use <- gene_marker_lists
all_markers_filtered <- unique(unlist(marker_lists_to_use))
assay_to_use <- "SCT"
all_markers_present <- intersect(all_markers_filtered, rownames(thisObj[[assay_to_use]]))

avg_expression <- AverageExpression(
  thisObj,
  features = all_markers_present,
  group.by = "ident",
  assay = assay_to_use,
  slot = "data" # Or "scale.data" if you prefer Z-scores directly
)

# Extract the matrix (handle potential list output from AverageExpression)
if (is.list(avg_expression)) {
  avg_expression_matrix <- avg_expression[[assay_to_use]]
} else {
  avg_expression_matrix <- avg_expression
}

# Scale the data (Z-score across clusters for each gene) for visualization
scaled_avg_expression <- t(scale(t(avg_expression_matrix)))
scaled_avg_expression[is.nan(scaled_avg_expression)] <- 0 # Handle cases with no variance

# --- Generate Heatmap using pheatmap ---
pheatmap_plot <- pheatmap(
  mat = scaled_avg_expression,
  scale = "none", # Data is already scaled
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE, # Set to FALSE if too many genes
  show_colnames = TRUE,
  fontsize_row = 8, # Adjust font size as needed
  fontsize_col = 8,
  main = "Average Scaled Marker Expression per Cluster"
  # You can customize colors using the 'color' argument, e.g.:
  # color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
)

























































































gene_marker_lists <- list(
  # 1. Neurons (Excitatory: Granule + UBCs)
  Pan_Neuronal = c("Neurod1", "Neurod2", "Neurod6", "Mecp2", "Rest"), # Key neuronal development and identity TFs
  Granule = c("Pax6", "Zeb1", "Neurod1"), # Pax6 and Neurod1 are crucial for granule cell development
  
  # 2. Interneurons (Inhibitory)
  Neuron_Inhibitory = c("Dlx2", "Dlx5", "Ascl1", "Hes1"), # Dlx family are master regulators of GABAergic neurons
  
  # 3. Purkinje Cells
  Purkinje = c("Rora", "Pou4f2", "Lmx1b", "Tbx2"), # Rora is a classic and essential Purkinje cell TF
  
  # 4. Glia (Astrocytes, Bergmann Glia)
  Astrocyte_Glia = c("Sox9", "Nfia", "Nfib", "Stat3", "Hes1", "Hesx1"), # Sox9 and Nfi family are key glial regulators
  
  # 5. Fibroblasts / Mesenchymal
  Fibroblast_VLMC = c("Foxc2", "Twist1", "Twist2", "Tcf7l2"),
  
  # 6. Microglia
  Microglia = c("Spi1", "Irf8", "Irf5", "Irf9", "Runx1", "Stat6"), # Spi1 (PU.1) and the IRF family are master regulators of microglia
  
  # 7. Oligodendrocytes & OPCs
  Oligodendrocyte_Lineage = c("Olig1", "Sox10", "Tcf7l2", "Nkx6-1"), # Olig1/2 and Sox10 are master regulators of this lineage
  
  # 8. Endothelial & Vascular Cells (Endothelial, Pericyte, VSMC)
  Vascular = c("Erg", "Fli1", "Foxc2", "Sox17", "Etv6") # ERG and FLI1 are key endothelial TFs
)

all_markers <- VariableFeatures(thisObj)
pca_loadings <- Loadings(thisObj[["pca"]])
markers_in_pca <- intersect(all_markers, rownames(pca_loadings))

marker_loadings <- pca_loadings[markers_in_pca, ]

top_pcs_per_marker <- apply(marker_loadings, 1, function(gene_row) {
  # Get absolute loadings, sort descending, take top N (e.g., top 3)
  sorted_loadings <- sort(abs(gene_row), decreasing = TRUE)
  names(head(sorted_loadings, 3)) # Return names of the top 3 PCs
})
print("Top 3 PCs for each marker (based on absolute loading):")
print(top_pcs_per_marker)

top_markers_per_pc <- apply(marker_loadings, 2, function(pc_col) {
  sorted_loadings <- sort(abs(pc_col), decreasing = TRUE)
  names(head(sorted_loadings, 5)) # Return names of top 5 markers for this PC
})
print("Top 5 marker genes for each PC (from your list):")
print(top_markers_per_pc[1:min(10, ncol(marker_loadings))]) # Adjust range as needed

pcs_to_analyze <- 1:50
pcs_to_analyze <- intersect(pcs_to_analyze, 1:ncol(pca_loadings))
if (length(pcs_to_analyze) == 0) {
  stop("No valid PCs selected for analysis.")
}
loadings_subset <- pca_loadings[, pcs_to_analyze, drop = FALSE]
enrichment_scores <- list()

for (list_name in names(gene_marker_lists)) {
  gene_list <- gene_marker_lists[[list_name]]
  
  # Find which genes from the list are present in the PCA loadings
  markers_in_pca <- intersect(gene_list, rownames(loadings_subset))
  
  if (length(markers_in_pca) > 0) {
    # Get loadings for these markers for the selected PCs
    marker_loadings_subset <- abs(loadings_subset[markers_in_pca, , drop = FALSE])
    
    # Calculate sum of absolute loadings for this list across the PCs
    # colSums calculates the sum for each column (PC)
    scores <- colSums(marker_loadings_subset)
    enrichment_scores[[list_name]] <- scores
  } else {
    # Handle case where no markers from the list are found in PCA
    enrichment_scores[[list_name]] <- rep(0, length(pcs_to_analyze))
    names(enrichment_scores[[list_name]]) <- colnames(loadings_subset)
    warning(paste("No markers from list", list_name, "found in PCA loadings."))
  }
}

enrichment_matrix <- do.call(rbind, enrichment_scores)
colnames(enrichment_matrix) <- paste0("PC_", pcs_to_analyze) # Ensure column names match PCs
print(enrichment_matrix)

top_pc_per_list <- apply(enrichment_matrix, 1, function(row) {
  colnames(enrichment_matrix)[which.max(row)]
})

print("Top PC for each gene list (based on sum of absolute loadings):")
print(top_pc_per_list)
plot.new()
PC_Heat_ColScale <- pheatmap(enrichment_matrix,
                             scale = "column",
                             main = "Relative Gene List Enrichment per PC (Scaled by Column)",
                             cluster_rows = TRUE,
                             cluster_cols = FALSE)


file_name <- paste0(local_WD, "PC_Heat_ColScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
PC_Heat_ColScale
dev.off() 

# Scale by row (emphasizes top PC(s) per list)
pheatmap(enrichment_matrix,
         scale = "row",
         main = "Relative PC Enrichment per Gene List (Scaled by Row)",
         cluster_rows = TRUE,
         cluster_cols = TRUE)
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**UMAP*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# dims_to_use <- c(1, 2, 3, 5, 6, 7, 11, 12, 13, 14, 16, 17, 19, 22, 29) #original list, PC selection from gene set enrichment across PCs
# dims_to_use <- c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 17, 19, 23, 27, 29, 30) #new list, PC selection from gene set enrichment across PCs
# dims_to_use <- c(1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13,14,15,16, 17, 18, 19,20,21,22, 23, 24,25,26) #attempt 3 only using the cerebellum gene list for PCA reduction features
# dims_to_use <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 18, 19, 20, 21, 22, 25, 28, 29, 30) #attempt 4 only using the NEW cerebellum gene list for PCA reduction features
# dims_to_use <- c(1, 2, 3, 5, 6, 10, 11, 12, 14, 15, 16, 24, 25, 27, 29, 30) #attempt 6 using new CB gene list but still using all features list cross-referenced with gene lists

# winner, maybe remove 31
dims_to_use <- c(15, 34, 11, 30, 27, 25, 33, 31, 14, 5, 12, 9, 21, 3, 1 ) #Analzyed 1:50 instead of 1:30

dims_to_use <- c(1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13,14,15,16, 17, 18, 19,20,21,22, 23, 24,25,26)
thisObj <- FindNeighbors(thisObj, dims = dims_to_use,
                         k.param = 15, assay = "SCT")
thisObj <- FindClusters(thisObj, algorithm = 3, reduction = "pca", resolution = 0.75)
thisObj <- RunUMAP(thisObj, dims = dims_to_use,  assay = "SCT", reduction = "pca", min.dist = 0.0001, spread = 3)

umap_p <- DimPlot(thisObj, reduction = "umap", raster = TRUE, label = TRUE, pt.size = 0.5) + ggtitle("UMAP")
print(umap_p)

ggsave(filename = paste0(local_WD, "umap_p", ".png"), plot = umap_p, dpi = 900, units = "in", width = 8, height = 9)

FeaturePlot(thisObj, reduction = "umap", features = "Pcp2")


saveRDS(thisObj, file.path(HDD_WD, "COSMONAUT_CB_Seurat_postPCA_newcountcutoff.RDS"))

thisObj <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_postPCA_newcountcutoff.RDS"))
# Plot spatial map from UMAP clusters

thisObj@meta.data$x_slide_micron <- thisObj@meta.data$x_slide_mm * 1000
thisObj@meta.data$y_slide_micron <- thisObj@meta.data$y_slide_mm * 1000


# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"

plotSpatialSample <- function(seurat_obj,
                              sample_id,
                              sample_col,
                              coord_x_micron = "x_slide_micron",
                              coord_y_micron = "y_slide_micron",
                              color_by = "ident",
                              custom_palette = NULL,
                              pt_size = 1,
                              alpha = 0.6,
                              plot_title = NULL,
                              # "ggplot" or "plotly" (interactive plot)
                              return_type = "ggplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (return_type == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' needed for return_type='plotly'. Please install it.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Sample column '", sample_col, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_x_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("X coordinate column '", coord_x_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_y_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Y coordinate column '", coord_y_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!(color_by == "ident" || color_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Coloring column/variable '", color_by, "' not found in Seurat object metadata or is not 'ident'."), call. = FALSE)
  }
  if (!sample_id %in% unique(seurat_obj@meta.data[[sample_col]])) {
    stop(paste("Sample ID '", sample_id, "' not found in column '", sample_col, "'."), call. = FALSE)
  }
  
  print(paste("Subsetting data for sample:", sample_id))
  subset_data <- subset(seurat_obj, subset = !!sym(sample_col) == sample_id)
  
  if (nrow(subset_data@meta.data) == 0) {
    stop(paste("No cells found after subsetting for sample ID '", sample_id, "'."), call. = FALSE)
  }
  
  plot_df <- subset_data@meta.data
  if (color_by == "ident") {
    plot_df$color_values <- as.character(Idents(subset_data))
    color_lab <- "Cluster ID"
  } else {
    plot_df$color_values <- as.character(plot_df[[color_by]])
    color_lab <- color_by
  }
  plot_df$color_values <- factor(plot_df$color_values) # Ensure factor for levels
  
  x_range <- range(plot_df[[coord_x_micron]], na.rm = TRUE)
  y_range <- range(plot_df[[coord_y_micron]], na.rm = TRUE)
  buffer <- 10
  xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
  ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Sample:", sample_id, "- Colored by:", color_lab)
  }
  
  # --- Determine Color Palette ---
  color_levels <- levels(plot_df$color_values)
  num_colors <- length(color_levels)
  
  if (is.null(custom_palette)) {
    # Generate default palette if none provided
    if (num_colors <= 12) {
      color_palette_values <- hcl.colors(num_colors, palette = "Set 3")
    } else {
      color_palette_values <- scales::hue_pal()(num_colors)
    }
    # Name palette for consistent mapping
    names(color_palette_values) <- color_levels
  } else {
    # Use the provided custom palette
    color_palette_values <- custom_palette
    # Optional: Add checks here if needed (e.g., ensure names match levels)
    # but adhering to 'no extra checks' request for now.
  }
  
  
  print("Generating ggplot...")
  p <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]], color = color_values)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_palette_values) + # Use the chosen palette
    labs(title = plot_title,
         x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
    theme_minimal() +
    coord_fixed() +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
  
  if (return_type == "plotly") {
    print("Generating plotly object...")
    plot_df$tooltip_text <- paste("Cell:", rownames(plot_df), "<br>", color_lab, ":", plot_df$color_values)
    p_plotly <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]],
                                    color = color_values, text=tooltip_text)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_manual(values = color_palette_values) + # Use the chosen palette
      labs(title = plot_title, x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
      theme_minimal() + coord_fixed() + xlim(xlim_vals) + ylim(ylim_vals) +
      theme( panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.box = "horizontal") +
      guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
    
    suppressWarnings({
      plotly_obj <- plotly::ggplotly(p_plotly, tooltip = "text")
    })
    return(plotly_obj)
    
  } else {
    print("Returning ggplot object.")
    return(p)
  }
}


YM04_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj,
                  sample_id = "YM04",
                  sample_col = "mus",
                  coord_x_micron = "x_slide_micron",
                  coord_y_micron = "y_slide_micron",
                  color_by = "ident",
                  custom_palette = NULL,
                  pt_size = 1,
                  alpha = 0.6,
                  plot_title = NULL,
                  # Return type "ggplot" or "plotly" (interactive)
                  return_type = "ggplot")

print(YM04_spatial_ggplot)

#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Cell Typing- Azumith Test*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Ref Data: Kozareva et al. (2021) Cerebellum Atlas (Broad Single Cell Portal)
# Ref Data link: https://singlecell.broadinstitute.org/single_cell/study/SCP795/a-transcriptomic-atlas-of-the-mouse-cerebellum

# Download and prepare reference dataset
# CB_ref_SeuratObj <- readRDS(file.path("C:/Users/m322332/Downloads/","cb_annotated_object.RDS"))
# CB_ref_SeuratObj <- UpdateSeuratObject(CB_ref_SeuratObj) #Update from v2 => v3
# options(future.globals.maxSize = 30 * 1024^3) # Increase memory allocation for reference file, reset after running!
# set conser.memory = TRUE since scREF is too large to compute after multiple different attempts; 4/17/2025
# CB_ref_SeuratObj <- SCTransform(CB_ref_SeuratObj, vst.flavor='v2', conserve.memory = TRUE, verbose = TRUE) #Calculate SCT to match query data
# CB_ref_SeuratObj <- RunPCA(CB_ref_SeuratObj, assay = "SCT", verbose = TRUE) #Process PCA for SCT
# saveRDS(CB_ref_SeuratObj, file.path(HDD_WD, "CB_ref_SeuratObj_SCTPCA"))
CB_ref_SeuratObj <- readRDS(file.path(HDD_WD, "CB_ref_SeuratObj_SCTPCA"))

# Calculating reference features to add additional genes to CB gene list
# 4.17.2025 Note
# I am iterating through the process determining appropriate genes and PCs and gene features
# that best differentiate cell types and denoise data and ID incorrectly segmented samples
# with the intent to run a quick ReSeg on problem samples or exclude
reference_markers <- FindAllMarkers(
  CB_ref_SeuratObj,
  assay = "SCT",
  only.pos = TRUE,            # Only find positive markers for each identity
  min.pct = 0.1,              # Detect in at least 10% of cells in the group
  logfc.threshold = 0.25,     # Minimum log2 fold change
  max.cells.per.ident = 1000, # Reduce samples per measure to limit memory issues, slows calculation
  verbose = TRUE
)

# Combine/clean markers list
reference_markers <- reference_markers %>%
  mutate(original_cluster = stringr::str_remove(cluster, "_\\d+$"))
head(reference_markers)

# saveRDS(reference_markers, file.path(HDD_WD, "CB_ref_SeuratObj_reference_genes_unfiltered"))
reference_markers <- readRDS(file.path(HDD_WD, "CB_ref_SeuratObj_reference_genes_unfiltered"))

reference_markers_top30 <- reference_markers %>%
  group_by(original_cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC, with_ties = FALSE) %>% # Get top 30 per cluster
  ungroup() 

print(reference_markers_top30, n = Inf)




print(head(CB_ref_SeuratObj@meta.data))
colnames(CB_ref_SeuratObj)
ref_annotation_col <- "cluster"

transfer_anchors <- FindTransferAnchors(
  reference = CB_ref_SeuratObj,
  query = thisObj,
  normalization.method = "SCT", # Use "SCT"
  reference.reduction = "pca",  # Assumes PCA is computed on reference
  reduction = "rpca",           # Use reciprocal PCA for robustness
  dims = 1:30                   # Adjust based on PCs used in your 'thisObj' analysis
  # May need additional parameters depending on object structure/assays
)

thisObj <- TransferData(
  anchorset = transfer_anchors,
  refdata = CB_ref_SeuratObj@meta.data[, ref_annotation_col], # Get the specific annotation column
  query = thisObj,
  weight.reduction = thisObj[['pca']], # Or use query PCA reduction if RPCA wasn't used above
  dims = 1:30 # Usually inherits from anchorset, but can specify
)

head(thisObj,20)

options(future.globals.maxSize = 3 * 1024^3) #reset memory allocation; may be unnecessary, just being safe

DimPlot(thisObj, group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Predicted Cell Types (Reference Mapping)")

hist(thisObj$predicted.id.score, main="Distribution of Max Prediction Scores")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Cell Typing*
#Now we run a similarity matrix of clusters x gene sets for initial cell type ID
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster_ident_name <- as.character(Idents(thisObj)[1]) # Get the name of one cluster ID to find column name later
print(paste("Using clusters from identity:", levels(Idents(thisObj))))

# Use the filtered marker list created previously
marker_lists_to_use <- gene_marker_lists # Or use cerebellum_marker_lists if preferred

# Get all unique genes from the lists
all_markers_filtered <- unique(unlist(marker_lists_to_use))

# Ensure markers exist in the assay to be used (e.g., SCT)
all_markers_present <- intersect(all_markers_filtered, rownames(thisObj[["SCT"]]))
if (length(all_markers_present) == 0) {
  stop("None of the marker genes found in the SCT assay of the Seurat object.")
}
print(paste("Analyzing enrichment for", length(all_markers_present), "marker genes present in the SCT assay."))


# 2. Calculate Average Expression per Cluster
# Using the SCT assay data slot is recommended since SCT was used for normalization/PCA
avg_expression <- AverageExpression(
  thisObj,
  features = all_markers_present,
  group.by = "ident", # Use the active cluster identity
  assay = "SCT",      
  slot = "data"       # Use normalized data
)

if (is.list(avg_expression)) {
  avg_expression_matrix <- avg_expression[["SCT"]] # Maybe "SCT" wasn't found?
} else {
  avg_expression_matrix <- avg_expression
}
# Error occurs here if avg_expression_matrix wasn't successfully created above
scaled_avg_expression <- t(scale(t(avg_expression_matrix)))
# Replace NaN values that can result from scaling genes with 0 variance across clusters
scaled_avg_expression[is.nan(scaled_avg_expression)] <- 0

cluster_enrichment_scores <- list()

for (list_name in names(marker_lists_to_use)) {
  gene_list <- marker_lists_to_use[[list_name]]
  markers_in_avg_exp <- intersect(gene_list, rownames(scaled_avg_expression))
  
  if (length(markers_in_avg_exp) > 0) {
    # Get scaled average expression for these markers
    list_scaled_avg_exp <- scaled_avg_expression[markers_in_avg_exp, , drop = FALSE]
    # Calculate mean scaled score per cluster for this list
    scores <- colMeans(list_scaled_avg_exp, na.rm = TRUE)
    cluster_enrichment_scores[[list_name]] <- scores
  } else {
    # Handle case where no markers from the list were found/calculated
    num_clusters <- ncol(scaled_avg_expression)
    cluster_enrichment_scores[[list_name]] <- rep(0, num_clusters)
    names(cluster_enrichment_scores[[list_name]]) <- colnames(scaled_avg_expression)
    warning(paste("No markers from list", list_name, "found in average expression matrix."))
  }
}

# 4. Combine into a Matrix
cluster_enrichment_matrix <- do.call(rbind, cluster_enrichment_scores)
# Column names should already be cluster IDs from AverageExpression

# 5. Visualize with pheatmap (Similar logic as with PCs)

print("Generating Heatmap: Relative Gene List Enrichment per Cluster (Scaled by Column)")
UMAP_Cluster_EnrichHeat_ColScale <- pheatmap(cluster_enrichment_matrix,
         scale = "column", # Which list(s) are prominent within each cluster?
         main = "Relative Gene List Enrichment per Cluster (Scaled by Column)",
         cluster_rows = TRUE, # Cluster rows (gene lists) by similarity
         cluster_cols = FALSE) # Cluster columns (clusters) by similarity

file_name <- paste0(local_WD, "UMAP_Cluster_EnrichHeat_ColScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
UMAP_Cluster_EnrichHeat_ColScale
dev.off() 

print("Generating Heatmap: Relative Cluster Enrichment per Gene List (Scaled by Row)")
UMAP_Cluster_EnrichHeat_RowScale <- pheatmap(cluster_enrichment_matrix,
         scale = "row", # Which cluster(s) are prominent for each list?
         main = "Relative Cluster Enrichment per Gene List (Scaled by Row)",
         cluster_rows = TRUE,
         cluster_cols = FALSE)

file_name <- paste0(local_WD, "UMAP_Cluster_EnrichHeat_RowScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
UMAP_Cluster_EnrichHeat_RowScale
dev.off() 



# Gene wise heatmap
marker_lists_to_use <- gene_marker_lists
all_markers_filtered <- unique(unlist(marker_lists_to_use))
assay_to_use <- "SCT"
all_markers_present <- intersect(all_markers_filtered, rownames(thisObj[[assay_to_use]]))

avg_expression <- AverageExpression(
  thisObj,
  features = all_markers_present,
  group.by = "ident",
  assay = assay_to_use,
  slot = "data" # Or "scale.data" if you prefer Z-scores directly
)

# Extract the matrix (handle potential list output from AverageExpression)
if (is.list(avg_expression)) {
  avg_expression_matrix <- avg_expression[[assay_to_use]]
} else {
  avg_expression_matrix <- avg_expression
}

# Scale the data (Z-score across clusters for each gene) for visualization
scaled_avg_expression <- t(scale(t(avg_expression_matrix)))
scaled_avg_expression[is.nan(scaled_avg_expression)] <- 0 # Handle cases with no variance

# --- Generate Heatmap using pheatmap ---
pheatmap_plot <- pheatmap(
  mat = scaled_avg_expression,
  scale = "none", # Data is already scaled
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE, # Set to FALSE if too many genes
  show_colnames = TRUE,
  fontsize_row = 5, # Adjust font size as needed
  fontsize_col = 8,
  main = "Average Scaled Marker Expression per Cluster"
  # You can customize colors using the 'color' argument, e.g.:
  # color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100)
)

# pheatmap usually draws directly, but storing in object allows saving etc.
# print(pheatmap_plot)

#'[We will also check the highest expressing markers and dotplots of cell type markers more closely]
all_markers <- FindAllMarkers(thisObj,
                              assay = "SCT",
                              only.pos = TRUE,    # Keep only positive markers
                              min.pct = 0.1,      # min detection rate
                              logfc.threshold = 0.25) # log fold change threshold

significant_markers <- all_markers %>%
  filter(p_val_adj <= 0.05)

# Get top 5 significant markers per cluster based on avg_log2FC
top_genes <- significant_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

print("Top 5 significant genes per cluster:")
print(head(top_genes))

top_genes <- top_genes %>% arrange(cluster, desc(avg_log2FC))

top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))

gene_plot <- ggplot(top_genes, aes(x = factor(cluster), y = gene, fill = avg_log2FC, size = pct.1)) + # Map avg_log2FC to fill
  geom_point(shape = 21, color = "black", stroke = 0.2) + # Use shape 21, set outline color and thinness (stroke)
  scale_fill_viridis_c(name = "Avg Log2 FC") + # Change to scale_fill_*
  scale_size_continuous(name = "Percent\nExpressed") +
  labs(title = "Top 5 Significant Genes per Cluster",
       x = "Cluster", # This will become y-axis after flip
       y = "Gene") + # This will become x-axis after flip
  theme_minimal() +
  coord_flip() + # Flip coordinates
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), # Now controls gene labels
        axis.text.y = element_text(size = 8), # Now controls cluster labels
        plot.title = element_text(hjust = 0.5)) # Center title

print(gene_plot)




plotClusterMarkers <- function(seurat_obj,
                               marker_results,
                               curated_marker_list,
                               target_cluster_id,
                               n_top = 25,
                               assay = "SCT",
                               pval_cutoff = 0.05,
                               group_by_col = NULL) {
  
  # --- Input Checks (Minimal) ---
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed. Please install it.", call. = FALSE)
  }
  if (!is.data.frame(marker_results) || !all(c("cluster", "gene", "avg_log2FC", "pct.1", "p_val_adj") %in% colnames(marker_results))) {
    stop("'marker_results' must be a data frame from FindAllMarkers with required columns.", call. = FALSE)
  }
  if (!is.list(curated_marker_list)) {
    stop("'curated_marker_list' must be a named list.", call. = FALSE)
  }
  
  # --- Set Active Identity ---
  orig_idents <- Idents(seurat_obj) # Store original
  if (!is.null(group_by_col)) {
    if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Specified group_by_col '", group_by_col, "' not found."), call. = FALSE)
    }
    Idents(seurat_obj) <- seurat_obj[[group_by_col]]
  }
  # Ensure target cluster exists
  if (!as.character(target_cluster_id) %in% levels(Idents(seurat_obj))) {
    warning(paste("Target cluster ID '", target_cluster_id, "' not found in active identities. Restoring original idents."), call. = FALSE)
    Idents(seurat_obj) <- orig_idents
    return(NULL) # Or stop()
  }
  
  
  # --- Data Preparation ---
  
  # Filter significant markers
  significant_markers <- marker_results %>%
    filter(p_val_adj <= pval_cutoff)
  
  # Get top N significant markers for the target cluster
  top_markers_df <- significant_markers %>%
    filter(cluster == target_cluster_id) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = n_top)
  
  if (nrow(top_markers_df) == 0) {
    warning(paste("No significant markers passed filters for cluster", target_cluster_id, ". Plots might be empty or inaccurate. Restoring original idents."), call. = FALSE)
    Idents(seurat_obj) <- orig_idents
    return(NULL) # Return NULL if no markers
  }
  
  # Get the curated list of marker genes
  curated_markers_all <- unique(unlist(curated_marker_list))
  curated_markers_present <- intersect(curated_markers_all, rownames(seurat_obj[[assay]])) # Check specified assay
  
  # --- Create Left Plot (Top N DEGs) ---
  top_markers_df$gene <- factor(top_markers_df$gene, levels = rev(top_markers_df$gene))
  
  plot_left <- ggplot(top_markers_df, aes(x = avg_log2FC, y = gene, color = avg_log2FC, size = pct.1)) +
    geom_point() +
    scale_color_viridis_c(name = "Avg Log2 FC") + # Viridis default
    scale_size_continuous(name = "Percent\nExpressed", limits = c(0, 1)) +
    labs(title = paste("Top", n_top, "Markers"), x = "Avg Log2 FC", y = "Gene") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_blank(),
          legend.position = "right")
  
  # --- Create Right Plot (Curated Markers) - CORRECTED ---
  plot_right <- NULL # Initialize
  if(length(curated_markers_present) > 0) {
    
    # --- Calculate DotPlot data across ALL clusters FIRST ---
    # Ensure original idents (or specified group_by column) are active
    if (!is.null(group_by_col)) {
      Idents(seurat_obj) <- seurat_obj[[group_by_col]]
    } else {
      # If group_by_col is NULL, ensure the original idents stored at function start are active
      # This assumes orig_idents holds the full cluster set
      Idents(seurat_obj) <- orig_idents
    }
    print(paste("Calculating DotPlot stats across identities:", paste(levels(Idents(seurat_obj)), collapse=", ")))
    
    dotplot_data_all <- tryCatch({
      # Calculate for all idents
      DotPlot(seurat_obj, features = curated_markers_present, assay = assay)$data
    }, error = function(e) {
      warning("Could not generate DotPlot data for curated markers: ", e$message, call. = FALSE)
      return(NULL)
    })
    # --- Restore original idents AFTER calculation (if changed) ---
    # Idents(seurat_obj) <- orig_idents # Optional if only using default idents initially
    
    if (!is.null(dotplot_data_all)) {
      # --- NOW filter for the target cluster ---
      dotplot_data_filt <- dotplot_data_all %>%
        filter(id == target_cluster_id) %>%
        # Add missing genes with 0 expression for plotting consistency
        tidyr::complete(features.plot = factor(curated_markers_present, levels = curated_markers_present),
                        id = factor(target_cluster_id),
                        # Fill with sensible values: 0 avg.exp, 0 pct.exp, and minimum observed scaled value across all clusters
                        fill = list(avg.exp = 0, pct.exp = 0, avg.exp.scaled = min(dotplot_data_all$avg.exp.scaled, na.rm=TRUE)-1))
      
      if (nrow(dotplot_data_filt) > 0) {
        # Order curated genes alphabetically for plotting
        dotplot_data_filt$features.plot <- factor(dotplot_data_filt$features.plot, levels = rev(sort(curated_markers_present)))
        
        # Generate the plot using correctly scaled values
        plot_right <- ggplot(dotplot_data_filt, aes(x = avg.exp.scaled, y = features.plot, color = avg.exp.scaled, size = pct.exp / 100)) +
          geom_point() +
          scale_color_viridis_c(option = "magma", name = "Avg Exp\nScaled") + # Using magma scale
          scale_size_continuous(name = "Percent\nExpressed", limits = c(0, 1)) +
          labs(title = "Curated Markers", x = "Avg Expression (Scaled)", y = "Gene") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 10),
                axis.text.y = element_text(size = 8),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_blank(),
                legend.position = "right")
      } else {
        warning(paste("No data remained after filtering DotPlot results for target cluster", target_cluster_id), call.=FALSE)
        plot_right <- ggplot() + theme_void() + ggtitle("Curated Markers\n(Filtering Error)") + theme(plot.title = element_text(hjust = 0.5, size = 10))
      }
    } else {
      # Create empty plot if dotplot data failed initially
      plot_right <- ggplot() + theme_void() + ggtitle("Curated Markers\n(Data Error)") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
  } else {
    warning("No curated marker genes found in the specified assay.", call. = FALSE)
    # Create empty plot if no markers
    plot_right <- ggplot() + theme_void() + ggtitle("Curated Markers\n(None Found)") + theme(plot.title = element_text(hjust = 0.5, size = 10))
  }

  Idents(seurat_obj) <- orig_idents
  
  # --- Combine Plots ---
  # Ensure plot_right exists before combining
  if (!is.null(plot_right)) {
    combined_plot <- plot_left + plot_right + plot_layout(guides = 'collect') & ggtitle(paste("Cluster:", target_cluster_id)) & theme(plot.title = element_text(hjust = 0.5))
  } else {
    # If plot_right failed, just show plot_left with title
    combined_plot <- plot_left + ggtitle(paste("Cluster:", target_cluster_id, "(Top Markers Only)")) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(combined_plot)
}


Cluster_Dotplots <- list()
for (i in 0:(length(levels(Idents(thisObj)))-1)) {
  
  result <- plotClusterMarkers(seurat_obj = thisObj,marker_results = significant_markers,curated_marker_list = gene_marker_lists,n_top = 25,assay = "SCT",
                                target_cluster_id = i)
  
  # Store the result in the list, naming it (e.g., Cluster_1, Cluster_2)
  Cluster_Dotplots[[paste0("Cluster_", i)]] <- result
  file_name <- paste0(local_WD, "Dotplot_Cluster_", i, ".png")
  ggsave(filename = file_name, plot = result,  dpi = 600)
}

print(Cluster_Dotplots$Cluster_0)
print(Cluster_Dotplots$Cluster_1)
print(Cluster_Dotplots$Cluster_2)
print(Cluster_Dotplots$Cluster_3)
print(Cluster_Dotplots$Cluster_4)
print(Cluster_Dotplots$Cluster_5)
print(Cluster_Dotplots$Cluster_6)
print(Cluster_Dotplots$Cluster_7)
print(Cluster_Dotplots$Cluster_8)
print(Cluster_Dotplots$Cluster_9)
print(Cluster_Dotplots$Cluster_10)
print(Cluster_Dotplots$Cluster_11)
print(Cluster_Dotplots$Cluster_12)
print(Cluster_Dotplots$Cluster_13)
print(Cluster_Dotplots$Cluster_14)
print(Cluster_Dotplots$Cluster_15)
print(Cluster_Dotplots$Cluster_16)
print(Cluster_Dotplots$Cluster_17)
print(Cluster_Dotplots$Cluster_18)
print(Cluster_Dotplots$Cluster_19)
print(Cluster_Dotplots$Cluster_20)
print(Cluster_Dotplots$Cluster_21)
print(Cluster_Dotplots$Cluster_22)

#'[Preliminary Cell Type IDs]
# Cluster_Labels <- c(
#   "ODC", # 0
#   "Mix_GL_1", # 1
#   "MLI_2", # 2 #MLI inhibitory of PKJ type 2
#   "Mix_GL_2", # 3
#   "Mix_GL_3", # 4
#   "Mix_GL_4", # 5
#   "Mix_GL_5", # 6
#   "Mix_GL_6", # 7 (ODC contaminated)
#   "Astrocyte_1", # 8
#   "Bergmann_Glia", # 9
#   "MLI_1", # 10 #Type 1 excitatory MLI <= double check location
#   "Mix_GL_7", # 11
#   "Mix_GL_8", # 12
#   "Mix_GL_9", # 13
#   "EC", # 14
#   "PKJ", # 15
#   "Mix_GL_10", # 16
#   "Mix_GL_11", # 17
#   "Mix_GL_12", # 18 (ODC contaminated)
#   "Vasc", # 19
#   "Mix_GL_11", # 20
#   "Mix_GL_12", # 21
#   "Mix_GL_13", # 22
#   "Mix_GL_14", # 23
#   "PL_Mix", # 24 #purkinje layer mix of interneuron, PKJ, and Bergmann_glia
#   "Astrocyte_2", # 25, possible inflammatory profile, check age
#   "Microglia", # 26
#   "MLI_3", # 27, possible MLI_1 check location
#   "Mix_GL_15", # 28
#   "OPC", # 29
#   "Choroid_Plexus" # 30
# )

Cluster_Labels <- c(
  "ODC", # 0
  "Mix_GL", # 1
  "MLI_2", # 2 #MLI inhibitory of PKJ type 2
  "Mix_GL", # 3
  "Mix_GL", # 4
  "Mix_GL", # 5
  "Mix_GL", # 6
  "Mix_GL", # 7 (ODC contaminated)
  "Astrocyte_1", # 8
  "Bergmann_Glia", # 9
  "Golgi", # 10 #Possible golgi based on location and inhbibitory expression
  "Mix_GL", # 11
  "Mix_GL", # 12
  "Mix_GL", # 13
  "EC", # 14
  "PKJ", # 15
  "Mix_GL", # 16
  "Mix_GL", # 17
  "Mix_GL", # 18 (ODC contaminated)
  "Vasc", # 19
  "Mix_GL", # 20
  "Mix_GL", # 21
  "Mix_GL", # 22
  "Mix_GL", # 23
  "PL_Mix", # 24 #purkinje layer mix of interneuron, PKJ, and Bergmann_glia
  "Astrocyte_2", # 25, possible inflammatory profile, check age
  "Microglia", # 26
  "DCN", # 27, confirmed DCN based on location
  "Mix_GL", # 28
  "OPC", # 29
  "Choroid_Plexus" # 30
)

thisObj.T <- thisObj
names(Cluster_Labels) <- levels(thisObj.T)
thisObj.T <- RenameIdents(thisObj.T, Cluster_Labels)
DimPlot(thisObj.T, reduction = "umap", group.by = "age", raster = FALSE, label = TRUE, pt.size = 0.1) + NoLegend()


# saveRDS(thisObj.T, file.path(HDD_WD, "COSMONAUT_CB_Seurat_GlobalCellTyped.RDS"))

thisObj.T <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_GlobalCellTyped.RDS"))
#'[Sub-cluster and check flagged clusters. They showed mixed expression and]
#'[may be segmentation/noise]



custom_col = list(
  "Astrocyte_1" = "orchid",
  "Astrocyte_2" = "orchid4",
  "Bergmann_Glia" = "orange",
  "EC" = "violet",
  "Golgi" = "blue",
  "Microglia" = "black",
  "Mix_GL" = "lightblue",
  "Golgi" = "green",
  "MLI_2" = "purple",
  "DCN" = "darkgreen",
  "ODC" = "yellow2",
  "OPC" = "tan",
  "PKJ" = "red",
  "Vasc" = "grey",
  "Choroid_Plexus" = "pink",
  "PL_Mix" = "lightgrey"
)

umap_p <- DimPlot(thisObj.T,
                  reduction = "umap",
                  raster = FALSE,
                  label = FALSE,
                  pt.size = 0.5,
                  repel = TRUE) + # Add repel = TRUE here
  scale_color_manual(values = custom_col) +
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5)) # Center title if desired

print(umap_p)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped", ".png")
ggsave(filename = file_name, plot = umap_p, dpi = 900, units = "in", width = 10, height = 11)


umap_age <- DimPlot(thisObj.T,
                          group.by = "age",
                          reduction = "umap",
                          raster = FALSE,
                          label = FALSE,
                          pt.size = 0.5,
                          repel = TRUE) +
  scale_color_manual(values = c("YNG" = "lightblue", "OLD" = "pink2")) +
  ggtitle("UMAP | Age") +
  theme(plot.title = element_text(hjust = 0.5)) # Center title if desired

print(umap_age)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped_Age", ".png")
ggsave(filename = file_name, plot = umap_age, dpi = 900, units = "in", width = 10, height = 11)



thisObj.T$age <- factor(thisObj.T$age, levels = c("YNG", "OLD"))
umap_age_facet <- DimPlot(thisObj.T,
                  group.by = "age",
                  split.by =  "age",
                  reduction = "umap",
                  raster = TRUE,
                  label = FALSE,
                  pt.size = 0.5,
                  repel = TRUE) + 
  scale_color_manual(values = c("YNG" = "deepskyblue", "OLD" = "red")) +
  ggtitle("UMAP | Age") +
  theme(plot.title = element_text(hjust = 0.5)) # Center title if desired

print(umap_age_facet)

file_name <- paste0(local_WD, "UMAP_CB_Celltyped_Age_facet", ".png")
ggsave(filename = file_name, plot = umap_age_facet, dpi = 900, units = "in", width = 20, height = 11)




# Plot spatial map from UMAP clusters


# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"

plotSpatialSample <- function(seurat_obj,
                              sample_id,
                              sample_col,
                              coord_x_micron = "x_slide_micron",
                              coord_y_micron = "y_slide_micron",
                              color_by = "ident",
                              custom_palette = NULL,
                              pt_size = 1,
                              alpha = 0.6,
                              plot_title = NULL,
                              # "ggplot" or "plotly" (interactive plot)
                              return_type = "ggplot") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (return_type == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' needed for return_type='plotly'. Please install it.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Sample column '", sample_col, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_x_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("X coordinate column '", coord_x_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!coord_y_micron %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Y coordinate column '", coord_y_micron, "' not found in Seurat object metadata."), call. = FALSE)
  }
  if (!(color_by == "ident" || color_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Coloring column/variable '", color_by, "' not found in Seurat object metadata or is not 'ident'."), call. = FALSE)
  }
  if (!sample_id %in% unique(seurat_obj@meta.data[[sample_col]])) {
    stop(paste("Sample ID '", sample_id, "' not found in column '", sample_col, "'."), call. = FALSE)
  }
  
  print(paste("Subsetting data for sample:", sample_id))
  subset_data <- subset(seurat_obj, subset = !!sym(sample_col) == sample_id)
  
  if (nrow(subset_data@meta.data) == 0) {
    stop(paste("No cells found after subsetting for sample ID '", sample_id, "'."), call. = FALSE)
  }
  
  plot_df <- subset_data@meta.data
  if (color_by == "ident") {
    plot_df$color_values <- as.character(Idents(subset_data))
    color_lab <- "Cluster ID"
  } else {
    plot_df$color_values <- as.character(plot_df[[color_by]])
    color_lab <- color_by
  }
  plot_df$color_values <- factor(plot_df$color_values) # Ensure factor for levels
  
  x_range <- range(plot_df[[coord_x_micron]], na.rm = TRUE)
  y_range <- range(plot_df[[coord_y_micron]], na.rm = TRUE)
  buffer <- 10
  xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
  ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Sample:", sample_id, "- Colored by:", color_lab)
  }
  
  # --- Determine Color Palette ---
  color_levels <- levels(plot_df$color_values)
  num_colors <- length(color_levels)
  
  if (is.null(custom_palette)) {
    # Generate default palette if none provided
    if (num_colors <= 12) {
      color_palette_values <- hcl.colors(num_colors, palette = "Set 3")
    } else {
      color_palette_values <- scales::hue_pal()(num_colors)
    }
    # Name palette for consistent mapping
    names(color_palette_values) <- color_levels
  } else {
    # Use the provided custom palette
    color_palette_values <- custom_palette
    # Optional: Add checks here if needed (e.g., ensure names match levels)
    # but adhering to 'no extra checks' request for now.
  }
  
  
  print("Generating ggplot...")
  p <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]], color = color_values)) +
    geom_point(size = pt_size, alpha = alpha) +
    scale_color_manual(values = color_palette_values) + # Use the chosen palette
    labs(title = plot_title,
         x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
    theme_minimal() +
    coord_fixed() +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
  
  if (return_type == "plotly") {
    print("Generating plotly object...")
    plot_df$tooltip_text <- paste("Cell:", rownames(plot_df), "<br>", color_lab, ":", plot_df$color_values)
    p_plotly <- ggplot(plot_df, aes(x = .data[[coord_x_micron]], y = .data[[coord_y_micron]],
                                    color = color_values, text=tooltip_text)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_manual(values = color_palette_values) + # Use the chosen palette
      labs(title = plot_title, x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = color_lab) +
      theme_minimal() + coord_fixed() + xlim(xlim_vals) + ylim(ylim_vals) +
      theme( panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.box = "horizontal") +
      guides(color = guide_legend(override.aes = list(size=3), ncol = min(num_colors, 6)))
    
    suppressWarnings({
      plotly_obj <- plotly::ggplotly(p_plotly, tooltip = "text")
    })
    return(plotly_obj)
    
  } else {
    print("Returning ggplot object.")
    return(p)
  }
}


thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000

custom_col = list(
  "Astrocyte_1" = "orchid",
  "Astrocyte_2" = "orchid4",
  "Bergmann_Glia" = "orange",
  "EC" = "violet",
  "Golgi" = "blue",
  "Microglia" = "black",
  "Mix_GL" = "lightblue",
  "Golgi" = "green",
  "MLI_2" = "purple",
  "DCN" = "darkgreen",
  "ODC" = "yellow2",
  "OPC" = "tan",
  "PKJ" = "red",
  "Vasc" = "grey",
  "Choroid_Plexus" = "white"
)



YM04_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj.T,
                                         sample_id = "YM04",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "ident",
                                         custom_palette = custom_col,
                                         pt_size = 1,
                                         alpha = 0.6,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(YM04_spatial_ggplot)


file_name <- paste0(local_WD, "YM04_spatial_ggplot", ".png")
ggsave(filename = file_name, plot = YM04_spatial_ggplot, dpi = 900, units = "in", width = 10, height = 11)



O243_spatial_ggplot <- plotSpatialSample(seurat_obj = thisObj.T,
                                         sample_id = "243",
                                         sample_col = "mus",
                                         coord_x_micron = "x_slide_micron",
                                         coord_y_micron = "y_slide_micron",
                                         color_by = "ident",
                                         custom_palette = custom_col,
                                         pt_size = 1,
                                         alpha = 0.6,
                                         plot_title = NULL,
                                         # Return type "ggplot" or "plotly" (interactive)
                                         return_type = "ggplot")

print(O243_spatial_ggplot)


file_name <- paste0(local_WD, "O243_spatial_ggplot", ".png")
ggsave(filename = file_name, plot = O243_spatial_ggplot, dpi = 900, units = "in", width = 10, height = 11)


# Quick cell count violin plots
meta_data <- thisObj.T@meta.data
meta_data$ClusterID <- Idents(thisObj.T)

sample_col <- "mus"
age_col <- "age"
sex_col <- "sex" # Define the sex column name
cluster_col <- "ClusterID"

sample_totals <- meta_data %>%
  group_by(!!sym(sample_col)) %>%
  summarise(TotalCellsInSample = n(), .groups = 'drop')

cluster_counts_per_sample <- meta_data %>%
  group_by(!!sym(sample_col), !!sym(cluster_col), !!sym(age_col), !!sym(sex_col)) %>% # Added sex_col
  summarise(CellCount = n(), .groups = 'drop')

plot_data <- cluster_counts_per_sample %>%
  left_join(sample_totals, by = sample_col) %>%
  mutate(Proportion = CellCount / TotalCellsInSample) %>%
  mutate(
    ClusterID = factor(!!sym(cluster_col)),
    Age = factor(!!sym(age_col), levels = c("YNG", "OLD")), # Keep specific levels if needed
    Sex = factor(!!sym(sex_col)) # Add Sex as factor
  )


# --- 2. Generate Violin Plots ---

# Option 1: All clusters on one plot (might be crowded)
plot_all <- ggplot(plot_data, aes(x = ClusterID, y = Proportion, fill = factor(Age, levels = c("YNG", "OLD")))) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = factor(Age, levels = c("YNG", "OLD"))), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 0.5, alpha = 0.6) +
  scale_fill_manual(name = "Age", values = c("YNG" = "lightblue", "OLD" = "salmon"), labels = c("YNG", "OLD")) +
  scale_color_manual(name = "Age", values = c("YNG" = "darkblue", "OLD" = "darkred"), labels = c("YNG", "OLD")) +
  labs(
    title = "Normalized Cell Type Proportions by Age Group",
    x = "Cluster ID",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_all)


# Option 2: Facet by cluster (often clearer)
plot_faceted <- ggplot(plot_data, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  # geom_boxplot(width=0.2, outlier.shape = NA, alpha=0.5) + # Optional boxplot
  scale_fill_manual(values = c("OLD" = "salmon", "YNG" = "lightblue")) + # Adjust colors as needed
  facet_wrap(~ ClusterID, scales = "free_y", ncol = 6) + # Adjust ncol as needed
  labs(
    title = "Normalized Cell Type Proportions by Age Group (Faceted by Cluster)",
    x = "Age Group",
    y = "Proportion of Cells per Sample"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot_faceted)


# Microglia only
plot_mg <- ggplot(plot_data_mg, aes(x = Age, y = Proportion, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(aes(shape = Sex), width = 0.1, size = 3.0, alpha = 0.8) +
  scale_fill_manual(values = c("OLD" = "salmon", "YNG" = "lightblue")) + # Hide redundant fill legend
  scale_shape_manual(values = c("F" = 16, "M" = 17)) + # Define shapes: 16=circle, 17=triangle. Adjust 'F'/'M' as needed.
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    title = "Microglia",
    x = NULL, # Remove x-axis title
    y = "Proportion of Cells per Sample",
    shape = "Sex" # Keep shape legend title
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14,face = "bold"), # Center title, increase size
    axis.title.y = element_text(size = 12,face = "bold"), # Increase Y title size
    axis.text.x = element_text(face = "bold", size = 12), # Bold and increase X text size
    axis.title.x = element_blank(), # Also remove x title space via theme
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.border = element_blank(), # Remove panel border
    axis.line = element_line(colour = "black"), # Add axis lines
    legend.position = "bottom" # Keep legend visible for shape
  )

print(plot_mg)


file_name <- paste0(local_WD, "plot_CB_Microglia_proportions", ".png")
ggsave(filename = file_name, plot = plot_mg, dpi = 900)




# --- Generate Stacked Bar Plot ---

# Ensure ClusterID is treated as a factor for discrete colors
plot_data$ClusterID <- factor(plot_data$ClusterID)
# Ensure Age is factor with desired order
plot_data$Age <- factor(plot_data$Age, levels = c("YNG", "OLD"))
# Ensure SampleID is factor if you want discrete x-axis positions
plot_data$SampleID <- factor(plot_data[[sample_col]]) # Use sample_col variable defined before


custom_col = list(
  "Astrocyte_1" = "orchid",
  "Astrocyte_2" = "orchid4",
  "Bergmann_Glia" = "orange",
  "EC" = "violet",
  "Golgi" = "blue",
  "Microglia" = "black",
  "Mix_GL" = "lightblue",
  "Golgi" = "green",
  "MLI_2" = "purple",
  "DCN" = "darkgreen",
  "ODC" = "yellow2",
  "OPC" = "tan",
  "PKJ" = "red",
  "Vasc" = "grey",
  "Choroid_Plexus" = "pink",
  "PL_Mix" = "lightgrey"
)

plot_data$Age <- factor(plot_data$Age, levels = c("YNG", "OLD"))
stacked_bar_plot <- ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = ClusterID)) +
  geom_col(position = "stack") +
  facet_grid(~Age, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_col) + # Use your custom color list here
  labs(
    title = "Cell Type Composition per Sample by Age Group",
    x = "Sample ID",
    y = "Percentage of Cells",
    fill = "Cluster ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing.x = unit(1, "lines")
  )

print(stacked_bar_plot)


# --- Optional: Facet by Sex as well ---
stacked_bar_plot_sex <- ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = ClusterID)) +
  geom_col(position = "stack") +
  facet_grid(Sex ~ Age, scales = "free_x", space = "free_x") + # Separate rows for Sex, columns for Age
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Cell Type Composition per Sample by Age and Sex",
    x = "Sample ID",
    y = "Percentage of Cells",
    fill = "Cluster ID"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing = unit(1, "lines")
  )

print(stacked_bar_plot_sex)






















































































#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Sen-GSEA*
# install.packages("UCell")
library(UCell)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
my_geneset = c("Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2",
                            "Bmp6", "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24", "Ccl26", "Ccl3",
                            "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55", "Cd9", "Csf1", "Csf2",
                            "Csf2rb", "Ctnnb1", "Ctsb", "Cxcl1", "Cxcl10", "Cxcl12", "Cxcl16",
                            "Cxcl2", "Cxcl3", "Cxcr2", "Dkk1", "Edn1", "Egf", "Egfr", "Ereg",
                            "Esm1", "Ets2", "Fas", "Fgf1", "Fgf2", "Fgf7", "Gdf15", "Gem",
                            "Gmfg", "Hgf", "Hmgb1", "Icam1", "Icam5", "Igf1", "Igfbp1",
                            "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7",
                            "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", "Il6",
                            "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", "Itpka", "Jun", "Lama3",
                            "Lamin1", "Lif", "Mif", "Mmp10", "Mmp12", "Mmp13", "Mmp2", "Mmp3",
                            "Mmp9", "Nap1l4", "Nrg1", "Pecam1", "Plat", "Plau", "Plaur", "Psap",
                            "Ptger2", "Ptges", "Ptgs2", "Scd2", "Selplg", "Serpinb2", "Serpine1",
                            "Sod1", "Spp1", "Stat3", "Timp2", "Tnf", "Tnfrsf1a", "Tnfsf10",
                            "Vcam1", "Vegfa", "Vegfc", "Wnt16", "Wnt2", "Wnt4")

my_geneset = c("Apoe", "App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Ccne1", "Cd26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
    "Cdk4", "Cdk8", "Cdk9", "Cdkn1a", "Cdkn1b", "Cdkn1c",
    "Cdkn2a", "Cdkn2b", "Cdkn2c","Cdkn2d",
    "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
    "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
    "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
    "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
    "Lmna", "Lpl", "Malat1", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
    "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
    "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
    "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vim", "Vtn", "Xcl1")

# senBrainList
#"Apoe",, "Vim", "Malat1"
my_geneset = c("App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                 "Ccne1", "CD26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
                 "Cdk4", "Cdk8", "Cdk9", "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                 "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19",
                 "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
                 "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
                 "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
                 "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
                 "Lmna", "Lpl", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
                 "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
                 "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
                 "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vtn", "Xcl1")
markers <- list(my_geneset = my_geneset)

# Calculate UCell scores for the custom gene set using the SCT assay
thisObj.T <- AddModuleScore_UCell(
  obj = thisObj.T,
  features = markers,
  assay = "SCT", 
  name = "_UCell" 
)

signature.names <- paste0(names(markers), "_UCell")

# Generate the UMAP with GSEA scores overlaid, faceted by age
gsea_umap_plot <- FeaturePlot(
  object = thisObj.T,
  features = signature.names,
  reduction = "umap",
  split.by = "age", # Facet by the 'age' metadata column
  raster = FALSE,
  keep.scale = "all",
  pt.size = 0.5,
  label = TRUE,
  repel = TRUE
)  +
  labs(color = "Enrichment Score") + # Add a label to the color bar
  # Use patchwork to apply a single title to the faceted plots
  patchwork::plot_annotation(
    title = "SenBrain Scores | UCell Module Scores",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )&
  theme(legend.position = "bottom")&
  scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .5)

print(gsea_umap_plot)





original_scores_col <- signature.names[1] 

# Define the name for your new, scaled metadata column
scaled_scores_col <- paste0(original_scores_col, "_scaled")

# 1. Extract the original scores from the Seurat object metadata
original_scores <- thisObj.T[[original_scores_col, drop = TRUE]]

# 2. Apply min-max scaling to transform data to a 0-1 range
scaled_scores <- (original_scores - min(original_scores)) / (max(original_scores) - min(original_scores))

# 3. Add the new 0-1 scaled scores back to the object's metadata
thisObj.T[[scaled_scores_col]] <- scaled_scores





gsea_umap_plot <- FeaturePlot(
  object = thisObj.T,
  features = scaled_scores_col,
  reduction = "umap",
  split.by = "age", # Facet by the 'age' metadata column
  raster = FALSE,
  keep.scale = "all",
  pt.size = 0.5,
  label = TRUE,
  repel = TRUE
)  +
  labs(color = "Enrichment Score") +
  patchwork::plot_annotation(
    title = "SenBrain Scores | Scaled UCell Module Scores",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )&
  theme(legend.position = "bottom")&
  scale_colour_gradient2(low = "grey", mid = "grey", high = "red", midpoint = .55)

print(gsea_umap_plot)




file_name <- paste0(local_WD, "SenScore_ageFacet_gsea_umap_plot", ".png")
ggsave(filename = file_name, plot = gsea_umap_plot, dpi = 900)











VlnPlot(
  object = thisObj.T,
  features = signature.names, # The confidence score we just calculated
  # group.by = "Cluster_Labels",        # The column with your cluster/cell type definitions
  split.by = "age",              # The metadata column to split the violins by
  pt.size = 0,                   # Hides the individual points on the plot
  combine = TRUE                 # Combines the plots into a single ggplot object
) +
  labs( # Use ggplot2 to add better labels
    title = "Cell Type Annotation Confidence by Age",
    y = "UCell Confidence Score",
    x = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
saveRDS(thisObj, paste0(WD,"thisObj_postPCA_v2_11.12.2024"))

# Global UMAP & Cell Typing####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**UMAP NOTES*
#>>dims = dimensions used for embedding the points into the UMAP
#less dims can tighten clustering and remove variance
#but too few dims can ignore/remove important differentiating features
#>>min.dist = minimum distance between data points, lower values can make clusters tighter
#but too high can eliminate important cluster distinctions.
#whereas too low can overly separate data points
#>>spread = distance between points within a cluster
#>>n.neighbors = the "zoom" of the embedding.
#i.e. lower neighbors => highlights local relationships between values, more specific and may reveal more refined sub-clusters but can also over-spread
#recommend 10-20 neighbors for local clusters
# recommend 30-50 for global
# higher neighbors => highlights more global relationships but may obscure more specialized expression within cell types/clusters
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# thisObj <- readRDS(paste0(WD,"thisObj_postPCA_v2_11.12.2024"))
thisObj.T = thisObj #Transfer to Test (.T) so we don't have to keep re-loading the saved object for testing


top_genes <- FindAllMarkers(thisObj.T, only.pos = TRUE) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
head(top_genes)

# Ensure that 'top_genes' is sorted by 'cluster' before creating the factor
top_genes <- top_genes[order(top_genes$cluster), ]

# Create a factor for 'gene' based on its occurrence in 'top_genes' ordered by 'cluster'
top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))

# Create the dot plot of the top 5 genes for each cluster
gene_plot <- ggplot(top_genes, aes(x = factor(cluster), y = gene, color = avg_log2FC, size = pct.1)) +
  geom_point() +
  scale_color_viridis_c() +  # Using viridis color scale for better visibility
  labs(title = "Top 5 Genes per Cluster",
       x = "Cluster",
       y = "Gene",
       color = "Log2 Fold Change",
       size = "Percentage Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(gene_plot)























# >>>>Select dims based on elbow plot per Seurat Vignette<<<<
# Gating variance is important to avoid including noise into cluster segmentation
thisObj.T <-FindNeighbors(thisObj.T, dims = 1:20)
# Select an initial starting resolution
thisObj.T <-FindClusters(thisObj.T, algorithm = 3, resolution = 0.1)
summary(Idents(thisObj.T))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Turn on parallelization to speed up optimization steps
plan("multicore", workers = 12) #>>workers=threads on computer
plan()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#** >>>>1. Initial UMAP Parameter Optimization<<<<**
#Run initial UPO to get rough idea of an optimized UMAP
source("Seurat_UMAP_Optimization.R")
run_flag <- FALSE #added check to avoid accidentally running UPO
if (run_flag) {
UPO_results <- Seurat_UMAP_Optimization(thisObj.T, dims = 1:50, neighbors_range = seq(5, 20, by = 2), min_dist_range = seq(0.001, 0.01, by = 0.001), spread_range = seq(0.2, 1, by = 0.1), downsample_pct = 0.1)
} else {
  print("run_flag is FALSE, skipping the code block")
}

thisObj.T <- RunUMAP(thisObj.T, n.neighbors = 20, dims = 1:50, min.dist = 0.008, spread = 0.5) #values from UPO_results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#** >>>>2. UMAP Cluster k Optimization Fxn<<<<*
# Function to find the optimal number of clusters based on initial UMAP optimization
# Run after initial UMAP run, this function calculates the ideal number of clusters
# based on silhoutte scoring
source("Seurat_UMAP_Cluster_Optimization")
if (run_flag) {
UCO_results <- Seurat_UMAP_Cluster_Optimization(thisObj.T, reduction = "umap", k_range = 2:10, downsample_frac = 0.1)
} else {
  print("run_flag is FALSE, skipping the code block")
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#** >>>>3. UMAP Resolution Optimization Fxn<<<<*
# Function to find the optimal resolution for a target cluster count
# Use after UCO to estimate ideal resolution based on ideal cluster #
source("Seurat_UMAP_Resolution_Optimization.R")
if (run_flag) {
URO_results <- Seurat_UMAP_Resolution_Optimization(thisObj.T, reduction = "umap", optimal_k = 8)
} else {
  print("run_flag is FALSE, skipping the code block")
}
thisObj.T <- URO_results$seurat_obj #push optimized clustering into seurat
#seurat object now contains optimized clustering based on optimized # of clusters
#determined from silhouette scoring
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#** >>>>4. Final UMAP Parameter Optimization<<<<*
#re-run UMAP Parameter Optimization now with the optimized cluster resolution+cluster_k
source("Seurat_UMAP_Optimization.R")
if (run_flag) {
UPO_results <- Seurat_UMAP_Optimization(thisObj.T, dims = 1:50, neighbors_range = seq(5, 50, by = 5), min_dist_range = seq(0.1, 1, by = 0.1), spread_range = seq(0.5, 2.0, by = 0.5), downsample_pct = 0.10)
} else {
  print("run_flag is FALSE, skipping the code block")
}

# Final RunUMAP
thisObj.T <- RunUMAP(thisObj.T, dims = 1:50, min.dist = UPO_results$min.dist, spread = UPO_results$spread, n.neighbors = UPO_results$n.neighbors)

#** >>>> 5. Save Umap<<<<*
saveRDS(thisObj.T, paste0(WD,"seurat_thisObj.T_postOptimizationUMAP_V5_11.14.2024"))

thisObj.T <-readRDS(paste0(WD,"seurat_thisObj.T_postOptimizationUMAP_V5_11.14.2024"))

# Turn off parallelization since it can throw errors for certain fxns
plan("sequential")

#**End Global UMAP Optimization*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# UMAP QC
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
thisObj.T <-readRDS(paste0(WD,"seurat_thisObj.T_postOptimizationUMAP_V5_11.14.2024"))
# Silhouette Calculation and Plotter Fxn
source("SilCal.R")
SilResults <- SilCal(thisObj.T, downsample_pct = 0.5) # 1 = no downsampling

install.packages("dimRed")
library(dimRed)
# Assume 'data' is the original data and 'embedding' is the low-dimensional embedding
umap_embedding <- Embeddings(thisObj.T, reduction = "umap")
trust <- trustworthiness(data, embedding, k = 15)

DimPlot(thisObj.T, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) + ggtitle("UMAP")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Global Top Genes
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Extract top 5 genes
top_genes <- FindAllMarkers(thisObj.T, only.pos = TRUE) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
head(top_genes)

# Ensure that 'top_genes' is sorted by 'cluster' before creating the factor
top_genes <- top_genes[order(top_genes$cluster), ]

# Create a factor for 'gene' based on its occurrence in 'top_genes' ordered by 'cluster'
top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))

# Create the dot plot of the top 5 genes for each cluster
gene_plot <- ggplot(top_genes, aes(x = factor(cluster), y = gene, color = avg_log2FC, size = pct.1)) +
  geom_point() +
  scale_color_viridis_c() +  # Using viridis color scale for better visibility
  labs(title = "Top 5 Genes per Cluster",
       x = "Cluster",
       y = "Gene",
       color = "Log2 Fold Change",
       size = "Percentage Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(gene_plot)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Global Cell Typing
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Define gene lists
gene_lists <- list(
  mgList = c('C1qa','C1qb','C1qc','Ctss','Csf1r','Cst3','Cx3cr1','Hexb','Itgam','P2ry12','Selplg','Tmem119'),
  astroList = c('Gja1','Gpr37l1','Plpp3','Mfge8','Gfap','Slc1a3','Bcan'),
  gabaergicList = c('Gad1','Gad2','Npy','Slc6a1','Slc32a1','Sst'),
  neuronList = c('Adcy1','Atp1b1','Camk2a','Dnm1','Epha6','Gria1', 'Grin2a','Grin2b','Nell2','Nrgn','Snap25','Syp','Syt1','Ywhag'),
  oligoList = c('Plp1','Mog','Mag','Mobp','Mbp','Cryab','Fa2h'),
  opcList = c('Cntn1','Cspg5','Gpr17','Pdgfra', 'Ptprz1','Serpine2'),
  vascList = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg'),
  vlmcList = c('Apod','Col1a1','Dcn','Il33','Ptgds','Ptn'),
  epenList = c('Dbi','Foxj1','Nnat','Tuba1a','Pltp','Vim'),
  periList = c('Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn'),
  nstemnblastList = c('Ctnna2','H3f3b','Nnat','Ptprs','Tcf4','Tuba1a','Tubb5'),
  chorplexList = c('Ttr','Htr2c','Lgals3bp','Cab39l','Igfbp7','Ppp1r1b'),
  vascList2 = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg','Apod','Col1a1','Dcn','Il33','Ptgds','Ptn',
                'Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn')
)

# Loop through each gene list and create a dot plot
for (list_name in names(gene_lists)) {
  genes <- gene_lists[[list_name]]
  
  # Generate and display the dot plot
  dot_plot <- DotPlot(thisObj.T, features = genes, cols = c('grey','black')) +
    ggtitle(paste("Dot Plot for", list_name)) +
    theme_minimal()
  
  print(dot_plot)  # Display each plot
}

#0 and 1 contain neurons
#1 has some gabaergic neurons
#2 is oligos
#3 is vascular and contains pericytes
#4 astrocyte like cells?
#5, astrocytes
#6 microglia
#7 OPCs
#8 has chorplex
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Global Cell type Labeling UMAP Clusters
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Define new labels for each cluster
new_cluster_labels <- c(
  "(0) Neurons", 
  "(1) Interneurons",
  "(2) Oligodendrocytes",
  "(3) Vascular/Pericyte",
  "(4) Purkinje",
  "(5) Astrocytes",
  "(6) MG",
  "(7) OPCs",
  "(8) Choroid Plexus"
)

names(new_cluster_labels) <- levels(thisObj.T)
thisObj.T <- RenameIdents(thisObj.T, new_cluster_labels)
DimPlot(thisObj.T, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1) + NoLegend()

saveRDS(thisObj.T, paste0(WD, "seurat_thisObj.T_UMAP_V5_cellTyped_11.14.2024"))

thisObj.T <- readRDS(paste0(local_WD, "seurat_thisObj.T_UMAP_V5_cellTyped_11.14.2024"))
raw_Obj <- thisObj.T # adding as a save space to avoid re-loading....
# SubCluster UMAP & Cell-Typing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Local Cluster Cell Typing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data
thisObj.T <- readRDS(paste0(local_WD, "seurat_thisObj.T_UMAP_V5_cellTyped_11.14.2024"))


# Select the global cluster label to then sub-cluster
source("saveRDSObject.R")
source("SubClusterCellTyping_Singles.R")

# Install for a faster clustering algorithim
install.packages('devtools')
devtools::install_github('immunogenomics/presto')

#0000000000000000000000000000000000000000000000000000000000000000000000000000000
# (0) Neuron subcluster
#0000000000000000000000000000000000000000000000000000000000000000000000000000000
cluster_label = "(0) Neurons"
SubC_0_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.05,
  k.param = 25,
  spread = 1,
  min.dist = 0.001
)
# Dims used: 15
top_gene_list <- SubC_0_Data$all_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  summarise(all_genes = list(gene)) %>%
  pull(all_genes)

print(top_gene_list)

# Load libraries
library(dplyr)
library(ggplot2)

my_genes <- list(c("Zic1", "Neurod1", "Gabra6", "Pcp4"))
SubC_0_Data$all_genes[SubC_0_Data$all_genes$gene %in% my_genes, ]

DotPlot(SubC_0_Data$subset_obj, features = "Pcp4", cols = c('grey','black'))
# 
# Dims used: 15
new_cluster_labels_0 <- c(
  "(1) Granule Cells",
  "(2) UBCs"
)

# Rename Subcluster clusters
names(new_cluster_labels_0) <- levels(SubC_0_Data$subset_obj)
subset_obj_C0 <- RenameIdents(SubC_0_Data$subset_obj, new_cluster_labels_0)
saveRDSObject(obj = SubC_0_Data, file_label = "SubC_0_Data", directory = local_WD) 

DimPlot(subset_obj_C0, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)

source("SubClusterCellTyping_Singles.R")
# Sub-UBC to distinguish excitatory vs inhibitory
cluster_label = "(2) UBCs"
SubC_0_UBCsub_Data <- SubClusterCellTyping_Singles(
  obj = subset_obj_C0,
  cluster_label = cluster_label,
  resolution = 0.3,
  k.param = 10,
  spread = 3,
  min.dist = 0.0001
)
top_gene_list <- SubC_0_UBCsub_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)
print(top_gene_list)





#1111111111111111111111111111111111111111111111111111111111111111111111111111111
# (1) Interneuron subcluster
#1111111111111111111111111111111111111111111111111111111111111111111111111111111
cluster_label = "(1) Interneurons"
SubC_1_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 20,
  spread = 1,
  min.dist = 0.1
)
# new_cluster_labels_1 <- c(
#   "(0) Inhibitory Golgi",
#   "(1) Excitatory Golgi",
#   "(2) Glial-associated Golgi",
#   "(3) Basket/Stellate"
# )

top_gene_list <- SubC_1_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

#Dims used: 13
new_cluster_labels_1 <- c(
  "Interneuron_1",
  "Interneuron_2",
  "Interneuron_3",
  "Interneuron_4"
)
# Rename Subcluster clusters
names(new_cluster_labels_1) <- levels(SubC_1_Data$subset_obj)
subset_obj_C1 <- RenameIdents(SubC_1_Data$subset_obj, new_cluster_labels_1)
saveRDSObject(obj = SubC_1_Data, file_label = "SubC_1_Data", directory = local_WD) 

DimPlot(subset_obj_C1, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)

#2222222222222222222222222222222222222222222222222222222222222222222222222222222
#(2) Oligo subcluster
#2222222222222222222222222222222222222222222222222222222222222222222222222222222
cluster_label = "(2) Oligodendrocytes"
SubC_2_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 20,
  spread = 1,
  min.dist = 0.1
)
# new_cluster_labels_2 <- c(
#   "(0) Myelinating Oligodendrocytes",
#   "(1) Bergmann glia-like cells",
#   "(2) Synaptogenic Oligodendrocytes"
# )

top_gene_list <- SubC_2_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

# Dims used?
new_cluster_labels_2 <- c(
  "OLG_1",
  "OLG_2",
  "OLG_3"
)
# Rename Subcluster clusters
names(new_cluster_labels_2) <- levels(SubC_2_Data$subset_obj)
subset_obj_C2 <- RenameIdents(SubC_2_Data$subset_obj, new_cluster_labels_2)
saveRDSObject(obj = SubC_2_Data, file_label = "SubC_2_Data", directory = local_WD) 

DimPlot(subset_obj_C2, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)

#3333333333333333333333333333333333333333333333333333333333333333333333333333333
#(3) Vascular subcluster
#3333333333333333333333333333333333333333333333333333333333333333333333333333333
cluster_label = "(3) Vascular/Pericyte"
SubC_3_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 20,
  spread = 1,
  min.dist = 0.1
)
# new_cluster_labels_3 <- c(
#   "(0) BBB-ECs",
#   "(1) Neuron-associated Pericytes",
#   "(2) Fibroblast-like Pericytes",
#   "(3) Synapse-associated Vascular Cells",
#   "(4) vSMCs"
# )
top_gene_list <- SubC_3_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

# Dims used: 15
new_cluster_labels_3 <- c(
  "VSMCs/PCs_1",
  "VSMCs/PCs_2",
  "VSMCs/PCs_3",
  "VSMCs/PCs_4",
  "VSMCs/PCs_5"
)
# Rename Subcluster clusters
names(new_cluster_labels_3) <- levels(SubC_3_Data$subset_obj)
subset_obj_C3 <- RenameIdents(SubC_3_Data$subset_obj, new_cluster_labels_3)
saveRDSObject(obj = SubC_3_Data, file_label = "SubC_3_Data", directory = local_WD) 

DimPlot(subset_obj_C3, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#4444444444444444444444444444444444444444444444444444444444444444444444444444444
#(4) Purkinje subcluster
#4444444444444444444444444444444444444444444444444444444444444444444444444444444
cluster_label = "(4) Purkinje"
SubC_4_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 20,
  spread = 1,
  min.dist = 0.1
)
# new_cluster_labels_4 <- c(
#   "(0) Neurotrophic Sub-type",
#   "(1) Canonical Purkinje",
#   "(2) Myelin-interactive",
#   "(3) Excitatory Sub-type",
#   "(4) Oligodendrocyte Associated"
# )
top_gene_list <- SubC_4_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

# Dims used: 15
new_cluster_labels_4 <- c(
  "Purkinje_1",
  "Purkinje_2",
  "Purkinje_3",
  "Purkinje_4"
)
# Rename Subcluster clusters
names(new_cluster_labels_4) <- levels(SubC_4_Data$subset_obj)
subset_obj_C4 <- RenameIdents(SubC_4_Data$subset_obj, new_cluster_labels_4)
saveRDSObject(obj = SubC_4_Data, file_label = "SubC_4_Data", directory = local_WD) 

DimPlot(subset_obj_C4, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#5555555555555555555555555555555555555555555555555555555555555555555555555555555
#(5) Astrocyte subcluster
#5555555555555555555555555555555555555555555555555555555555555555555555555555555
cluster_label = "(5) Astrocytes"
SubC_5_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 20,
  spread = 1,
  min.dist = 0.1
)
# new_cluster_labels_5 <- c(
#   "(0) Neuroprotective Sub-type",
#   "(1) Inflammatory",
#   "(2) Neuron-associated",
#   "(3) CSF-associated"
# )
top_gene_list <- SubC_5_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

#Dims used: 
new_cluster_labels_5 <- c(
  "Astrocyte_1",
  "Astrocyte_2",
  "Astrocyte_3",
  "Astrocyte_4"
)
# Rename Subcluster clusters
names(new_cluster_labels_5) <- levels(SubC_5_Data$subset_obj)
subset_obj_C5 <- RenameIdents(SubC_5_Data$subset_obj, new_cluster_labels_5)
saveRDSObject(obj = SubC_5_Data, file_label = "SubC_5_Data", directory = local_WD) 

DimPlot(subset_obj_C5, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#6666666666666666666666666666666666666666666666666666666666666666666666666666666
#(6) MG subcluster
#6666666666666666666666666666666666666666666666666666666666666666666666666666666
cluster_label = "(6) MG"
SubC_6_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.6,
  k.param = 5,
  spread = 3,
  min.dist = 0.001
)

# new_cluster_labels_6 <- c(
#   "(0) Myelin/Phagocytic-associated",
#   "(1) Neuron-associated",
#   "(2) Astrocyte-associated"
# )

# Dims used: 10
top_gene_list <- SubC_6_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

new_cluster_labels_6 <- c(
  "MG_1_Lgals3+",
  "MG_2",
  "MG_3",
  "MG_4",
  "MG_5",
  "MG_6",
  "MG_7",
  "MG_8",
  "MG_9",
  "MG_10"
)
# Rename Subcluster clusters
names(new_cluster_labels_6) <- levels(SubC_6_Data$subset_obj)
subset_obj_C6 <- RenameIdents(SubC_6_Data$subset_obj, new_cluster_labels_6)
saveRDSObject(obj = SubC_6_Data, file_label = "SubC_6_Data", directory = local_WD) 

DimPlot(subset_obj_C6, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#7777777777777777777777777777777777777777777777777777777777777777777777777777777
# (7) OPC subcluster
#7777777777777777777777777777777777777777777777777777777777777777777777777777777
cluster_label = "(7) OPCs"
SubC_7_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 10,
  spread = 2,
  min.dist = 0.01
)
top_gene_list <- SubC_7_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

# new_cluster_labels_7 <- c(
#   "(0) Pre-synaptic",
#   "(1) Immune-associated",
#   "(2) Neuron-associated",
#   "(3) Committed Pre-myelinating"
# )
#Dims used:15
new_cluster_labels_7 <- c(
  "OPC_1",
  "OPC_2",
  "OPC_3",
  "OPC_4"
)
# Rename Subcluster clusters
names(new_cluster_labels_7) <- levels(SubC_7_Data$subset_obj)
subset_obj_C7 <- RenameIdents(SubC_7_Data$subset_obj, new_cluster_labels_7)
saveRDSObject(obj = SubC_7_Data, file_label = "SubC_7_Data", directory = local_WD) 

DimPlot(subset_obj_C7, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#8888888888888888888888888888888888888888888888888888888888888888888888888888888
# (8) Choroid Plexus subcluster
#8888888888888888888888888888888888888888888888888888888888888888888888888888888
cluster_label = "(8) Choroid Plexus"
SubC_8_Data <- SubClusterCellTyping_Singles(
  obj = thisObj.T,
  cluster_label = cluster_label,
  resolution = 0.1,
  k.param = 5,
  spread = 2,
  min.dist = 0.01
)

top_gene_list <- SubC_8_Data$top_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(top_genes = list(gene)) %>%
  pull(top_genes)

print(top_gene_list)

# new_cluster_labels_8 <- c(
#   "(0) Homeostatic",
#   "(1) Neurotransmitter Signaling",
#   "(2) Immune-responsive"
# )

#Dims used: 15
new_cluster_labels_8 <- c(
  "CP_1",
  "CP_2",
  "CP_3"
)

# Rename Subcluster clusters
names(new_cluster_labels_8) <- levels(SubC_8_Data$subset_obj)
subset_obj_C8 <- RenameIdents(SubC_8_Data$subset_obj, new_cluster_labels_8)
saveRDSObject(obj = SubC_8_Data, file_label = "SubC_8_Data", directory = local_WD) 

DimPlot(subset_obj_C8, reduction = "umap", raster = FALSE, label = TRUE, pt.size = 0.1)+ggtitle(paste0("Sub-Cluster UMAP | ", cluster_label))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

subset_objs_list <- list(
  "Neurons"       = subset_obj_C0,  # Neurons
  "Interneurons"  = subset_obj_C1,  # Interneurons
  "OLG"           = subset_obj_C2,  # OLG
  "VSCMCs+PCs"    = subset_obj_C3,  # VSMCs/PCs
  "Purkinje"      = subset_obj_C4,  # Purkinje
  "Astrocyte"     = subset_obj_C5,  # Astrocyte
  "MG"            = subset_obj_C6,  # MG
  "OPC"           = subset_obj_C7,  # OPC
  "PC"            = subset_obj_C8   # CP
)

# Function to efficiently update sub-cluster labels in the global Seurat object
# Re-combines the subset clusters back into the main object for easier saving/data access
source("update_subcluster_labels.R")

thisObj.T <- update_subcluster_labels_iterative(thisObj.T, subset_objs_list)
unique(thisObj.T@meta.data$subtype)
saveRDSObject(obj = thisObj.T, file_label = "Cerebellum_CosmonautObj_SubClusterLabeled_2.24.2025", directory = local_WD) 

#**End Sub-Cluster Cell-Typing*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Analyzing Global & Sub-Clusters ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Load objects and sub objects*
# Notes: thisObj.T contains the correct subtype labels, but there is a weird
# error when re-subsetting the main object by the Idents(). The subtype column  
# is returned empty for some reason to get around this, I just load the 
# subobject directly, but in the future it would be easier to just subset
# directly from  the main seurat object (thisObj)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisObj <- readRDS(paste0(local_WD, "Cerebellum_CosmonautObj_SubClusterLabeled_2.24.2025.rds"))
thisObj.T <- thisObj

SubC_0_Data <- readRDS(paste0(local_WD, "SubC_0_Data.rds"))
SubC_1_Data <- readRDS(paste0(local_WD, "SubC_1_Data.rds"))
SubC_2_Data <- readRDS(paste0(local_WD, "SubC_2_Data.rds"))
SubC_3_Data <- readRDS(paste0(local_WD, "SubC_3_Data.rds"))
SubC_4_Data <- readRDS(paste0(local_WD, "SubC_4_Data.rds"))
SubC_5_Data <- readRDS(paste0(local_WD, "SubC_5_Data.rds"))
SubC_6_Data <- readRDS(paste0(local_WD, "SubC_6_Data.rds"))
SubC_7_Data <- readRDS(paste0(local_WD, "SubC_7_Data.rds"))
SubC_8_Data <- readRDS(paste0(local_WD, "SubC_8_Datards"))


subset_objs_list <- list(
  "Neurons"       = SubC_0_Data$subset_obj,  # Neurons
  "Interneurons"  = SubC_1_Data$subset_obj,  # Interneurons
  "OLG"           = SubC_2_Data$subset_obj,  # OLG
  "VSCMCs+PCs"    = SubC_3_Data$subset_obj,  # VSMCs/PCs
  "Purkinje"      = SubC_4_Data$subset_obj,  # Purkinje
  "Astrocyte"     = SubC_5_Data$subset_obj,  # Astrocyte
  "MG"            = SubC_6_Data$subset_obj,  # MG
  "OPC"           = SubC_7_Data$subset_obj,  # OPC
  "PC"            = SubC_8_Data$subset_obj   # CP
)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Generate basic global UMAP
p <- DimPlot(thisObj.T, reduction = "umap", raster = TRUE, label = TRUE, repel = TRUE) +
  ggtitle("Global UMAP") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p)
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Generate subtype global UMAP
n_colors <- length(unique(thisObj.T@meta.data$subtype))
if(n_colors <= 12){
  distinct_colors <- brewer.pal(n_colors, "Set1")
} else {
  distinct_colors <- colorRampPalette(brewer.pal(12, "Set1"))(n_colors)
}

p <- DimPlot(thisObj.T, reduction = "umap", group.by = "subtype", 
             raster = TRUE, label = TRUE, repel = TRUE) +
  scale_color_manual(values = distinct_colors) +
  ggtitle("Global UMAP | Grouped by Subtype") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p)
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Global UMAP by Age
grouping <- 'age'
p <- DimPlot(thisObj.T, reduction = "umap",
             group.by = grouping, 
             label = FALSE,
             raster = TRUE,
             pt.size = 1
) +ggtitle(paste("Global UMAP | Group by ", grouping))
print(p)
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Generate sub-cluster UMAPs by Age (or any other variable)
cell_type = "Neurons"     
cell_type = "Interneurons" 
cell_type = "OLG"          
cell_type = "VSCMCs+PCs"  
cell_type = "Purkinje"    
cell_type = "Astrocyte"  
cell_type = "MG"   
cell_type = "OPC" 
cell_type = "PC"   

sub_obj <- subset_objs_list[[cell_type]]
sub_obj <- SubC_6_Data$subset_obj
p_1a <- DimPlot(sub_obj, reduction = "umap",
             label = FALSE,
             raster = TRUE,
             pt.size = 2
) +
  ggtitle(paste("Sub-cluster UMAP:", cell_type))
# print(p_1a)
p_1b <- DimPlot(sub_obj, reduction = "umap",
                group.by = grouping,
                label = FALSE,
                raster = TRUE,
                pt.size = 2
) +
  scale_color_manual(values = c("gold", "purple")) +  # Custom colors
  ggtitle(paste("Sub-cluster UMAP:", cell_type, "Group by ", grouping))

# print(p_1b)
p_1a + p_1b


# BATCH Method
grouping <- 'age'
for(cell_type in names(subset_objs_list)) {
  sub_obj <- subset_objs_list[[cell_type]]
  print(paste("Processing", cell_type))
  
  # Example plotting code; adjust as needed for your analysis:
  p <- DimPlot(sub_obj, reduction = "umap",
               group.by = grouping, 
               label = FALSE,
               raster = TRUE,
               pt.size = 2
               ) +
    ggtitle(paste("Sub-cluster UMAP:", cell_type, "Group by ", grouping))
  print(p)
  
}
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Check if a specific gene is within our object
gene_to_check <- "Lpl"  # Change to your gene of interest

gene_exists <- gene_to_check %in% rownames(sub_obj)
if (gene_exists) {
  cat(gene_to_check, "is present in the Seurat object.\n")
} else {
  cat(gene_to_check, "is NOT found in the Seurat object.\n")
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FeaturePlot(subset_objs_list, features = "Lgals3", cols = c("gray", "yellow", "red"))
FeaturePlot(subset_objs_list$MG, features = c("Lgals3", "Lpl"), blend = TRUE, split.by = "age")

FeaturePlot(subset_obj, features = "Average_GAL3enrichedModule_Expr", pt.size = 1, shape.by = "age") +
  ggtitle("Average Expression of Genes in GAL3enrichedModule")+
  scale_color_viridis(option = "viridis", name = "Expression")+
  scale_shape_manual(values = c("YNG" = 1, "OLD" = 15)) + # Customize shapes
  theme_minimal()

FeatureScatter(subset_objs_list$MG, feature1 = "Lgals3", feature2 = "age") +
  ggtitle("Lgals3 Expression vs Age")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Single Sample Spatial Plots + Plotly ####
# Convert coordinates from mm to microns
thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000


# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# SUBSET DATA
# Subset the Seurat object by slide and sample
subset_data <- subset(thisObj.T, subset = mus == sample_name)
# subset_data <- RenameIdents(subset_data, new_cluster_labels)
# Get the min and max values of the x and y coordinates
x_range <- range(subset_data@meta.data$x_slide_micron)
y_range <- range(subset_data@meta.data$y_slide_micron)

# Define the plot limits with a small buffer to avoid clipping points
buffer <- 10
xlim_vals <- c(min(x_range) - buffer, max(x_range) + buffer)
ylim_vals <- c(min(y_range) - buffer, max(y_range) + buffer)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# GLOBAL CELL-TYPE SPATIAL PLOT
# Create the spatial plot w/ global cell-types
ggplot(subset_data@meta.data, aes(x = x_slide_micron, y = y_slide_micron, color = Idents(subset_data))) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = paste("Sample:", sample_name),
       x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = "Cell Type") +
  theme_minimal() +
  coord_fixed() +
  xlim(xlim_vals) +
  ylim(ylim_vals)+
  theme(
  panel.grid = element_blank(),
  axis.text = element_blank())    # Hide axis ticks
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# LOCAL SUB-TYPE SPATIAL PLOT
# Create the spatial plot w/ local cell-types
p <- ggplot(subset_data@meta.data, aes(x = x_slide_micron, y = y_slide_micron, color = subset_data$subtype)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(title = paste("Sample:", sample_name),
       x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = "Cell Type") +
  theme_minimal() +
  coord_fixed() +
  xlim(xlim_vals) +
  ylim(ylim_vals)+
  theme(
    panel.grid = element_blank(),   
    axis.text = element_blank())+ 
    theme(legend.position = "bottom")
print(p)

p_interactive <- ggplotly(p, tooltip = "text")
p_interactive
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# CREATE A SPATIAL PLOT W/ HIGHLIGHTED SPECIFIC CELL SUBTYPES
# Function to highlight a specific subtype
# Fxn can be modified to highlight other features too in the spatial map
library(crayon)
cat(green("View available subtypes:\n"))
unique(subset_data@meta.data$subtype)

highlight_names <- c("Purkinje_4", "Purkinje_3", "Purkinje_2", "Purkinje_1")
highlight_colors <- c("blue", "green", "orange", "red")

highlight_subtype <- function(data, highlight_names, highlight_colors = "red", default_color = "gray") {
  # If a single color is provided for multiple highlight names, replicate it.
  if (length(highlight_colors) == 1 && length(highlight_names) > 1) {
    highlight_colors <- rep(highlight_colors, length(highlight_names))
  }
  
  # Check that the length of highlight_colors matches highlight_names
  if (length(highlight_colors) != length(highlight_names)) {
    stop("Length of highlight_colors must be either 1 or equal to the number of highlight_names.")
  }
  
  # Create a mapping from highlight names to their colors
  color_mapping <- setNames(highlight_colors, highlight_names)
  
  # Set the default color for all cells
  data$highlight_color <- default_color
  
  # Assign the corresponding highlight color for cells whose subtype matches any highlight name
  for (name in highlight_names) {
    data$highlight_color[data$subtype == name] <- color_mapping[[name]]
  }
  
  return(data)
}
subset_data@meta.data <- highlight_subtype(subset_data@meta.data, highlight_names, highlight_colors)

ggplot(subset_data@meta.data, aes(x = x_slide_micron, y = y_slide_micron, color = highlight_color)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_identity() +  # Use the colors as specified, without additional scaling
  labs(
    title = paste("Sample:", sample_name, "- Highlighted:", paste(highlight_names, collapse = ", ")),
    x = "X Coordinate (microns)",
    y = "Y Coordinate (microns)"
  ) +
  theme_minimal() +
  coord_fixed() +
  xlim(xlim_vals) +
  ylim(ylim_vals) +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.text = element_blank()     # Hide axis ticks
  )

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# HIGHLIGHT LOCALIZATION OF SPECIFIED GENES
# Gene of interest
QueryList = "Apoe" #Must be exact name match for gene name
# Extract gene expression data for the gene of interest
expression_data <- FetchData(subset_data, vars = QueryList)

# Add the expression data as a new column in the metadata
subset_data@meta.data$expression_level <- expression_data[[QueryList]]

# Create the plot with color gradient based on the expression levels
ggplot(subset_data@meta.data, aes(x = x_slide_micron, y = y_slide_micron, color = expression_level)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_gradient(low = "gray", high = "red") +  # Heatmap gradient from low to high expression
  labs(title = paste("Sample:", sample_name, "- Expression of", QueryList),
       x = "X Coordinate (microns)", y = "Y Coordinate (microns)", color = paste(QueryList, "Expression")) +
  theme_minimal() +
  coord_fixed() +
  xlim(xlim_vals) +
  ylim(ylim_vals) +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.text = element_blank()    # Hide axis ticks
  )

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Banksy Main ####
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Banksy w/ Clustering by only known cell markers*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Install missing packages if necessary
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("Raphael-LS/Banksy")
}
if (!requireNamespace("scuttle", quietly = TRUE)) install.packages("scuttle")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

# Load all necessary libraries
library(Banksy)          
library(SummarizedExperiment)  
library(SpatialExperiment)  
library(scuttle)        
library(scater)        
library(cowplot)        
library(ggplot2)        

thisObj <- readRDS(file.path(project_WD, "COSMONAUT_CB_ONLY-Seurat_PostQC_v3.RDS"))
thisObj.T <- thisObj

# Define unique slide offsets to prevent overlapping coordinates
slide_offsets <- list(
  "305" = c(0, 0),
  "YM04" = c(5000, 0),
  "304" = c(10000, 0),
  "YM03" = c(15000, 0),
  "YM02" = c(20000, 0),
  "302" = c(0, 5000),
  "YM01" = c(5000, 5000),
  "301" = c(10000, 5000),
  "YF04" = c(15000, 5000),
  "243" = c(20000, 5000),
  "YF03" = c(0, 10000),
  "239" = c(5000, 10000),
  "YF02" = c(10000, 10000),
  "237" = c(15000, 10000),
  "YF01" = c(20000, 10000),
  "236" = c(0, 15000)
)

# Convert mm to microns for consistent spatial scaling
thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000

# Apply unique offsets per slide/sample
thisObj.T@meta.data$x_slide_global <- thisObj.T@meta.data$x_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][1])

thisObj.T@meta.data$y_slide_global <- thisObj.T@meta.data$y_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][2])

# Verify transformed coordinates
head(thisObj.T@meta.data[, c("mus", "x_slide_global", "y_slide_global")])

# Extract corrected spatial coordinates
locs <- as.matrix(thisObj.T@meta.data[c("x_slide_global", "y_slide_global")])

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
gene_lists <- list(
  mgList = c('C1qa','C1qb','C1qc','Ctss','Csf1r','Cst3','Cx3cr1','Hexb','Itgam','P2ry12','Selplg','Tmem119'),
  astroList = c('Gja1','Gpr37l1','Plpp3','Mfge8','Gfap','Slc1a3','Bcan'),
  gabaergicList = c('Gad1','Gad2','Npy','Slc6a1','Slc32a1','Sst'),
  neuronList = c('Adcy1','Atp1b1','Camk2a','Dnm1','Epha6','Gria1', 'Grin2a','Grin2b','Nell2','Nrgn','Snap25','Syp','Syt1','Ywhag'),
  oligoList = c('Plp1','Mog','Mag','Mobp','Mbp','Cryab','Fa2h'),
  opcList = c('Cntn1','Cspg5','Gpr17','Pdgfra', 'Ptprz1','Serpine2'),
  vascList = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg','Apod','Col1a1','Dcn','Il33','Ptgds','Ptn',
               'Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn'), #combined vascList and vascList2
  epenList = c('Dbi','Foxj1','Nnat','Tuba1a','Pltp','Vim'),
  periList = c('Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn'),
  nstemnblastList = c('Ctnna2','H3f3b','Nnat','Ptprs','Tcf4','Tuba1a','Tubb5'),
  chorplexList = c('Ttr','Htr2c','Lgals3bp','Cab39l','Igfbp7','Ppp1r1b'),
  GAL3enrichedModule = c('Apoe','Cd9','Clec7a','Ctsb','Ctsd','Fth1','Gpnmb','Itgax',
                         'Lyz1','Lilrb4a','Mmp12','Spp1','Tmsb4x','Tyrobp'),
  my_marker_genes = c(
    "Rbfox3", "Snap25", "Syp", "Stmn2",  # Pan-Neuronal
    "Calb2", "Grm1", "Plcb4",  # UBC-specific
    "Slc17a7", "Grm4", "Nhlh1", "Ebf2", "Cbln3",  # Granule Cell Alternative Markers
    "Gad1", "Gad2", "Slc6a1", "Htr3a", "Nrgn",  # Generic GABAergic
    "Pvalb", "Crh", "Grm2", "Npy", "Htr2a",  # Specific for Golgi, Stellate, Basket Interneurons
    "Calb1", "Pcp2", "Grid2", "Itpr1", "Rora",  # Purkinje Cells
    "Gfap", "S100b", "Aqp4", "Slc1a3", "Sox9",  # Bergmann Glia
    "Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum",  # Fibroblasts
    "P2ry12", "Tmem119", "Cx3cr1", "C1qa", "Trem2",  # Microglia
    "Mbp", "Plp1", "Mog", "Mag",  # Oligodendrocytes
    "Pdgfra", "Olig2", "Gpr17", "Tcf7l2",  # Oligodendrocyte Precursor Cells (OPCs)
    "Cldn5", "Pecam1", "Flt1", "Esam",  # Endothelial Cells
    "Gfap", "Aldh1l1", "Slc1a3", "Aqp4", "S100b",  # Astrocytes
    "Pdgfrb", "Rgs5", "Des", "Acta2",  # Pericytes
    "Acta2", "Tagln", "Notch3",  # Vascular Smooth Muscle Cells (VSMCs)
    "Ttr", "Cldn5", "Otx2", "Kcne2", "Aqp1" # Choroid Plexus
  )
)
gcm <- thisObj.T@assays$RNA@counts
all_genes_to_keep <- unique(unlist(gene_lists))

# Get the gene names present in the gcm matrix
genes_in_gcm <- rownames(gcm)

# Find the intersection of genes to keep and genes in gcm
genes_to_filter <- intersect(all_genes_to_keep, genes_in_gcm)

# Filter the gcm matrix to only include the genes in genes_to_filter
gcm_filtered <- gcm[genes_to_filter, ]

keep_cells <- colSums(gcm_filtered) > 10 # or some other threshold
gcm_filtered <- gcm_filtered[, keep_cells]
locs <- locs[keep_cells,]
thisObj.T@meta.data <- thisObj.T@meta.data[keep_cells,]

# Create SpatialExperiment object with spatial data
se <- SpatialExperiment(
  assay = list(counts = gcm_filtered), 
  spatialCoords = locs,
  colData = DataFrame(thisObj.T@meta.data)  
)

# Store slide/sample identity explicitly
colData(se)$slide_id <- thisObj.T@meta.data$mus

# Define parameters
lambda <- c(0,0.1,0.2)
res <- 0.25
k_geom <- c(15, 30)    

# Normalize counts
se <- computeLibraryFactors(se)
assay(se, aname) <- normalizeCounts(se, log = FALSE)

# Compute Banksy features per slide
se_list <- lapply(unique(colData(se)$slide_id), function(slide) {
  se_subset <- se[, colData(se)$slide_id == slide]  
  se_subset <- Banksy::computeBanksy(se_subset, assay_name = "counts" , compute_agf = TRUE, k_geom = k_geom)
  return(se_subset)
})

# Merge processed slides
se <- do.call(cbind, se_list)

# Ensure reproducibility
set.seed(1000)

# Run Banksy PCA
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)

# Run Banksy UMAP
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

# Perform clustering
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = res)

# Extract spatial coordinates
coords <- spatialCoords(se)
colData(se)$sdimx <- coords[, 1]
colData(se)$sdimy <- coords[, 2]

# Get unique sample names
unique_samples <- unique(colData(se)$slide_id)

# View all dim reductions
names(reducedDims(se))

# Identify Banksy clustering column dynamically
clustering_columns <- grep("^clust", colnames(colData(se)), value = TRUE)
print(clustering_columns)  # Debug: View available cluster labels

if (length(clustering_columns) == 0) {
  stop("Error: No clustering columns found in colData(se). Ensure Banksy clustering completed successfully.")
}

# Select first clustering column (most relevant)
clust_column <- clustering_columns[5]  
print(paste("Using clustering column:", clust_column))  

# Identify Banksy's UMAP embedding dynamically
umap_columns <- grep("^UMAP", names(reducedDims(se)), value = TRUE)
umap_columns

if (length(umap_columns) == 0) {
  stop("Error: No UMAP embedding found in reducedDims(se). Ensure Banksy UMAP completed successfully.")
}

umap_dimred <- umap_columns[1]  # Select first available UMAP embedding
umap_dimred

# Banksy Full UMAP plot
plot_umap <- plotReducedDim(se, 
                            dimred = umap_dimred,
                            point_size = 0.1,
                            colour_by = clust_column) +
  ggtitle("Banksy UMAP") +
  theme_minimal()+
  theme(legend.text = element_text(size = 12))

print(plot_umap)

file_name <- paste0(WD, "banksy_se_object.rds")
saveRDS(se, file = file_name)
# BankObj <- readRDS(file_name)
# se <- BankObj
# Loop through each sample and generate spatial + UMAP plots
for (sample in unique_samples) {
  se_subset <- se[, colData(se)$slide_id == sample]
  
  # Spatial reconstruction plot
  plot_spatial <- plotColData(se_subset,
                              x = "sdimx", y = "sdimy",
                              point_size = 0.6, colour_by = clust_column) +
    ggtitle(paste("Spatial Reconstruction -", sample)) +
    theme_minimal()
  
  # Banksy UMAP plot
  plot_umap <- plotReducedDim(se_subset,
                              dimred = umap_dimred,
                              colour_by = clust_column) +
    ggtitle(paste("Banksy UMAP -", sample)) +
    theme_minimal()
  
  # Arrange plots side by side
  print(plot_grid(plot_spatial, plot_umap, ncol = 2))
}


# se <- readRDS(paste0(WD, "banksy_se_object.rds"))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**Banksy & Harmony Testing Main*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Install missing packages if necessary
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("Raphael-LS/Banksy")
}
if (!requireNamespace("scuttle", quietly = TRUE)) install.packages("scuttle")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

# Load all necessary libraries
library(Banksy)          
library(SummarizedExperiment)  
library(SpatialExperiment)  
library(scuttle)        
library(scater)        
library(cowplot)        
library(ggplot2)        

# Define unique slide offsets to prevent overlapping coordinates
slide_offsets <- list(
  "305" = c(0, 0),
  "YM04" = c(5000, 0),
  "304" = c(10000, 0),
  "YM03" = c(15000, 0),
  "YM02" = c(20000, 0),
  "302" = c(0, 5000),
  "YM01" = c(5000, 5000),
  "301" = c(10000, 5000),
  "YF04" = c(15000, 5000),
  "243" = c(20000, 5000),
  "YF03" = c(0, 10000),
  "239" = c(5000, 10000),
  "YF02" = c(10000, 10000),
  "237" = c(15000, 10000),
  "YF01" = c(20000, 10000),
  "236" = c(0, 15000)
)

# Convert mm to microns for consistent spatial scaling
thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000

# Apply unique offsets per slide/sample
thisObj.T@meta.data$x_slide_global <- thisObj.T@meta.data$x_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][1])

thisObj.T@meta.data$y_slide_global <- thisObj.T@meta.data$y_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][2])

# Verify transformed coordinates
head(thisObj.T@meta.data[, c("mus", "x_slide_global", "y_slide_global")])

# Extract corrected spatial coordinates
locs <- as.matrix(thisObj.T@meta.data[c("x_slide_global", "y_slide_global")])

# Extract gene expression count matrix
gcm <- thisObj.T@assays$RNA@counts

# Create SpatialExperiment object with spatial data
se <- SpatialExperiment(
  assay = list(counts = gcm), 
  spatialCoords = locs,
  colData = DataFrame(thisObj.T@meta.data)  
)

# Store slide/sample identity explicitly
colData(se)$slide_id <- thisObj.T@meta.data$mus

# Define parameters
lambda <- 0
res = 0.25
k_geom <- c(15, 30)    
aname <- "normcounts"
# Normalize counts
se <- computeLibraryFactors(se)
assay(se, aname) <- normalizeCounts(se, log = TRUE)

# Compute Banksy features per slide
se_list <- lapply(unique(colData(se)$slide_id), function(slide) {
  se_subset <- se[, colData(se)$slide_id == slide]  
  se_subset <- Banksy::computeBanksy(se_subset, assay_name = "counts" , compute_agf = TRUE, k_geom = k_geom)
  return(se_subset)
})

# Merge processed slides
se <- do.call(cbind, se_list)

# Ensure reproducibility
set.seed(1000)

# Run Banksy PCA
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)

# Run Banksy UMAP
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

# Perform clustering
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = res)

# Extract spatial coordinates
coords <- spatialCoords(se)
colData(se)$sdimx <- coords[, 1]
colData(se)$sdimy <- coords[, 2]

# Get unique sample names
unique_samples <- unique(colData(se)$slide_id)

# View all dim reductions
names(reducedDims(se))

# Identify Banksy clustering column dynamically
clustering_columns <- grep("^clust", colnames(colData(se)), value = TRUE)
print(clustering_columns)  # Debug: View available cluster labels

if (length(clustering_columns) == 0) {
  stop("Error: No clustering columns found in colData(se). Ensure Banksy clustering completed successfully.")
}

# Select first clustering column (most relevant)
clust_column <- clustering_columns[2]  
print(paste("Using clustering column:", clust_column))  

# Identify Banksy's UMAP embedding dynamically
umap_columns <- grep("^UMAP", names(reducedDims(se)), value = TRUE)
umap_columns

if (length(umap_columns) == 0) {
  stop("Error: No UMAP embedding found in reducedDims(se). Ensure Banksy UMAP completed successfully.")
}

umap_dimred <- umap_columns[2]  # Select first available UMAP embedding
umap_dimred

# Banksy Full UMAP plot
plot_umap <- plotReducedDim(se, 
                            dimred = umap_dimred,
                            point_size = 0.1,
                            colour_by = clust_column) +
  ggtitle("Banksy UMAP") +
  theme_minimal()+
  theme(legend.text = element_text(size = 12))

print(plot_umap)

file_name <- paste0(WD, "banksy_se_object.rds")
saveRDS(se, file = file_name)
# BankObj <- readRDS(file_name)
# se <- BankObj
# Loop through each sample and generate spatial + UMAP plots
for (sample in unique_samples) {
  se_subset <- se[, colData(se)$slide_id == sample]
  
  # Spatial reconstruction plot
  plot_spatial <- plotColData(se_subset,
                              x = "sdimx", y = "sdimy",
                              point_size = 0.6, colour_by = clust_column) +
    ggtitle(paste("Spatial Reconstruction -", sample)) +
    theme_minimal()
  
  # Banksy UMAP plot
  plot_umap <- plotReducedDim(se_subset,
                              dimred = umap_dimred,
                              colour_by = clust_column) +
    ggtitle(paste("Banksy UMAP -", sample)) +
    theme_minimal()
  
  # Arrange plots side by side
  print(plot_grid(plot_spatial, plot_umap, ncol = 2))
}


# se <- readRDS(paste0(WD, "banksy_se_object.rds"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Harmony*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Required Libraries
# Load Required Libraries
library(harmony)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(Banksy)

# Step 1: Transfer Metadata from Seurat to SpatialExperiment (`se`)
# Extract Seurat metadata (Ensure row names match column names in `se`)
seurat_metadata <- thisObj.T@meta.data[, c("age", "sex", "slide_ID_numeric")]

# Ensure row names match `se` sample names
rownames(seurat_metadata) <- colnames(thisObj.T)

# Check if all sample names match
if (!all(rownames(seurat_metadata) %in% colnames(se))) {
  stop("Error: Some sample names in Seurat metadata do not match the `se` object!")
}

# Add metadata to `se`
colData(se) <- cbind(colData(se), seurat_metadata[colnames(se), ])  # Ensure correct ordering

# Step 2: Extract Banksy PCA embeddings
pca_columns <- grep("^PCA", names(reducedDims(se)), value = TRUE)
if (length(pca_columns) == 0) stop("Error: No PCA embeddings found in reducedDims(se). Ensure Banksy PCA has been run.")

pca_dimred <- reducedDims(se)[[pca_columns[1]]]  # Select first available PCA embedding

# Step 3: Define multiple covariates for batch correction
batch_variables <- c("slide_ID_numeric", "sex", "age")  # Match exact names from metadata

# Step 4: Run Harmony to correct for multiple covariates
harmony_embeddings <- HarmonyMatrix(
  data_mat = pca_dimred,  # Banksy PCA embeddings
  meta_data = colData(se),  # Now contains `age`, `sex`, `slide_ID_numeric`
  vars_use = batch_variables,  # Specify multiple covariates
  do_pca = FALSE  # Already PCA-transformed
)

# Step 5: Store Harmony embeddings back into `se`
reducedDims(se)$harmony <- harmony_embeddings

# Step 6: Re-run UMAP using Harmony-corrected embeddings
se <- Banksy::runBanksyUMAP(se, dimred = "harmony", use_agf = TRUE, lambda = lambda)

# Step 7: Re-run Clustering using Harmony embeddings
se <- Banksy::clusterBanksy(se, dimred = "harmony", use_agf = TRUE, lambda = lambda, resolution = res)

# Step 8: Visualize the updated UMAP
umap_columns <- grep("^UMAP", names(reducedDims(se)), value = TRUE)
if (length(umap_columns) == 0) stop("Error: No UMAP embedding found after Harmony adjustment.")
print(umap_columns)

umap_dimred <- umap_columns[2]  # Select first available UMAP embedding
clustering_columns <- grep("^clust", colnames(colData(se)), value = TRUE)
clustering_columns

clust_column <- clustering_columns[4] 
clust_column
# Generate updated UMAP plot
plot_umap_harmony <- plotReducedDim(se, 
                                    dimred = umap_dimred,
                                    point_size = 0.1,
                                    colour_by = clust_column) +
  ggtitle("Banksy + Harmony (Multiple Covariates) UMAP") +
  theme_minimal()

print(plot_umap_harmony)


for (sample in unique_samples) {
  se_subset <- se[, colData(se)$slide_id == sample]
  
  # Spatial reconstruction plot
  plot_spatial <- plotColData(se_subset,
                              x = "sdimx", y = "sdimy",
                              point_size = 0.6, colour_by = clust_column) +
    ggtitle(paste("Spatial Reconstruction -", sample)) +
    theme_minimal()
  
  # Banksy UMAP plot
  plot_umap <- plotReducedDim(se_subset,
                              dimred = umap_dimred,
                              colour_by = clust_column) +
    ggtitle(paste("Banksy UMAP -", sample)) +
    theme_minimal()
  
  # Arrange plots side by side
  print(plot_grid(plot_spatial, plot_umap, ncol = 2))
}


file_name <- paste0(project_WD, "banksyHarmonized_se_object.rds")
saveRDS(se, file = file_name)
se <- readRDS(file_name)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Asessing harmonization improvements
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(cluster)

names(reducedDims(se))
# Silhoutte Scores between pre- and post- harmony
# Extract UMAP and cluster labels
umap_pre <- reducedDims(se)$UMAP_M1_lam0  # Pre-Harmony UMAP
umap_post <- reducedDims(se)$UMAP_harmony  # Post-Harmony UMAP

clusters_pre <- colData(se)$clust_M1_lam0_k50_res0.25
clusters_post <- colData(se)$clust_harmony_k50_res0.25

# Set a max sample size (e.g., 10,000 points) to avoid memory issues
sample_size <- min(50000, length(clusters_pre))  # Limit to 10K cells if dataset is larger

# Randomly select subset of cells
set.seed(42)
sample_indices <- sample(seq_along(clusters_pre), sample_size)

# Compute silhouette score on subset
sil_pre <- silhouette(as.numeric(as.factor(clusters_pre[sample_indices])), 
                      dist(umap_pre[sample_indices, ]))

sil_post <- silhouette(as.numeric(as.factor(clusters_post[sample_indices])), 
                       dist(umap_post[sample_indices, ]))

# Compare silhouette scores
boxplot(sil_pre[, 3], sil_post[, 3], names = c("Pre-Harmony", "Post-Harmony"),
        main = "Silhouette Scores (Subset)", ylab = "Silhouette Score")

cat("Mean Silhouette Score (Subset) - Pre-Harmony:", mean(sil_pre[, 3]), "\n")
cat("Mean Silhouette Score (Subset) - Post-Harmony:", mean(sil_post[, 3]), "\n")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Gene marker distribution check
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Convert Banksy objects to Seurat for marker gene analysis
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)


counts_mat <- assay(se, "counts")  # Modify to "normcounts" if needed
metadata <- as.data.frame(colData(se))
head(metadata)

clustering_columns <- grep("^clust", colnames(metadata), value = TRUE)
clustering_columns

# Step 4: Create Seurat object with extracted data
seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = metadata)

# Step 5: Assign first available clustering column to `Idents()`
if (length(clustering_columns) > 0) {
  cluster_id <- clustering_columns[1]  # Use the first clustering column
  Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_id]]
}
head(Idents(seurat_obj))

# Step 6: Extract and Add Banksy PCA Embeddings (Pre and Post-Harmony)
# Extract all four reduced dimensions from `se`
pca_columns <- reducedDims(se)

# Iterate through all four and store them in the Seurat object
for (pca_name in names(pca_columns)) {
  # Clean the key to ensure it meets Seurat's naming conventions
  clean_key <- gsub("[^a-zA-Z0-9]", "", pca_name)  # Remove non-alphanumeric characters
  clean_key <- paste0(clean_key, "_")  # Append underscore to match Seurat key format
  
  # Store each PCA/UMAP reduction in the Seurat object
  seurat_obj[[pca_name]] <- CreateDimReducObject(
    embeddings = pca_columns[[pca_name]],  # Extract matrix correctly
    key = clean_key,
    assay = DefaultAssay(seurat_obj)
  )
}

# Final Step: Print Confirmation
cat("Successfully converted `se` back to Seurat object.\n")
print(seurat_obj)

# Identify variable features using Seurat's method
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000, verbose = TRUE)
cat("Number of variable features in original seurat obj:", length(VariableFeatures(thisObj.T)), "\n")
cat("Number of variable features after correction:", length(VariableFeatures(seurat_obj)), "\n")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Plotting different UMAPs
# View reductions available in the seurat obj

names(seurat_obj@reductions)
names(seurat_obj@meta.data)

umap_pre_df <- data.frame(
  Embeddings(seurat_obj, reduction = "UMAP_M1_lam0"), 
  Cluster = as.factor(seurat_obj$clust_M1_lam0_k50_res0.25)
)

umap_post_df <- data.frame(
  Embeddings(seurat_obj, reduction = "UMAP_harmony"), 
  Cluster = as.factor(seurat_obj$clust_harmony_k50_res0.25)
)

# Rename columns to match expected UMAP coordinate names
colnames(umap_pre_df)[1:2] <- c("UMAP1", "UMAP2")
colnames(umap_post_df)[1:2] <- c("UMAP1", "UMAP2")

color_palette_1 <- brewer.pal(min(length(unique(umap_pre_df$Cluster)), 12), "Set3")  
color_palette_2 <- brewer.pal(min(length(unique(umap_post_df$Cluster)), 12), "Set3")  
# Plot Pre-Harmony UMAP
p1 <- ggplot(umap_pre_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.6) +  # Slightly larger points
  scale_color_manual(values = color_palette_1) +  # Apply Set1 colors
  ggtitle("Pre-Harmony UMAP") +
  theme_minimal() +
  theme(legend.text = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 10))) 

# Plot Post-Harmony UMAP
p2 <- ggplot(umap_post_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = color_palette_2) +
  ggtitle("Post-Harmony UMAP") +
  theme_minimal() +
  theme(legend.text = element_text(size = 12))+
  guides(color = guide_legend(override.aes = list(size = 10))) 

# Display side-by-side
p1 + p2
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Find markers for pre- and post- harmony Banksy embeddings

library(presto)
# devtools::install_github('immunogenomics/presto')

se_counts <- assay(se, "counts")
seurat_counts <- GetAssayData(seurat_obj, layer = "counts")
identical_counts <- identical(se_counts, seurat_counts)
cat("Count matrices identical:", identical_counts, "\n")

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- JoinLayers(seurat_obj)

# Step 1: Ensure the cluster identities are correctly set
Idents(seurat_obj) <- seurat_obj$clust_M1_lam0_k50_res0.25  # Pre-Harmony clusters


# Step 2: Run FindAllMarkers() for Pre-Harmony Clusters
cat("Running FindAllMarkers() for Pre-Harmony Clusters...\n")
pre_harmony_markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = FALSE
)

# Step 3: Run FindAllMarkers() for Post-Harmony Clusters
cat("Running FindAllMarkers() for Post-Harmony Clusters...\n")
Idents(seurat_obj) <- seurat_obj$clust_harmony_k50_res0.25  # Post-Harmony clusters

post_harmony_markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = FALSE
)

# Step 5: Print top markers for comparison
cat("\nTop 10 Pre-Harmony Markers:\n")
print(head(pre_harmony_markers, 10))

cat("\nTop 10 Post-Harmony Markers:\n")
print(head(post_harmony_markers, 10))


top_genes <- pre_harmony_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
head(top_genes)

# Ensure that 'top_genes' is sorted by 'cluster' before creating the factor
top_genes <- top_genes[order(top_genes$cluster), ]

# Create a factor for 'gene' based on its occurrence in 'top_genes' ordered by 'cluster'
top_genes$gene <- factor(top_genes$gene, levels = unique(top_genes$gene))

# Create the dot plot of the top 5 genes for each cluster
gene_plot <- ggplot(top_genes, aes(x = factor(cluster), y = gene, color = avg_log2FC, size = pct.1)) +
  geom_point() +
  scale_color_viridis_c() +  # Using viridis color scale for better visibility
  labs(title = "Top 5 Genes per Cluster",
       x = "Cluster",
       y = "Gene",
       color = "Log2 Fold Change",
       size = "Percentage Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 7))

print(gene_plot)

gene_lists <- list(
  mgList = c('C1qa','C1qb','C1qc','Ctss','Csf1r','Cst3','Cx3cr1','Hexb','Itgam','P2ry12','Selplg','Tmem119'),
  astroList = c('Gja1','Gpr37l1','Plpp3','Mfge8','Gfap','Slc1a3','Bcan'),
  gabaergicList = c('Gad1','Gad2','Npy','Slc6a1','Slc32a1','Sst'),
  neuronList = c('Adcy1','Atp1b1','Camk2a','Dnm1','Epha6','Gria1', 'Grin2a','Grin2b','Nell2','Nrgn','Snap25','Syp','Syt1','Ywhag'),
  oligoList = c('Plp1','Mog','Mag','Mobp','Mbp','Cryab','Fa2h'),
  opcList = c('Cntn1','Cspg5','Gpr17','Pdgfra', 'Ptprz1','Serpine2'),
  vascList = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg'),
  vlmcList = c('Apod','Col1a1','Dcn','Il33','Ptgds','Ptn'),
  epenList = c('Dbi','Foxj1','Nnat','Tuba1a','Pltp','Vim'),
  periList = c('Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn'),
  nstemnblastList = c('Ctnna2','H3f3b','Nnat','Ptprs','Tcf4','Tuba1a','Tubb5'),
  chorplexList = c('Ttr','Htr2c','Lgals3bp','Cab39l','Igfbp7','Ppp1r1b'),
  vascList2 = c('Cldn5','Vtn','Rgs5','Flt1','Itm2a','Vim','Pecam1','Bsg','Apod','Col1a1','Dcn','Il33','Ptgds','Ptn',
                'Acta2','Crip1','Higd1b','Flt1','Myl9','Myl6','Mfge8','Rgs5','Vtn'),
  GAL3enrichedModule = c('Apoe','Cd9','Clec7a','Ctsb','Ctsd','Fth1','Gpnmb','Itgax',
                         'Lyz1','Lilrb4a','Mmp12','Spp1','Tmsb4x','Tyrobp')
)

gene_lists <- list(
  # 1. Neurons (Excitatory: Granule + UBCs)
  Pan_Neurons = c("Rbfox3", "Snap25", "Syp", "Stmn2"),  # Pan-Neuronal
  UBC = c("Calb2", "Grm1", "Plcb4"),  # UBC-specific
  Granule = c("Slc17a7", "Grm4", "Nhlh1", "Ebf2", "Cbln3"),  # Granule Cell Alternative Markers
  
  # 2. Interneurons (Inhibitory: Golgi, Stellate, Basket, Lugaro, DCN GABAergic)
  Interneurons = c("Gad1", "Gad2", "Slc6a1", "Htr3a", "Nrgn",  # Generic GABAergic
  "Pvalb", "Crh", "Grm2", "Npy", "Htr2a"),  # Specific for Golgi, Stellate, Basket Interneurons
  
  # 3. Purkinje Cells
  Purkinje = c("Calb1", "Pcp2", "Grid2", "Itpr1", "Rora"),
  
  # 4. Bergmann Glia
  BergmannGlia = c("Gfap", "S100b", "Aqp4", "Slc1a3", "Sox9"),
  
  # 5. Fibroblasts
  Fibroblasts = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum"),
  
  # 6. Microglia
  MG = c("P2ry12", "Tmem119", "Cx3cr1", "C1qa", "Trem2"),
  
  # 7. Oligodendrocytes
  Olig = c("Mbp", "Plp1", "Mog", "Mag"),
  
  # 8. Oligodendrocyte Precursor Cells (OPCs)
  OPC = c("Pdgfra", "Olig2", "Gpr17", "Tcf7l2"),
  
  # 9. Endothelial Cells
  EC = c("Cldn5", "Pecam1", "Flt1", "Esam"),
  
  # 10. Astrocytes
  Astro  = c("Gfap", "Aldh1l1", "Slc1a3", "Aqp4", "S100b"),
  
  # 11. Pericytes
  Peri = c("Pdgfrb", "Rgs5", "Des", "Acta2"),
  
  # 12. Vascular Smooth Muscle Cells (VSMCs)
  VSCMCs = c("Acta2", "Tagln", "Notch3"),
  
  # 13. Choroid Plexus
  ChoroidPlexus = c("Ttr", "Cldn5", "Otx2", "Kcne2", "Aqp1")
)

# Loop through each gene list and create a dot plot
for (list_name in names(gene_lists)) {
  genes <- gene_lists[[list_name]]
  
  # Generate and display the dot plot
  dot_plot <- DotPlot(seurat_obj, features = genes, cols = c('grey','black')) +
    ggtitle(paste("Dot Plot for", list_name)) +
    theme_minimal()
  
  print(dot_plot)  # Display each plot
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Spacia Data Prep & Run from R wrapper ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Move obj to general obj, writing to eventually wrap into a FXN
obj <- thisObj.T
assay <- "RNA"
cell_labels <- "subtype" #or Idents(obj)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj@meta.data$x_slide_micron <- obj@meta.data$x_slide_mm * 1000
obj@meta.data$y_slide_micron <- obj@meta.data$y_slide_mm * 1000
# 1. Extract normalized expression data (from the "data" slot)
#    Note: Seurat stores data as a genes x cells matrix, so we transpose it.
expr_mat <- GetAssayData(
  object = obj,
  assay = assay,
  layer = "data"
)
expr_mat <- as.data.frame(t(expr_mat))
expr_mat <- cbind(cell = rownames(expr_mat), expr_mat)

# Write expression matrix to file (cells as rows, genes as columns)
write.table(expr_mat, file = "counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2. Prepare metadata file
#    Adjust the column names below if your spatial coordinates are stored differently.
meta <- obj@meta.data
meta_spacia <- meta[, c("x_slide_micron", "y_slide_micron", cell_labels)]

# Add cell names as a column
meta_spacia <- cbind(cell = rownames(meta), meta_spacia)

# Rename columns to match Spacia's expectations: cell, X, Y, cell_type
colnames(meta_spacia) <- c("cell", "X", "Y", "cell_type")

# Write metadata to file
write.table(meta_spacia, file = paste0(local_WD,"spacia_metadata.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
# Run Spacia from R wrapper
#$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
# Required R packages for spacia
install.packages(c("coda", "ggmcmc", "Rcpp", "RcppArmadillo", "rjson", "diptest", "filelock", "data.table", "jsonlite", "optparse", "patchwork", "ggplot2"))





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Banksy Testing ####

# Running Banksy tutorial vignette
library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)

install.packages("scuttle")

thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000
# Sample names
# Select one to generate spatial plot
sample_name <- "305"
sample_name <- "YM04"
sample_name <- "304"
sample_name <- "YM03"
sample_name <- "YM02"
sample_name <- "302"
sample_name <- "YM01"
sample_name <- "301"
sample_name <- "YF04"
sample_name <- "243"
sample_name <- "YF03"
sample_name <- "239"
sample_name <- "YF02"
sample_name <- "237"
sample_name <- "YF01"
sample_name <- "236"
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Subset the Seurat object by slide and sample
subset_data <- subset(thisObj.T, subset = mus == sample_name)

head(subset_data@meta.data[c("x_slide_micron", "y_slide_micron")])
# Extract spatial coordinates from your Seurat object's metadata
locs <- as.matrix(subset_data@meta.data[c("x_slide_micron", "y_slide_micron")])

# Extract the count matrix from the RNA assay (adjust if your counts are stored elsewhere)
gcm <- subset_data@assays$RNA@counts

# Load the required package for creating SpatialExperiment objects
library(SpatialExperiment)

# Create a SpatialExperiment object using the count matrix and spatial coordinates
se <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)

# Perform quality control based on total counts
qcstats <- perCellQCMetrics(se)
thres <- quantile(qcstats$total, c(0.05, 0.98))
keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
se <- se[, keep]

# Normalize to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)


#%%%%%%%%%%%%%%
# Set parameters
lambda <- c(0, 0.2)
k_geom <- c(15, 30)

# Run BANKSY to compute the neighborhood-augmented matrix.
# 'compute_agf=TRUE' computes both the weighted neighborhood mean (ℳ) and the azimuthal Gabor filter (𝒢)
se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
# Run PCA, UMAP and clustering on the BANKSY matrix
set.seed(1000)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)

# Optionally, relabel clusters to minimize differences across clustering runs
se <- Banksy::connectClusters(se)
# e.g. clust_M1_lam0.2_k50_res1.2 --> clust_M1_lam0_k50_res1.2

# Ensure the spatial coordinates are available as "sdimx" and "sdimy" in colData.
# (If your SpatialExperiment object was built using 'locs', these may not be named yet.)
coords <- spatialCoords(se)
colData(se)$sdimx <- coords[, 1]
colData(se)$sdimy <- coords[, 2]

# Extract the clustering labels from colData (labels starting with "clust")
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
# Bind the spatial coordinates (if not already present)
colData(se) <- cbind(colData(se), coords)

# Visualize the non-spatial clustering (lambda = 0) versus BANKSY clustering (lambda = 0.2)
plot_nsp <- plotColData(se,
                        x = "sdimx", y = "sdimy",
                        point_size = 0.6, colour_by = cnames[1]
)
plot_bank <- plotColData(se,
                         x = "sdimx", y = "sdimy",
                         point_size = 0.6, colour_by = cnames[2]
)

# Arrange the two plots side by side
library(cowplot)
plot_grid(plot_nsp + coord_equal(), plot_bank + coord_equal(), ncol = 2)

plot_grid(
  plot_nsp + facet_wrap(~colour_by),
  plot_bank + facet_wrap(~colour_by),
  ncol = 2
)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Banksy Testing w/ sample x,y grid buffering
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Install missing packages if necessary
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("Raphael-LS/Banksy")
}
if (!requireNamespace("scuttle", quietly = TRUE)) install.packages("scuttle")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

# Load all necessary libraries
library(Banksy)          
library(SummarizedExperiment)  
library(SpatialExperiment)  
library(scuttle)        
library(scater)        
library(cowplot)        
library(ggplot2)        

# Define unique slide offsets to prevent overlapping coordinates
slide_offsets <- list(
  "305" = c(0, 0),
  "YM04" = c(5000, 0),
  "304" = c(10000, 0),
  "YM03" = c(15000, 0),
  "YM02" = c(20000, 0),
  "302" = c(0, 5000),
  "YM01" = c(5000, 5000),
  "301" = c(10000, 5000),
  "YF04" = c(15000, 5000),
  "243" = c(20000, 5000),
  "YF03" = c(0, 10000),
  "239" = c(5000, 10000),
  "YF02" = c(10000, 10000),
  "237" = c(15000, 10000),
  "YF01" = c(20000, 10000),
  "236" = c(0, 15000)
)

# Convert mm to microns for consistent spatial scaling
thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000

# Apply unique offsets per slide/sample
thisObj.T@meta.data$x_slide_global <- thisObj.T@meta.data$x_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][1])

thisObj.T@meta.data$y_slide_global <- thisObj.T@meta.data$y_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][2])

# Verify transformed coordinates
head(thisObj.T@meta.data[, c("mus", "x_slide_global", "y_slide_global")])

# Extract corrected spatial coordinates
locs <- as.matrix(thisObj.T@meta.data[c("x_slide_global", "y_slide_global")])

# Extract gene expression count matrix
gcm <- thisObj.T@assays$RNA@counts

# Create SpatialExperiment object with spatial data
se <- SpatialExperiment(
  assay = list(counts = gcm), 
  spatialCoords = locs,
  colData = DataFrame(thisObj.T@meta.data)  
)

# Store slide/sample identity explicitly
colData(se)$slide_id <- thisObj.T@meta.data$mus

# Define parameters
lambda <- 0
res = 0.25
k_geom <- c(15, 30)    

# Normalize counts
se <- computeLibraryFactors(se)
assay(se, aname) <- normalizeCounts(se, log = FALSE)

# Compute Banksy features per slide
se_list <- lapply(unique(colData(se)$slide_id), function(slide) {
  se_subset <- se[, colData(se)$slide_id == slide]  
  se_subset <- Banksy::computeBanksy(se_subset, assay_name = "counts" , compute_agf = TRUE, k_geom = k_geom)
  return(se_subset)
})

# Merge processed slides
se <- do.call(cbind, se_list)

# Ensure reproducibility
set.seed(1000)

# Run Banksy PCA
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)

# Run Banksy UMAP
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)

# Perform clustering
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = res)

# Extract spatial coordinates
coords <- spatialCoords(se)
colData(se)$sdimx <- coords[, 1]
colData(se)$sdimy <- coords[, 2]

# Get unique sample names
unique_samples <- unique(colData(se)$slide_id)

# Identify Banksy clustering column dynamically
clustering_columns <- grep("^clust", colnames(colData(se)), value = TRUE)
print(clustering_columns)  # Debug: View available cluster labels

if (length(clustering_columns) == 0) {
  stop("Error: No clustering columns found in colData(se). Ensure Banksy clustering completed successfully.")
}

# Select first clustering column (most relevant)
clust_column <- clustering_columns[3]  
print(paste("Using clustering column:", clust_column))  

# Identify Banksy's UMAP embedding dynamically
umap_columns <- grep("^UMAP", names(reducedDims(se)), value = TRUE)


if (length(umap_columns) == 0) {
  stop("Error: No UMAP embedding found in reducedDims(se). Ensure Banksy UMAP completed successfully.")
}

umap_dimred <- umap_columns[1]  # Select first available UMAP embedding


# Banksy Full UMAP plot
plot_umap <- plotReducedDim(se, 
                            dimred = umap_dimred,
                            point_size = 0.1,
                            colour_by = clust_column) +
  ggtitle("Banksy UMAP") +
  theme_minimal()

print(plot_umap)



# Loop through each sample and generate spatial + UMAP plots
for (sample in unique_samples) {
  se_subset <- se[, colData(se)$slide_id == sample]

  # Spatial reconstruction plot
  plot_spatial <- plotColData(se_subset,
                              x = "sdimx", y = "sdimy",
                              point_size = 0.6, colour_by = clust_column) +
    ggtitle(paste("Spatial Reconstruction -", sample)) +
    theme_minimal()

  # Banksy UMAP plot
  plot_umap <- plotReducedDim(se_subset,
                              dimred = umap_dimred,
                              colour_by = clust_column) +
    ggtitle(paste("Banksy UMAP -", sample)) +
    theme_minimal()

  # Arrange plots side by side
  print(plot_grid(plot_spatial, plot_umap, ncol = 2))
}
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Iterative testing of banksy parameters
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("Raphael-LS/Banksy")
}
if (!requireNamespace("scuttle", quietly = TRUE)) install.packages("scuttle")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

# Load all necessary libraries
library(Banksy)          
library(SummarizedExperiment)  
library(SpatialExperiment)  
library(scuttle)        
library(scater)        
library(cowplot)        
library(ggplot2)        

# Define unique slide offsets to prevent overlapping coordinates
slide_offsets <- list(
  "305" = c(0, 0),
  "YM04" = c(5000, 0),
  "304" = c(10000, 0),
  "YM03" = c(15000, 0),
  "YM02" = c(20000, 0),
  "302" = c(0, 5000),
  "YM01" = c(5000, 5000),
  "301" = c(10000, 5000),
  "YF04" = c(15000, 5000),
  "243" = c(20000, 5000),
  "YF03" = c(0, 10000),
  "239" = c(5000, 10000),
  "YF02" = c(10000, 10000),
  "237" = c(15000, 10000),
  "YF01" = c(20000, 10000),
  "236" = c(0, 15000)
)

# Convert mm to microns for consistent spatial scaling
thisObj.T@meta.data$x_slide_micron <- thisObj.T@meta.data$x_slide_mm * 1000
thisObj.T@meta.data$y_slide_micron <- thisObj.T@meta.data$y_slide_mm * 1000

# Apply unique offsets per slide/sample
thisObj.T@meta.data$x_slide_global <- thisObj.T@meta.data$x_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][1])

thisObj.T@meta.data$y_slide_global <- thisObj.T@meta.data$y_slide_micron + 
  sapply(thisObj.T@meta.data$mus, function(m) slide_offsets[[m]][2])

# Verify transformed coordinates
head(thisObj.T@meta.data[, c("mus", "x_slide_global", "y_slide_global")])

# Extract corrected spatial coordinates
locs <- as.matrix(thisObj.T@meta.data[c("x_slide_global", "y_slide_global")])

# Extract gene expression count matrix
gcm <- thisObj.T@assays$RNA@counts

# Create SpatialExperiment object with spatial data
# se <- SpatialExperiment(
#   assay = list(counts = gcm), 
#   spatialCoords = locs,
#   colData = DataFrame(thisObj.T@meta.data)  
# )

#%%%%%%Subsetting Data for iterative testing%%%%%%%%%%%%%%%%%%%%%%%%
# Define the proportion of the dataset to sample
subset_fraction <- 0.5  # Adjust as needed (e.g., 10% of total cells)

# Ensure we don't sample more than available cells
num_total_cells <- ncol(gcm)
num_subset_cells <- round(num_total_cells * subset_fraction)

# Randomly select a subset of column names (cells)
set.seed(12)  # Ensure reproducibility
selected_cells <- sample(colnames(gcm), num_subset_cells)

# Subset `gcm` (gene expression matrix)
gcm_subset <- gcm[, selected_cells]

# Subset `colData` (metadata) to match `gcm`
colData_subset <- thisObj.T@meta.data[selected_cells, ]

# Subset `locs` (spatial coordinates) to match the selected cells
locs_subset <- as.matrix(colData_subset[c("x_slide_global", "y_slide_global")])

# Create the SpatialExperiment object for the subset
se_subset.T <- SpatialExperiment(
  assay = list(counts = gcm_subset), 
  spatialCoords = locs_subset,
  colData = DataFrame(colData_subset)  
)


# Store slide/sample identity explicitly
colData(se_subset.T)$slide_id <- colData_subset$mus

# Define parameters
aname <- "counts"  
lambda <- 0.2 # Recommended by Banksy package for CosMX
res <- 0.4 #seq(0.1, 0.8, by = 0.1)

# Define possible k_geom values
k_geom_single <- c(5,10) #seq(5, 40, by = 5)  # Range of k_geom values
# Generate all possible pairs of k_geom values (ensuring k1 < k2)
k_geom_pairs <- combn(k_geom_single, 2, simplify = FALSE)

# Print generated pairs
print(k_geom_pairs)

# Store results
banksy_results <- list()

for (k_geom in k_geom_pairs) {
  k_geom <- unlist(k_geom)  # Ensure numeric vector
  
  for (r in res) {  # Loop over different resolution values
    cat("Testing k_geom =", k_geom, "Resolution =", r, "\n")
    
    # Normalize counts
    se_norm <- computeLibraryFactors(se_subset.T)
    assay(se_norm, aname) <- normalizeCounts(se_norm, log = FALSE)
    
    # Compute Banksy features
    se_list <- lapply(unique(colData(se_norm)$slide_id), function(slide) {
      se_subset <- se_norm[, colData(se_norm)$slide_id == slide]  
      se_subset <- Banksy::computeBanksy(se_subset, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
      return(se_subset)
    })
    
    # Merge processed slides
    se_processed <- do.call(cbind, se_list)
    
    # Run PCA, UMAP, and Clustering
    se_processed <- Banksy::runBanksyPCA(se_processed, use_agf = TRUE, lambda = lambda)
    se_processed <- Banksy::runBanksyUMAP(se_processed, use_agf = TRUE, lambda = lambda)
    se_processed <- Banksy::clusterBanksy(se_processed, use_agf = TRUE, lambda = lambda, resolution = r)
    
    # Store results
    banksy_results[[paste0("k", k_geom[1], "_", k_geom[2], "_res", r)]] <- se_processed
    
    # Ensure spatial coordinates exist
    if (!is.null(spatialCoords(se_processed))) {
      coords <- spatialCoords(se_processed)
      colData(se_processed)$sdimx <- coords[, 1]
      colData(se_processed)$sdimy <- coords[, 2]
    } else {
      warning("Spatial coordinates missing in se_processed!")
    }
    
    # Detect clustering columns
    clustering_columns <- grep("^clust", colnames(colData(se_processed)), value = TRUE)
    if (length(clustering_columns) == 0) next
    
    # Detect UMAP embedding
    umap_columns <- grep("^UMAP", names(reducedDims(se_processed)), value = TRUE)
    if (length(umap_columns) == 0) next
    umap_dimred <- umap_columns[1]
    
    # Generate UMAP plots
    for (clust_column in clustering_columns) {
      plot_umap <- plotReducedDim(se_processed, dimred = umap_dimred, point_size = 0.1, colour_by = clust_column) +
        ggtitle(paste0("Banksy UMAP: ", clust_column)) + theme_minimal()
      print(plot_umap)
      file_name <- paste0("C:/Users/m322332/Documents/UMAP_Banksy_test_",clustering_columns, ".png")
      ggsave(filename = file_name, plot = plot_umap, dpi = 600)
    }
  }
}
# Res = 0.4 seems best for Banksy
# 

for (result_name in names(banksy_results)) {
  cat("Generating UMAP plot for:", result_name, "\n")
  
  # Extract the Banksy SpatialExperiment object
  se_obj <- banksy_results[[result_name]]
  
  # Identify Banksy's UMAP embedding dynamically
  umap_columns <- grep("^UMAP", names(reducedDims(se_obj)), value = TRUE)
  
  # Ensure at least one UMAP embedding exists
  if (length(umap_columns) == 0) {
    warning(paste("No UMAP embedding found for:", result_name))
    next
  }
  
  # Select the first available UMAP embedding
  umap_dimred <- umap_columns[1]
  cat("Using UMAP embedding:", umap_dimred, "\n")
  
  # Identify the correct clustering column dynamically
  clustering_columns <- grep("^clust", colnames(colData(se_obj)), value = TRUE)
  
  # Ensure at least one clustering column exists
  if (length(clustering_columns) == 0) {
    warning(paste("No clustering column found for:", result_name))
    next
  }
  
  # Select the first available clustering column
  clust_column <- clustering_columns[1]
  cat("Using clustering column:", clust_column, "\n")
  
  # Generate UMAP plot for the current clustering column
  plot_umap <- plotReducedDim(se_obj, 
                              dimred = umap_dimred,
                              point_size = 0.1,
                              colour_by = clust_column) +
    ggtitle(paste0("Banksy UMAP: ", result_name)) +
    theme_minimal()
  
  # Print the UMAP plot
  print(plot_umap)

}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Spacer ####

# Cluster Expression Plots ####
QueryList = c("Lgals3", "Glb1")
DotPlot(thisObj.T, features = QueryList) + 
  RotatedAxis()

DoHeatmap(subset(thisObj.T, downsample = 100),
          group.by= "age",
          features = senBrainList,
          size = 3)+
  ggtitle("senBrainList")

RidgePlot(thisObj.T, features = QueryList, ncol = 2)
VlnPlot(thisObj.T, features = QueryList, group.by = "subtype", split.by = 'age',
        split.plot = TRUE, alpha = 0.1)


# Define outlier removal function
remove_outliers <- function(expression_values, threshold = 2) {
  # Identify outliers using Z-scores
  z_scores <- scale(expression_values)
  non_outliers <- which(abs(z_scores) <= threshold)
  return(expression_values[non_outliers])
}

# Example gene lists
gene_lists <- list(
  senBrainList = c("Apoe", "App", "Bcl2", "Bcl2l1", "C4b", "Casp3", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                   "Ccne1", "CD26", "Cd274", "Cd36", "Cd74", "Cd9", "Cdc6", "Cdc27",
                   "Cdk4", "Cdk8", "Cdk9", "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                   "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19",
                   "Cenpb", "Col3a1", "Clu", "Csf1", "Csf1r", "Csf3r", "Cxcl14", "Cxcl2", "Elmo1", "Fxyd5", "Gadd45a",
                   "Gdf15", "Glb1", "Gpnmb", "H2-Aa", "H2-Ab1", "Hmgb1", "Ifitm3", "Ifng",
                   "Il10", "Il1a", "Il1b", "Il23a", "Il23r", "Il33", "Il6", "Irf5",
                   "Irf7", "Klre1", "Lamb1", "Lamp2", "Lgals3", "Lgals3bp",
                   "Lmna", "Lpl", "Malat1", "Mapk1", "Mapk8", "Mapk10", "Mapk14", "Meox1", "Mmp12", "Mmp14",
                   "Mmp2", "Mmp3", "Mmp9", "Nfkb1", "Nol3", "Notch3", "Pak1", "Pak2", "Pak3", "Pdcd1", "Plaur", "Pten",
                   "Rbl2", "Serpina3n", "Serpinb1a", "Serpine1", "Serpine2", "Sirt1", "Smad2", "Smad3", "Spp1",
                   "Tgfb1", "Tgfb2", "Tgfbr1", "Tmem173", "Tnf", "Tp53bp1", "Trp53", "Vim", "Vtn", "Xcl1"),
  cellCycleDomain = c("Cdk4", "Cdk8","Cdk9", "Cdc6", "Cdc27",
                      "Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
                      "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19", 
                      "Tp53bp1", "Trp53", "Jun","Rb1","Rbl2"),
  #2&3
  nuclearStressDomain = c("Cenpb","Gadd45a","Hmgb1","Irf7","Mapk14","Lmna","Nfkb1",
                          "Pak1","Pak2","Pak3","Phc2","Phc3","Smad2","Smad3","Tp53bp1","Trp53","Xcl1"),
  #4
  scapDomain = c("Akt1","Akt3","Bcl2", "Bcl2l1","Clu","Casp3","Pik3ca", "Pik3cb", "Pik3r1", "Serpina3n", "Serpinb1a","Serpine1","Serpine2"),
  #5
  saspDomain = c("Ccl2","Ccl3","Ccl4","Ccl5","Csf1","Cd74","Cd81","Cxcl2","Cxcl14",
                 "Il1a","Il1b","Il6","Il10","Il23a","Il33","Gdf15","Mmp12", "Mmp14", "Mmp2", "Mmp3", "Mmp9",
                 "Spp1","Tgfb1","Tgfb2","Tnf","Ifng","Irf5","Irf7"),
  #6&7 Lysomal stress, metabolic adaptations, mitochondrial ROS
  organelleStressDomain = c("Cyba","Fth1","Glb1","Hmgb1","Hspa1a","Mapk1", "Mapk8", "Mapk10", "Mapk14","Lamp2","Lgals3","Lgals3bp","Sirt1","Sod1"),
  #8 Cell surface markers
  cellSurfaceDomain = c("Cd9","CD26","Cd274","Csf1r","Gpnmb","H2-Aa","H2-Ab1","Ifitm3","Igfbp7","Il23r","Notch3","Plaur","Tgfbr1","Tmem173","Vim"),
  #0 Aged brain upregulated
  #ageSelectiveDomain = c("Apoe","App","Clu","Malat1","Lpl"),
  
  nanostringSenQuiet = c("Ago3","Akt1","Akt3","Anapc16","Atm","Atr","Btrc","Cacna1d",
                         "Calm1","Calm2","Calm3","Cdc27","Cdk6","Fbxw11","Foxo1","Foxo3",
                         "H3f3b","Hipk2","Igfbp7","Il6","Itpr1","Itpr2","Jun","Map2k1","Map2k4","Map2k5",
                         "Map2k6","Map3k5","Map4k4","Mapk1","Mapk10","Mapk14","Mapk8",
                         "Mcu","Nfkb1","Phc2","Phc3","Pik3ca","Pik3cb","Pik3r1","Ppp3ca","Ppp3cb","Ppp3cc",
                         "Pten","Raf1","Rb1","Rps6ka2","Rps6ka3","Slc25a4","Smad2","Smad3","Stat3",
                         "Terf2ip","Tfdp2","Tgfb1","Tgfb2","Tgfbr1","Tnik","Tnrc6a","Tnrc6b","Tnrc6c",
                         "Trp53","Trpm7","Ube2e1"),
  
  damList = sort(c('Apoe', 'B2m', 'C1qa','Ctsb', 'Ctsd', 'Ctss', 'Csf1','Csf1r', 'Cst3', 'Cst7', 'Cx3cr1', "Cxcl14",
                   'Lgals3', 'Ftl1', 'Fth1', 'Clec7a', 'Cd74', 'H2-Aa','H2-Ab1', 'Itgax', 'Lyz1/2','Tgfbr1','Tmem173','Trem2','Tyrobp','Spp1')),
  #mgList = c('C1qa','C1qb','C1qc','Ctss','Csf1r','Cst3','Cx3cr1','Hexb','Itgam','P2ry12','Selplg','Tmem119'),
  
  cdknList = c("Cdkn1a/p21", "Cdkn1b/p27kip1", "Cdkn1c/p57kip2",
               "Cdkn2a/p16", "Cdkn2b/p15ink4b", "Cdkn2c/p18ink4c","Cdkn2d/p19"),
  
  #Cd9 low, Ctsd low
  macList = c('C2','C3','C4b','Crip1','Cd74','Lyz1/2','H2-Ab1','H2-Aa','Hsp90ab1',
              'Ldha','Pltp','Bsg','Hmgb1','H3f3b',
              'Ifitm3','Fxyd5','Vim','Ctsd','Cd9','Cd63','Cd81'),
  
  #Cd81 low, Ptprc and Tgfb1 too prevalent to include
  immuneinfiltrateList = c('Ccl5', 'Cd3d', 'Cd3e', 'Cd3g', 'Cd44','Cd46',
                           'Cd74', 'Cd274', 'Crip1', 'Fau','H2-Aa', 'H2-Ab1',
                           'Ifitm3', 'Il6', 'Il23a', 'Il33', 'Klre1', 'Ltbp1',
                           'Ms4a4a', 'Nfkb1', 'Nkg7', 'Pdcd1', 'Tmsb10')
  
)
# DOT PLOT GENE LISTS AVERAGE OF EACH CLUSTER
# Initialize a data frame to store results
average_expression_df <- data.frame()

# Loop through each gene list
for (gene_list_name in names(gene_lists)) {
  genes <- gene_lists[[gene_list_name]]
  # Filter out genes not in the dataset
  genes <- intersect(genes, rownames(seurat_obj))
  
  # Extract expression matrix for the gene list
  expression_data <- FetchData(seurat_obj, vars = genes)
  
  # Calculate average expression for each cluster
  cluster_means <- tapply(rowMeans(expression_data, na.rm = TRUE), Idents(seurat_obj), mean, na.rm = TRUE)
  
  # Store results
  average_expression_df <- rbind(
    average_expression_df,
    data.frame(
      GeneList = gene_list_name,
      Cluster = names(cluster_means),
      AverageExpression = cluster_means
    )
  )
}

# Prepare data for DotPlot
average_expression_df$Cluster <- factor(average_expression_df$Cluster, levels = unique(Idents(seurat_obj)))
average_expression_df$GeneList <- factor(average_expression_df$GeneList)

# Create DotPlot
ggplot(average_expression_df, aes(x = GeneList, y = Cluster, size = AverageExpression, color = AverageExpression)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(
    x = "Gene List",
    y = "Cluster",
    size = "Average Expression",
    color = "Expression Level"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# DOT PLOT AVERAGE SPLIT BY AGE
# Initialize a data frame to store results
average_expression_df <- data.frame()

# Check if "Age" metadata exists
if (!"age" %in% colnames(seurat_obj@meta.data)) {
  stop("The Seurat object does not have an 'Age' metadata column.")
}

# Loop through each age group and each gene list
for (age_group in c("OLD", "YNG")) {
  # Subset the Seurat object by age group
  age_subset <- subset(seurat_obj, subset = age == age_group)
  
  for (gene_list_name in names(gene_lists)) {
    genes <- gene_lists[[gene_list_name]]
    # Filter out genes not in the dataset
    genes <- intersect(genes, rownames(age_subset))
    
    # Extract expression matrix for the gene list
    expression_data <- FetchData(age_subset, vars = genes)
    
    # Calculate average expression for each cluster
    cluster_means <- tapply(rowMeans(expression_data, na.rm = TRUE), Idents(age_subset), mean, na.rm = TRUE)
    
    # Store results
    average_expression_df <- rbind(
      average_expression_df,
      data.frame(
        GeneList = gene_list_name,
        Cluster = names(cluster_means),
        AverageExpression = cluster_means,
        Age = age_group
      )
    )
  }
}

# Prepare data for DotPlot
average_expression_df$Cluster <- factor(average_expression_df$Cluster, levels = unique(Idents(seurat_obj)))
average_expression_df$GeneList <- factor(average_expression_df$GeneList)
average_expression_df$Age <- factor(average_expression_df$Age, levels = c("YNG", "OLD"))

ggplot(average_expression_df, aes(x = GeneList, y = Cluster, size = AverageExpression, color = AverageExpression)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(
    x = "Gene List",
    y = "Cluster",
    size = "Average Expression",
    color = "Expression Level"
  ) +
  facet_wrap(~ Age) +  # Separate plots for OLD and YNG
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#DELTA AGE DOTPLOT
# Initialize a data frame to store results
average_expression_df <- data.frame()

# Loop through each age group and each gene list
for (age_group in c("OLD", "YNG")) {
  # Subset the Seurat object by age group
  age_subset <- subset(seurat_obj, subset = age == age_group)
  
  for (gene_list_name in names(gene_lists)) {
    genes <- gene_lists[[gene_list_name]]
    # Filter out genes not in the dataset
    genes <- intersect(genes, rownames(age_subset))
    
    # Extract expression matrix for the gene list
    expression_data <- FetchData(age_subset, vars = genes)
    
    # Calculate average expression for each cluster
    cluster_means <- tapply(rowMeans(expression_data, na.rm = TRUE), Idents(age_subset), mean, na.rm = TRUE)
    
    # Store results
    average_expression_df <- rbind(
      average_expression_df,
      data.frame(
        GeneList = gene_list_name,
        Cluster = names(cluster_means),
        AverageExpression = cluster_means,
        Age = age_group
      )
    )
  }
}

# Split into OLD and YNG subsets
old_df <- subset(average_expression_df, Age == "OLD")
yng_df <- subset(average_expression_df, Age == "YNG")

# Merge OLD and YNG data on GeneList and Cluster
merged_df <- merge(old_df, yng_df, by = c("GeneList", "Cluster"), suffixes = c("_OLD", "_YNG"))

# Calculate the difference
merged_df$Difference <- merged_df$AverageExpression_OLD - merged_df$AverageExpression_YNG

# Prepare data for DotPlot
merged_df$Cluster <- factor(merged_df$Cluster, levels = unique(Idents(seurat_obj)))
merged_df$GeneList <- factor(merged_df$GeneList)

ggplot(merged_df, aes(x = GeneList, y = Cluster, size = abs(Difference), color = Difference)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Gene List",
    y = "Cluster",
    size = "Magnitude of Difference",
    color = "OLD - YNG"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Dimplots ####
DimPlot(thisObj.T,
        group.by = "age",
        label = FALSE,
        reduction = "umap",
        raster = FALSE,
        pt.size = 0.1) +
  ggtitle("UMAP")

# Cell Count/Proportion Plots ####
seurat_obj <- thisObj.T
# Count the number of cells per cluster
cell_counts <- data.frame(
  Cluster = Idents(seurat_obj),
  Count = as.numeric(table(Idents(seurat_obj)))
)

ggplot(cell_counts, aes(x = Cluster, y = Count)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Number of Cells",
    title = "Number of Cells per Cluster"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




cell_counts <- seurat_obj@meta.data %>%
  group_by(Cluster = Idents(seurat_obj), age, sex) %>%
  summarise(Count = n(), .groups = "drop")

ggplot(cell_counts, aes(x = Cluster, y = Count, fill = interaction(age, sex))) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Number of Cells",
    fill = "Age and Sex",
    title = "Number of Cells per Cluster by Age and Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Group the data by Cluster, Age, and Sex, count the cells, and calculate proportions
cell_counts <- seurat_obj@meta.data %>%
  group_by(Cluster = Idents(seurat_obj), age, sex) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count))  # Normalize counts per cluster

# Create the normalized stacked column plot
ggplot(cell_counts, aes(x = Cluster, y = Proportion, fill = interaction(age, sex))) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Proportion of Cells",
    fill = "Age and Sex",
    title = "Proportion of Cells per Cluster by Age and Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Group by subtype too
# Group the data by Cluster, Age, and Sex, count the cells, and calculate proportions
cell_counts <- seurat_obj@meta.data %>%
  group_by(Cluster = seurat_obj$subtype, age, sex) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count))  # Normalize counts per cluster

# Create the normalized stacked column plot
ggplot(cell_counts, aes(x = Cluster, y = Proportion, fill = interaction(age, sex))) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Proportion of Cells",
    fill = "Age and Sex",
    title = "Proportion of Cells per Sub-Cluster by Age and Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# CellTalkDB####
thisObj.T <- readRDS("C:/Users/WA30480/Documents/seurat_thisObj.T_UMAP_V5_cellTyped_11.14.2024")

# download the source package of scsrctdb-1.0.tar.gz and install it
# ensure the right directory for scsrctdb-1.0.tar.gz
# ensure the dependency packages have been installed.('BiocManager','circlize','limma','igraph','gplots','grDevices','edgeR','SIMLR','data.table','pheatmap','stats','Rtsne','graphics','stringr','foreach','multtest','scran')
"C:/Users/WA30480/Downloads/scsrctdb-1.0.tar.gz"
install.packages(pkgs = "C:/Users/WA30480/Downloads/scsrctdb-1.0.tar.gz",repos = NULL, type = "source")

new_cluster_labels <- c(
  "Neurons", 
  "Interneurons",
  "Oligodendrocytes",
  "VascularandPericyte",
  "Purkinje",
  "Astrocytes",
  "MG",
  "OPCs",
  "Choroid Plexus"
)

# Extract normalized expression matrix
data <- as.matrix(thisObj.T@assays$RNA@data)  # Or the relevant assay
write.table(data, 
            file = "seurat_rna_matrix.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA)

file <- read.table(
  file = "seurat_rna_matrix.txt", 
  sep = "\t", 
  header = TRUE, 
  row.names = 1
)

# Extract cluster assignments
cluster <- as.vector(Idents(thisObj.T))
genes <- rownames(data)

# Ensure that Idents are set correctly as a character vector
Idents(thisObj.T) <- as.character(Idents(thisObj.T))

# Initialize the list to store DE results
DE_results <- list()

# Loop over each unique cluster identity in the Seurat object
for (cluster_id in levels(thisObj.T)) {
  # Print cluster_id to monitor progress
  message(paste("Processing cluster:", cluster_id))
  
  # Perform DE analysis for each cluster
  DE_results[[paste0("cluster_", cluster_id)]] <- FindMarkers(
    thisObj.T, ident.1 = cluster_id, test.use = "wilcox"
  )
  
  # Save DE results as text files for each cluster
  write.table(DE_results[[paste0("cluster_", cluster_id)]], 
              file = paste0("cluster-analysis/table_dge_cluster_", cluster_id, ".txt"), 
              sep = "\t", quote = FALSE, row.names = TRUE)
}
data <- data_prepare(file="seurat_rna_matrix.txt")
data <- data[rowSums(data) > 10, ]  # Filter low-expression genes

clust <- clustering(data = data, n.cluster = 9, n = 10, method = "simlr",write = FALSE,pdf=FALSE)
clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = clust$cluster, write = FALSE)

signal <- cell_signaling(
  data = data,
  genes = genes,
  cluster = cluster_numeric,
  c.names = unique(Idents(thisObj.T)),
  write = FALSE,
  species = "mus musculus",
  verbose=TRUE
)

inter.net <- inter_network(
  data = data,
  signal = signal, 
  genes = genes, 
  cluster = cluster_numeric,
  c.names = unique(Idents(thisObj.T)),
  write = FALSE
)






file <- "example_data.txt"
data <- data_prepare(file = file)
clust <- clustering(data = data, n.cluster = 4, n = 10, method = "simlr",write = FALSE,pdf=FALSE)
clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = clust$cluster, write = FALSE)
signal <- cell_signaling(data = data, genes = rownames(data), cluster = clust$cluster, write = FALSE)



# Trash 2####
# >>>>START ANALYSIS - PCA<<<< ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**SC Transform for noisy data*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisObj <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_preALL_newcountcutoff.RDS"))

# SCT transform offers a noise stabilizing version (V2) and a built in variable feature
# calculation. We do not/should not call FindVariableFeatures; would be redundant
# and is not appropriate for noisy data
# Increase memory allocation to handle the size of the transform
options(future.globals.maxSize = 3 * 1024^3)
thisObj <-SCTransform(thisObj,
                      assay = "RNA",
                      vst.flavor = "v2",
                      variable.features.n = 3000, #Find variable features, do not call FindVariableFeatures() if you use this
                      verbose = TRUE)


print("SCTransform potentially completed. Checking variable features...")
print(paste("Number of variable features:", length(VariableFeatures(thisObj[["SCT"]]))))
top20 <- head(VariableFeatures(thisObj[["SCT"]]), 20)
print("Top 20 Variable Features (SCT):")
print(top20)
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Plotting variable features and labeling top genes*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot1 <- VariableFeaturePlot(thisObj)
# plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
# 
# # Customize the plots (e.g., setting the legend position and axis limits)
# plot2 <- plot2 + theme(legend.position = "bottom")
# plot2_modified <- plot2 + theme(legend.position = "bottom") + ylim(0, 30) + 
#   theme(plot.title = element_text(hjust = 0.5))  # Center the title
# 
# # Use grid.arrange to combine title and plots
# grid.arrange(plot2, plot2_modified, top = textGrob(expression(bold("Cerebellum Variable Features")), rot = 0,
# gp = gpar(col = "black", fontsize = 20)), ncol=2)
#**^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**SC Cerebellum Reference dataset for marker comparison*
# Trying to improve gene list markers to better separate clusters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Download and prepare reference dataset
# CB_ref_SeuratObj <- readRDS(file.path("C:/Users/m322332/Downloads/","cb_annotated_object.RDS"))
# CB_ref_SeuratObj <- UpdateSeuratObject(CB_ref_SeuratObj) #Update from v2 => v3
# options(future.globals.maxSize = 30 * 1024^3) # Increase memory allocation for reference file, reset after running!
# set conser.memory = TRUE since scREF is too large to compute after multiple different attempts; 4/17/2025
# CB_ref_SeuratObj <- SCTransform(CB_ref_SeuratObj, vst.flavor='v2', conserve.memory = TRUE, verbose = TRUE) #Calculate SCT to match query data
# CB_ref_SeuratObj <- RunPCA(CB_ref_SeuratObj, assay = "SCT", verbose = TRUE) #Process PCA for SCT
# saveRDS(CB_ref_SeuratObj, file.path(HDD_WD, "CB_ref_SeuratObj_SCTPCA"))
CB_ref_SeuratObj <- readRDS(file.path(HDD_WD, "CB_ref_SeuratObj_SCTPCA"))

# Calculating reference features to add additional genes to CB gene list
# 4.17.2025 Note
# I am iterating through the process determining appropriate genes and PCs and gene features
# that best differentiate cell types and denoise data and ID incorrectly segmented samples
# with the intent to run a quick ReSeg on problem samples or exclude
# reference_markers <- FindAllMarkers(
#   CB_ref_SeuratObj,
#   assay = "SCT",
#   only.pos = TRUE,            # Only find positive markers for each identity
#   min.pct = 0.1,              # Detect in at least 10% of cells in the group
#   logfc.threshold = 0.25,     # Minimum log2 fold change
#   max.cells.per.ident = 1000, # Reduce samples per measure to limit memory issues, slows calculation
#   verbose = TRUE
# )

# Combine/clean markers list
reference_markers <- reference_markers %>%
  mutate(original_cluster = stringr::str_remove(cluster, "_\\d+$"))
head(reference_markers)

# saveRDS(reference_markers, file.path(HDD_WD, "CB_ref_SeuratObj_reference_genes_unfiltered"))
reference_markers <- readRDS(file.path(HDD_WD, "CB_ref_SeuratObj_reference_genes_unfiltered"))

reference_markers_top30 <- reference_markers %>%
  group_by(original_cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC, with_ties = FALSE) %>% # Get top 30 per cluster
  ungroup() 

print(reference_markers_top30, n = Inf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Cerebellum Markers*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gene_marker_lists <- list(
  # 1. Neurons (Excitatory: Granule + UBCs)
  Pan_Neuronal = c("Rbfox3", "Snap25", "Syp", "Stmn2"),  # Pan-Neuronal in list
  UBC = c("Calb2", "Grm1", "Plcb4"),            # UBC-specific
  Granule = c("Slc17a7"), # Granule Cell Markers
  
  # 2. Interneurons (Inhibitory: Golgi, Stellate, Basket, Lugaro, DCN GABAergic)
  Neuron_Inhibitory = c("Gad1", "Gad2", "Slc6a1", "Nrgn", "Sst", "Slc32a1"), # Generic GABAergic in list
  
  Interneuron = c("Pvalb", "Npy"), # Specific Interneurons
  
  # 3. Purkinje Cells
  Purkinje = c("Calb1", "Pcp2", "Grid2", "Itpr1", "Rora","Car8"),
  
  # 4. Bergmann Glia
  Bergmann_Glia = c("Sox9"),
  
  # 5. Fibroblasts (similar to VLMCs)
  Fibroblast_VLMC = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum","Apod"), # Note Pdgfra overlap with OPC
  
  # 6. Microglia
  Microglia = c("P2ry12", "Tmem119", "Cx3cr1", "C1qa", "Trem2","Ctss"),
  
  # 7. Oligodendrocytes
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag","Myrf"),
  
  # 8. Oligodendrocyte Precursor Cells (OPCs)
  OPC = c("Pdgfra", "Olig2", "Gpr17", "Tcf7l2"),
  
  # 9. Endothelial Cells
  Endothelial = c("Cldn5", "Pecam1", "Flt1", "Esam"),
  
  # 10. Astrocytes (General)
  Astrocyte = c("Gfap", "Aldh1l1", "Slc1a3", "Aqp4", "S100b"), # Note overlap with Bergmann
  
  # 11. Pericytes
  Pericyte = c("Pdgfrb", "Rgs5", "Des", "Acta2"),
  
  # 12. Vascular Smooth Muscle Cells (VSMCs)
  VSMC = c("Acta2", "Tagln", "Notch3"),
  
  # 13. Choroid Plexus
  Choroid_Plexus = c("Ttr", "Cldn5", "Otx2", "Kcne2", "Aqp1") # Note Cldn5 overlap with Endothelial
)

gene_marker_lists <- lapply(gene_marker_lists, unique)

print("Newly defined cerebellum marker list structure:")
print(str(gene_marker_lists, max.level = 1))

all_markers <- unique(unlist(gene_marker_lists))

# Compare our list against the SC CB reference data top genes
top_ref_genes <- unique(reference_markers_top30$gene)

genes_in_thisObj <- rownames(thisObj[["SCT"]])

common_top_ref_genes <- intersect(top_ref_genes, genes_in_thisObj)
reference_markers_detected <- reference_markers_top30 %>%
  filter(gene %in% common_top_ref_genes)

head(reference_markers_detected)





# Now includes top genes from reference dataset to make our list more robust
gene_marker_lists_v2 <- list(
  # 1. Neurons (Excitatory: Granule + UBCs)
  Pan_Neuronal = c("Rbfox3", "Snap25", "Syp", "Stmn2"),  # Pan-Neuronal in list
  UBC = c("Calb2", "Grm1", "Plcb4"),            # UBC-specific
  Granule = c("Slc17a7"), # Granule Cell Markers
  
  # 2. Interneurons (Inhibitory: Golgi, Stellate, Basket, Lugaro, DCN GABAergic)
  Pan_Inhibitory = c("Gad1", "Gad2", "Slc6a1", "Nrgn", "Sst", "Slc32a1"), # Generic GABAergic in list
  
  Interneuron = c("Pvalb", "Npy","Gcgr"), # Specific Interneurons
  
  Golgi = c("Gcgr","Cacng3","Gdnf","Gabrg3","Chrm2"), #golgi specific
  
  Granule = c("Chrm3","Fos","Pcsk5","Dcc","Atp2b4"), #granule specific
  
  MLI = c("Grm8","Ret","Penk","Itga8"), # Molecular Layer Interneuron Specific (MLI-1/2 combined from reference markers)
  
  # 3. Purkinje Cells
  Purkinje = c("Calb1", "Pcp2", "Grid2", "Itpr1", "Rora","Car8", "Htr2a","Ngf","Oprd1"),
  
  # 4. Bergmann Glia
  Bergmann_Glia = c("Sox9","Tnc","Gria1","Slc1a3"), # Note overlap with astrocytes
  
  # 5. Fibroblasts (similar to VLMCs)
  Fibroblast_VLMC = c("Col1a1", "Col1a2", "Pdgfra", "Dcn", "Lum","Apod"), # Note Pdgfra overlap with OPC
  
  # 6. Microglia
  Microglia = c("P2ry12", "Tmem119", "Cx3cr1", "C1qa", "Trem2","Ctss","Selplg","Trem2", "Csf1r","Ptprc"),
  
  Macrophage = c("Ms4a4a","Msr1"),
  
  # 7. Oligodendrocytes
  Oligodendrocyte = c("Mbp", "Plp1", "Mog", "Mag","Myrf","Il33","Serpinb1a","Gpr17"), #note overlap with OPC
  
  # 8. Oligodendrocyte Precursor Cells (OPCs)
  OPC = c("Pdgfra", "Olig2", "Tcf7l2"), # note overlap with ODC
  
  # 9. Endothelial Cells
  Pan_Endothelial = c("Cldn5", "Pecam1", "Flt1", "Esam","Ace2","Higd1b","Ednra","Flt1","Cldn5"),
  
  # 10. Astrocytes (General)
  Astrocyte = c("Gfap", "Aldh1l1", "Slc1a3", "Aqp4", "S100b","Slc6a11","Agt","Cd44"), # Note overlap with Bergmann
  
  # 11. Pericytes
  Pericyte = c("Pdgfrb", "Rgs5", "Des", "Acta2"),
  
  # 12. Vascular Smooth Muscle Cells (VSMCs)
  VSMC = c("Acta2", "Tagln", "Notch3"),
  
  # 13. Choroid Plexus
  Choroid_Plexus = c("Ttr", "Cldn5", "Otx2", "Kcne2", "Aqp1","Htr2c","Crhr2") # Note Cldn5 overlap with Endothelial
)

all_markers_v2 <- unique(unlist(gene_marker_lists_v2))
gene_marker_lists <- gene_marker_lists_v2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#**Linear Dimensional Reduction PCA*
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisObj<- RunPCA(thisObj,
                 assay = "SCT",
                 features = all_markers_v2, #VariableFeatures(thisObj[["SCT"]]),#all_markers_v2, #all_markers,
                 npcs = 50, 
                 verbose = TRUE)


# SCT is the better assay option since it regresses out noisy genes allowing better
# differentiation between expressions.
saveRDS(thisObj, file.path(HDD_WD, "COSMONAUT_CB_Seurat_postPCA_newcountcutoff.RDS"))

thisObj <- readRDS(file.path(HDD_WD, "COSMONAUT_CB_Seurat_postPCA_newcountcutoff.RDS"))
# Check PCA results
print(ElbowPlot(thisObj, ndims = 50)) #Elbow plot
print(VizDimLoadings(thisObj, dims = 1:4, ncol = 2, nfeatures = 20)) #View top loadings

# Heatmaps of top loading genes
DimHeatmap(thisObj, dims = 1:1, cells = 500, balanced = TRUE)
DimHeatmap(thisObj, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(thisObj, dims = 16:30, cells = 500, balanced = TRUE)

ScreePlot <- function(seurat_object, dims, annotate_thresholds) {
  # Calculate variance of each PC
  variance <- (Stdev(seurat_object, reduction = 'pca'))^2
  
  # Calculate percentage of variance explained for each PC
  pc_var_explained <- variance / sum(variance) * 100
  
  # Calculate cumulative variance
  cumvar <- cumsum(pc_var_explained)
  
  # Find the number of PCs needed to explain 80%, 90%, and 95% of variance
  n_80 <- which(cumvar >= 80)[1]
  n_90 <- which(cumvar >= 90)[1]
  n_95 <- which(cumvar >= 95)[1]
  
  # Print results
  cat("Number of PCs to explain at least:\n")
  cat("80% variance:", n_80, "\n")
  cat("90% variance:", n_90, "\n")
  cat("95% variance:", n_95, "\n")
  
  # Create a data frame for plotting
  Scree.df <- data.frame(
    PC = factor(1:dims),  # PCs as a factor for better plotting
    Variance = pc_var_explained[1:dims]  # First n_components values
  )
  
  # Start creating the plot
  p <- ggplot(Scree.df, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "black") +
    geom_point(color = "red", size = 3) +
    labs(title = "Variance Explained by Principal Components",
         x = "Principal Components",
         y = "Percentage of Variance") +
    theme_minimal()
  
  # Conditionally add annotations if annotate_thresholds is TRUE
  if (annotate_thresholds) {
    p <- p +
      # Add vertical lines for 80%, 90%, and 95% thresholds
      geom_vline(xintercept = n_80, linetype = "dashed", color = "yellow") +
      geom_vline(xintercept = n_90, linetype = "dashed", color = "orange") +
      geom_vline(xintercept = n_95, linetype = "dashed", color = "red") +
      # Annotate the plot with labels for each threshold
      annotate("text", x = n_80, y = max(pc_var_explained[1:dims]), 
               label = paste("80% (PC", n_80, ")"), color = "black", hjust = -0.1) +
      annotate("text", x = n_90, y = max(pc_var_explained[1:dims]) * 0.9, 
               label = paste("90% (PC", n_90, ")"), color = "black", hjust = -0.1) +
      annotate("text", x = n_95, y = max(pc_var_explained[1:dims]) * 0.8, 
               label = paste("95% (PC", n_95, ")"), color = "black", hjust = -0.1)
  }
  
  # Output variance thresholds
  ScreeResults <- list(
    ScreeResults = Scree.df,
    n_80 = n_80,
    n_90 = n_90,
    n_95 = n_95)
  
  
  # Print the plot
  print(p)
  return(ScreeResults)
}
ScreeResults <- ScreePlot(thisObj, dims = 50, annotate_thresholds = TRUE)

# check loadings of our researched gene markers

pca_loadings <- Loadings(thisObj[["pca"]])
markers_in_pca <- intersect(all_markers, rownames(pca_loadings))

length(VariableFeatures(thisObj[["SCT"]])) # num of all variable features
length(all_markers) # number of our gene markers
length(markers_in_pca) # number of our markers x feature list

marker_loadings <- pca_loadings[markers_in_pca, ]

top_pcs_per_marker <- apply(marker_loadings, 1, function(gene_row) {
  # Get absolute loadings, sort descending, take top N (e.g., top 3)
  sorted_loadings <- sort(abs(gene_row), decreasing = TRUE)
  names(head(sorted_loadings, 3)) # Return names of the top 3 PCs
})
print("Top 3 PCs for each marker (based on absolute loading):")
print(top_pcs_per_marker)

top_markers_per_pc <- apply(marker_loadings, 2, function(pc_col) {
  sorted_loadings <- sort(abs(pc_col), decreasing = TRUE)
  names(head(sorted_loadings, 5)) # Return names of top 5 markers for this PC
})
print("Top 5 marker genes for each PC (from your list):")
print(top_markers_per_pc[1:min(10, ncol(marker_loadings))]) # Adjust range as needed

#'[NOTE 4.11.2025: PCA inspection does NOT show clear differentiation suggesting
#'[contamination and poor segmentation. To improve SNR, we will inspect
#'[each PC/loadings and similarity of known cell markers against each PC using
#'[an assymetric similarity matrix to decide which PCs to exclude.]]
pcs_to_analyze <- 1:50
pcs_to_analyze <- intersect(pcs_to_analyze, 1:ncol(pca_loadings))
if (length(pcs_to_analyze) == 0) {
  stop("No valid PCs selected for analysis.")
}
loadings_subset <- pca_loadings[, pcs_to_analyze, drop = FALSE]
enrichment_scores <- list()

for (list_name in names(gene_marker_lists)) {
  gene_list <- gene_marker_lists[[list_name]]
  
  # Find which genes from the list are present in the PCA loadings
  markers_in_pca <- intersect(gene_list, rownames(loadings_subset))
  
  if (length(markers_in_pca) > 0) {
    # Get loadings for these markers for the selected PCs
    marker_loadings_subset <- abs(loadings_subset[markers_in_pca, , drop = FALSE])
    
    # Calculate sum of absolute loadings for this list across the PCs
    # colSums calculates the sum for each column (PC)
    scores <- colSums(marker_loadings_subset)
    enrichment_scores[[list_name]] <- scores
  } else {
    # Handle case where no markers from the list are found in PCA
    enrichment_scores[[list_name]] <- rep(0, length(pcs_to_analyze))
    names(enrichment_scores[[list_name]]) <- colnames(loadings_subset)
    warning(paste("No markers from list", list_name, "found in PCA loadings."))
  }
}

enrichment_matrix <- do.call(rbind, enrichment_scores)
colnames(enrichment_matrix) <- paste0("PC_", pcs_to_analyze) # Ensure column names match PCs
print(enrichment_matrix)

top_pc_per_list <- apply(enrichment_matrix, 1, function(row) {
  colnames(enrichment_matrix)[which.max(row)]
})

print("Top PC for each gene list (based on sum of absolute loadings):")
print(top_pc_per_list)
plot.new()
PC_Heat_ColScale <- pheatmap(enrichment_matrix,
                             scale = "column",
                             main = "Relative Gene List Enrichment per PC (Scaled by Column)",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE)


file_name <- paste0(local_WD, "PC_Heat_ColScale", ".png")
png(file = file_name, width = 10, height = 5, units = "in", res = 900)
PC_Heat_ColScale
dev.off() 

# Scale by row (emphasizes top PC(s) per list)
pheatmap(enrichment_matrix,
         scale = "row",
         main = "Relative PC Enrichment per Gene List (Scaled by Row)",
         cluster_rows = TRUE,
         cluster_cols = TRUE)


# Ridge plot version
library(ggplot2)
library(ggridges)
library(tidyr)
library(dplyr)
library(tibble) # For rownames_to_column


# 1. Convert matrix to long data frame suitable for ggplot
enrichment_long <- enrichment_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_list") %>%
  pivot_longer(
    cols = -gene_list,
    names_to = "component",
    values_to = "enrichment_score"
  ) %>%
  # Extract numeric component number for plotting on x-axis
  mutate(
    component_num = as.numeric(gsub("PC_", "", component))
  )

# Ensure gene_list is treated as a factor for ordered plotting
enrichment_long$gene_list <- factor(enrichment_long$gene_list, levels = rev(rownames(enrichment_matrix)))

# 2. Create the ridge plot
ridge_plot <- ggplot(enrichment_long, aes(x = component_num, y = gene_list, height = enrichment_score, group = gene_list, fill = gene_list)) +
  geom_ridgeline(aes(y = as.numeric(factor(gene_list)) - 0.5), scale = 0.9, alpha = 0.7) +
  scale_y_continuous(
    breaks = (1:length(levels(factor(enrichment_long$gene_list)))) - 0.5,
    labels = levels(factor(enrichment_long$gene_list)),
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(
    breaks = sort(unique(enrichment_long$component_num)),
    labels = sort(unique(enrichment_long$component_num))
  ) +
  labs(
    title = "Gene List Enrichment Score Across Components",
    x = "Principal Component",
    y = "Gene List"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8) # Rotate labels
  )

print(ridge_plot)

# PCs with high scores for a SINGLE list are likely good biological PCs.
# PCs with high scores for MULTIPLE lists simultaneously (especially oligoList + neuronList,
# or other known contaminating pairs) are candidates for EXCLUSION.
# PCs with low scores across ALL lists might be noisy or represent biology not
# captured by your marker lists, requiring inspection of the actual top loading genes.

#'[Added additional metrics to decide which PC's to exclude; 4.18.2025]
# Sum of squares of gene sets to each PC, tells us how "noisy" a PC is relative to each gene set

enrichment_scores_sq <- list()
gene_marker_lists <- gene_marker_lists_v2

for (list_name in names(gene_marker_lists)) {
  gene_list <- gene_marker_lists[[list_name]]
  markers_in_pca <- intersect(gene_list, rownames(loadings_subset))
  
  if (length(markers_in_pca) > 0) {
    # Get squared loadings for these markers
    marker_loadings_subset_sq <- (loadings_subset[markers_in_pca, , drop = FALSE])^2
    # Calculate sum of squared loadings for this list across the PCs
    scores <- colSums(marker_loadings_subset_sq)
    enrichment_scores_sq[[list_name]] <- scores
  } else {
    # Handle case where no markers from the list are found in PCA
    enrichment_scores_sq[[list_name]] <- rep(0, length(pcs_to_analyze))
    names(enrichment_scores_sq[[list_name]]) <- colnames(loadings_subset)
    # warning(paste("No markers from list", list_name, "found in PCA loadings.")) # Optional warning
  }
}

enrichment_matrix_sq <- do.call(rbind, enrichment_scores_sq)

# Calculate total score per PC
total_score_per_pc <- colSums(enrichment_matrix_sq)
# Avoid division by zero if a PC has zero total score from markers
total_score_per_pc[total_score_per_pc == 0] <- 1e-9

# Calculate specificity matrix (each cell = Score(List, PC) / TotalScore(PC))
# Use sweep to divide each row of enrichment_matrix_sq by total_score_per_pc
specificity_matrix <- sweep(enrichment_matrix_sq, 2, total_score_per_pc, FUN = "/")

# Visualize Specificity Matrix
pheatmap(specificity_matrix,
         scale = "none", # Already scaled as proportions
         main = "PC Specificity Score per Gene List",
         cluster_rows = FALSE,
         cluster_cols = FALSE, # Cluster PCs based on specificity profiles
         fontsize_row = 8,
         fontsize_col = 8)


#'[4.18.2025; PCs are still non-specific especially for neuron cell types]
#'[Issue is expected since neurons/granule layer is the most contaminated]

#'[Attemping ICA instead to try to deconvolute noise]
#'
thisObj<- RunICA(thisObj,
                 assay = "SCT",
                 features = VariableFeatures(thisObj[["SCT"]]),#all_markers_v2, #all_markers,
                 nics = 50, 
                 verbose = TRUE)

ics_to_analyze <- 1:50
ica_loadings <- Loadings(thisObj[["ica"]])
ics_to_analyze <- intersect(ics_to_analyze, 1:ncol(ica_loadings))
loadings_subset <- ica_loadings[, ics_to_analyze, drop = FALSE]
enrichment_scores <- list()

for (list_name in names(gene_marker_lists)) {
  gene_list <- gene_marker_lists[[list_name]]
  markers_in_ica <- intersect(gene_list, rownames(loadings_subset))
  if (length(markers_in_ica) > 0) {
    marker_loadings_subset <- abs(loadings_subset[markers_in_ica, , drop = FALSE])
    scores <- colSums(marker_loadings_subset)
    enrichment_scores[[list_name]] <- scores
  } else {
    enrichment_scores[[list_name]] <- rep(0, length(ics_to_analyze))
    names(enrichment_scores[[list_name]]) <- colnames(loadings_subset)
  }
}

enrichment_matrix <- do.call(rbind, enrichment_scores)
colnames(enrichment_matrix) <- paste0("IC_", ics_to_analyze)

top_ic_per_list <- apply(enrichment_matrix, 1, function(row) {
  colnames(enrichment_matrix)[which.max(row)]
})

library("pheatmap")
pheatmap(enrichment_matrix,
         scale = "column",
         main = "Relative Gene List Enrichment per IC (Scaled by Column)",
         cluster_rows = TRUE,
         cluster_cols = FALSE)

pheatmap(enrichment_matrix,
         scale = "row",
         main = "Relative IC Enrichment per Gene List (Scaled by Row)",
         cluster_rows = TRUE,
         cluster_cols = FALSE)
#'[Attempting ICA reduction for UMAP]

ic2use <- c(1,2,3,4,5,7,8,10,11,14,15,17,23,30,33,38,47) # 23 maybe
ic2use <- c(14, 12, 15, 17, 29, 49, 40, 36, 19, 11, 25, 4, 10, 5, 2, 3 )

thisObj <- FindNeighbors(thisObj, dims = ic2use, reduction = "ica",  assay = "SCT")
thisObj <- FindClusters(thisObj, algorithm = 3, reduction = "ica", resolution = 0.75)
thisObj <- RunUMAP(thisObj, dims = dims_to_use,  assay = "SCT", reduction = "ica", min.dist = 0.0001, spread = 4)

DimPlot(thisObj, reduction = "umap", raster = TRUE, label = TRUE, pt.size = 0.5) + ggtitle("UMAP")


#'[Conclusion 4.18.2025, ICA ineffictive for this dataset, sticking with PCA]


















































# Trash####
#neigbors 25
#min.dist 0.03
#spread 0.5
# thisObj.T <- RunUMAP(thisObj.T, dims = 1:50, min.dist = 0.05, spread = 0.5, n.neighbors = 50)
# DimPlot(thisObj.T)
# 
# #dims 13 w/ new range possible
# #mindist 0.015
# 
# #new new, nNeighbor=19
# 
# thisObj.T2 <- RunUMAP(thisObj.T, n.neighbors = 20, dims = 1:50, min.dist = 0.008, spread = 0.5)
# DimPlot(thisObj.T2, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP with Cell Type Labels")