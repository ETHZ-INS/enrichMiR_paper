# Functions and objects for the benchmark -- Tomás Germade

TPs <- c( 
  let.7a = "GAGGUAG", lsy.6 = "UUUGUAU", miR.1 = "GGAAUGU", miR.124 = "AAGGCAC", 
  miR.137 = "UAUUGCU", miR.139 = "CUACAGU", miR.143 = "GAGAUGA", 
  miR.144 = "ACAGUAU", miR.153 = "UGCAUAG", miR.155 = "UAAUGCU", 
  miR.182 = "UUGGCAA", miR.199a = "CCAGUGU", miR.204 = "UCCCUUU", 
  miR.205 = "CCUUCAU", miR.216b = "AAUCUCU", miR.223 = "GUCAGUU", 
  miR.7 = "GGAAGAC", miR.122 = "GGAGUGU", miR.133 = "UGGUCCC", 
  miR.138 = "GCUGGUG", miR.145 = "UCCAGUU", miR.184 = "GGACGGA", 
  miR.190a = "GAUAUGU", miR.200b = "AAUACUG", miR.216a = "AAUCUCA", 
  miR.217 = "ACUGCAU", miR.219a = "GAUUGUC", miR.375 = "UUGUUCG", 
  miR.451a = "AACCGUU", "DKOvWT"="GCUACAU", "DKO"="GCUACAU", 
  "218DKOvWT"="UGUGCUU", "218DKO"="UGUGCUU",
  "miR.138vNeg"="GCUGGUG", "miR.499vNeg"="UAAGACU", 
  "miR.138"="GCUGGUG", "miR.499"="UAAGACU")


#################################################
#' TSperm
#'
#' @param TS TargetScan DataFrame of miRNA tx targets
#' @param TP string of miRNA seed (true-positive), eg. "ACAUCAC"
#' @param genes string vector of gene names
#' @param props vector of proportions, e.g. c(.2,.3,.4,.5)
#'
#' @return TS DataFrame with permuted miRNA targets
#' 
TSperm <- function(TS, TP, genes, props){
  # get TargetScan data for each treatment miRNA family
  TS.part <- TS[TS$set==TP,]
  lapply(props, FUN=function(p){
    sn <- floor(p*nrow(TS.part))
    genes <- genes[!genes %in% TS.part$feature]
    i <- sample(seq_len(nrow(TS.part)), sn)
    j <- sample(seq_len(length(genes)), sn)
    TS.part$feature <- as.character(TS.part$feature)
    TS.part$feature[i] <- as.character(genes[j])
    TS <- rbind(TS.part, TS[TS$set!=TP,])
    return(TS)
  })
}

#################################################
#' getDF
#'
#' @param dea.names 
#' @param e.list 
#' @param TP 
#'
#' @return
#' @export
#'
#' @examples
getDF <- function(dea.names, e.list, TP){ 
  ### summarize results
  #### generate the benchmarking scores
  BM.list <- lapply(setNames(dea.names,dea.names), FUN=function(x){
    lapply(e.list[[x]], FUN=function(y){
      doBenchmark(y, TP[x])
    })
  })
  #### generate a results df for plotting)
  BM.list2 <- lapply(BM.list, FUN=function(x)
                        as.data.frame(dplyr::bind_rows(x, .id = "prop.rep")))
  BM.df <- dplyr::bind_rows(BM.list2, .id="treatment")
  BM.df$prop <- unlist(lapply(strsplit(BM.df$prop.rep, "[.]"), FUN=function(x) x[1]))
  BM.df$prop <- factor(BM.df$prop, levels=unique(BM.df$prop))
  
  return(BM.df)
}

#################################################
#' doBenchmark
#'
#' @param res enrichMiR test results (e object)
#' @param TP character vector: contains the seeds (families) for a miRNA treatment
#'
#' @return a dataframe containing scores for each enrichMiR test
#'
doBenchmark <- function(res, TP){
  res <- lapply(res, FUN=function(x){
    x <- x[order(x$FDR),]
    if("family" %in% colnames(x)) row.names(x) <- x$family
    x$truth <- row.names(x) %in% TP
    x$FDR[is.na(x$FDR)] <- 1
    x
  })
  data.frame( method=names(res),
              sets = sapply(res, FUN = function(x) {x$sets[1]}),
              detPPV = sapply(res, FUN=function(x) 1/which(x$truth)[1] ),
              FP.atFDR05 = sapply(res, FUN=function(x) sum(!x$truth & x$FDR<0.05)),
              log10QDiff = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                tp1-fp1
              }),
              log10QrelIncrease = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                (tp1-fp1)/(min(tp1,fp1))
              }),
              TP.atFDR05 = sapply(res, FUN=function(x) sum(x$truth[x$FDR<0.05]))
  )
  
}



#################################################
#' getDF
#'
#' @param dea.names 
#' @param e.list 
#' @param TP 
#'
#' @return
#' @export
#'
#' @examples
getDF2 <- function(dea.names, e.list, TP){ 
  ### summarize results
  #### generate the benchmarking scores
  BM.list <- lapply(setNames(dea.names,dea.names), FUN=function(x){
    lapply(e.list[[x]], FUN=function(y){
      doBenchmark2(y, TP[x])
    })
  })
  #### generate a results df for plotting)
  BM.list2 <- lapply(BM.list, FUN=function(x)
    as.data.frame(dplyr::bind_rows(x, .id = "prop.rep")))
  BM.df <- dplyr::bind_rows(BM.list2, .id="treatment")
  BM.df$prop <- unlist(lapply(strsplit(BM.df$prop.rep, "[.]"), FUN=function(x) x[1]))
  BM.df$prop <- factor(BM.df$prop, levels=unique(BM.df$prop))
  
  return(BM.df)
}

#################################################
#' doBenchmark
#'
#' @param res enrichMiR test results (e object)
#' @param TP character vector: contains the seeds (families) for a miRNA treatment
#'
#' @return a dataframe containing scores for each enrichMiR test
#'
doBenchmark2 <- function(res, TP){
  res <- lapply(res, FUN=function(x){
    x <- x[order(x$FDR),]
    if("family" %in% colnames(x)) row.names(x) <- x$family
    x$truth <- row.names(x) %in% TP
    x$FDR[is.na(x$FDR)] <- 1
    x
  })
  data.frame( method=names(res),
              detPPV = sapply(res, FUN=function(x) 1/which(x$truth)[1] ),
              FP.atFDR05 = sapply(res, FUN=function(x) sum(!x$truth & x$FDR<0.05)),
              log10QDiff = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                tp1-fp1
              }),
              log10QrelIncrease = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                (tp1-fp1)/(min(tp1,fp1))
              }),
              TP.atFDR05 = sapply(res, FUN=function(x) sum(x$truth[x$FDR<0.05]))
  )
  
}




enrichBenchmark <- function( params, tests=c("siteoverlap", "areamir"),
                             props=c(.2,.35,.5), nrep=3, cores=6, seed=1234){
  names(i) <- i <- 1:nrep
  if(is.null(names(props))) names(props) <- round(100*props)
  props <- unlist(lapply( props, FUN=function(x) lapply(i, FUN=function(y) x)))
  
  lapply(params, FUN=function(p){
    p$name <- gsub("\\.rds","",p$ds)
    if(!is.null(p$miRNAs)) p$name <- paste0(p$name,".mirexp")
    print(paste(p$name, "using", p$ts))
    RD <- rowData(readRDS(file.path("data", p$ds)))
    dea.df.names <- colnames(RD[,grepl("DEA",colnames(RD))])
    dea.names <- gsub("DEA\\.","", dea.df.names)
    dea.names <- gsub("spliced\\.","", dea.names)
    dea.names <- gsub("^\\.","", dea.names)
    w <- which(dea.names %in% names(TPs))
    if(length(w)==0) stop("Unknown DEAs (", paste(dea.names, collapse=", "), ")")
    dea.names <- dea.names[w]
    dea.df.names <- dea.df.names[w]
    #### get list of all DEAs
    dea.list <- lapply(RD[,dea.df.names,drop=FALSE], FUN=function(x){
      if(sum(grepl("^ENSG|^ENSRNOG|^ENSMUSG", row.names(x)))>100){
        x <- x[order(x$FDR, x$PValue),]
        g <- gsub("^ENS[^\\.]+\\.","",row.names(x))
        x <- x[!duplicated(g),]
        row.names(x) <- g[!duplicated(g)]
      }
      x
    })
    names(dea.list) <- dea.names
    #### specify expressed miRNA (celltype-specific)
    mirexpr <- p$miRNAs
    #### load TS object
    TS <- readRDS( file.path("data",p$ts) )
    TS$sites <- round(TS$sites)
    #### get true positive miRNA treatments
    TP <- TPs[dea.names]
    ### enrichMiR
    
    e.list <- bplapply(dea.list, BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                        FUN=function(dea){
       set.seed(seed)
       dea <- as.data.frame(dea)
       eres <- lapply(c(list(TS),TSperm(TS, TP, row.names(dea), props)), FUN=function(ts){
         e <- enrichMiR(DEA=dea, TS=TS, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)
         lapply(e@res, FUN=function(x){
             x$features <- NULL
             x$sets <- p$ts
             x
         })
       })
       names(eres) <- c("original",names(props))
       eres
    })
    attr(e.list, "TP") <- TP
    saveRDS(e.list, file=paste0("./sets/",gsub("DEA\\.SE","enrichMiR_sets",p$name),"_",p$ts,".rds"))
    ### summarize enrichMiR results in a dataframe
    bench <- getDF(names(TP), e.list, TP)
    saveRDS(bench, file=paste0("./sets/",gsub("DEA\\.SE","benchmark_sets",p$name),"_",p$ts,".rds"))
    NULL
  })
}


enrichBenchmark2 <- function( params, tests=c("siteoverlap"),
                             props=c(.2,.35,.5), nrep=3, cores=6, seed=1234){
  names(i) <- i <- 1:nrep
  if(is.null(names(props))) names(props) <- round(100*props)
  props <- unlist(lapply( props, FUN=function(x) lapply(i, FUN=function(y) x)))
  
  lapply(params, FUN=function(p){
    p$name <- gsub("\\.rds","",p$ds)
    if(!is.null(p$miRNAs)) p$name <- paste0(p$name,".mirexp")
    print(paste(p$name, "using", p$ts))
    RD <- rowData(readRDS(file.path("data", p$ds)))
    dea.df.names <- colnames(RD[,grepl("DEA",colnames(RD))])
    dea.names <- gsub("DEA\\.","", dea.df.names)
    dea.names <- gsub("spliced\\.","", dea.names)
    dea.names <- gsub("^\\.","", dea.names)
    w <- which(dea.names %in% names(TPs))
    if(length(w)==0) stop("Unknown DEAs (", paste(dea.names, collapse=", "), ")")
    dea.names <- dea.names[w]
    dea.df.names <- dea.df.names[w]
    #### get list of all DEAs
    dea.list <- lapply(RD[,dea.df.names,drop=FALSE], FUN=function(x){
      if(sum(grepl("^ENSG|^ENSRNOG|^ENSMUSG", row.names(x)))>100){
        x <- x[order(x$FDR, x$PValue),]
        g <- gsub("^ENS[^\\.]+\\.","",row.names(x))
        x <- x[!duplicated(g),]
        row.names(x) <- g[!duplicated(g)]
      }
      x
    })
    names(dea.list) <- dea.names
    #### specify expressed miRNA (celltype-specific)
    mirexpr <- p$miRNAs
    #### load TS object
    TS <- readRDS( file.path("data",p$ts) )
    TS$sites <- round(TS$sites)
    #### get true positive miRNA treatments
    TP <- TPs[dea.names]
    ### enrichMiR
    
    e.list <- bplapply(setNames(dea.names,dea.names), BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                       FUN=function(dea.name){
                         set.seed(seed)
                         dea <- as.data.frame(dea.list[[dea.name]])
                         tp <- TP[dea.name]
                         eres <- lapply(c(list(TS),TSperm(TS, tp, row.names(dea), props)), FUN=function(ts){
                           e <- enrichMiR(DEA=dea, TS=ts, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)
                           lapply(e@res, FUN=function(x){
                             x$features <- NULL
                             x
                           })
                         })
                         names(eres) <- c("original",names(props))
                         eres
                       })
    
    e.list2 <- bplapply(setNames(dea.names,dea.names), BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                        FUN=function(dea.name){
                          set.seed(seed)
                          dea <- as.data.frame(dea.list[[dea.name]])
                          tp <- TP[dea.name]
                          back <- setdiff(unique(TS$feature), row.names(dea))
                          dea2 <- data.frame(row.names=back, logFC=rep(0,length(back)), 
                                             PValue=rep(0.95,length(back)), FDR=rep(1,length(back)))
                          dea <- rbind(dea[,colnames(dea2)], dea2)
                          eres2 <- lapply(c(list(TS),TSperm(TS, tp, row.names(dea), props)), FUN=function(ts){
                            e <- testEnrichment(x=dea2, sets=ts, sets.properties=mirexpr, tests=tests)
                            lapply(e@res, FUN=function(x){
                              x$features <- NULL
                              x
                            })
                          })
                          names(eres2) <- c("original",names(props))
                          names(eres2) <- paste0(names(eres2),"_wo_bg")
                          eres2
                        })
    e.listc <- Map(c, e.list, e.list2)
    attr(e.listc, "TP") <- TP
    saveRDS(e.listc, file=paste0("./back/",gsub("DEA\\.SE","enrichMiR_back",p$name),".rds"))
    ### summarize enrichMiR results in a dataframe
    bench <- getDF2(names(TP), e.listc, TP)
    saveRDS(bench, file=paste0("./back/",gsub("DEA\\.SE","benchmark_back",p$name),".rds"))
    NULL
  })
}

# miRNAs (mirexpr vector generation)
getMirexpr <- function(tissue){
  library(microRNAome)
  data("microRNAome")
  counts <- rowMedians(cpm(calcNormFactors(DGEList(assay(microRNAome)[,microRNAome$cell_tissue==tissue]))))
  names(counts) <- rownames(microRNAome)
  return(counts[counts>=10])
}

allMirExprs <- function(){
  library(edgeR)
  mousemirna <- readRDS("data/miRNA_celltype_mouseBrain_GSE30286.SE.rds")
  mouse <- rowMedians(cpm(calcNormFactors(DGEList(assay(mousemirna)[,mousemirna$type!="Cerebellum"]))))
  names(mouse) <- rownames(mousemirna)
  mouse <- mouse[mouse>=10]
  mouse.df <- data.frame(expression=mouse, id=sub("\\*","",names(mouse)), row.names=names(mouse))
  mouse.df <- aggregate(mouse.df$expression, by=list(mouse.df$id), FUN="max")
  mouse.df <- data.frame(expression=rep(mouse.df$x,3),
                         id=c(as.character(mouse.df$Group.1),
                              sapply(mouse.df$Group.1, function(x) paste0(x,"-3p")),
                              sapply(mouse.df$Group.1, function(x) paste0(x,"-5p"))) )
  mouse <- structure(mouse.df$expression, names=as.character(mouse.df$id))
  
  ratmirna <- read.csv("data/DIV20_Neuron_HC_miRNA.csv")
  list(
    hela=c(getMirexpr("hela"), structure(rep(1000,16), names=c(
    "hsa-let-7a-5p","hsa-miR-1-3p","hsa-miR-124-3p.1","hsa-miR-137","hsa-miR-139-5p",
    "hsa-miR-143-3p","hsa-miR-144-3p","hsa-miR-153-3p","hsa-miR-155-5p","hsa-miR-182-5p",
    "hsa-miR-199a-5p","hsa-miR-204-5p","hsa-miR-205-5p","hsa-miR-216b-5p","hsa-miR-223-3p",
    "hsa-miR-7-5p")) ),
    hek=c(getMirexpr("hek293"), structure(rep(1000,12), names=c(
    "hsa-miR-122-5p","hsa-miR-133a-3p.1","hsa-miR-138-5p","hsa-miR-145-5p",
    "hsa-miR-184","hsa-miR-190a-5p","hsa-miR-200b-3p","hsa-miR-216a-5p",
    "hsa-miR-217","hsa-miR-219a-5p","hsa-miR-375","hsa-miR-451a")) ),
    rat=setNames(2^(ratmirna$logCPM), ratmirna$Gene.name),
    mouse=mouse
  )
}



aggregateBenchmarkResults <- function(files=NULL){
  if(is.null(files)) files <- list.files(path = "./sets",pattern="benchmark_sets",full.names = TRUE)
  names(files) <- gsub("\\.benchmark_sets","",files)
  names(files) <- gsub("\\.rds","",names(files))
  e <- lapply(files, FUN=readRDS)
  e2 <- dplyr::bind_rows(e, .id="dataset")
  e <- aggregate(e2[,c("detPPV","FP.atFDR05","TP.atFDR05")],
                 by=e2[,c("dataset","treatment","prop","method","sets")], na.rm=TRUE, FUN=mean)
  e <- e[e$treatment!="lsy.6" & e$dataset!="jeong",]
  e$mirexp <- grepl("mirexp", e$dataset)
  e$dataset <- gsub("\\.mirexp","",e$dataset)
  e$dataset <- gsub("\\./sets/","",e$dataset)
  e$dataset <- sapply(strsplit(e$dataset,"_"),`[`,1)
  e$rank <- 1/e$detPPV
  e$FDR <- e$FP.atFDR05/(e$FP.atFDR05+e$TP.atFDR05)
  e$FDR[is.nan(e$FDR)] <- 0
  e
}

aggregateBenchmarkResults2 <- function(files=NULL){
  if(is.null(files)) files <- list.files(path = "./back",pattern="benchmark_back",full.names = TRUE)
  names(files) <- gsub("\\.benchmark_back","",files)
  names(files) <- gsub("\\.rds","",names(files))
  e <- lapply(files, FUN=readRDS)
  e2 <- dplyr::bind_rows(e, .id="dataset")
  e <- aggregate(e2[,c("detPPV","FP.atFDR05","TP.atFDR05")],
                 by=e2[,c("dataset","treatment","prop.rep","method")], na.rm=TRUE, FUN=mean)
  e <- e[e$treatment!="lsy.6" & e$dataset!="jeong",]
  e$mirexp <- grepl("mirexp", e$dataset)
  e$dataset <- gsub("\\.mirexp","",e$dataset)
  e$exonic <- grepl("exon", e$dataset)
  e$dataset <- gsub("\\.exon","",e$dataset)
  e$dataset <- gsub("\\./back/","",e$dataset)
  e$background <- !grepl("_wo_bg",e$prop.rep)
  e$rank <- 1/e$detPPV
  e$FDR <- e$FP.atFDR05/(e$FP.atFDR05+e$TP.atFDR05)
  e$FDR[is.nan(e$FDR)] <- 0
  e
}
