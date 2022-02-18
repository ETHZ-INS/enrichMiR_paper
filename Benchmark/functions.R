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

TPs2 <- c(
  miR.122 = "hsa-miR-122-5p", miR.133 = "hsa-miR-133a-3p", miR.138 = "hsa-miR-138-5p",
  miR.145 = "hsa-miR-145-5p", miR.184 = "hsa-miR-184", miR.190a = "hsa-miR-190b-5p", 
  miR.200b = "hsa-miR-200b-3p", miR.216a = "hsa-miR-216a-5p", miR.217 = "hsa-miR-217-5p", 
  miR.219a = "hsa-miR-219a-5p", miR.375 = "hsa-miR-375-3p", miR.451a = "hsa-miR-451a",
  let.7a = "hsa-let-7a-5p", miR.1 = "hsa-miR-1-3p", miR.124 = "hsa-miR-124-3p", 
  miR.137 = "hsa-miR-137-3p", miR.139 = "hsa-miR-139-5p", miR.143 = "hsa-miR-143-3p",
  miR.144 = "hsa-miR-144-3p", miR.153 = "hsa-miR-153-3p", miR.155 = "hsa-miR-155-5p", 
  miR.182 = "hsa-miR-182-5p", miR.204 = "hsa-miR-204-5p", miR.205 = "hsa-miR-205-5p",
  miR.216b = "hsa-miR-216b-5p", miR.223 = "hsa-miR-223-3p", miR.7 = "hsa-miR-7-5p",
  "218DKOvWT"="mmu-miR-218-5p")


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
  TS.part <- TS[(TS$set %in% TP),]
  print(dim(TS.part))
  lapply(props, FUN=function(p){
    sn <- floor(p*nrow(TS.part))
    genes <- genes[!genes %in% TS.part$feature]
    i <- sample(seq_len(nrow(TS.part)), sn)
    j <- sample(seq_len(length(genes)), sn)
    TS.part$feature <- as.character(TS.part$feature)
    TS.part$feature[i] <- as.character(genes[j])
    TS2 <- rbind(TS.part, TS[!(TS$set %in% TP),])
    metadata(TS2) <- metadata(TS)
    return(TS2)
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
      doBenchmark(y, TP[[x]])
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


enrichBenchmark <- function( params, tests=c("siteoverlap", "areamir"), prefix="",
                             props=c(.2,.35,.5), nrep=3, cores=6, seed=1234){
  names(i) <- i <- 1:nrep
  if(is.null(names(props))) names(props) <- round(100*props)
  props <- unlist(lapply( props, FUN=function(x) lapply(i, FUN=function(y) x)))
  
  lapply(params, FUN=function(p){
    p$name <- gsub("\\.rds","",p$ds)
    p$suffix <- gsub("\\.rds$","",basename(p$ts))
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
    
    fams <- metadata(TS)$families

    ### enrichMiR

    e.list <- bplapply(setNames(dea.names,dea.names),
                       BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                        FUN=function(dea.name){
       set.seed(seed)
       dea <- as.data.frame(dea.list[[dea.name]])
       tp <- TP[dea.name]
       if(!is.null(mirexpr)){
         #### ensure that the TP is among the expressed miRNAs (for transfections)
         tpe <- names(fams)[fams %in% tp]
         tpe <- setNames(rep(1000,length(tpe)),tpe)
         mirexpr <- c(mirexpr, tpe[setdiff(names(tpe),names(mirexpr))])
       }
       eres <- lapply(c(list(original=TS),TSperm(TS, tp, row.names(dea), props)), FUN=function(ts){
         #e <- enrichMiR(DEA=dea, TS=ts, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)
         e <- testEnrichment(dea, sets=ts, sets.properties=mirexpr, tests=tests)
         lapply(e@res, FUN=function(x){
             x$features <- NULL
             x
         })
       })
       attr(eres, "TP") <- tp
       eres
    })
    saveRDS(e.list, file=paste0(prefix,gsub("DEA\\.SE","enrichMiR",p$name),p$suffix,".rds"))
    ### summarize enrichMiR results in a dataframe
    bench <- getDF(names(TP), e.list, TP)
    saveRDS(bench, file=paste0(prefix,gsub("DEA\\.SE","benchmark",p$name),p$suffix,".rds"))
    NULL
  })
}


# miRNAs (mirexpr vector generation)
getMirexpr <- function(tissue){
  library(microRNAome)
  data("microRNAome")
  counts <- rowMedians(cpm(calcNormFactors(DGEList(assay(microRNAome)[,microRNAome$cell_tissue==tissue]))))
  names(counts) <- rownames(microRNAome)
  counts <- counts[counts>0]
  counts[counts>=median(counts)]
}

allMirExprs <- function(){
  library(edgeR)
  mousemirna <- readRDS("data/miRNA_celltype_mouseBrain_GSE30286.SE.rds")
  mouse <- rowMedians(cpm(calcNormFactors(DGEList(assay(mousemirna)[,mousemirna$type!="Cerebellum"]))))
  names(mouse) <- rownames(mousemirna)
  mouse <- mouse[mouse>0]
  mouse <- mouse[mouse>=median(mouse)]
  mouse.df <- data.frame(expression=mouse, id=sub("\\*","",names(mouse)), row.names=names(mouse))
  mouse.df <- aggregate(mouse.df$expression, by=list(mouse.df$id), FUN="max")
  mouse.df <- data.frame(expression=rep(mouse.df$x,3),
                         id=c(as.character(mouse.df$Group.1),
                              sapply(mouse.df$Group.1, function(x) paste0(x,"-3p")),
                              sapply(mouse.df$Group.1, function(x) paste0(x,"-5p"))) )
  mouse <- structure(mouse.df$expression, names=as.character(mouse.df$id))
  
  ratmirna <- read.csv("data/DIV20_Neuron_HC_miRNA.csv")
  list(
    hela=getMirexpr("hela"),
    hek=getMirexpr("hek293"),
    rat=setNames(2^(ratmirna$logCPM), ratmirna$Gene.name),
    mouse=mouse
  )
}



aggregateBenchmarkResults <- function(files=NULL){
  if(is.null(files)) files <- list.files(pattern="benchmark")
  if(all(unlist(lapply(files,is.character)))){
    names(files) <- gsub("\\.benchmark","",files)
    names(files) <- gsub("\\.rds","",names(files))
    e <- lapply(files, FUN=readRDS)
  }else{
    e <- files
  }
  e2 <- dplyr::bind_rows(e, .id="dataset")
  e <- aggregate(e2[,c("detPPV","FP.atFDR05","TP.atFDR05")],
                 by=e2[,c("dataset","treatment","prop","method")], na.rm=TRUE, FUN=mean)
  e <- e[e$treatment!="lsy.6" & e$dataset!="jeong",]
  e$mirexp <- grepl("mirexp", e$dataset)
  e$dataset <- gsub("\\.mirexp","",e$dataset)
  e$exonic <- grepl("exon", e$dataset)
  e$dataset <- gsub("\\.exon","",e$dataset)
  e$rank <- 1/e$detPPV
  e$FDR <- e$FP.atFDR05/(e$FP.atFDR05+e$TP.atFDR05)
  e$FDR[is.nan(e$FDR)] <- 0
  e
}


addTestCombinations <- function(x){
  x <- addTestCombination(x)
  x <- addTestCombination(x, combName="comb.geoMean", agfn=function(x){
    x[which(x==0)] <- 10^-300
    10^mean(log10(x),na.rm=TRUE)
  })
  x <- addTestCombination(x, agfn=function(x) median(x,na.rm=TRUE), combName="comb.medianP")
  x <- addTestCombination(x, agfn=function(x){
    x[which(x==0)] <- 10^-300
    x[which(x==1)] <- 0.99
    tryCatch(as.numeric(metap::sumz(x)$p), error=function(e) median(x,na.rm=TRUE))
  }, combName="comb.sumz")
  x
}
addTestCombination <- function(x, agfn=function(x) aggregation::fisher(x[!is.na(x)]), combName="comb.fisher", siguse="pvalue"){
  x1 <- merge(x$siteoverlap.down[,c("enrichment",siguse)],
              x$siteoverlap.up[,c("enrichment",siguse)],
              by="row.names",all=TRUE)
  x1$wh <- apply(x1[,grep(siguse,colnames(x1))],1,FUN=which.min)
  x2 <- x$siteoverlap.down[x1$Row.names[which(x1$wh==1)],c("enrichment",siguse)]
  x2$enrichment <- -1*x2$enrichment
  x2 <- rbind(x2,x$siteoverlap.up[x1$Row.names[which(x1$wh==2)],c("enrichment",siguse)])
  m <- merge(x$areamir[,c("enrichment",siguse)], x2, by="row.names")
  m <- merge(m, x$lmadd[,c("combined.coef",ifelse(siguse=="pvalue","combined.pvalue",siguse))], 
             by.x="Row.names", by.y="row.names", all=TRUE)
  row.names(m) <- m$Row.names
  m <- m[,-1]
  m <- data.frame(row.names=row.names(m),
                  enrichment=rowMeans(m[,grep(siguse,colnames(m),invert=TRUE)]),
                  pvalue=apply(m[,grep(siguse,colnames(m))],1,FUN=agfn))
  if(siguse=="FDR"){
    m$FDR <- m$pvalue
  }else{
    m$FDR <- p.adjust(m$pvalue)
  }
  x[[combName]] <- m[order(m$pvalue),]
  x
}


bmhm <- function(e, value.var="rank", by=NULL, doSort=FALSE, ..., 
                 show_row_names=FALSE, show_row_dend=FALSE,
                 cluster_columns=(is.null(by) && !doSort)){
  hlm <- list()
  cols <- viridis::inferno(10)
  if(value.var=="rank"){
    breaks <- c(1,2,3,4,5,10)
    hlm <- list(at=sqrt(breaks), break_dist=rep(1,length(breaks)-1), title="rank\nof true\nmiRNA",
                labels=c("top",breaks[2:(length(breaks)-1)],paste0(rev(breaks)[1],"+")))
    cols <- circlize::colorRamp2(sqrt(breaks), colors = viridis::viridis(6))
  }
  doSplit <- by
  if(!is(by, "formula")){
    if(is.null(by)){
      by <- dataset+treatment~method
    }else{
      by <- as.formula(paste0("dataset+treatment~method+",by))
    }
  }
  m <- reshape2::dcast(e, by, fun.aggregate = mean, value.var=value.var)
  row.names(m) <- paste(m$dataset,m$treatment,sep=":")
  ds <- c(ratPolyA="rat", "amin"="Amin et al.", "bartel.hek"="HEK", "bartel.hela"="HeLa")[as.character(m$dataset)]
  m <- as.matrix(m[,-1:-2])
  if(doSort) m <- m[,order(-matrixStats::colMedians(m, na.rm=TRUE)-sqrt(colMeans(m, na.rm=TRUE)))]
  column_split <- NULL
  cll <- colnames(m)
  if(!is.null(doSplit)){
    if(is.logical(e[[doSplit]])){
      column_split <- gsub("_TRUE|_FALSE","",colnames(m))
    }else{
      column_split <- gsub(paste(paste0("_",unique(e[[doSplit]])),collapse="|"),"",colnames(m))
    }
    cll <- gsub(".+_","",colnames(m))
  }
  Heatmap(m, column_split=column_split, col=cols, row_split=ds,
          show_row_names=show_row_names, show_row_dend=show_row_dend, heatmap_legend_param=hlm, 
          name=value.var, ..., row_title_rot = 0, cluster_columns=cluster_columns, 
          row_labels=gsub("vNeg","",gsub(".+:","",row.names(m))), row_names_gp=gpar(fontsize=9),
          column_labels=cll, column_title_gp=gpar(fontsize=11),
          column_names_gp=gpar(fontsize=11))
}


renameDatasets <- function(e, removeUnpublished=TRUE){
  ds <- c(ratPolyA="rat", "amin"="Amin et al.", "bartel.hek"="HEK", "bartel.hela"="HeLa")
  stopifnot(all(unique(e$dataset) %in% names(ds)))
  e$dataset <- ds[as.character(e$dataset)]
  if(removeUnpublished) e <- e[e$dataset!="rat",]
  e
}
