#' Run enrichment analysis
#' This function calculate results of enrichment.
#'
#' @param membership membership data.frame with columns \code{names} and
#'                   \code{membership}
#' @param annotation annotation data.frame with columns \code{names},
#'                   \code{termID} and \code{termName} for each object-term
#'                   pair.
#' @param printTwoSided print two-sided p-value
#' @param usePrintAlt print also alternative side, i.e. depletion
#' @param usePrintID print annotation ID, else annotation description
#'
#' @return
#' @export
#'
#' @examples
run_enrichment<-function(membership,annotation,gname,printTwoSided=TRUE,
                         usePrintAlt=TRUE,
                         usePrintID=TRUE){
    ## check input data.structure
    if(!is.data.frame(membership)){
        stop('Membership have to be a data.frame.\n')
    }
    if(any(is.na(match(c('names','membership'),names(membership))))){
        stop("Membership should contains columns 'names' and 'membership'.\n")
    }
    if(!is.data.frame(annotation)){
        stop('Annotation have to be a data.frame.\n')
    }
    if(any(is.na(match(c('names','termID','termName'),names(annotation))))){
        stop("Membership should contains columns 'names', ",
             "'termID' and 'termName'.\n")
    }
    anno<-annotation[,c('termID','termName','names')]
    ## load clustering and annotation data
    rEnrich::load(x=membership,anno=anno)

    ## run enrichment analysis on loaded data
    rEnrich::run()

    ## get enrichment values
    res = rEnrich::getResults(printTwoSided=as.integer(printTwoSided),
                              usePrintAlt=as.integer(usePrintAlt),
                              usePrintID=as.integer(usePrintID))
    df<-makeDF(res,unique(anno[,1]),gname,usePrintAlt)
    return(df)

}

#' Convert result matrix to data.frame.
#'
#' Creates data.frame with the following columns:
#' * Fn annotation term ID
#' * C cluster/group ID
#' * Cn size of the cluster
#' * N size of the universe
#' * Mu actual number of term found in cluster
#' * OR
#' * CIl lower CI boundary
#' * CIu upper CI boundary
#' * Pv p.value
#' * Ap adjusted p.value
#' * PvALT p-value for the alternative side, i.e. depletion
#' if \code{\link{usePrintAlt}} is set to \code{TRUE}.
#' * ApALT adjusted p.value for the alternative side, i.e. depletion
#' if \code{\link{usePrintAlt}} is set to \code{TRUE}.
#' * Fe
#' * E(Mu) expected number of term found in cluster/group by chance
#' * F total number of term instances in the universe
#' * Alg name of grouping algorithm/approach
#' * Fc
#'
#' @param tt result matrix from \code{\link{getResults}}
#' @param terms list of unique annotation terms
#' @param gname name of grouping
#' @param usePrintAlt
#'
#' @return data.frame with the following columns:
#'         "Fn","C","Cn","N","Mu","OR","CIl","CIu","Pv","Ap","PvALT","ApALT",
#'         "Fe","E(Mu)","F","Alg","Fc"
makeDF<-function(tt,terms,gname,usePrintAlt){
    indx <- match(terms,colnames(tt))
    fn   <- colnames(tt)[indx]
    FN   <- length(fn)
    CN   <- length(tt[,1])
    DF <- matrix("",nrow=(FN*CN), ncol=17)
    colnames(DF) <- c("Fn","C","Cn","N","Mu","OR","CIl","CIu","Pv","Ap",
                      "PvALT","ApALT","Fe","E(Mu)","F","Alg","Fc")
    DF[,1]  = rep(fn,CN) # "Fn"
    DF[,4]  = rep(N,(FN*CN)) # "N"
    DF[,15] = rep(as.vector(unlist(tt[1,indx])),CN) # "F"
    DF[,16] = rep(gname,(FN*CN)) # "Alg"
    temp1  <- c()
    temp2  <- c()
    temp3  <- c()
    temp4  <- c()
    temp5  <- c()
    temp6  <- c()
    temp7  <- c()
    temp8  <- c()
    temp9  <- c()
    temp10 <- c()
    temp11 <- c()
    temp12 <- c()
    for( i in 1:length(tt[,1]) ){
        cc    <- rep(tt[i,1],FN)
        temp1 <- c(temp1, cc)
        cn    <- rep(tt[i,2],FN)
        temp2 <- c(temp2, cn)
        ov = as.vector(unlist(tt[i,grepl("actual",colnames(tt))]))
        temp3 <- c(temp3, ov)
        ep = as.vector(unlist(tt[i,grepl("expected",colnames(tt))]))
        temp12 <- c(temp12, ep)
        or = as.vector(unlist(tt[i,grepl("OR",colnames(tt))]))
        temp4 <- c(temp4, or)
        ci = as.vector(unlist(tt[i,grepl("CI",colnames(tt))]))
        ci = gsub("\\[","",ci)
        ci = gsub("\\]","",ci)
        temp5 <- c(temp5,as.vector(sapply(ci, getCI, indx=1)))
        temp6 <- c(temp6,as.vector(sapply(ci, getCI, indx=2)))
        pv = as.vector(unlist(tt[i,grepl("p.value", colnames(tt)) &
                                     !grepl("p_value", colnames(tt)) &
                                     !grepl("X.p_value",colnames(tt)) &
                                     !grepl("ALT",colnames(tt))]))
        temp7 <- c(temp7, pv)
        pa = as.vector(unlist(tt[i,grepl("adjusted",colnames(tt)) &
                                     !grepl("ALT",colnames(tt))]))
        temp8 <- c(temp8, pa)
        palt = as.vector(unlist(tt[i,grepl("p.value.ALT.",colnames(tt))]))
        temp9 <- c(temp9, palt)
        alt = as.vector(unlist(tt[i,grepl("adjusted.ALT.",colnames(tt))]))
        temp10 <- c(temp10, alt)
        fc = (as.numeric(ov)/as.numeric(cn)) / (as.numeric(cn)/as.numeric(N))
        temp11 <- c(temp11, fc)
        if( max(fc) > FcMAX ){
            FcMAX = max(fc)
        }
        #rm(cc,cn,ov,ep,or,ci,pv,pa,palt,alt,fc)
    }
    DF[,2]  = temp1 # "C"
    DF[,3]  = temp2 # "Cn"
    DF[,5]  = temp3 # "Mu"
    DF[,6]  = temp4 # "OR"
    DF[,7]  = temp5 # "CIl"
    DF[,8]  = temp6 # "CIu"
    DF[,9]  = temp7 # "Pv"
    DF[,10] = temp8 # "Ap"
    DF[,13] = temp11 # "Fe"
    DF[,14] = temp12 # "E(Mu)"
    DF[,17] = (as.numeric(DF[,5])/as.numeric(DF[,15])) /
        (as.numeric(DF[,3])/as.numeric(DF[,4]) ) # "Fc"
    if( max(as.numeric(DF[,17])) > FeMAX ){
        FeMAX = max(as.numeric(DF[,17]))
    }
    if(usePrintAlt){
        DF[,11] = temp9 # "PvALT"
        DF[,12] = temp10 # "ApALT"
    }else{
        DF<-DF[,-c(11,12)]
    }
    rm(cc,cn,ov,or,ci,pv,pa,palt,alt,fc,ep)
    rm(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12)
    DF = as.data.frame(DF)
    stridx<-match(c("Fn",'Alg'),names(df))
    DF[,-stridx]<-lapply(DF[,-stridx], function(.x){as.numeric(.x)})
    DF[,stridx]<-lapply(DF[,stridx], function(.x){factor(.x)})
    return(DF)
}
