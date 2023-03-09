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
run_enrichment<-function(membership,annotation,printTwoSided=1,
                         usePrintAlt=1,
                         usePrintID=1){
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
    res = rEnrich::getResults(printTwoSided=1,
                              usePrintAlt=1,
                              usePrintID=1)
    return(res)

}
