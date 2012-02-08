#' Plot individually estimated parameters
#' 
#' Description of \code{parplot}
#' 
#' 
#' @aliases parplot
#' @param object DESCRIPTION
#' @param type DESCRIPTION
#' @param group DESCRIPTION
#' @param breaks DESCRIPTION
#' @param ... DESCRIPTION
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords OLS, growth
#' @export parplot
#' @examples 
#' 
#' \dontrun{
#' data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
#' mod <- OLScurve(~ time, data = data)	
#' parplot(mod)
#' 
#' }
parplot <- function(object, ...){
	UseMethod('parplot')
}

#' @S3method parplot OLScurve
#' @rdname parplot 
#' @method parplot OLScurve 
parplot.OLScurve <- function(object, type = 'hist', group = NULL, 
	breaks = NULL, ...)
{
	pars <- object$pars
	npars <- ncol(pars)
	f <- deparse(object$formula[[3]])
	tmp <- unlist(strsplit(f,"\\+"))
	Names <- c("(intercept-", tmp)
	if(any(grep("I\\(",Names))){
		Names <- gsub("I\\(","",Names)		
		Names <- gsub("\\)\\)","-",Names)
		Names <- gsub("\\)","",Names)		
	}
	Names <- gsub(" ","",Names)
	Names <- gsub("-",")",Names)	
	pars2 <- data.frame(pars)
	colnames(pars2) <- paste('X',1:npars,sep='')		
	forms <- paste("~ X",1:npars,sep='')	
	if(type == 'splom'){
		pars <- data.frame(pars)
		colnames(pars) <- Names
		if(is.null(group)) print(splom(~pars, data = pars, main = 'Growth Parameters'))
		else {
			pars$group <- as.factor(na.omit(data.frame(object$orgdata,group))$group)
			splom(~pars|group, data = pars, main = 'Growth Parameters')		
		}
	}
	devAskNewPage(ask=TRUE)
	if(type == 'hist'){
		for(i in 1:npars){      
			form <- as.formula(forms[i])
			print(histogram(form,pars2,xlab = Names[i],breaks = breaks,main = 'Parameter Distributions'))			
		}
	}
	if(type == 'boxplot'){
		for(i in 1:npars){
			form <- as.formula(forms[i])
			print(bwplot(form,pars2,xlab = Names[i],main = 'Parameter Distributions'))
		}
	}
}
