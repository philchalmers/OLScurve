#' Plot distribution of parameters
#' 
#' A plotting function for displaying the distribution of the OLS parameter
#' estimates. 
#' 
#' 
#' @aliases parplot
#' @param object an object of class \code{OLScurve}
#' @param type type of plot to display; can be \code{'hist'}, \code{'boxplot'}, or \code{'splom'}
#'    for a histogram, boxplot, or scatter plot matrix
#' @param group a \code{factor} grouping variable used to parition the results
#' @param breaks number of breaks to be used in plotting the histogram
#' @param ... additional arguments to be passed
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
