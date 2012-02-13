#' Plot individually estimated parameters
#' 
#' A plotting function for displaying the individuals trajectories and their 
#' modelled functional form. Useful for detecting aberrant individual trajectories.
#' 
#' 
#' @aliases subjplot
#' @param object an object of class \code{OLScurve}
#' @param group a \code{factor} grouping variable used to parition the results
#' @param layout a variable to be passed to \code{xyplot} to adjust the graphical layout
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords OLS, growth
#' @export subjplot
#' @examples 
#' 
#' \dontrun{
#' data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
#' mod <- OLScurve(~ time, data = data)	
#' subjplot(mod)
#' 
#' 
#' }
subjplot <- function(object, ...){
	UseMethod('subjplot')
}

#' @S3method subjplot OLScurve
#' @rdname subjplot 
#' @method subjplot OLScurve 
subjplot.OLScurve <- function(object, group = NULL, layout = NULL, ...)
{
	data <- object$data
	N <- nrow(data)
	fn <- fn1 <- object$formula
	id.o <- id <- as.numeric(rownames(data))
	data <- data.frame(data)	
	if(is.null(layout)) layout <- c(ceiling(log(N)),ceiling(log(N)))
	
	plotOLScurve <- function(data, fn, group = NULL, layout = NULL) {
		devAskNewPage(ask=TRUE)
		data <- data.frame(data) 
		if(is.null(data$id)) data$id.o<-1:nrow(data) 
			else data$id.o <- data$id 
		if(!is.null(group)) 
			data$group <- group
		if(is.null(data$group)) { 
			data <- data[order(data$id.o),] 
			plot.fn <- y ~ time
		} else { 
			data <- data[order(data$group,data$id.o),]
			data$id.o <- paste(data$group,data$id.o) 
			plot.fn <- y ~ time | group
		}
		ys <- colnames(data)[!(colnames(data)=="id")&
			!(colnames(data)=="group")&
			!(colnames(data)=="id.o")]
		datalg <- reshape(data, idvar="id",
					varying = list(ys),
					v.names = c("y"), 
					times = c(1:length(ys)),
					direction="long")
		ch <- as.character(fn)
		ch[2] <- gsub("x","y",ch[2],fixed=TRUE)
		ch[3] <- gsub("time","x",ch[3],fixed=TRUE)
		fn1 <- paste(ch[2],ch[1],ch[3])

		###### PLOT INDIVIDUAL PARTICIPANTS ######
		mypanel = function(x, y, ...){
			panel.xyplot(x, y, ...)
			fn2 <- as.formula(fn1)
			fm <- lm(fn2)
			panel.lines(x,predict(fm))
		}
		  
		subjectPlots <- xyplot( y ~ time | factor(id), datalg, groups = id, layout = layout,
			xlab = "time",
			ylab = "",
			main = 'Subject plots',
			panel = mypanel)
		  
		print(subjectPlots)	  				
	}
	plotOLScurve(data, fn, group, layout)
}
