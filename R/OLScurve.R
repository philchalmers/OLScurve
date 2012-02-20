#' Ordinary least squares growth curve trajectories
#' 
#' The \code{OLScurve} provides a simple way of specifying ordinary least squares (OLS)
#' growth curve models in R. Individual OLS trajectories are fit to each case and an 
#' OLScurve object is returned which can be passed to several graphical and summary
#' function within the package. 
#' 
#' As Bollen and Curran (2006) note, there are a variety of advantages 
#' to using the case-by-case approach for estimating trajectory parameters. First of all, OLS 
#' estimation is intuitively appealing, making it a good pedagogical tool for introducing how to model
#' trajectories, and illuminates many essential conditions and assumptions necessary for LCMs.
#' Second, prediction of the parameters for individual trajectory estimates are calculated for
#' each case in the sample, which can lead to several diagnostics by statistical and graphical
#' means. Also, summary statistics can be computed for these estimates (which can also be
#' graphically portrayed) and if need be these estimates can be analyzed further by other
#' statistical frameworks.
#'
#' Unfortunately there are also several limitation to OLS estimation for LCMs, namely:
#' overall tests of fit are not readily available, the structure of the error variances must be
#' unrealistically constrained to estimate a pooled standard error, the latent factors cannot be
#' regressed without error on other exogenous or time-varying variables, and analytic significance 
#' tests are often not readily available (Bollen & Curran, 2006). However, OLS estimation
#' may still be useful in the preliminary stages of latent curve modeling for (a) selecting appropriate 
#' functional forms of growth, (b) examining unconditional population homogeneity, (c)
#' observing whether the relationship between growth factors are linear, and for (d) detecting
#' influential outliers (Carrig et al., 2004).
#' 
#' 
#' @aliases OLScurve
#' @param formula a \code{formula} specifying how the functional form of \code{time} should be 
#'    coded. The keyword \code{time} must be used, and additive components such as powers, square 
#'    roots, and exponentials can be included
#' @param data a data frame in the wide (one subject per row) format containting only the time
#'    related variables. Can be of class \code{matrix} or \code{data.frame} 
#' @param time the relative spacing between time points. Default is equal spacing
#' @param x an \code{OLScurve} object
#' @param group a \code{factor} grouping variable used to parition the results
#' @param SE logical; print a list containing the standard errors?
#' @param digits number of digits to round
#' @param sep logical; should the plots be seperated?
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords OLS, growth
#' @export OLScurve
#' @seealso \code{\link{parplot}}, \code{\link{subjplot}}
#' @references 
#'
#' Bollen, K. A. & Curran, P. J. (2006). \emph{Latent Curve Models: A Structural Equation Perspective}. John Wiley & Sons.
#'
#' Carrig, M. M., Wirth, R. J., & Curran, P. J. (2004). A SAS Macro for Estimating and Visualizing Individual Growth Curves. \emph{Structural Equation Modeling, 11}, 132-149.
#'
#'
#' @examples 
#' 
#' \dontrun{
#' ##linear
#' data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
#' mod1 <- OLScurve(~ time, data = data)	
#' mod1
#' plot(mod1)
#'
#' ##quadratic
#' data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2)  
#' mod2 <- OLScurve(~ time + I(time^2), data = data)    
#' mod2
#' plot(mod2)
#' 
#' ##sqrt
#' data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(5:1))   
#' mod3 <- OLScurve(~ sqrt(time), data = data)    
#' mod3
#' plot(mod3)
#' 
#' ##exponential
#' data <- t(t(matrix(rnorm(1000,0,5),200)) + exp(0:4))  
#' mod4 <- OLScurve(~ sqrt(time), data = data)    
#' mod4
#' plot(mod4)
#'
#' ##combination
#' data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(1:5))  
#' mod5 <- OLScurve(~ time + sqrt(time), data = data)    
#' mod5
#' plot(mod5)
#'
#' ##two group analysis with linear trajectories
#' data1 <- t(t(matrix(rnorm(500),100)) + 1:5) 
#' data2 <- t(t(matrix(rnorm(500),100)) + 9:5)
#' data <- rbind(data1,data2) 
#' group <- c(rep('male',100),rep('female',100)) 
#'
#' mod <- OLScurve(~ time, data)
#' print(mod,group)
#' plot(mod,group)
#' }
OLScurve <- function(formula, data, time = 0:(ncol(data)-1), ...){
	call <- match.call()
	ch <- as.character(formula)
	if(length(ch) == 2) ch <- c(ch[1], "y", ch[2])
		else ch[2] <- "y"
	formula <- paste(ch[2],ch[1],ch[3])
	N <- nrow(data)
	J <- ncol(data)
	formula <- as.formula(formula)
	ind <- 1:N
	rownames(data) <- as.character(ind)
	orgdata <- data
	data <- na.omit(as.matrix(data))
	y <- data[1,]
	npars <- length(lm(formula)$coef)	
	pars <- matrix(0,nrow(data),npars)
	res <- pred <- data	
	for(i in 1:nrow(data)){
		y <- as.numeric(data[i,])
		mod <- lm(formula)
		pars[i,] <- mod$coef
		pred[i,] <- predict(mod)
		res[i, ] <- residuals(mod)
	}	
	mod <- list(pars = pars, pred = pred, res=res, data = data.frame(data), orgdata = orgdata, 
		formula = formula, omitted = nrow(data)/nrow(orgdata), time = time, call = call)
	class(mod) <- "OLScurve"
	mod
}

#' @S3method print OLScurve
#' @rdname OLScurve 
#' @method print OLScurve 
print.OLScurve <- function(x, group = NULL, SE = TRUE, digits = 3, ...){
	data <- x$data
	orgdata <- x$orgdata	
	J <- ncol(data)
	N <- nrow(data)	
	pars <- x$pars
	npars <- ncol(pars)
	lowerind <- matrix(FALSE,npars,npars)
	time <- x$time	
	varEi <- rowSums((x$res^2))/(length(time) - 2)
	varE <- mean(varEi)	
	for(i in 1:npars)
		for(j in 1:npars)
			if(i < j) lowerind[j,i] <- TRUE
	if(!is.null(group)){
		group <- as.vector(group)
		dat <- cbind(group,x$orgdata)
		group <- na.omit(dat)[,1]
		u <- unique(group)
		loops <- length(u)		
	} else {
		loops <- 1
		u <- 'fulldata'
		group <- rep('fulldata',N)
	}
	f <- deparse(x$formula[[3]])
	tmp <- unlist(strsplit(f,"\\+"))
	Names <- c("(intercept-", tmp)
	if(any(grep("I\\(",Names))){
		Names <- gsub("I\\(","",Names)		
		Names <- gsub("\\)\\)","-",Names)
		Names <- gsub("\\)","",Names)		
	}
	Names <- gsub(" ","",Names)
	Names <- gsub("-",")",Names)	
	Meanslist <- Covlist <- SElist <- list()	
	for(i in 1:loops){
		Means <- colMeans(pars[group == u[i],])
		Ntmp <- nrow(pars[group == u[i],])
		MeanSEs	<- sqrt(colSums(((pars[group == u[i],] - Means)^2)/(Ntmp-1))/Ntmp)
		varEtmp <- mean(varEi[group == u[i]])
		Cov <- cov(pars[group == u[i],])		
		Cor <- cor(pars[group == u[i],])
		Cov[lowerind] <- Cor[lowerind]
		names(MeanSEs) <- names(Means) <- colnames(Cov) <- rownames(Cov) <- Names				
		Meanslist[[i]] <- round(Means,digits)
		Covlist[[i]] <- round(Cov, digits)		
		SElist[[i]] <- round(MeanSEs, digits)
	}
	names(Meanslist) <- names(Covlist) <- names(SElist) <- u
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
		"\n", sep = "")
	cat("Note: ", (1-x$omitted)*100,
		"% ommited cases. \n", sep='')	
	cat("Pooled standard error = ",varE,"\n\n")	
	cat("MEANS:\n\n")	
	print(Meanslist)	
	if(SE){		
		cat("Standard Errors for Means:\n\n")	
		print(SElist)	
	}
	cat("\nCOVARIANCE (correlations on lower off-diagonal):\n\n")
	print(Covlist)
	invisible(list(Meanslist,Covlist,SElist))
}

#' @S3method plot OLScurve
#' @rdname OLScurve
#' @method plot OLScurve 
plot.OLScurve <- function(x, group = NULL, sep = FALSE, ...){
	data <- x$data
	N <- nrow(data)
	fn <- fn1 <- x$formula
	id.o <- id <- as.numeric(rownames(data))
	data <- data.frame(data)		
	
	plotOLScurve <- function(data, fn, group = NULL, layout = NULL) {		
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
			fn2 <- as.formula(fn1)
			fm <- lm(fn2)
			panel.lines(x,predict(fm)) 
		}
		groupPlots <-xyplot( plot.fn , datalg, groups = factor(id), 
			layout = layout,
			xlab = "Time",
			ylab = "Predicted",
			main = 'Group Plots',
			panel = panel.superpose,
			panel.groups = mypanel)
		print(groupPlots)				
	}
	plotOLScurve(data, fn, group, layout = NULL)
}
