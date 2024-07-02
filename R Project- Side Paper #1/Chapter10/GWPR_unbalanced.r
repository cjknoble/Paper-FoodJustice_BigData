# This is an attempt for using unbalanced panel estimation technique for GWPR calibration
# The general idea goes like this: if we weigh the panel data using geographic effect (spatial effect) across all time periods simultaneously, it is likely that there is going to be overfitting problem espeically when there are longer periods of time. Even if the temporal period is short, we still find when applying the spatially varying coefficients of spatial filtering approach to panel data, that the results are almost the same as if there are no spatial effects - because the individual effects somehow absorbed the spatial effects (practically the same, anyways for geographic units as individuals)
# If that's the case, then why shall we weigh all the observations simultaneously? It probably will be better if we weigh each year's observations separately


# There is one very important note: the coords parameter MUST be created when the units are sorted alphabetically ALPHABETICALLY!
# This is because when the data was transformed through pdata.frame, the arrangement will be automatically alphabetical!
# This rule shall apply to all following functions.

# This version is now workable - Need to add more to make sure that it works with both balanced and unbalanced panles.
# This is called version 1, finished on December, 7, 2019, at 2:14 am
# Author: Danlin Yu
#
##########################################################################################################################################################################################################

gwpr.unbalanced <- function(formula, data = list(), coords, bandwidth = vector(mode = "numeric", length = nyears), nearestnb = vector(mode = "numeric", length = nyears), gweight=gwr.adaptive, 
		effect = c("individual","time","twoways"),
    model = c("within","random","ht","between","pooling","fd"),
    random.method = c("swar","walhus","amemiya","nerlove"),
    inst.method = c("bvk","baltagi"),
	index = NULL, criterion = "AIC") #the name and year index have to be given in the argument list, makes things easier
    														#name and year shall be two column vector that follows exactly the order as how the
    														#data is organized in the datalist, that is, stacked over time.
    														#plm requires the data to be stacked over cross-sectional units, that can be done
    														#via pdata.frame() function
{
	this.call <- match.call()
	# Check to see if the data is a pdata.frame. If not, convert it to be one.
	# Check to see if the data is a pdata.frame. If not, convert it to one:
	# The coordinates parameter, coords, must be obtained when the spatial units are ordered alphabetically!!! This is very important!!!
	require(plm)
	if (class(data)[1] != "pdata.frame") {
		# Convert the data to plm recognizable format to extract individual years:
		anal.data <- pdata.frame(data, index = index)		
	}
	else {
		anal.data <- data
	}	

	# get the index values - number of units and number of time periods.
	
	nunits <- dim(unique(index(anal.data)[1]))[1]
	nyears <- dim(unique(index(anal.data)[2]))[1]
	

	# Now, after the data is converted to pdata.frame, per the data structure, the observations are now organized as (each column):
	# Unit(1-year1), Unit(1-year2), Unit(1-year3) ...
	# Such arrangement makes the arrangement of the annual spatial weights a bit awkward.
	# The annual spatial weights is produced in an nunits * nunits square matrix, there will be nyears such square matrices.
	# Now, a wide format spatial weight matrix must be produced to weigh both the dependent and independent variables.
	# The wide format, however, is not simply concatenate the nyears spatial weight matrices together, instead
	# the format must follow (for each row) Unit(1-year1), Unit(1-year2), Unit(1-year3) ... format so that the later muliplying these weights to the observations
	# makes sense. To arrange the wide format spatial weight matrix in such way, we will need to use matrix index trick when obtaining the weights in later segement.
	# Basically, a wide format weight matrix is defined first with nunits rows and nunits * nyears columns.
	# Second, when assigning weights to this wide weight matrix, the column index will be first stored in a vector that starts from the first spaial units in the each year 
	# for (i in 1:nyears) and ends in the last spatial units of the wide matrix (nunits * nyears), but with an interval of nyears 
	# (so it can be expressed as vc <- seq(i, nunits * nyears, by = nyears))
	# Third, this vc index will be used to assign the obtained weights.
		
	mt <- terms(formula, data = anal.data)
	mf <- lm(formula, data = anal.data, method="model.frame")
	lm <- lm(formula, data = anal.data, model=FALSE)
	if (missing(coords))
		stop("Observation coordinates have to be given")
	if (missing(bandwidth) && missing (nearestnb)) 
	stop("Either a bandwidth or a nearest neighbor value has to be given")
	if ((unique(bandwidth < 0)) || (unique(nearestnb < 0)))
	stop("Invalid bandwidth or nearest neighbor")
	xdiff <- diff(range(coords[,1]))
	ydiff <- diff(range(coords[,2]))
	if (min(bandwidth) > min(xdiff,ydiff)) stop("Invalid bandwidth")
	require(stats)
	y <- model.response(mf, "numeric")
	x <- as.data.frame(model.matrix(mt, mf)[,-1]) #I don't need the intercept to be in the x, so an index of -1 would do the trick, as.data.frame keeps x to be a data frame
	colnames(x) <- names(mf)[-1] # the name for the dependent variable is removed.
	n <- nunits #we can't use NROW(x) because x will contain repetitive rows as the panel is stacked over the temproal dimension
	m <- NCOL(x) #+ 1 # the intercept column needs to be added back for the gwpr.b/gwpr.pt/etc. - for some reason, plm does not report intercept anymore, so the 1 does not need to be added back.
	dist2 <- (as.matrix(dist(coords)))^2
	w <- matrix(data = 0, nrow = nunits, ncol = nunits * nyears)
	# To get either the fixed or adaptive kernal
	# Here, however, the weights shall be obtained for each year instead for the entire panel
	# This section needs to be re-written: instead of construct an n by n weight matrix
	# The weight matrix must be a matrix of dimension nt by n (nt rows, n columns)
	# This way, when weighting the x and y, I don't need the rep() function. I will just use the w to weight x and y.
	
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		for (i in 1:nyears) {
			vc <- seq(i, nyears * nunits, by = nyears)
			w[,vc] <- gweight(dist2, bandwidth[i])
			}
	}
	else if (identical(gweight, gwr.adaptive))
	{
		if (min(nearestnb) <= m) stop ("Invalid neareast neighbor")
		for (i in 1:nyears){
			vc <- seq(i, nyears * nunits, by = nyears)
			w[,vc] <- gweight(dist2,nearestnb[i])
			# bandwidth <- gweight(dist2,nearestnb)$bw adaptive bandwidth is not needed since it changes for every location.
		}
	}
	

	# Because the name and year indexes are given explicitly in the argument list, we can combine these two to the extended weight matrix:
	# This and the entire reshaping operation are not necessary anymore - I simply changed the assignment of weights from "to rows" to "to columns"
	# I must be too dizzy to realize this is the easier solution with less complilcation
	# w <- cbind(index(anal.data), w)
	# The weight matrix must be pretreated here to make sure that it can be used to weigh the observations:
	# The idea is that the first line and all the following "first lines" of the above stacked w weight matrix need to become a line of the new weight matrix, which has a dimension of n by nt
	# where w is of dimension of nt by n. The new matrix is not a transpose of w, since w is not symatric.
	
	# The reshape function from the base stat package do this nicely. The trick is that the name and year parameter must be specified!
	# After reshape, it will produce a data.frame object, this needes further treatment. The first column (name column) must be removed
	# Then the dataframe must be converted to a matrix (so that w.wide[i, ] makes sense):
	
	# w.wide <- reshape(w, idvar = names(index(anal.data))[1], timevar = names(index(anal.data))[2], direction = "wide")
	# w.wide <- w.wide[, -1] # Remove the first column
	# w.wide <- as.matrix(w.wide)
	
	# The containers holding the critical values.
	
	gwpr.b <- as.data.frame(matrix(nrow=n, ncol=m))
	gwpr.pv<- as.data.frame(matrix(nrow=n, ncol=m))
	gwpr.psuedot<- as.data.frame(matrix(nrow=n, ncol=m))
	gwpr.r.squared <- as.data.frame(matrix(nrow = n, ncol = 2))
	colnames(gwpr.r.squared) <- c("rsq", "adjrsq")
	gwpr.pooltest<-numeric(n) #the pooltest for each location's coefficients
	
	# For each spatial units, I will hold the values for the dependent variable and independent variables that participate in the local unbalanced panel regression:
	# This information can be directly extracted from the plm object's model element.
	
	indvmd <- list()
		
	# hold the demeaned residuals:
	resd <- double(nunits*nyears)
	# hold the demeaned values
	#data.y.x<-pdata.frame(as.data.frame(cbind(name,year,y,x)),index=c(colnames(name),colnames(year)))
	#data.y.x.noindex<-data.y.x[,-1][,-1]
	#data.ybar.xbar<-matrix(nrow=n,ncol=m+1)
	#	for (i in 1:n)
	#		{
	#		o = (i-1)*nyears + 1
	#		d = i*nyears
	#		data.ybar.xbar[i,] <- colMeans(data.y.x.noindex[o:d,])
	#		}

	#gwpr.R2 <- numeric(n) #not much use for now
	#gwpr.resid <- numeric(n) #not much use for now
	#yiybar <- (y - mean(y)) #not much use for now
	
	colnames(gwpr.b) <- colnames(as.data.frame(model.matrix(mt,mf)))[-1] #remove the "Intercept", same below.
	colnames(gwpr.pv) <- colnames(as.data.frame(model.matrix(mt,mf)))[-1]
	colnames(gwpr.psuedot) <- colnames(as.data.frame(model.matrix(mt,mf)))[-1]
	
	#lhat <- matrix(nrow=n, ncol=n)
	#require(plm)
	
	require(dplyr)	# The dplyr package contains the filter function that can filter out the zeros.
		#The x and y need to be weighted separately for each year to construct the dependent and independent variables
		#This is because the weights are applied to each year simultaneously

		#Let's assume that the years will be specified in the parameter list
		#The plm requires data to be stacked over cross-sectional units, the pdata.frame function is able to create a plm compatible file format
		
		#Each row of the weight matrix can be inflated via repeate the row for the number of years
		#This will make a inflated row matrix and it can be used to multiply the panel data to get a geographically weighted panel data
		#And this geographically weighted panel data will particiapte for panel regression using plm.
		#The coefficients generated from the plm will be the varying coefficients of that particular geographic location
		#Repeating a specific row of the weight matrix can be down by function rep(w[i,],nyears)
		
		# dep<-numeric(length(y))
		# colnames(dep)<-colnames(y)
		# indp<-matrix(nrow=NROW(x),ncol=NCOL(x))
		# colnames(indp)<-colnames(x)
		
	for (i in 1:n) {
		# Because the weight matrix is not symmetric, I cannot use the following weight assignment expressions, neither can I use the transpose of the weight matrix to do that:
		# Instead, the weight matrix must be pre-created before these two expressions can happen:
		
		dep <- w[i,] * y
		indp <- w[i,] * x

		# Get the name from the data. This is a bit tricky:
		# The names are stored in the index of the anal.data, the first argument. It takes a few steps:
		
		# First, extract the unique values in the first argument of the index as characters:
		
		tmpname <- lapply(unique(index(anal.data)[1]), as.character)
		
		# This produces a list (tmpname is a list), it needs to be unlisted and converted to vector:
		
		tmpname <- as.vector (unlist(tmpname))
		
		# Assign the first name to the first element of tmpname:
		
		name.i <- tmpname[i]
		
		fm <- formula
		# Now the data:
		# First, transfer the data to plm compatible data

		tempdata <-pdata.frame(as.data.frame(cbind(index(anal.data), dep, indp)), index= index)
				
		# then get rid of all the zero components, this is necessary for pooltest later on
		# 
		plmdata <- filter(tempdata, tempdata[,-(1:2)][1] !=0)
		# colnames(plmdata[,-(1:2)])<-colnames(mf)
		colnames(plmdata)[which(names(plmdata) == "dep")] <- colnames(mf)[1] # Changing the column names is a real pain here.
		#tmp.i <- try(plm(fm, plmdata, effect = effect, model = model, random.method = random.method, inst.method = inst.method,
		#			index= names(index(anal.data))), silent = TRUE) # in case there is an error
		plm.i <- plm(fm, plmdata, effect = effect, model = model, random.method = random.method, inst.method = inst.method,
					index= names(index(anal.data)))			
		#if (class(tmp.i) != "try-error") { 
		#	plm.i <- tmp.i 
			indvmd[[i]] <- plm.i$model		
			gwpr.b[i,] <- coefficients(plm.i)
			rownames(gwpr.b)[i] <- as.vector(name.i)
			gwpr.psuedot[i,] <- summary(plm.i)$coefficients[,3]
			rownames(gwpr.psuedot)[i] <- as.vector(name.i)
			gwpr.pv[i,] <- summary(plm.i)$coefficients[,4]
			rownames(gwpr.pv)[i] <- as.vector(name.i)
			gwpr.r.squared[i,] <- summary(plm.i)$r.squared
			rownames(gwpr.r.squared)[i] <- as.vector(name.i)
			
		for (j in seq(from = 1, to=dim(plmdata)[1])) {
				if (attr(plm.i$residuals, "index")[j, 1] == name.i) {
					resd[((i-1)*nyears+1):(i*nyears)] <- plm.i$residuals[j:(j+nyears-1)]
					break
					}
				}			
			
		#	}# using only the name as the index to estimate unbalanced panel
		#else {
		#	plm.i = NULL
		#	indvmd[[i]] <- NULL		
		#	#gwpr.b[i,] <- NULL
		#	rownames(gwpr.b)[i] <- as.vector(name.i)
		#	#gwpr.psuedot[i,] <- NULL
		#	rownames(gwpr.psuedot)[i] <- as.vector(name.i)
		#	#gwpr.pv[i,] <- NULL
		#	rownames(gwpr.pv)[i] <- as.vector(name.i)
		#	#gwpr.r.squared[i,] <- NULL
		#	rownames(gwpr.r.squared)[i] <- as.vector(name.i)			
		#	}
			
					
		
		# store the individual location's participating dependent and independent variables:
		


		# resd[i] <- data.ybar.xbar[i,1] - (t(gwpr.b[i,]) %*% data.ybar.xbar[i, 2:ncol(data.ybar.xbar)])
		# ybar.i <- t(b.i)%*%x[((i-1)*nyears+1):(i*nyears),]
		# resd[((i-1)*nyears+1):(i*nyears)] <- y[((i-1)*nyears+1):(i*nyears)] - ybar.i
		# resd[((i-1)*nyears+1):(i*nyears)] <- resd[((i-1)*nyears+1):(i*nyears)] - mean(resd[((i-1)*nyears+1):(i*nyears)])
		# The following loop finds the geographic units that has the same name as the ith geographic units that the current local calibration is centered on, using the attr(plm.i$residuals, "index")[j, 1] == name.i to make the judgement call:
		# The following codes works fine for the current location because there will always be the nyears of observations of the current location for the local calibration:
		#
		################# TODO
		# Need to take care of the NULL plm.i here:
		

		}
		
		# sp <- summary(plm.i)		
		u.hat <- resd # extract residuals
		df <- cbind(as.vector(u.hat), attr(u.hat, "index"))
		names(df) <- c("resid", "Country", "Time")
		n.N = length(u.hat) # number of data
		s.sq  <- log(sum(u.hat^2)/(n.N)) # local sigmas
		if (effect == "individual"){
		n = nunits
		np = m+n+1 # update number of parameters
		}

		if (effect == "time"){
		T = nyears
		np = m+T+1 # update number of parameters
		}

		if (effect == "twoways"){
		n = nunits
		T = nyears
		np = m+n+T # update number of parameters
		}
		aic <- round(       2*np  +  n.N * (  log(2*pi) + s.sq  + 1 ),1)
		bic <- round(log(n.N)*np  +  n.N * (  log(2*pi) + s.sq  + 1 ),1)
	
		if (criterion == "AIC") { gwpr.ic = aic}
			else { gwpr.ic = bic }
		


		
		#
		# If the number of time is less than the number of parameters, then the pooltest won't work
		# So a check is needed:

		# No pooltest is needed for now.
#		if (nyears > ncol(x)){
#		gwpr.pooltest[i] <- pooltest(fm,plmdata,model="within")$p.value
#		}
		#R2 is hard to define yet, the following will be commented out
		#gwpr.R2[i] <- 1 - ((t(ei) %*% diag(w[i,]) %*% ei) / 
		#    (t(yiybar) %*% diag(w[i,]) %*% yiybar))
		#I don't think the hat matrix is well defined in the GWPR yet, so the following section will be commented out temporarily
		#Z <- matrix(0, m, m) 
		#Z[upper.tri(Z, diag=TRUE)] <- lm.i$qr$qr[upper.tri(lm.i$qr$qr, diag=TRUE)]
		#Z[lower.tri(Z)] <- lm.i$qr$qr[upper.tri(lm.i$qr$qr)]

		#In case the inversion cannot be completed, that row of the
		#hat matrix will be assigned zeros
		#inv.Z <- try(chol2inv(Z))
		#if (class(inv.Z) != "try-error"){
		#	lhat[i,] <- x[i,] %*% inv.Z %*% t(x) %*% diag(w[i,])
		#}
		#else {
		#	lhat[i,] <- 0
		#}
			
	#This section calculates the local standard errors, localse, for
	#local psuedo-t scores, localpt.
	#Reference: The GWR book page 54-56, formula 2.14
	#However, in panel scenario, these statistics need further thoughts before they can be calcuated.
	#So they are temporarily commented out
	
	#localse <- matrix(0,nrow = nrow(coords), ncol = ncol(x))
	#localpt <- matrix(0,nrow = nrow(coords), ncol = ncol(x))
	#for (i in 1:n) {
	#	CC <- (solve(t(x)%*%diag(w[i,])%*%x))%*%t(x)%*%diag(w[i,])
	#	locvar <- as.numeric(sigma2)*(CC%*%t(CC))
	#	localse[i,] <- sqrt(diag(locvar))
	#	localpt[i,] <- gwpr.b[i,]/localse[i,]
	#}
	#colnames(localpt) <- colnames (x)
	#colnames(localse) <- colnames (x)

	z <- list(this.call=this.call, lm=lm, x=x, y=y, bandwidth=bandwidth, w = w, 
	gwpr.b=gwpr.b, localpt = gwpr.psuedot, localpv=gwpr.pv, gwpr.pooltest=gwpr.pooltest, name = unique(index(anal.data)[1]), timeperiod = unique(index(anal.data)[2]), gwpr.ic = gwpr.ic, gwpr.residual = resd, gwpr.r.sqr = gwpr.r.squared, indvmd = indvmd)
	class(z) <- "gwpr"
	invisible(z)
}

#The following function is for bandwidth calculation
#bisquare and gauss:

gwr.bisquare <- function(dist2, d) {
	d2 <- d^2
	w <- ifelse(dist2 > d2, 0, (1 - (dist2/d2))^2)
	w
}

gwr.gauss <- function(dist2, bandwidth) {
	w <- exp((-dist2)/(bandwidth^2))
	w
}

#This one calculates the varying weights
#The calculated w is no longer symmetric

gwr.adaptive <- function (dist2,N)
{
	n <- nrow(dist2)
	d2 <- numeric(n)
	w <- matrix(0, nrow = n, ncol = n)
	for (i in 1:n)
	{
		d2[i] <- sort(dist2[i,])[N]
		w[i,] <- ifelse(dist2[i,] > d2[i], 0, (1-(dist2[i,]/d2[i]))^2)
	}
	bw <- sqrt(d2)
	result <- list(w = w, bw = bw)
#	invisible(result)
	w <- result$w
	w
}


# The bandwidth selection uses the same gwr.cv /gwr.cv.f and gwr.aic / gwr.aic.f functions as in the cross-sectional gwr
# This is because:
# If we use the same weight matrix over the years, the assumption is that the spatial process that generates and influences the data will remain unchanged over the years.
# This is a rather strong assumption and likely not going to hold for the panel setting.
# In addition, it is possible that if we assume the spatial process does not change over time, the spatial effect might be absorbed by the individual effect for the panel model.
# To prevent the individual effect to be inseparable from the spatial effects, we separate the spatial effects out for each cross-section (each year), estimate the optimal bandwidth
# or nearest neighbor for each year, and use that individual optimal bandwiths or nearest neighbors for each year to weigh the yearly data (for corresponding years).
# At each location, the weighting will very likely produce unbalanced panel, but it can be estimated by only specify the individual index in the plm package (as has been done above):



gwpr.bdwt <- function(formula, data = list(), coords, gweight = gwr.gauss, index = NULL, fullsearch = FALSE, out = FALSE, method = "cv") {
	# There is one very important note: the coords parameter MUST be created when the units are sorted alphabetically ALPHABETICALLY!
	# This is because when the data was transformed through pdata.frame, the arrangement will be automatically alphabetical!
	# Check to see if the data is already a panel data frame (pdata.frame):
	# Here the coords argument must be ordered alphabetically!!!
	require(plm)
	if (class(data)[1] != "pdata.frame") {
		# Convert the data to plm recognizable format to extract individual years:
		anal.data <- pdata.frame(data, index = index)
	}
	else {
		anal.data <- data
	}
	# We need filter() from tidyverse package for the following data extraction:
	require(tidyverse)
	#data.for.anal <- list()
	timelen <- dim(unique(index(anal.data)[2]))[1]
	bdwt <- vector(mode = "numeric", length = timelen) # index function will produce the indexes of the panel data, the second column must be the year (temporal) column.
	# Using as.name() can remove the quotation marks from the colnames (names) functions' results:
	# To extract the Year index as numeric values is a bit tricky:
	# The second column of the index(anal.data) is going to be a "pindex" and "data frame" object
	# Each row of that object is a factor
	# To convert the factor to numeric value, the as.numeric(level(factor)) function has to be used:
	# timeindex <- as.numeric(levels((unique(anal.data)[2])[1,]))
	# I do not need the timeindex anymore, it seems the as.numeric will return the levels to generic orders starting from 1.
	# What I really need to know is the timelen, this will be sufficient.
	for (i in 1:timelen) {
		data.for.anal <- filter(anal.data, as.numeric(index(anal.data)[2][,1]) == i)
		if (method == "cv") {
			bdwt[i] <- gwr.cv(formula = formula, data = data.for.anal, coords = coords, gweight = gweight, fullsearch = fullsearch, out = out)
		}
		else {
			bdwt[i] <- gwr.aic(formula = formula, data = data.for.anal, coords = coords, gweight = gweight, fullsearch = fullsearch, out = out)
		}
	}

	bdwt
	invisible(bdwt)
}

gwr.cv <- function(formula, data = list(), coords, gweight=gwr.gauss, fullsearch = FALSE, out=FALSE)
#argument out controls the printout of the gwr.cv.f's score.
{
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, method="model.frame")
	if (missing(coords))
		stop("Observation coordinates have to be given")
	require(stats)
	dist2 <- (as.matrix(dist(coords)))^2
	y <- model.response(mf, "numeric")
	x <- model.matrix(mt, mf)
	if (NROW(x) != NROW(dist2))
		stop("Input data and coordinates have different dimensions")

	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		xdiff <- diff(range(coords[,1]))
		ydiff <- diff(range(coords[,2]))
		difmin <- min(xdiff, ydiff)
		beta1 <- difmin/1000
		beta2 <- difmin
		opt <- optimize(gwr.cv.f, lower=beta1, upper=beta2, maximum=F,
			y=y, x=x, dist2=dist2, gweight=gweight, fullsearch=fullsearch, out=out)
		bdwt <- opt$minimum
		invisible(bdwt)
	}
	else if (identical(gweight, gwr.adaptive))
	{
		beta1 <- ncol(x) + 1
		beta2 <- length(y) - 1
		if (!fullsearch)
		{
			opt <- optimize(gwr.cv.f, lower=beta1, upper=beta2, maximum=F,
				y=y, x=x, dist2=dist2, gweight=gweight, fullsearch=fullsearch, out=out)
			kthnb <- opt$minimum
			invisible(kthnb)
		}
		else
		{
			#to store the full search scores
			fullsc <- numeric (beta2 - beta1 + 1)
			#use a matrix to store the full search scores and corresponding nearest neighbors
			fullscmt <- matrix (nrow = beta2 - beta1 + 1, ncol = 2)
			for (i in beta1:beta2)
			{
				fullsc[i-beta1+1] <- gwr.cv.f(i,y=y,x=x,dist2=dist2,gweight=gweight,
					fullsearch=fullsearch, out=out)
				fullscmt[i-beta1+1,] <- cbind(i,fullsc[i-beta1+1])
			}
			kthnb <- which.min(fullsc) + beta1 - 1
			#If using fullsearch option, I would like the full search sorce be recorded as well
			z <- list (kthnb = kthnb, fullscmt = fullscmt)
			invisible(z)
		}
	}
}

gwr.cv.f <- function(BorN, y, x, dist2, gweight, fullsearch=FALSE, out=FALSE)
{
	n <- NROW(x)
	m <- NCOL(x)
	cv <- numeric(n)
	
	#Determine the weights according to weighting schemes
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		w <- gweight(dist2, BorN)
	}
	else if (identical(gweight, gwr.adaptive))
	{
		w <- gweight(dist2, BorN)
	}


# Now we have the weight matrix w, it is time to double weight the weight matrix:
		
	w <- apply (w, 2, function(col) col / sum(col))		
				
	for (i in 1:n) {
		xx <- x[i,]
		w.i <- w[i,]
		w.i[i] <- 0
		lm.i <- lm.wfit(y=y, x=x, w=w.i)
		b <- coefficients(lm.i)
		cv[i] <- y[i] - (t(b) %*% xx)
	}
	score <- sqrt(sum(t(cv) %*% cv) / n)
	
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		cat("Bandwidth:", BorN, "CV score:", score, "\n")
	}
	else if (identical (gweight, gwr.adaptive))
	{
		if (!fullsearch && out) {
			cat("Nearest neighbor:", BorN, "CV score:", score, "\n")
		}
	}
	score
}

#These two functions are going to select the best bandwidth or nearest neighbor
#through minimize the AIC score

gwr.aic <- function (formula, data = list(), coords, gweight=gwr.gauss, fullsearch = FALSE, out=FALSE)
{
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, method="model.frame")
	if (missing(coords))
		stop("Observation coordinates have to be given")
	require(stats)
	dist2 <- (as.matrix(dist(coords)))^2
	y <- model.response(mf, "numeric")
	x <- model.matrix(mt, mf)
	if (NROW(x) != NROW(dist2))
		stop("Input data and coordinates have different dimensions")
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		xdiff <- diff(range(coords[,1]))
		ydiff <- diff(range(coords[,2]))
		difmin <- min(xdiff, ydiff)
		beta1 <- difmin/1000
		beta2 <- difmin
		opt <- optimize(gwr.aic.f, lower=beta1, upper=beta2, maximum=F,
			y=y, x=x, dist2=dist2, gweight=gweight, fullsearch=fullsearch, out=out)
		bdwt <- opt$minimum
		invisible(bdwt)
	}
	else if (identical(gweight, gwr.adaptive))
	{
		#beta1 <- ncol(x) + 1
		#beta2 <- length(y) - 1
		#According to Prof. Brunsdon's codes, a search between ncol(x)+1
		#to length(y)-1 is not necessary - as been proved in practice as well
		#the suggested search range in his codes starts from 30/7 and stops
		#at 10*30 = 300 (this might be very suspicious, especially for very large
		#dataset, stops at 300 may not be a good idea, I need more reference for
		#selecting a working and flexible range. For testing purpose, here I will
		#use this range.
		#beta1 <- 30/7
		#beta2 <- 300
		#I want to try something else
		beta1 <- max(10, length(y)/20)
		beta2 <- length(y) - 1
		
		if (!fullsearch) {
			opt <- optimize(gwr.aic.f, lower=beta1, upper=beta2, maximum=F,
				y=y, x=x, dist2=dist2, gweight=gweight, fullsearch=fullsearch, out=out)
			kthnb <- opt$minimum
			invisible(kthnb)
		}
		else 
		{
			#to store the full search scores
			fullsc <- numeric (beta2 - beta1 + 1)
			#use a matrix to store the full search scores and corresponding nearest neighbors
			fullscmt <- matrix (nrow = beta2 - beta1 + 1, ncol = 2)
			for (i in beta1:beta2) {
				fullsc[i-beta1+1] <- gwr.aic.f(i,y=y,x=x,dist2=dist2, gweight=gweight, 
					fullsearch=fullsearch, out=out)
				fullscmt[i-beta1+1,] <- cbind(i,fullsc[i-beta1+1])
			}
			kthnb <- which.min(fullsc) + beta1 - 1
			#If using fullsearch option, I would like the full search sorce be recorded as well
			z <- list (kthnb = kthnb, fullscmt = fullscmt)
			invisible(z)
		}		
	}
}

gwr.aic.f <- function (BorN, y, x, dist2, gweight, fullsearch = FALSE, out=FALSE)
{
	n <- NROW(x)
	m <- NCOL(x)
	AIC <- numeric(n)
	iden <- diag(1, nrow = n, ncol = n)
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare))
	{
		w <- gweight(dist2, BorN)
	}
	else if (identical(gweight, gwr.adaptive))
	{
		w <- gweight(dist2, BorN)
	}

	# Now we have the weight matrix w, it is time to double weight the weight matrix:
		
	w <- apply (w, 2, function(col) col / sum(col))

	#Get the hat matrix and its trace to calculating AIC score
	
		
	lhat <- matrix(nrow=n, ncol=n)	
	for (i in 1:n) {
		lm.i <- lm.wfit(y=y, x=x, w=w[i,])
		Z <- matrix(0, m, m)
		Z[upper.tri(Z, diag=TRUE)] <- lm.i$qr$qr[upper.tri(lm.i$qr$qr,
			diag=TRUE)]
		Z[lower.tri(Z)] <- lm.i$qr$qr[upper.tri(lm.i$qr$qr)]

		inv.Z <- try(chol2inv(Z), silent = TRUE)
		if (class(inv.Z)[1] != "try-error"){
			lhat[i,] <- x[i,] %*% inv.Z %*% t(x) %*% diag(w[i,])
		}
		else {
			lhat[i,] <- 0
		}
		#inv.Z <- chol2inv(Z)
		#lhat[i,] <- x[i,] %*% inv.Z %*% t(x) %*% diag(w[i,])
	}
	
	hattrace <- sum(diag(lhat))
	
	#Get the sigma2 for AIC calculation
	
	rss <- t(y)%*%t((iden-lhat))%*%(iden-lhat)%*%y
	sigma2 <- rss/n
	
	AIC <- n*log(sigma2) + n*log(2*3.14) + n * (n + hattrace)/(n - 2 - hattrace)
	
	if (identical(gweight,gwr.gauss) || identical(gweight, gwr.bisquare)){
		cat("Bandwidth:", BorN, "Akaike Information Criterion:", AIC, "\n")
	}
	else {
		if (!fullsearch && out) {
			cat("Nearest neighbor:", BorN, "Akaike Information Criterion:", AIC, "\n")
		}
	}
	AIC
}





# To-do:
# Bisquare doesn't work -  Error in uniqval[as.character(effect), , drop = F] : 
#  incorrect number of dimensions 
#
# Need further investigation
#
gwr.bisquare <- function(dist2, d) {
	d2 <- d^2
	w <- ifelse(dist2 > d2, 0, (1 - (dist2/d2))^2)
	w
}