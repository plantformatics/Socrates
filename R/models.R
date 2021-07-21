###################################################################################################
###################################################################################################
###################################################################################################
#' Accessibility residuals from unconstrained quasibinomial logistic regression
#'
#' Run unregularized logistic regression. Output is unconstrained Pearson's residuals
#' for each cell and ACR.
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix colSums
#' @importFrom CMF matrix_to_triplets
#' @import doSNOW
#' @import iterators
#' @import itertools
#'
#' @param obj, list containing a slot named 'counts' for the dgCMatrix object with binary cell
#' (columns) x peak (rows) accessibility scores, and a slot named 'meta' containing a data.frame
#' with barcode meta data information.
#' @param r.variates model formula specifying variable to regress out. Should be specified as
#' a character string. Defaults to "~log10nSites". Changing this parameter is not recommended.
#' @param variates model formula as character string to regress out variables from the initial
#' residuals using QR residuals in a second round of regression. Ideal for specifying
#' batch effects, cell-cycle and other technical sources of variation. Defaults to NULL.
#' @param nthreads numeric, number of threads to run logistic model in parallel. Depends on
#' doSNOW. Defaults to 1.
#' @param chunks numeric, number of peaks to run for each daughter process. Defaults to 256.
#' @param type character, specify the residual type to extract from the model. Defaults to
#' "pearson".
#' @param link character, specify the link-function type for logistic regression. Possible
#' parameters are "logit" and "probit". Defaults to "logit".
#' @param center.resid logical, whether to zero-center residuals. Defaults to TRUE.
#' @param scale.resid logical, whther to standardize residuals. Defaults to FALSE.
#' @param make.sparse logical, whether or not to set negative values to 0, and reduce memory usage.
#' Setting this parameter to TRUE has negligible effects on downstream results. Setting make.sparse
#' to TRUE will override center.resid and scale.resid arguments, setting them to FALSE. Defaults to
#' FALSE.
#' @param verbose logical. Defaults to FALSE.
#' @param slotName character, specify the slot name for saving residuals. Useful for saving
#' multiple normalization steps. Note, make sure to update the slotName argument for
#' downstream functions. Defaults to "residuals".
#'
#' @rdname logisticModel
#' @export
#'
logisticModel <- function(obj,
                          r.variates="~log10nSites",
                          variates=NULL,
                          nthreads=1,
                          chunks=256,
                          type="pearson",
                          link="logit",
                          center.resid=TRUE,
                          scale.resid=FALSE,
                          make.sparse=FALSE,
                          verbose=FALSE,
                          slotName="residuals"){

    # verbose
    if(verbose){message(" - running logistic regression ...")}

    # set-up data
    x <- obj$counts
    y <- obj$meta
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    yo <- y
    form <- as.formula(paste0("z",r.variates))
    if(verbose){message("   * formula: ", form)}

    # use method
    do.func <- function(form, df){

        # run GLM
        suppressWarnings(mod <- glm(form,
                                    data=df,
                                    family=quasibinomial(link = link)))

        # return residuals
        return(residuals.glm(mod, type=type))
    }

    # parallel
    if(nthreads == 1){

        # transpose
        x <- t(x)

        # begin iterate
        its <- 0
        res <- lapply(seq(1:ncol(x)), function(i, y, form){

            # verbose about progress
            its <<- its + 1
            if(its %% 1000 == 0){if(verbose){message("   * estimated residual devaince for ", its, " peaks ...")}}

            # create DF
            df <- cbind(y, x[,i])
            colnames(df) <- c(colnames(y), "z")

            # run GLM
            do.func(form, df)

        }, y=y, form=form)

        # reform res
        res <- as.matrix(do.call(rbind, res))
        rownames(res) <- colnames(x)
        colnames(res) <- rownames(x)

    }else{

        # control looping/chunks
        peaknames <- rownames(x)

        # run logistic regression for each peak independently - in parallel
        res <- mclapply(seq(1:nrow(x)), function(j){

            # create DF
            df <- cbind(y, x[j,])
            colnames(df) <- c(colnames(y), "z")

            # run GLM
            do.func(form, df)

        }, mc.cores=nthreads)

        # return results
        res <- as.matrix(do.call(rbind, res))
        colnames(res) <- colnames(x)
        rownames(res) <- rownames(x)

    }


    ###########################################################################################
    # variables to regress from residuals -----------------------------------------------------
    ###########################################################################################
    if(!is.null(variates)){

        # verbose
        if(verbose){message(" - removing confounding variables ...")}

        # model
        form.nr <- as.formula(paste0("z",variates))
        if(verbose){message("   * formula: ", form.nr)}

        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)

        # control looping/chunks
        chunks <- 1000
        peaknames <- rownames(x)
        niter <- nrow(x)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }

        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it)
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }

        # calc first QR
        y1 <- yo[colnames(res),]
        regression.mat <- cbind(y1, res[1,])
        colnames(regression.mat) <- c(colnames(y1),"z")
        qr <- lm(form.nr, data = regression.mat, qr = TRUE)$qr
        rm(regression.mat)
        r.names <- rownames(res)

        # run logistic regression for each peak independently - in parallel
        res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts) %dopar% {

            # run chunks
            its <- c(n$i:(n$i+n$m-1))
            dfout <- lapply(its, function(j){

                # run GLM
                qr.resid(qr = qr, y = res[j,])

            })

            # reformat
            dfout <- do.call(rbind, dfout)
            rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
            return(dfout)

        }

        # close connections
        close(pb)
        stopCluster(cl)

        # clean
        res <- res2
        rm(res2)
        res <- res[r.names,]

    }

    # return
    obj[[slotName]] <- res
    obj$norm_method <- "logisticModel"
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Regularized quasibinomial logistic regression
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix colSums
#' @import doSNOW
#' @import glmnet
#' @import iterators
#' @import itertools
#'
#' @param obj, list containing a slot named 'counts' for the dgCMatrix object with binary cell
#' (columns) x peak (rows) accessibility scores, and a slot named 'meta' containing a data.frame
#' with barcode meta data information.
#' @param r.variates model formula specifying variable to regress out. Should be specified as
#' a character string. Defaults to "~log10nSites". Changing this parameter is not recommended.
#' @param variates model formula as character string to regress out variables from the initial
#' residuals using QR residuals in a second round of regression. Ideal for specifying
#' batch effects, cell-cycle and other technical sources of variation. Defaults to NULL.
#' @param subpeaks numeric, number of ACRs to select for regularization parameter estimates.
#' Defaults to 5000.
#' @param bins numeric, number of bins to split ACRs into. Defaults to 256.
#' @param bw_adjust numeric, sets the bandwidth for kernel regression during parameter
#' regularization. Defaults to 10.
#' @param type character, specify the residual type to extract from the model. Possible choices
#' inlucde: "pearson", "deviance", "response", and "working". Defaults to "pearson".
#' @param link character, specify the link-function type for logistic regression. Possible
#' parameters are "logit" and "probit". Defaults to "logit".
#' @param nthreads numeric, number of threads to run logistic model in parallel. Depends on
#' doSNOW. Defaults to 1.
#' @param weights sets weights to glmnet regression to remove technical effects.
#' @param method character, type of secondary regression to use for removing additional technical
#' effects. Choices include "elasticNet" for glmnet based regularized regression or "lm" for
#' linear regression. Option has no effect when variates is NULL. Using glmnet instead of "lm" is
#' still experimental and time consuming. Final accessibility residuals should be carefully
#' analyzed if using "elasticNet". Defaults to "lm".
#' @param alpha numeric, sets the alpha parameter for glmnet. 0 for LASSO like regression, 1 for
#' ridge regression. Default 0.5 for elastic net regression.
#' @param center.resid logical, whether to zero-center residuals. Defaults to FALSE.
#' @param scale.resid logical, whther to standardize residuals. Defaults to FALSE.
#' @param make.sparse logical, whether or not to set negative values to 0, and reduce memory usage.
#' Setting this parameter to TRUE has negligible effects on downstream results. Setting make.sparse
#' to TRUE (default) will override center.resid and scale.resid arguments, setting them to FALSE. make.sparse
#' is mutally exclusive with the argument 'variates'. Defaults to TRUE
#' @param verbose logical. Defaults to FALSE.
#' @param slotName character, specify the slot name for saving residuals. Useful for saving
#' multiple normalization steps. Note, make sure to update the slotName argument for
#' downstream functions. Defaults to "residuals".
#'
#' @rdname regModel
#' @export
#'
regModel <- function(obj,
                     r.variates='~log10nSites',
                     variates=NULL,
                     subpeaks=5000,
                     bins=256,
                     bw_adjust=10,
                     type="pearson",
                     link="logit",
                     nthreads=1,
                     weights=NULL,
                     method="lm",
                     alpha=0.5,
                     center.resid=F,
                     scale.resid=F,
                     make.sparse=T,
                     verbose=FALSE,
                     slotName="residuals"){

    # fix type if necessary
    if(type!="deviance" & type !="pearson" & type!="response" & type !="working"){
        message(" - *type* incorrect, setting residual to default (pearson) ...")
        type <- "pearson"
    }

    # functions
    is_outlier          <- function(y, x, th=10) {
        bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
        eps <- .Machine$double.eps * 10
        breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
        breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
        score1 <- robust_scale_binned(y, x, breaks1)
        score2 <- robust_scale_binned(y, x, breaks2)
        return(pmin(abs(score1), abs(score2)) > th)
    }
    robust_scale_binned <- function(y, x, breaks) {
        bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
        tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
        score <- rep(0, length(x))
        o <- order(bins)
        if (inherits(x = tmp$x, what = 'list')) {
            score[o] <- unlist(tmp$x)
        } else {
            score[o] <- as.numeric(t(tmp$x))
        }
        return(score)
    }
    deviance_residual   <- function(y, mu, wt){
        d.res <- sqrt(pmax((binomial()$dev.resid)(y, mu, wt),0))
        d.res <- ifelse(y > mu, d.res, -d.res)
        return(d.res)
    }
    pearson_residual    <- function(y, mu, theta){(y-mu)/sqrt((mu*(1-mu))*theta)}
    robust_scale        <- function(x){return((x - median(x)) / (mad(x) + .Machine$double.eps))}
    ChunkPoints         <- function(dsize, csize) {
        return(vapply(
            X = 1L:ceiling(x = dsize / csize),
            FUN = function(i) {
                return(c(
                    start = (csize * (i - 1L)) + 1L,
                    end = min(csize * i, dsize)
                ))
            },
            FUN.VALUE = numeric(length = 2L)
        ))
    }
    invprobit           <- function(x){
        thresh <- -qnorm(.Machine$double.eps)
        x <- pmin(pmax(x, -thresh), thresh)
        pnorm(x)
    }

    # check if make.sparse
    if(make.sparse){
        variates <- NULL
    }

    ###############################################################################################
    # start ---------------------------------------------------------------------------------------
    ###############################################################################################
    if(verbose){message(" - regularizing logistic model parameters ...")}

    # set-up data
    y <- obj$meta
    x <- obj$counts
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    form <- as.formula(paste0("z",r.variates))
    if(verbose){message("   * formula: ", form)}

    #select subsample of peaks
    set.seed(1111)
    con <- 1
    if(verbose){message("   * estimating geometric mean ...")}
    log_peak_mean <- log(Matrix::rowMeans(x))
    if(verbose){message("   * density sampling on peak space ...")}
    log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
    sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean)$y + .Machine$double.eps)
    peak.sub <- sample(x=rownames(x), size=subpeaks, prob=sample_prob)
    x.sub <- x[peak.sub,]
    log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))

    # fit model to subset
    if(verbose){message("   * fitting parameters to ",nrow(x.sub)," peaks ...")}
    pars <- lapply(seq(1:nrow(x.sub)), function(j){
        df <- cbind(y, x.sub[j,])
        colnames(df) <- c(colnames(y), "z")
	    df$z <- as.factor(df$z)
        suppressWarnings(mod <- glm(form,
                                    data=df,
                                    family=quasibinomial(link = link)))
        theta <- summary(mod)$dispersion
        names(theta) <- "theta"
        return(c(theta, mod$coefficients))
    })

    # merge
    pars <- do.call(rbind, pars)
    rownames(pars) <- rownames(x.sub)


    ###############################################################################################
    ### regularize --------------------------------------------------------------------------------
    ###############################################################################################

    # outliers
    if(verbose){message("   * finding outliers ...")}
    peaks <- names(log_peak_mean)
    outliers <- apply(pars, 2, function(k) is_outlier(k, log_peak_mean_sub))
    outliers <- apply(outliers, 1, any)
    if (sum(outliers) > 0){
        if(verbose){message("   * found ", sum(outliers), " outliers ...")}
        pars <- pars[!outliers, ]
        peaks_step1 <- rownames(pars)
        log_peak_mean_sub <- log_peak_mean_sub[!outliers]
    }

    # select bw
    bw <- bw.SJ(log_peak_mean) * bw_adjust

    # parameter for predictions
    x_points <- pmax(log_peak_mean, min(log_peak_mean))
    x_points <- pmin(x_points, max(log_peak_mean))

    # take results from step 1 and fit/predict parameters to all genes
    if(verbose){message("   * regularizing all coefficients ...")}
    o <- order(x_points)
    pars_fit <- matrix(NA_real_, length(peaks), ncol(pars),
                       dimnames = list(peaks, colnames(pars)))

    # global fit / regularization for all coefficients
    for (i in 1:ncol(pars)){
        pars_fit[o, i] <- ksmooth(x = log_peak_mean_sub, y = pars[, i],
                                  x.points = x_points, bandwidth = bw, kernel='normal')$y
        #suppressWarnings(fit <- smooth.spline(log_peak_mean_sub, pars[,i], cv=T))
        #pars_fit[o, i] <- predict(fit, x=x_points[o])$y
    }


    ###############################################################################################
    # fit data ------------------------------------------------------------------------------------
    ###############################################################################################

    # affine transform
    if(verbose){message("   * estimating residuals on the full data set ...")}
    regressor_data <- model.matrix(as.formula(r.variates), data=y)

    # split peaks into bins
    bin_ind <- ceiling(x = 1:length(x = peaks) / bins)
    max_bin <- max(bin_ind)

    # prepare residual  matrix
    #res <- matrix(NA_real_, length(peaks), nrow(regressor_data),
    #              dimnames = list(peaks, rownames(regressor_data)))

    # iterate
    #if(verbose){
    #    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    #}
    res <- mclapply((1:max_bin), function(i){
        if(verbose){if((i %% 10 == 0)){message("   * finished collecting residuals for ", i, " of ", max_bin," peaks/bins ...")}}
        peaks_bin <- peaks[bin_ind == i]
        if(link=="logit"){
            mu <- exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)) /
                (1+exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }else if(link == "probit"){
            mu <- invprobit(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data))
        }else if(link =="cloglog"){
            mu <- 1-exp(-exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }
        y <- as.matrix(x[peaks_bin, , drop=FALSE])
        if(type=="deviance"){
            dat <- deviance_residual(y, mu, 1)
            #res[peaks_bin, ] <- deviance_residual(y, mu, 1)
        }else if(type=="pearson"){
            dat <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
            #res[peaks_bin,] <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }else if(type=="response"){
            dat <- (y - mu)
            #res[peaks_bin,] <- (y - mu)
        }else{
            dat <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
            #res[peaks_bin,] <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }
        spm <- matrix_to_triplets(dat)
        colnames(spm) <- c("i","j","x")
        spm <- as.data.frame(spm)
        spm$i <- peaks_bin[spm$i]
        spm$j <- rownames(regressor_data)[spm$j]
        spm <- subset(spm, spm$x > 0)
        return(spm)
    }, mc.cores=nthreads)
        # if(is.null(variates)){
        #     if(make.sparse==F){
        #         if(center.resid==T & scale.resid==F){
        #             res[peaks_bin,] <- t(apply(res[peaks_bin,], 1, function(pp){
        #                 pp - mean(pp, na.rm=T)
        #             }))
        #         }
        #         if(scale.resid==T){
        #             res[peaks_bin,] <- t(apply(res[peaks_bin,], 1, function(pp){
        #                 (pp - mean(pp, na.rm=T))/sd(pp, na.rm=T)
        #             }))
        #         }
        #     }
        # }

    #if(verbose){
    #    close(pb)
    #}

    # make sparse
    res <- do.call(rbind, res)
    res$i <- factor(res$i)
    res$j <- factor(res$j)
    res <- sparseMatrix(i=as.numeric(res$i),
                        j=as.numeric(res$j),
                        x=as.numeric(res$x),
                        dimnames=list(levels(res$i), levels(res$j)))
    res <- res[peaks, rownames(regressor_data)]
    res@x[res@x < 0] <- 0
    res <- drop0(res, tol=0)

    # remove na
    res[is.na(res)] <- 0
    res <- res[Matrix::rowSums(res == 0) != ncol(res),]

    # report intial residuals
    res.range <- range(res)
    if(verbose){message("   * residual range: ",res.range[1], " - ", res.range[2])}


    ###############################################################################################
    # variables to regress from residuals ---------------------------------------------------------
    ###############################################################################################

    if(!is.null(variates) && method=="lm"){

        # verbose
        if(verbose){message(" - removing effects from confounding variables ...")}

        # model
        form.nr <- as.formula(paste0("z",variates))
        if(verbose){message("   * formula: ", form.nr)}

        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)

        # control looping/chunks
        chunks <- 1000
        peaknames <- rownames(x)
        niter <- nrow(x)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }

        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it)
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }

        # calc first QR
        y1 <- obj$meta[colnames(res),]
        regression.mat <- cbind(y1, res[1,])
        colnames(regression.mat) <- c(colnames(y1),"z")
        qr <- lm(form.nr, data = regression.mat, qr = TRUE)$qr
        rm(regression.mat)
        r.names <- rownames(res)

        # run logistic regression for each peak independently - in parallel
        res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts) %dopar% {

            # run chunks
            its <- c(n$i:(n$i+n$m-1))
            dfout <- lapply(its, function(j){

                # run GLM
                val <- qr.resid(qr = qr, y = res[j,])
                return(val)

            })

            # reformat
            dfout <- do.call(rbind, dfout)
            rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]

            # center/scale
            if(center.resid==T & scale.resid==F){
                dfout <- t(apply(dfout, 1, function(pp){
                    pp - mean(pp, na.rm=T)
                }))
            }
            if(scale.resid==T){
                dfout <- t(apply(dfout, 1, function(pp){
                    (pp - mean(pp, na.rm=T))/sd(pp, na.rm=T)
                }))
            }

            # return
            return(dfout)

        }

        # close connections
        close(pb)
        stopCluster(cl)

        # clean
        res <- res2
        rm(res2)
        res <- res[r.names,]

    }else if(!is.null(variates) & method=="elasticNet"){

        # verbose
        if(verbose){message(" - removing effects from confounding variables with elasticNet ...")}

        # prep variates
        y1 <- obj$meta[colnames(res),]
        if(verbose){message("   * number of cells in meta data = ", nrow(y1))}
        ind.vars <- gsub("~","",variates)
        ind.vars <- gsub(" ","",ind.vars)
        ind.vars <- strsplit(ind.vars, "\\+")
        ind.vars <- do.call(c, ind.vars)
        n.vars <- c()
        f.vars <- c()
        for(i in ind.vars){
            if(is.numeric(y1[,i])){
                n.vars <- c(n.vars, i)
            }else if(is.factor(y1[,i])){
                f.vars <- c(f.vars, i)
            }else if(is.character(y1[,i])){
                f.vars <- c(f.vars, i)
                y1[,i] <- factor(y1[,i])
            }
        }
        sqr <- paste("I(",n.vars,"^0.5)",sep="")
        sqr <- paste(sqr, collapse=" + ")
        sq <- paste("I(",n.vars,"^2)",sep="")
        sq <- paste(sqr, collapse=" + ")
        if(length(n.vars) > 1){
            combs <- combn(ind.vars, 2)
            combs <- paste(combs[1,], combs[2,], sep=":")
            combs <- paste(combs, collapse=" + ")
            final <- as.formula(paste(variates, "+", combs, "+", sqr, "+", sq, sep=" "))
        }else if(is.null(n.vars) & length(f.vars) > 0){
            combs <- combn(ind.vars, 2)
            combs <- paste(combs[1,], combs[2,], sep=":")
            combs <- paste(combs, collapse=" + ")
            final <- as.formula(paste(variates,"+",combs,sep=" "))
        }else{
            final <- as.formula(variates)
            #final <- as.formula(paste(variates, "+", sqr, "+", sq, sep=" "))
        }

        # prep train/full
        if(length(f.vars) > 1){
            if(verbose){message("   * factor variables > 1")}
            xvars <- sparse.model.matrix(final, y1,
                                         contrasts.arg=lapply(y1[,colnames(y1) %in% f.vars], contrasts, contrasts=F))
        }else if(length(f.vars)==1){
            if(verbose){message("   * factor variables = 1")}
            f.list <- list(f.vars=contrasts(y1[,f.vars], contrasts=F))
            names(f.list) <- f.vars
            xvars <- sparse.model.matrix(final, y1, contrasts.arg=f.list)
        }else{
            if(verbose){message("   * factor variables = 0")}
            xvars <- sparse.model.matrix(final, y1)
        }

        # save peak ids
        r.names <- rownames(res)

        # get initial estimates of lambda
        if(verbose){message("   * getting regularized value of lambda from 1000 random peaks ...")}
        lambda <- c()
        for(i in 1:1000){
            rand <- sample(nrow(res), 1)
            yvar <- res[rand,]
            cv.lasso <- cv.glmnet(xvars, yvar, alpha=alpha, family="gaussian")
            lambda <- c(lambda, cv.lasso$lambda.min)
        }
        ave.lambda <- median(lambda[lambda > 0], na.rm=T)
        if(verbose){message("   * average lambda = ", ave.lambda)}

        # regularize glmnet function
        regularizedR  <- function(y, xvars, ave.lambda){

            # elastic-net
            fit <- glmnet(xvars, y, family="gaussian", alpha=alpha, lambda=ave.lambda)
            return(as.numeric(y-predict(fit, s=ave.lambda, newx=xvars, type="response")[,1]))
        }

        # parallel parameters
        cl <- makeSOCKcluster(nthreads)
        registerDoSNOW(cl)

        # control looping/chunks
        chunks <- bins
        peaknames <- rownames(res)
        niter <- nrow(res)
        tasks <- ceiling(niter/chunks)
        if(tasks < nthreads){
            message("Tasks (",tasks,") should not be less than number of CPUs (",nthreads,")")
        }

        # track vars
        pb <- txtProgressBar(max = tasks, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        package.list <- c("glmnet")

        # functions
        idivix <- function(n, chunkSize) {
            i <- 1
            it <- idiv(n, chunkSize=chunkSize)
            nextEl <- function() {
                m <- nextElem(it)
                value <- list(i=i, m=m)
                i <<- i + m
                value
            }
            obj <- list(nextElem=nextEl)
            class(obj) <- c('abstractiter', 'iter')
            obj
        }

        # run logistic regression for each peak independently - in parallel
        if(verbose){message("   * fitting elasticNet regressions and extracting residuals ...")}
        res2 <- foreach(n=idivix(niter, chunks), .combine='rbind', .inorder=F, .options.snow=opts,
                        .packages=package.list) %dopar% {

                            # run chunks
                            its <- c(n$i:(n$i+n$m-1))
                            dfout <- lapply(its, function(j){

                                # run GLMNET
                                return(regularizedR(res[j,], xvars, ave.lambda))

                            })

                            # reformat
                            dfout <- do.call(rbind, dfout)
                            rownames(dfout) <- rownames(res)[n$i:(n$i+n$m-1)]
                            colnames(dfout) <- colnames(res)

                            # center/scale
                            if(center.resid==T & scale.resid==F){
                                dfout <- t(apply(dfout, 1, function(pp){
                                    pp - mean(pp, na.rm=T)
                                }))
                            }
                            if(scale.resid==T){
                                dfout <- t(apply(dfout, 1, function(pp){
                                    (pp - mean(pp, na.rm=T))/sd(pp, na.rm=T)
                                }))
                            }

                            # return
                            return(dfout)

                        }

        # close connections
        close(pb)
        stopCluster(cl)

        # clean
        res <- res2
        rm(res2)
        res <- res[r.names,]

    }

    # attach residuals to object
    obj[[slotName]] <- res
    obj$norm_method <- "regModel"

    # return
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' Regularized quasibinomial logistic regression
#'
#' Requires a column titled 'library' in the meta slot of the input object.
#' @import Matrix
#' @import parallel
#'
#' @param obj, list containing a slot named 'counts' for the dgCMatrix object with binary cell
#' (columns) x peak (rows) accessibility scores, and a slot named 'meta' containing a data.frame
#' with barcode meta data information.
#' @param r.variates model formula specifying variable to regress out. Should be specified as
#' a character string. Defaults to "~log10nSites". Changing this parameter is not recommended.
#' @param subpeaks numeric, number of ACRs to select for regularization parameter estimates.
#' Defaults to 5000.
#' @param do.sample logical, whether or not to sub-sample cells for regularized regression.
#' Defaults to TRUE when the number of cells is greater 5000, FALSE otherwise. Set do.sample
#' to NULL to override subsampling.
#' @param subcells numeric, number of cells to select for regularization. Defaults to 1000.
#' @param propC numeric, proportion of cells to sample from each factor in meta data variable
#' 'library'. Requires 'library' is a column name in meta data with more than two factors.
#' Ultimately, the number of cells that are selected from each library is the greater of
#' propC * number of cells in library and subcells.
#' @param bins numeric, number of bins to split ACRs into. Defaults to 256.
#' @param bw_adjust numeric, sets the bandwidth for kernel regression during parameter
#' regularization.
#' @param type character, specify the residual type to extract from the model. Possible choices
#' inlucde: "pearson", "deviance", "response", and "working". Defaults to "pearson".
#' @param link character, specify the link-function type for logistic regression. Possible
#' parameters are "logit" and "probit". Defaults to "logit".
#' @param nthreads numeric, number of threads to run logistic model in parallel. Depends on
#' doSNOW. Defaults to 1.
#' @param center.resid logical, whether to zero-center residuals. Defaults to TRUE.
#' @param scale.resid logical, whther to standardize residuals. Defaults to FALSE.
#' @param make.sparse logical, whether or not to set negative values to 0, and reduce memory usage.
#' Setting this parameter to TRUE has negligible effects on downstream results. Setting make.sparse
#' to TRUE will override center.resid and scale.resid arguments, setting them to FALSE. Defaults to
#' FALSE.
#' @param verbose logical. Defaults to FALSE.
#' @param slotName character, specify the slot name for saving residuals. Useful for saving
#' multiple normalization steps. Note, if changing from default, make sure to provide
#' downstream functions with non-default slotName with useSlot. Defaults to "residuals".
#'
#' @rdname regModel2
#' @export
#'
regModel2 <- function(obj,
                      r.variates='~log10nSites',
                      do.sample=T,
                      subpeaks=5000,
                      subcells=1000,
                      propC=0.1,
                      cellNames=NULL,
                      bins=256,
                      bw_adjust=10,
                      type="pearson",
                      link="logit",
                      nthreads=1,
                      center.resid=T,
                      scale.resid=F,
                      make.sparse=F,
                      verbose=FALSE,
                      slotName="residuals"){

    # load necessary data
    x <- obj$counts
    y <- obj$meta

    # fix type if necessary
    if(type!="deviance" & type !="pearson" & type!="response" & type !="working"){
        message(" - *type* incorrect, setting residual to default (pearson) ...")
        type <- "pearson"
    }

    # hidden functions
    is_outlier          <- function(y, x, th = 10) {
        bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
        eps <- .Machine$double.eps * 10
        breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
        breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
        score1 <- robust_scale_binned(y, x, breaks1)
        score2 <- robust_scale_binned(y, x, breaks2)
        return(pmin(abs(score1), abs(score2)) > th)
    }
    robust_scale_binned <- function(y, x, breaks) {
        bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
        tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
        score <- rep(0, length(x))
        o <- order(bins)
        if (inherits(x = tmp$x, what = 'list')) {
            score[o] <- unlist(tmp$x)
        } else {
            score[o] <- as.numeric(t(tmp$x))
        }
        return(score)
    }
    deviance_residual   <- function(y, mu, wt){
        d.res <- sqrt(pmax((binomial()$dev.resid)(y, mu, wt),0))
        d.res <- ifelse(y > mu, d.res, -d.res)
        return(d.res)
    }
    pearson_residual    <- function(y, mu, theta){(y-mu)/sqrt((mu*(1-mu))*theta)}
    robust_scale        <- function(x){return((x - median(x)) / (mad(x) + .Machine$double.eps))}
    ChunkPoints         <- function(dsize, csize) {
        return(vapply(
            X = 1L:ceiling(x = dsize / csize),
            FUN = function(i) {
                return(c(
                    start = (csize * (i - 1L)) + 1L,
                    end = min(csize * i, dsize)
                ))
            },
            FUN.VALUE = numeric(length = 2L)
        ))
    }
    invprobit           <- function(x){
        thresh <- -qnorm(.Machine$double.eps)
        x <- pmin(pmax(x, -thresh), thresh)
        pnorm(x)
    }
    sampleCells         <- function(x, min.cells.per=1000, propLib=0.05){

        # number of cells
        x$library <- as.character(x$library)
        libs <- unique(x$library)
        out <- lapply(libs, function(z){

            # subset by library
            y <- x[x$library == z,]

            # num cells
            num.cells <- ceiling(nrow(y)*propLib)
            if(num.cells < min.cells.per){
                num.cells <- min.cells.per
            }

            # sample
            rownames(y)[sample(seq(1:nrow(y)), size=num.cells)]
        })
        out <- do.call(c, out)
        if(verbose){message("   * sampled a total of ",length(out)," cells ...")}
        return(out)

    }
    stddev              <- function(r.dev, bins=256, threads=1){

        # verbose
        if(verbose){message(" - standardizing deviance scores ...")}

        # set up bins
        bin_ind <- ceiling(x = 1:nrow(r.dev) / bins)
        max_bin <- max(bin_ind)
        ids <- rownames(r.dev)

        # run in bins
        dev <- lapply(seq(1:max_bin), function(x){
            peaks_bin <- rownames(r.dev)[bin_ind == x]
            t(as.matrix(scale(t(r.dev[peaks_bin,]))))
        })#, mc.cores=threads)
        rm(r.dev)

        # merge and return
        dev <- do.call(rbind, dev)
        dev <- dev[ids,]

        return(dev)
    }


    ###############################################################################################
    # start ---------------------------------------------------------------------------------------
    ###############################################################################################
    if(verbose){message(" - regularizing logistic model parameters ...")}

    # set-up data
    yo <- y
    x <- x[,colnames(x) %in% rownames(y)]
    y <- y[colnames(x),]
    form <- as.formula(paste0("z",r.variates))
    if(verbose){message("   * formula: ", form)}

    if(ncol(x) > 5000){
        do.sample <- T
    }
    if(is.null(do.sample)){
        do.sample <- F
    }

    # sample cells
    if(do.sample){
        if(verbose){message("   * sampling cells ...")}
        sub.cells <- sampleCells(y, min.cells.per=subcells, propLib=propC)
        x.sub1 <- x[,sub.cells]
        x.sub1 <- x.sub1[Matrix::rowSums(x.sub1)>0,]
        x.sub1 <- x.sub1[,Matrix::colSums(x.sub1)>0]
        y.sub1 <- y[colnames(x.sub1),]

        #select subsample of peaks
        set.seed(1111)
        con <- 1
        if(verbose){message("   * estimating geometric mean ...")}
        log_peak_mean <- log(Matrix::rowMeans(x))
        log_peak_mean_red <- log(Matrix::rowMeans(x.sub1))
        r.lpm <- range(log_peak_mean)
        log_peak_mean_red <- log_peak_mean_red[log_peak_mean_red > r.lpm[1] & log_peak_mean_red < r.lpm[2]]

        if(verbose){message("   * density sampling on peak space ...")}
        log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
        sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean_red)$y + .Machine$double.eps)
        peak.sub <- sample(x=names(log_peak_mean_red), size=subpeaks, prob=sample_prob)
        x.sub <- x.sub1[peak.sub,]
        log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))

        # fit model to subset
        if(verbose){message("   * fitting parameters to ",nrow(x.sub)," peaks ...")}
        pars <- lapply(seq(1:nrow(x.sub)), function(j){
            df <- cbind(y.sub1, x.sub[j,])
            colnames(df) <- c(colnames(y), "z")
            suppressWarnings(mod <- glm(form,
                                        data=df,
                                        family=quasibinomial(link = link)))
            theta <- summary(mod)$dispersion
            names(theta) <- "theta"
            return(c(theta, mod$coefficients))
        })

    }else{

        #select subsample of peaks
        set.seed(1111)
        con <- 1
        if(verbose){message("   * estimating geometric mean ...")}
        log_peak_mean <- log(Matrix::rowMeans(x))

        if(verbose){message("   * density sampling on peak space ...")}
        log_peak_dens <- density(x = log_peak_mean, bw = 'nrd', adjust = 1)
        sample_prob <- 1 / (approx(x=log_peak_dens$x, y=log_peak_dens$y, xout=log_peak_mean)$y + .Machine$double.eps)
        peak.sub <- sample(x=rownames(x), size=subpeaks, prob=sample_prob)
        x.sub <- x[peak.sub,]
        log_peak_mean_sub <- log(Matrix::rowMeans(x.sub))

        # fit model to subset
        if(verbose){message("   * fitting parameters to ",nrow(x.sub)," peaks ...")}
        pars <- lapply(seq(1:nrow(x.sub)), function(j){
            df <- cbind(y, x.sub[j,])
            colnames(df) <- c(colnames(y), "z")
            suppressWarnings(mod <- glm(form,
                                        data=df,
                                        family=quasibinomial(link = link)))
            theta <- summary(mod)$dispersion
            names(theta) <- "theta"
            return(c(theta, mod$coefficients))
        })

    }

    # merge
    pars <- do.call(rbind, pars)
    rownames(pars) <- rownames(x.sub)


    ###############################################################################################
    ### regularize --------------------------------------------------------------------------------
    ###############################################################################################

    # outliers
    if(verbose){message("   * finding outliers ...")}
    peaks <- names(log_peak_mean)
    outliers <- apply(pars, 2, function(k) is_outlier(k, log_peak_mean_sub))
    outliers <- apply(outliers, 1, any)
    if (sum(outliers) > 0){
        if(verbose){message("   * found ", sum(outliers), " outliers ...")}
        pars <- pars[!outliers, ]
        peaks_step1 <- rownames(pars)
        log_peak_mean_sub <- log_peak_mean_sub[!outliers]
    }

    # select bw
    bw <- bw.SJ(log_peak_mean) * bw_adjust

    # parameter for predictions
    x_points <- pmax(log_peak_mean, min(log_peak_mean))
    x_points <- pmin(x_points, max(log_peak_mean))

    # take results from step 1 and fit/predict parameters to all genes
    if(verbose){message("   * regularizing all coefficients ...")}
    o <- order(x_points)
    pars_fit <- matrix(NA_real_, length(peaks), ncol(pars),
                       dimnames = list(peaks, colnames(pars)))

    # global fit / regularization for all coefficients
    for (i in 1:ncol(pars)){
        pars_fit[o, i] <- ksmooth(x = log_peak_mean_sub, y = pars[, i],
                                  x.points = x_points, bandwidth = bw, kernel='normal')$y
    }


    ###############################################################################################
    # fit data ------------------------------------------------------------------------------------
    ###############################################################################################

    # affine transform
    if(verbose){message("   * estimating residuals on the full data set ...")}
    regressor_data <- model.matrix(as.formula(r.variates), data=y)

    # split peaks into bins
    bin_ind <- ceiling(x = 1:length(x = peaks) / bins)
    max_bin <- max(bin_ind)

    # iterate
    res <- mclapply(seq(1:max_bin), function(i){
        peaks_bin <- peaks[bin_ind == i]
        if(link=="logit"){
            mu <- exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)) /
                (1+exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }else if(link == "probit"){
            mu <- invprobit(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data))
        }else if(link =="cloglog"){
            mu <- 1-exp(-exp(tcrossprod(pars_fit[peaks_bin, -1, drop=FALSE], regressor_data)))
        }
        y <- as.matrix(x[peaks_bin, , drop=FALSE])
        if(type=="deviance"){
            out.1 <- deviance_residual(y, mu, 1)
        }else if(type=="pearson"){
            out.1 <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }else if(type=="response"){
            out.1 <- (y - mu)
        }else{
            out.1 <- pearson_residual(y, mu, pars_fit[peaks_bin,'theta'])
        }
        colnames(out.1) <- rownames(regressor_data)
        rownames(out.1) <- peaks_bin
        if(make.sparse==F){
            if(center.resid==T & scale.resid==F){
                out.1 <- t(apply(out.1, 1, function(pp){
                    pp - mean(pp, na.rm=T)
                }))
            }
            if(scale.resid==T){
                out.1 <- t(apply(out.1, 1, function(pp){
                    (pp - mean(pp, na.rm=T))/sd(pp, na.rm=T)
                }))
            }
        }
        return(out.1)
    }, mc.cores=nthreads)

    # merge
    res <- do.call(rbind, res)
    res <- res[peaks, ]

    # remove na
    res[is.na(res)] <- 0
    res <- res[rowSums(res)!=0,]

    # make sparse?
    if(make.sparse){
        res <- Matrix(res, sparse=T)
        res@x[res@x < 0] <- 0
        res <- drop0(res, tol=0)
    }

    # report intial residuals
    res.range <- range(res)
    if(verbose){message("   * residual range: ",res.range[1], " - ", res.range[2])}

    # return
    obj[[slotName]] <- res
    obj$norm_method <- "regModel2"
    return(obj)
}


###################################################################################################
###################################################################################################
###################################################################################################
#' TF-IDF normalization
#'
#' Run TF-IDF on binary cell x peak matrix. Returns normalized TF-IDF values in the residuals slot.
#' Code adapted from Andrew Hill (http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/).
#'
#' @import Matrix
#' @import DelayedArray
#'
#' @param obj list object containing dgCMatrix in 'counts' slot.
#' @param frequencies logical, whether to scale matrix by barcode peak sums. Defaults to TRUE.
#' @param log_scale_tf logical, whether to log1p transform the term frequency (TF).
#' Defaults to TRUE.
#' @param scale_factor numeric,
#' @param slotName character, specify the slot name for saving residuals. Useful for saving
#' multiple normalization steps. Note, make sure to update the slotName argument for
#' downstream functions. Defaults to "residuals".
#'
#' @rdname tfidf
#' @export
#'
tfidf <- function(obj,
                  frequencies=T,
                  log_scale_tf=T,
                  scale_factor=10000,
                  slotName="residuals"){

    # set bmat
    bmat <- obj$counts

    # hidden functions
    .safe_tfidf       <- function(tf, idf,  block_size=2000e6){
        result = tryCatch({
            result = tf * idf
            result
        }, error = function(e) {
            options(DelayedArray.block.size=block_size)
            DelayedArray:::set_verbose_block_processing(TRUE)

            tf = DelayedArray(tf)
            idf = as.matrix(idf)

            result = tf * idf
            result
        })
        return(result)
    }

    # Use either raw counts or divide by total counts in each cell
    if (frequencies) {
        # "term frequency" method
        tf = t(t(bmat) / Matrix::colSums(bmat))
    } else {
        # "raw count" method
        tf = bmat
    }

    # Either TF method can optionally be log scaled
    if (log_scale_tf) {
        if (frequencies) {
            tf@x = log1p(tf@x * scale_factor)
        } else {
            tf@x = log1p(tf@x * 1)
        }
    }

    # IDF w/ "inverse document frequency smooth" method
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))

    # TF-IDF
    tf_idf_counts = .safe_tfidf(tf, idf)
    rownames(tf_idf_counts) = rownames(bmat)
    colnames(tf_idf_counts) = colnames(bmat)
    obj[[slotName]] <- Matrix(tf_idf_counts, sparse=T)
    obj$norm_method <- "tfidf"

    # return
    return(obj)
}
