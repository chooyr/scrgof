# Modified version of secr::sim.detect
# Returns expected capture frequencies and capture histories simulated from 
# a fitted model.
sim.detect.new <- function (object, popnlist, maxperpoly = 100, renumber = FALSE) {

    dummycapthist <- function(capthist, pop, fillvalue = 1) {
        if (ms(capthist)) {
            output <- list()
            for (i in 1:nsession) output[[i]] <- dummycapthist(capthist[[i]], 
                pop = pop[i], fillvalue = fillvalue)
            class(output) <- c("capthist", "list")
            session(output) <- session(capthist)
            output
        } else {
            newdim <- dim(capthist)
            newdim[1] <- nrow(pop[[1]])
            output <- array(fillvalue, dim = newdim)
            output[, , newdim[3]] <- 1
            class(output) <- "capthist"
            traps(output) <- traps(capthist)
            session(output) <- session(capthist)
            covariates(output) <- covariates(pop[[1]])
            output
        }
    }
    if (!is.null(object$groups) & (object$details$param > 1)) 
        stop("simulation does not extend to groups when param>1")
    if ("telemetry" %in% unlist(detector(traps(object$capthist)))) 
        stop("telemetry models are not supported in simulate.secr")
    Markov <- any(c("B", "Bk", "K") %in% object$vars)
    btype <- which(c("b", "bk", "k") %in% tolower(object$vars))
    if (length(btype) > 1) 
        stop("all behavioural responses must be of same type in sim.detect")
    if (length(btype) == 0) 
        btype <- 0
    if (is.null(object$details$ignoreusage)) 
        object$details$ignoreusage <- FALSE
    if (is.null(object$details$miscparm)) 
        object$details$miscparm <- numeric(4)
    N <- sapply(popnlist, nrow)
    nocc <- if (ms(object)) {
        sapply(object$capthist, ncol)
    }else ncol(object$capthist)
    ndet <- if (ms(object)) { 
        sapply(traps(object$capthist), nrow)
    }else nrow(traps(object$capthist))
    sessionlevels <- session(object$capthist)
    nsession <- length(sessionlevels)
    sessmask <- object$mask
    if (!ms(sessmask)) 
        sessmask <- list(sessmask)
    grplevels <- secr:::group.levels(object$capthist, object$groups, 
        sep = ".")
    beta <- object$fit$par
    userd <- is.null(object$details$userdist)
    ncores <- setNumThreads()
    grain <- if (ncores == 1) {
        0
    }else 1
    if (length(unique(object$design0$PIA)) == 1) {
        design0 <- object$design0
        dim0 <- dim(design0$PIA)
        design0$PIA <- array(1, dim = c(dim0[1], max(N), dim0[3:5]))
        dummyCH <- NULL
    } else {
        dummyCH <- dummycapthist(object$capthist, popnlist, fillvalue = 1)
        design0 <- secr.design.MS(dummyCH, object$model, object$timecov, 
            object$sessioncov, object$groups, object$hcov, object$dframe, 
            naive = TRUE, CL = object$CL, ignoreusage = object$details$ignoreusage, 
            contrasts = object$details$contrasts)
    }
    for (i in 1:nsession) {
        if (is.null(object$groups)) {
            matchedgroup <- rep(1, N[i])
        }
        else {
            newgroup <- secr:::group.factor(popnlist[[i]], object$groups)
            matchedgroup <- (1:nlevels(newgroup))[as.numeric(newgroup)]
        }
        design0$PIA[i, 1:N[i], , , ] <- object$design0$PIA[i, 
            matchedgroup, , , ]
    }
    realparval0 <- secr:::makerealparameters(design0, beta, object$parindx, 
        object$link, object$fixed)
    if (btype > 0) {
        if (is.null(dummyCH)) 
            dummyCH <- secr:::dummycapthist(object$capthist, popnlist, 
                fillvalue = 1)
        design1 <- secr.design.MS(dummyCH, object$model, object$timecov, 
            object$sessioncov, object$groups, object$hcov, object$dframe, 
            ignoreusage = object$details$ignoreusage, contrasts = object$details$contrasts, 
            CL = object$CL)
        realparval1 <- secr:::makerealparameters(design1, beta, object$parindx, 
            object$link, object$fixed)
    } else {
        design1 <- design0
        realparval1 <- realparval0
    }
    if (!object$CL) {
        D <- secr:::getD(object$designD, beta, sessmask, object$parindx, 
            object$link, object$fixed, grplevels, sessionlevels, 
            parameter = "D")
    }
    NE <- secr:::getD(object$designNE, beta, sessmask, object$parindx, 
        object$link, object$fixed, grplevels, sessionlevels, 
        parameter = "noneuc")
    output <- list()
    for (sessnum in 1:nsession) {
        if (ms(object)) {
            s <- ncol(object$capthist[[sessnum]])
            session.traps <- traps(object$capthist)[[sessnum]]
            session.mask <- if (userd | (object$details$param %in% 
                c(2, 6))) {
                object$mask[[sessnum]]
           } else NULL
            Dtemp <- if (object$details$param %in% c(4:6)){ 
                predict(object)[[sessnum]]["D", "estimate"]
           } else NA
            nmash <- attr(object$capthist[[sessnum]], "n.mash")
        } else {
            s <- ncol(object$capthist)
            session.traps <- traps(object$capthist)
            session.mask <- if (userd | (object$details$param %in% 
                c(2, 6))) {
                object$mask
           } else NULL
            Dtemp <- if (object$details$param %in% c(4:6)) {
                predict(object)["D", "estimate"]
           } else NA
            nmash <- attr(object$capthist, "n.mash")
        }
        Xrealparval0 <- secr:::reparameterize(realparval0, object$detectfn, 
            object$details, session.mask, session.traps, Dtemp, 
            s)
        Xrealparval1 <- secr:::reparameterize(realparval1, object$detectfn, 
            object$details, session.mask, session.traps, Dtemp, 
            s)
        session.animals <- popnlist[[sessnum]]
        nmiscparm <- length(object$details$miscparm)
        if (nmiscparm > 0) {
            miscindx <- max(unlist(object$parindx)) + (1:nmiscparm)
            attr(session.mask, "miscparm") <- coef(object)[miscindx, 
                1]
        }
        dettype <- secr:::detectorcode(session.traps, MLonly = FALSE, 
            noccasions = nocc[sessnum])
        if (!all(dettype %in% c(-1, 0, 1, 2, 3, 4, 5, 6, 7))) 
            stop("detector type ", paste(detector(session.traps)), 
                " not implemented")
        binomN <- secr:::expandbinomN(object$details$binomN, dettype)
        if (all(detector(session.traps) %in% secr:::.localstuff$polydetectors)) {
            k <- c(table(polyID(session.traps)), 0)
            K <- length(k) - 1
        }else {
            k <- nrow(session.traps)
            K <- k
        }
        usge <- usage(session.traps)
        if (is.null(usge) | object$details$ignoreusage) 
            usge <- matrix(1, nrow = K, ncol = s)
        NR <- N[sessnum]
        if (all(detector(session.traps) %in% secr:::.localstuff$exclusivedetectors)) {
            maxdet <- NR * s
        } else {
            maxdet <- NR * s * K * maxperpoly
        }
        if ((object$detectfn == 12) || (object$detectfn == 13)) {
            object$details$miscparm[2:3] <- beta[max(unlist(object$parindx)) + 
                (1:2)]
        }
        knownclass <- secr:::getknownclass(session.animals, object$details$nmix,
            object$hcov)
        pmix <- secr:::getpmixall(design0$PIA, Xrealparval0)
        if (length(unique(dettype)) > 1 | length(unique(binomN)) > 
            1) 
            stop("simulation not yet updated for varying detector type")
        if (all(dettype %in% c(-1, 0, 1, 2, 5, 8))) {
            if (is.function(object$details$userdist)) {
                noneuc <- secr:::getmaskpar(!is.null(NE), NE, nrow(session.mask), 
                  sessnum, FALSE, NULL)
                density <- secr:::getmaskpar(!object$CL, D, nrow(session.mask), 
                  sessnum, object$details$unmash, nmash)
                distmat2 <- secr:::getuserdist(session.traps, session.animals, 
                  object$details$userdist, sessnum, noneuc[, 
                    1], density[, 1], object$details$miscparm, 
                  object$detectfn == 20)
            } else {
                distmat2 <- secr:::getdistmat2(session.traps, session.animals, 
                  object$details$userdist, object$detectfn == 
                    20)
            }
            gkhk0 <- secr:::makegkPointcpp(as.integer(object$detectfn), 
                as.integer(grain), as.integer(ncores), as.matrix(Xrealparval0), 
                as.matrix(distmat2), as.double(object$details$miscparm))
            gkhk <- secr:::makegkPointcpp(as.integer(object$detectfn), 
                as.integer(grain), as.integer(ncores), as.matrix(Xrealparval1), 
                as.matrix(distmat2), as.double(object$details$miscparm))
            if (any(dettype == 8)) {
                stop("sim.secr not ready for capped detectors")
            }
        }
				
				all_expected0 <- aperm(
					array(gkhk0$gk, dim=c(nrow(Xrealparval0), K, NR)),
					perm=c(1,3,2))
				all_expected <- aperm(
					array(gkhk$gk, dim=c(nrow(Xrealparval1), K, NR)),
					perm=c(1,3,2))
				expected <- array(0, dim=dim(design0$PIA)[2:5])
				
				# Cumulative captures to infer behavioural effects
				if(btype > 0){
					carray <- abind::abind(array(0, dim=c(NR, 1, K)),
						abind::abind(aperm(apply(object$capthist, c(1,3), cumsum),
								perm=c(2,1,3)), 
							array(0, dim=c(NR - nrow(object$capthist), s, K)), along=1), 
						along=2)[, -(s + 1), ] 
					if(Markov){
						# Find days with transient effect for each individual and trap
						Carray <- abind::abind(carray[,1:2,],
							(carray[, 3:s, ] - carray[, 3:s - 1, ]) > 0,
							along=2)
					} else { # Find days with learned effect for each individual and trap
						Carray <- carray > 0
					}
					if(btype == 1){
						# Apply behavioural effect to all traps following capture of individual
						Carray <- replicate(K, rowSums(Carray, dim=2) > 0)
					}
					if(btype == 3){
						# Apply behavioural effect to all individuals following capture at trap
						Carray <- aperm(replicate(NR, colSums(Carray) > 0), perm=c(3,1,2))
					}
					PIA1 <- design0$PIA * 2 - 1
					for(i in 1:dim(PIA1)[5]){
						PIA1[1, , , , i] <- PIA1[1, , , , i] + Carray
						for(j in 1:NR){
							PIA1[1, j, , , i] <- PIA1[1, j, , , i] * t(usge)
						}
					}
				}
					
						
				if(btype == 0){ # with no behavioural effect
					for(i in 1:nrow(Xrealparval0)){
						for(j in 1:s){
							for(k in 1:dim(design0$PIA)[5]){
								expected[ ,j,,k][design0$PIA[1, ,j, , k] == i] <- 
									all_expected0[i,,][design0$PIA[1, ,j, , k] == i]
							}
						}
					}
				} else {
					for(i in 1:nrow(Xrealparval1)){
						for(j in 1:s){
							for(k in 1:dim(PIA1)[5]){
								expected[ ,j,,k][PIA1[1, ,j, , k] == i] <- 
									all_expected[i,,][PIA1[1, ,j, , k] == i]
							}
						}
					}
					design1$PIA <- PIA1
				}
				
				if(dim(design0$PIA)[5] > 1){
					indcov <- covariates(object$capthist)[object$hcov]
					for(i in 2:dim(design0$PIA)[5]){
						expected[indcov %in% levels(indcov)[i], , , 1] <- 
							expected[indcov %in% levels(indcov)[i], , , i]
					}
					expected <- array(expected[,,,1], dim=c(dim(expected)[1:3], 1))
				}
				
				if(dim(design0$PIA)[3] == 1){
					expected[,1,,] <- expected[,1,,] * t(replicate(NR, usge, 
						simplify="matrix"))
				}
				
        if (all(dettype %in% c(-1, 0, 1, 2, 8))) {
            temp <- secr:::simdetectpointcpp(as.integer(dettype[1]), 
                as.integer(NR), as.integer(nrow(Xrealparval0)), 
                as.integer(nrow(Xrealparval1)), as.double(gkhk0$gk), 
                as.double(gkhk$gk), as.double(gkhk0$hk), as.double(gkhk$hk), 
                as.integer(design0$PIA[sessnum, 1:NR, 1:s, 1:K, 
                  ]), as.integer(design1$PIA[sessnum, 1:NR, 1:s, 
                  1:K, ]), as.integer(object$details$nmix), as.integer(knownclass), 
                as.double(pmix), as.matrix(usge), as.integer(btype), 
                as.integer(Markov), as.integer(binomN))
        } else if (all(dettype %in% c(3, 4, 6, 7))) {
            temp <- secr:::simdetectpolycpp(as.integer(dettype[1]), 
                as.integer(object$detectfn), as.integer(object$details$nmix), 
                as.integer(btype), as.integer(Markov), as.integer(k), 
                as.matrix(session.animals), as.matrix(session.traps), 
                as.matrix(Xrealparval0), as.matrix(Xrealparval1), 
                as.integer(design0$PIA[sessnum, 1:NR, 1:s, 1:K, 
                  ]), as.integer(design1$PIA[sessnum, 1:NR, 1:s, 
                  1:K, ]), as.integer(knownclass), as.double(pmix), 
                as.matrix(usge), as.integer(binomN), as.integer(maxperpoly))
            if ((temp$resultcode == 2) & (any(dettype %in% c(6, 
                7)))) {
                stop(">100 detections per animal per polygon per occasion")
            }
        } else if (all(dettype %in% c(5))) {
            if (is.null(object$details$cutval)) 
                stop("sim.detect for signal model requires object$details$cutval")
            temp <- simdetectsignalcpp(as.integer(dettype[1]), 
                as.integer(object$details$nmix), as.integer(object$detectfn), 
                as.integer(object$details$cutval), as.matrix(Xrealparval0), 
                as.integer(design0$PIA[sessnum, 1:NR, 1:s, 1:K, 
                  ]), as.integer(pmix), as.integer(knownclass), 
                as.matrix(session.animals), as.matrix(session.traps), 
                as.matrix(distmat2), as.matrix(usge), as.double(object$details$miscparm))
        }
        if (temp$resultcode != 0) {
            stop("simulated detection failed, code ", temp$resultcode)
        }
			
				
        w <- array(temp$value, dim = c(s, K, NR), dimnames = list(1:s, 
            NULL, NULL))
        w <- aperm(w, c(3, 1, 2))
        w <- w[apply(w, 1, sum) > 0, , , drop = FALSE]
        class(w) <- "capthist"
        traps(w) <- session.traps
        session(w) <- sessionlevels[sessnum]
        if (!is.null(covariates(popnlist))) {
            covariates(w) <- subset(covariates(popnlist[[sessnum]]), 
                subset = as.logical(temp$caught))
        }
        if (any(dettype %in% c(5, 12)) & (temp$n > 0)) {
            nd <- sum(abs(w))
            signal(w) <- temp$signal[1:nd]
            if ((object$detectfn == 12) || (object$detectfn == 
                13)) 
                noise(w) <- temp$signal[(nd + 1):(2 * nd)]
            attr(w, "cutval") <- object$details$cutval
        }  else {
            attr(w, "signalframe") <- NULL
            attr(w, "cutval") <- NULL
        }
        if (any(dettype %in% c(3, 4, 6, 7))) {
            nd <- sum(abs(w))
            if (nd > 0) {
                xymat <- temp$detectedXY[1:nd, 1:2, drop = FALSE]
           } else xymat <- matrix(nrow = 0, ncol = 2)
            detectedXY <- data.frame(xymat)
            names(detectedXY) <- c("x", "y")
            attr(w, "detectedXY") <- detectedXY
        }  else attr(w, "detectedXY") <- NULL
        if (renumber & (temp$n > 0)) {
            rownames(w) <- 1:temp$n
        }   else {
            rownames(w) <- rownames(session.animals)[as.logical(temp$caught)]
        }
        output[[sessnum]] <- w
    }
    if (nsession == 1) {
        output <- output[[1]]
		} else {
        names(output) <- sessionlevels
        class(output) <- c("capthist", "list")
    }
    
		return(list(expected=expected, capthist=output))
}