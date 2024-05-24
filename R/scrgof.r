# Function for testing GoF in SCR models
scrgof <- function(object, nsamp = 1000){
  
	require(spatstat)
	
  nind <- nrow(object$capthist) # no. of individuals
  xy <- object$mask[,c("x","y")] # mask pixels
  gof_obs_ik <- gof_obs_i <- gof_obs_j <- numeric(nsamp)
  gof_sim_ik <- gof_sim_i <- gof_sim_j <- numeric(nsamp)

  pars <- MASS::mvrnorm(n = nsamp,
                        mu = object$fit$estimate, 
                        Sigma = object$beta.vcv)  
  
  pb <- txtProgressBar(0, nsamp, style = 3)

  for(nS in 1:nsamp){
    
    tmp_object <- object
    tmp_object$fit$par <- pars[nS,] # use resampled parameters
    tmp_object$fit$estimate <- pars[nS,] # use resampled estimates
		
		# total individuals in population
    ntot <- ceiling(sum(exp(c(tmp_object$designD %*% 
			pars[nS, tmp_object$parindx$D])) *  attributes(tmp_object$mask)$area)) 
    newac <- matrix(NA, nrow = ntot, ncol = 2)
		
		if(class(predict(tmp_object)) == "list"){
			tmp_detpar <- lapply(predict(tmp_object), function(x){ 
				detpars <- as.list(x$estimate[-1])
				names(detpars) <- rownames(x)[-1]
				return(detpars)
			})
		} else {
			tmp_detpar <- list(detectpar(tmp_object))
		}

    pdobs <- fxi.secr(tmp_object)
		D <- predictDsurface(tmp_object, mask = tmp_object$mask)
    D <- covariates(D)$D.0
		pdunobs <- list()
		for(i in 1:length(tmp_detpar)){
			pd <- pdot(X=tmp_object$mask, traps=traps(tmp_object$capthist), 
				detectfn=tmp_object$detectfn, detectpar=tmp_detpar[[i]],
				noccasions=ncol(tmp_object$capthist))
			pdunobs[[i]] <- D * (1 - pd)
		}

    for(i in 1:nind){
		
			pd <- as.im(cbind(object$mask, pdobs[[i]]))
			tmpIndex <- rpoint(1, pd)				
      newac[i,] <- c(tmpIndex$x, tmpIndex$y)

    }
		
		if(length(tmp_detpar) > 1){
			indcov <- as.factor(unlist(covariates(object$capthist)[object$hcov]))
			tmp_indcov_assign <- rmultinom(ntot - nind, 1, 
				sapply(tmp_detpar, function(x) x$pmix))
			tmp_indcov <- rep(NA, ntot - nind)
			for(i in 1:length(tmp_detpar)){
				tmpIndex <- rpoint(sum(tmp_indcov_assign[i,]), 
					as.im(cbind(tmp_object$mask, pdunobs[[i]])))
				newac[((nind + 1):ntot)[tmp_indcov_assign[i,] != 0], ] <- 
					cbind(tmpIndex$x, tmpIndex$y)
				tmp_indcov[tmp_indcov_assign[i,] != 0] <- levels(indcov)[i]
			}
			tmp_indcov <- c(indcov, tmp_indcov)
		} else {
			tmpIndex <- rpoint(ntot - nind,	as.im(cbind(tmp_object$mask, pdunobs)))
			newac[(nind + 1):ntot, ] <- cbind(tmpIndex$x, tmpIndex$y)
		}

    tmp_ac <- data.frame(x = newac[,1], y = newac[,2])
    class(tmp_ac) <- c("popn", "data.frame")
    attr(tmp_ac, "boundingbox") <- attr(tmp_object$mask, "boundingbox")
		

    if(length(tmp_detpar) > 1){
      covariates(tmp_ac) <- data.frame(tmp_indcov)
			names(covariates(tmp_ac)) <- object$hcov
    }
		
    popnlist <- list(tmp_ac)
		test_stats <- sim.detect.new(tmp_object, popnlist)
		
    i_obs <- rownames(test_stats$capthist)
    i_notobs <- paste0(1:ntot)
    i_notobs <- i_notobs[!(i_notobs %in% i_obs)]
    newCH <- abind::abind(test_stats$capthist, 
                          array(0,c(ntot-dim(test_stats$capthist)[1],
														dim(test_stats$capthist)[2:3])),
                          along = 1)
    rownames(newCH) <- c(i_obs,i_notobs)
    newCH <- newCH[order(as.numeric(rownames(newCH))),,,drop=FALSE]


    ## Test: individual x traps
    EY <- apply(test_stats$expected, c(1,3), sum)
    OYobs <- rbind(apply(tmp_object$capthist,c(1,3),sum),
                         matrix(0,ntot-nind,nrow(traps(tmp_object$capthist))))
    OYsim <- apply(newCH,c(1,3),sum)
    gof_obs_ik[nS] <- FT_stat(OYobs, EY)
    gof_sim_ik[nS] <- FT_stat(OYsim, EY)

    ## Test: individual
    EY2 <- apply(test_stats$expected,1,sum)
    OYobs2 <- apply(OYobs,1,sum)
    OYsim2 <- apply(OYsim,1,sum)
    gof_obs_i[nS] <- FT_stat(OYobs2, EY2)
    gof_sim_i[nS] <- FT_stat(OYsim2, EY2)

    ## Test: traps
    EY3 <-  apply(test_stats$expected,3,sum)
    OYobs3 <- apply(OYobs,2,sum)
    OYsim3 <- apply(OYsim,2,sum)
    gof_obs_j[nS] <- FT_stat(OYobs3, EY3)
    gof_sim_j[nS] <- FT_stat(OYsim3, EY3)

    setTxtProgressBar(pb, nS)
  }

  close(pb)
  par(mfrow=c(2,2), oma=c(0,0,0,0))
  mm <- c(min(c(gof_sim_ik, gof_obs_ik)),max(c(gof_sim_ik, gof_obs_ik)))
  clr <- ifelse(gof_sim_ik > gof_obs_ik,adjustcolor("navy",0.5),adjustcolor("darkgreen",0.5))
  plot(gof_sim_ik, gof_obs_ik, pch=16, col=clr, las=1, asp=1, 
       xlim=c(mm[1],mm[2]), ylim=c(mm[1],mm[2]), 
       main = paste("pr(sim > obs):",round(mean(gof_sim_ik > gof_obs_ik),2)))
  abline(0, 1, col=2, lwd = 2)

  mm <- c(min(c(gof_sim_i, gof_obs_i)),max(c(gof_sim_i, gof_obs_i)))
  clr <- ifelse(gof_sim_i > gof_obs_i, adjustcolor("navy",0.5), adjustcolor("darkgreen",0.5))
  plot(gof_sim_i, gof_obs_i, pch=16, col=clr, las=1, asp=1, 
       xlim=c(mm[1],mm[2]), ylim=c(mm[1],mm[2]), 
       main = paste("pr(sim > obs):",round(mean(gof_sim_i > gof_obs_i),2)))
  abline(0, 1, col=2, lwd = 2)
  
  mm <- c(min(c(gof_sim_j, gof_obs_j)),max(c(gof_sim_j, gof_obs_j)))
  clr <- ifelse(gof_sim_j > gof_obs_j,adjustcolor("navy",0.5),adjustcolor("darkgreen",0.5))
  plot(gof_sim_j, gof_obs_j, pch=16, col=clr, las=1, asp=1, 
       xlim=c(mm[1],mm[2]), ylim=c(mm[1],mm[2]), 
       main = paste("pr(sim > obs):",round(mean(gof_sim_j > gof_obs_j),2)))
  abline(0, 1, col=2, lwd = 2)

  return(list(scrgof_pval = c(`FT-ind-trap`=mean(gof_sim_ik > gof_obs_ik),
                              `FT-individuals`=mean(gof_sim_i > gof_obs_i),
                              `FT-traps`=mean(gof_sim_j > gof_obs_j)),
              gof_ik = data.frame(simulated=gof_sim_ik, observed=gof_obs_ik),
              gof_i = data.frame(simulated=gof_sim_i, observed=gof_obs_i),
              gof_j = data.frame(simulated=gof_sim_j, observed=gof_obs_j)))
}


