
Brownie <- function(tag_data, fix_M = FALSE, M = 0.06, latency = 4, report_rate = 1,
                    fix_retain = TRUE, retain = 1) {
  get_TMB_model <- try(dyn.load("Brownie/brownie.dll"), silent = TRUE)
  if(is.character(get_TMB_model)) TMB::compile("Brownie/brownie.cpp")
  dyn.load("Brownie/brownie.dll")
  on.exit(dyn.unload("Brownie/brownie.dll"))
  
  TMB_data <- list(N_rel = tag_data$N_tag, N_NR = tag_data$N_tag - rowSums(tag_data$N_recap, na.rm = TRUE),
                   N_recap = tag_data$N_recap, n_group = length(tag_data$N_tag), n_y = ncol(tag_data$N_recap),
                   latency = latency)
  
  TMB_params <- list(log_M = log(M), log_F = log(rep(0.05, TMB_data$n_y)),
                     log_retain = log(retain), log_report_rate = log(report_rate), 
                     log_latent_report = rep(0.001, latency) %>% log())

  map <- list()
  if(fix_M) map$log_M <- factor(NA)
  if(fix_retain) map$log_retain <- factor(NA)
  map$log_report_rate <- factor(NA)
  #map$log_latent_report <- c(NA, 1:(latency - 1)) %>% factor()
  
  obj <- TMB::MakeADFun(data = TMB_data, parameters = TMB_params, map = map, DLL = "brownie", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  SD <- TMB::sdreport(obj)
  report <- obj$report(opt$par)
  return(list(obj = obj, opt = opt, SD = SD, report = report))
}

plot_F_Brownie <- function(Brownie, ylim = c(0, 0.1)) {
  Find <- rownames(Brownie$SD$cov.fixed) == "log_F"
  
  SD_logF <- Brownie$SD$cov.fixed[Find, Find] %>% diag() %>% sqrt()
  #SD_logF[1] <- 0.001
  FF_lower <- log(Brownie$report$F) - 2 * SD_logF
  FF_upper <- log(Brownie$report$F) + 2 * SD_logF
  plot_timeseries <- function(Year, obs, fit = NULL, obs_upper = NULL, obs_lower = NULL, fit_linewidth = 3,
                              fit_color = "red", label = "Fishing mortality") {
    plot(Year, obs, typ = 'o', ylab = label, ylim = ylim)
    arrows(Year, obs_lower, Year, obs_upper, length = 0.025, angle = 90, code = 3)
    if(!is.null(fit)) lines(Year, fit, lwd = fit_linewidth, col = fit_color)
    abline(h = 0, col = 'grey')
    invisible()
  }
  plot_timeseries(2009:2020, Brownie$report$F, obs_lower = exp(FF_lower), obs_upper = exp(FF_upper))
}

plot_Brownie_fit <- function(Brownie, diag = TRUE) {
  
  plot_tags <- function(Year_rel = 2009:2020, Year_recap = 2009:2020, obs, fit = NULL, N_rel, 
                        ylab = "Recaptures", fit_linewidth = 3, fit_color = "red") {
    obs[lower.tri(obs, !diag)] <- NA_real_
    if(!is.null(fit)) fit[lower.tri(fit, !diag)] <- NA_real_
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    par(mfcol = c(4, 3), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))
    ylim <- c(0, 1.1 * max(obs, fit, na.rm = TRUE))
    yaxp <- c(0, max(pretty(ylim, n = 4)), 4)
    las <- 1
    
    for(i in 1:length(Year_rel)) {
      yaxt <- ifelse(i %% 16 %in% c(1:4), "s", "n") # TRUE = first column
      xaxt <- ifelse(i < length(Year_rel) & i %% 4 %in% c(1:3), "n", "s") # TRUE = first three rows
      
      plot(Year_recap, obs[i, ], typ = "n", ylim = ylim, yaxp = yaxp, xaxt = xaxt, yaxt = yaxt, las = las)
      abline(h = 0, col = "grey")
      abline(v = Year_rel[i], lty = 3)
      lines(Year_recap, obs[i, ], typ = "o")
      if(!is.null(fit)) {
        if(sum(!is.na(fit[i, ])) == 1) {
          points(Year_recap, fit[i, ], pch = 16, lwd = fit_linewidth, col = fit_color)
        } else {
          lines(Year_recap, fit[i, ], lwd = fit_linewidth, col = fit_color)
        }
      }
      legend(ifelse(i == 1, "topright", "topleft"), legend = c(Year_rel[i], paste0("N released = ", N_rel[i])), 
             bty = "n", xjust = 1)
      
      if(i %% 16 == 1) {
        mtext("Year", side = 1, line = 3, outer = TRUE)
        mtext(ylab, side = 2, line = 3.5, outer = TRUE)
      }
    }
    return(invisible())
  }
  plot_tags(obs = Brownie$obj$env$data$N_recap, fit = Brownie$report$C_pred, N_rel = Brownie$obj$env$data$N_rel)
}

plot_Brownie_resid <- function(Brownie, bubble_adj = 5) {
  
  plot_resid <- function(Year_rel = 2009:2020, Year_recap = 2009:2020, obs, fit = NULL, N_rel, 
                         bubble_adj) {
    bubble_color = c("black", "white")
    
    obs[lower.tri(obs)] <- NA_real_
    if(!is.null(fit)) fit[lower.tri(fit)] <- NA_real_
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    obs_prob <- obs/N_rel
    fit_prob <- fit/N_rel
    resid <- N_rel * (obs_prob - fit_prob) / sqrt(N_rel * fit_prob * (1 - fit_prob))
    
    diameter_max <- bubble_adj / pmin(10, max(abs(resid), na.rm = TRUE))
    plot(NULL, NULL, typ = 'n', xlim = range(Year_recap), xlab = "Year recapture", 
         ylim = range(Year_rel), ylab = "Year release")
    abline(a = 0, b = 1, lty = 3)
    x_mat <- matrix(Year_recap, ncol = ncol(resid), nrow = nrow(resid), byrow = TRUE)
    y_mat <- matrix(Year_rel, ncol = ncol(resid), nrow = nrow(resid))
    isPositive <- resid > 0
    points(x_mat[!isPositive], y_mat[!isPositive], cex = pmin(0.5 * diameter_max * abs(resid[!isPositive]), diameter_max), pch = 21, bg = bubble_color[1])
    points(x_mat[isPositive], y_mat[isPositive], cex = pmin(0.5 * diameter_max * resid[isPositive], diameter_max), pch = 21, bg = bubble_color[2])
    legend("topleft", legend = c("<-10", "-1", "1", ">10"),
           pt.cex = c(diameter_max, 0.5 * diameter_max, 0.5 * diameter_max, diameter_max),
           pt.bg = rep(bubble_color, each = 2), pch = 21, horiz = TRUE)

    return(invisible())
  }
  plot_resid(obs = Brownie$obj$env$data$N_recap, fit = Brownie$report$C_pred, N_rel = Brownie$obj$env$data$N_rel,
             bubble_adj = bubble_adj)
}

