ddply(.data = subset(results, stat == "phi"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "zeta"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "logit phi"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "log zeta"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "mu"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "lambda"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "log mu"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
ddply(.data = subset(results, stat == "log lambda"), 
      .variables = .(K_intervals, k_priors, theta), summarize, bias = mean(bias),
      median_bias = mean(median_bias), variance = mean(var))
