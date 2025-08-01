half_p_curve_ev = function(X2){
  crit = qchisq(.025,1,lower.tail = FALSE)
  X2 = X2[X2>crit]
  pcurve_X2(X2,1,alpha_bound = .025)
}

both_p_curve_ev = function(X2, test = '2015'){
  crit = qchisq(.05,1,lower.tail = FALSE)
  X2 = X2[X2>crit]
  c(full = pcurve_X2(X2)[2,1,test], half = half_p_curve_ev(X2)[2,1,test])
}

construct_nonmono_data = function(X2,test,nudge=.00001,max=10,len=2048){
  crit025 = qchisq(.025,1,lower.tail = FALSE)
  crit05 = qchisq(.05,1,lower.tail = FALSE)
  
  tibble(
    x2 = seq(crit05+nudge,max,length.out = len),
    pvals = purrr::map_df(x2, ~both_p_curve_ev(c(X2,.),test=test)),
  ) |> 
    tidyr::unnest(cols=c(pvals)) |>
    add_row(
      x2 = crit05, full = if_else(test=='2015',1,NA)
    ) |>
    add_row(
      x2 = crit025, half = ifelse(test=='2015',1,NA)
    ) |>
    arrange(x2) |>
    tidyr::fill(full, half, .direction = 'up') |>
    mutate(
      sig = (half<=.1 & full<=.1) | half<=.05,
      sig_reason = case_when(
        (half<=.05 & full<=.1) ~ "Both",
        (half<=.1 & full<=.1) ~ "Half<.1,Full<.1",
        half<.05 ~ "Half<.05",
        TRUE ~ 'No rejection'
      )
    ) |>
    tidyr::pivot_longer(cols = c(full,half), names_to = 'pcurve', values_to = 'p') |>
    mutate(
      group = case_when(
        pcurve=='full' ~ 'full',
        pcurve=='half' & x2<crit025 ~ 'half0',
        pcurve=='half' & x2>=crit025 ~ 'half1',
        TRUE ~ NA_character_
      )
    )
}

construct_nonmono_data_2d = function(X2,test,nudge=.00001,max=10,len=512){
  crit025 = qchisq(.025,1,lower.tail = FALSE)
  crit05 = qchisq(.05,1,lower.tail = FALSE)
  
  expand.grid(
    x2_1 = seq(crit05+nudge,max,length.out = len),
    x2_2 = seq(crit05+nudge,max,length.out = len)
  ) |>
    mutate(
      pvals = furrr::future_map2_dfr(x2_1, x2_2, ~both_p_curve_ev(c(X2,.x,.y),test=test),.progress = interactive()),
      #pvals = purrr::map2_dfr(x2_1, x2_2, ~both_p_curve_ev(c(X2,.x,.y),test=test)),
    ) |> 
    tidyr::unnest(cols=c(pvals)) |>
    arrange(x2_1,x2_2) |>
    tidyr::fill(full, half, .direction = 'up') |>
    mutate(
      half_sig = half<=.05,
      both_sig = (half<=.1 & full<=.1),
      sig = half_sig | both_sig
    ) 
}