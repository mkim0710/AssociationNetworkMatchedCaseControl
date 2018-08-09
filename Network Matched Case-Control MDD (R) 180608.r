
#@@@@ Network =====

function.deparse.cat.width = function(x, width.cutoff = 500) {
  cat(deparse(x, width.cutoff = width.cutoff), sep = "\n")
}

function.binary2numeric = function(x) {
  if (is.character(x)) {
    x = as.factor(x)
  }
  if (length(levels(x)) == 2) {
    x = as.numeric(x)
  }
  if (is.logical(x)) {
    x = as.numeric(x)
  } 
  x
}
   
function.pairwise.data_frame = function(vars) {
  library(tidyverse)
  vars.outer = outer(vars, vars, function(x, y) paste(x, y, sep = "&"))
  tmp = data_frame(vars = vars.outer[lower.tri(vars.outer)])
  out = tmp %>% separate(col = "vars", into = c("var_j", "var_i"), sep = "&")
  out = out[c("var_i", "var_j")]
  out
}

#@ x1x2z.partial_correlation_dbl() revision 180523 ==========
x1x2z.partial_correlation_dbl = function(x1, x2, z, cor_method = "pearson", convert_binary2numeric = F, p.value = F) {  # revision 180523
    if (convert_binary2numeric == T) {
        x1.binary2numeric = function.binary2numeric(x1)
        x2.binary2numeric = function.binary2numeric(x2)
    } else if (length(levels(x1)) > 2) {
        stop("error - length(levels(x1)) > 2")
    } else if (length(levels(x2)) > 2) {
        stop("error - length(levels(x2)) > 2")
    } else {
        x1.binary2numeric = (x1)
        x2.binary2numeric = (x2)
    }

    resid1 = lm(x1.binary2numeric ~ z)$residuals
    resid2 = lm(x2.binary2numeric ~ z)$residuals
    
    if (p.value == F) {
        unadjusted_cor = cor(x1.binary2numeric, x2.binary2numeric, method = cor_method)
        partial_cor = cor(resid1, resid2, method = cor_method)
        out = list(unadjusted_cor = unadjusted_cor, partial_cor = partial_cor) %>% unlist
    } else if (p.value == T) {
        unadjusted_cor = cor.test(x1.binary2numeric, x2.binary2numeric, method = cor_method) %>% {.[c("estimate", "conf.int", "p.value", "statistic")]} 
        if (!is.null(unadjusted_cor$conf.int)) unadjusted_cor$conf.int = unadjusted_cor$conf.int %>% as.vector %>% setNames(c("LL", "UL"))
        names(unadjusted_cor) = paste0("unadjusted_cor.", names(unadjusted_cor))
        unadjusted_cor = unadjusted_cor %>% unlist 
        names(unadjusted_cor)[1] = "unadjusted_cor"
        partial_cor = cor.test(resid1, resid2, method = cor_method) %>% {.[c("estimate", "conf.int", "p.value", "statistic")]} 
        if (!is.null(partial_cor$conf.int)) partial_cor$conf.int = partial_cor$conf.int %>% as.vector %>% setNames(c("LL", "UL"))
        partial_cor = partial_cor %>% unlist  
        names(partial_cor) = paste0("partial_cor.", names(partial_cor))
        names(partial_cor)[1] = "partial_cor"
        out = c(unadjusted_cor, partial_cor)  
    }
    names(out) = paste0(names(out), "_", cor_method)
    out
}

 
#@ function.CaseCtrl.pairwise.partial_correlation ======
function.CaseCtrl.pairwise.partial_correlation = function( 
  input
  , varname4Case = "is.Case"
  , varnames4pairwise = 
      c("AcquiredHypothyroidism", "Dementia", "Anemia", "Asthma", 
        "Cataract", "ChronicKidney", "COPD", "Diabetes",
        "HeartFailure", "Hyperlipidemia", "Hypertension",
        "IschemicHeart", "Osteoporosis", "Arthritis", "StrokeTIA")  
  , varnames4adjustment = "strata"
  , p.value = T
  , convert_binary2numeric = F
) {
  print("* dim(input)")
  print(dim(input))
  out = list()
  input2 = input[,which( 
    colSums( input[input[[varname4Case]] == T, ] %>% map_df(as.numeric) ) != 0 &
    colSums( input[input[[varname4Case]] == F, ] %>% map_df(as.numeric) ) != 0
  )]
  input2[[varname4Case]] = input[[varname4Case]]
  print("* Dropping variables with all-zero values, either in the case or the control:")
  function.deparse.cat.width(setdiff(names(input), names(input2)))
  
  varnames4pairwise2 = varnames4pairwise[varnames4pairwise %in% names(input2)]
  print("* Final variables for pairwise computation:")
  function.deparse.cat.width(varnames4pairwise2)
  print("* dim(input2)")
  print(dim(input2))

  pairwise.data_frame = function.pairwise.data_frame(varnames4pairwise2)
  
  if(is.logical(input2[[varname4Case]])) {
    input.list = list()
    input.list$Case = input2[input2[[varname4Case]] == T, ]
    input.list$Ctrl = input2[input2[[varname4Case]] == F, ]
  } else {
    print("error - !is.logical(input[[varname4Case]])")
    stop()
  }
  list4map = 1:nrow(pairwise.data_frame)
  names(list4map) = pairwise.data_frame %>% unite(col = "unite", var_i, var_j, sep = " & ") %>% unlist %>% unname
  list4map = as.list(list4map)
  
  map.map.out = map(
    input.list
    , function(input.CaseCtrl) {
      map(list4map
        , function(i) {
          x1x2z.partial_correlation_dbl(
            x1 = input.CaseCtrl[[ pairwise.data_frame$var_i[i] ]]
            , x2 = input.CaseCtrl[[ pairwise.data_frame$var_j[i] ]]
            , z = input.CaseCtrl[[ varnames4adjustment ]]
            , cor_method = "pearson"
            , p.value = p.value
            , convert_binary2numeric = convert_binary2numeric)
        }
      )
    }
  )
  out = map.map.out %>% map(function(ls) rownames_to_column(as.data.frame(t(as.data.frame(ls)))) )
}



#@ dataset.list.byCase.pairwise.partial_cor.p =====
dataset.list.byCase.pairwise.partial_cor.p =
    dataset.list %>%
    map(function(df) {

        df %>% function.CaseCtrl.pairwise.partial_correlation(
            varname4Case = "is.Case"
            , varnames4pairwise =
				c("AcquiredHypothyroidism", "Dementia", "Anemia", "Asthma", 
				"Cataract", "ChronicKidney", "COPD", "Diabetes",
				"HeartFailure", "Hyperlipidemia", "Hypertension",
				"IschemicHeart", "Osteoporosis", "Arthritis", "StrokeTIA")  
            , varnames4adjustment = "strata"
            , p.value = T
            , convert_binary2numeric = T
        )
    })

	

#@@@@ Network Rewire =====
dataset.list.byCase.pairwise.partial_cor.p.join =
    dataset.list.byCase.pairwise.partial_cor.p %>% 
    map(function(ls) {
        df_Case = ls$Case
        names(df_Case) = names(df_Case) %>% gsub("unadjusted_cor", "r_Case", .) %>% gsub("partial_cor", "partial_Case", .) %>% gsub("_pearson", "", .)
        df_Ctrl = ls$Ctrl
        names(df_Ctrl) = names(df_Ctrl) %>% gsub("unadjusted_cor", "r_Ctrl", .) %>% gsub("partial_cor", "partial_Ctrl", .) %>% gsub("_pearson", "", .)
        
        out = full_join(df_Case, df_Ctrl)
        
        out = out %>% separate(rowname, into = c("var_i", "var_j"), sep = "\\.{3}")
        
        out$r_diff = out$r_Case - out$r_Ctrl
        out$partial_diff = out$partial_Case - out$partial_Ctrl

        out
    })

function.Fisher_rz_Transform = function(r) {
    .5 * log ( (1+r)/(1-r) )
}

function.rz_diff = function(r1, r2, n1, n2 = NULL) {
    if(is.null(n2)) {
        n2 = n1
    }
    rz1 = .5 * log ( (1+r1)/(1-r1) )
    rz2 = .5 * log ( (1+r2)/(1-r2) )
    rz_diff = (rz1 - rz2) / ( 1/(n1-3) + 1/(n2-3) )^.5
}



dataset.list.byCase.pairwise.partial_cor.p.join.rz = map2(
    dataset.list.byCase.pairwise.partial_cor.p.join
    , (dataset.list %>% map(nrow))
    , function(df, n) {
        df$rz_Case = function.Fisher_rz_Transform(df$r_Case)
        df$rz_Ctrl = function.Fisher_rz_Transform(df$r_Ctrl)
        df$rz_diff = function.rz_diff(df$r_Case, df$r_Ctrl, n1 = n/2, n2 = n/2)
        df$rz_Case.partial = function.Fisher_rz_Transform(df$partial_Case)
        df$rz_Ctrl.partial = function.Fisher_rz_Transform(df$partial_Ctrl)
        df$rz_diff.partial = function.rz_diff(df$partial_Case, df$partial_Ctrl, n1 = n/2, n2 = n/2)
        df
    }
)






#@@@@ Network Rewire (Re2) Shuffle Test =====
#@ x1x2z.partial_correlation_scalar() ==========
x1x2z.partial_correlation_scalar = function(x1, x2, z, cor_method = "pearson", convert_binary2numeric = F) {

    if (convert_binary2numeric == T) {
        x1.binary2numeric = function.binary2numeric(x1)
        x2.binary2numeric = function.binary2numeric(x2)
    } else if (length(levels(x1)) > 2) {
        stop("error - length(levels(x1)) > 2")
    } else if (length(levels(x2)) > 2) {
        stop("error - length(levels(x2)) > 2")
    } else {
        x1.binary2numeric = (x1)
        x2.binary2numeric = (x2)
    }

    resid1 = lm(x1.binary2numeric ~ z)$residuals
    resid2 = lm(x2.binary2numeric ~ z)$residuals
    
    partial_cor = cor(resid1, resid2, method = cor_method)
    partial_cor
}


function.MatchingPairID_isExposed_PERSON_ID_shuffle = function(MatchingPairID_isExposed_PERSON_ID, var_MatchingPairID = "MatchingPairID", var_PERSON_ID = "PERSON_ID", seed = NULL) {
    MatchingPairID_isExposed_PERSON_ID = MatchingPairID_isExposed_PERSON_ID[order(MatchingPairID_isExposed_PERSON_ID[[var_MatchingPairID]]), ]
    
    sample.vec = sample(0:1, nrow(MatchingPairID_isExposed_PERSON_ID)/2, replace = T)
    index4Exposed = ( 1 : (nrow(MatchingPairID_isExposed_PERSON_ID)/2) ) * 2 - sample.vec
    
    out = MatchingPairID_isExposed_PERSON_ID
    out$isExposed_shuffle = F
    out$isExposed_shuffle[index4Exposed] = T
    out
}




#@ arguments to pass ----
var_MatchingPairID = "MatchingPairID"
var_PERSON_ID = "PERSON_ID"
var_isExposed = "is.Case"
var_i = dataset.list.byCase.pairwise.data_frame$var_i
var_j = dataset.list.byCase.pairwise.data_frame$var_j
var_z = "strata"
iteration = 2000
set.seed(1)
dataset.list.pairwise.shuffle2k_partial_cor_byExposed = map(
    dataset.list
    , function(data) {
        replicate(
            iteration
            , simplify2array(pmap(list(i = var_i, j = var_j, z = var_z), .f = function(i, j, z) {
                out = list()
                data.shuffle = function.MatchingPairID_isExposed_PERSON_ID_shuffle(data)
                data.var_isExposed_T = data.shuffle[data.shuffle[[var_isExposed]] == T, ]
                data.var_isExposed_F = data.shuffle[data.shuffle[[var_isExposed]] == F, ]
                out$var_isExposed = x1x2z.partial_correlation_scalar(data.var_isExposed_T[[i]], data.var_isExposed_T[[j]], data.var_isExposed_T[[z]])
                out$var_isUnexposed = x1x2z.partial_correlation_scalar(data.var_isExposed_F[[i]], data.var_isExposed_F[[j]], data.var_isExposed_F[[z]])
                out$diff = out$var_isExposed - out$var_isUnexposed
                out
            }))
        )
    }
)








