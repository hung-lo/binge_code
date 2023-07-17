#LLM with actual opto data

# ---- activate packages ----
library(lme4)
library(ggplot2)
library(ggbeeswarm)
library(brms)
library(glue)
library(afex)
library(patchwork)
library(sjPlot)
library(sjstats)
library(lattice)
library(emmeans)
library(RColorBrewer)

## plotting setting


my_color_map <- c('#56b4e9','#e69f00','#009e73','#f0e442','#0072b2','#d55e00',
                  '#cc79a7','#999999','#000000')

## Set fig output folder
fig_output_folder <- '/Users/hunglo/Documents/inscopix_csv/fig_output/'


# loop through all cell types in the df
stats_output <-data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("estimate", "pvals","feeding_type", "Celltype"))))
df <- read.csv('/Users/hunglo/Documents/inscopix_csv/processed_csv/slow_binge_cohensd_window_all.csv',header=TRUE)

Celltype_list <- unique(unlist(df$Celltype))

for (i in 1:length(Celltype_list)){
  Celltype = Celltype_list[i]
  print(Celltype_list[i])
  print(nrow(df))
  df_sub <- subset(df,anosmic=='False' & fasted=='False' & window_size=='0_4' & Celltype == Celltype_list[i])
  print(nrow(df_sub))

    # # ---- replace actual mouse ids with M1,M2...for simplicity  ----
  # for(i in  seq_along(unique(df_sub$mouse_id))){
  #   df_sub$mouse_id[df_sub$mouse_id==unique(df_sub$mouse_id)[i]] <- glue('M{i}')
  # }

  # ---- Quick plotting of the data  ----
  
  # ## plotting for LED state vs milk consumption with individual mice as a column
  # ggplot(data = df_sub, aes(x=cohens_d,y=milk_consumption))+
  #   geom_point(alpha=0.75)+
  #   geom_smooth(method='lm',alpha=0.5)+
  #   facet_wrap(~feeding.type+Celltype,nrow = 2)+
  #   scale_color_manual(values=my_color_map)
  
  ## fitting LMM
  
  df_binge <- subset(df_sub, feeding_type=='binge')
  formula_inter <- milk_consumption ~ 1 +cohens_d+ (1|mouse_id)
  df_binge.model = lmer(formula_inter, data=df_binge,REML=FALSE)
  
  # check again if necessary for limiting data point here
  ## if limiting only the dates when there is binge data 
  # df_slow <-df_sub[is.element(df_sub$date, df_binge$date),] 
  # df_slow <- subset(df_slow, feeding_type=='slow')
  df_slow <- subset(df_sub, feeding_type=='slow')
  formula_inter <- milk_consumption ~ 1 +cohens_d+(1|mouse_id)
  df_slow.model = lmer(formula_inter, data=df_slow,REML=FALSE)
  
  ## get coef for cohen's d on milk consumption
  est_b <- df_binge.model@beta[2]
  p_b <- anova(df_binge.model)$"Pr(>F)"
  
  est_s <- df_slow.model@beta[2]
  p_s <- anova(df_slow.model)$"Pr(>F)"
  
  stats_output[nrow(stats_output)+1,] = c(est_b,p_b,"binge",Celltype)
  stats_output[nrow(stats_output)+1,] = c(est_s,p_s,"slow",Celltype)
  
}

# stats_output <- data.frame(estimate= c(est_b,est_s),pvals = c(p_b,p_s), feeding_type = c('binge','slow'),Celltype = Celltype)
stats_output

write.csv(stats_output, file = "/Users/hunglo/Documents/inscopix_csv/stats_output/LMM_cohensd_ensure_consumption.csv")



