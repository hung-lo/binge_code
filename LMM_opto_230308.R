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
library(forcats)

library(dplyr)



# ---- load csv file ----
df <- read.csv('/Users/hunglo/Documents/inscopix_csv/opto_csv/processed_csv/opto_data.csv',header=TRUE)
df <- df %>% group_by(mouse_id) %>% mutate(m = -mean(milk.consumption)) %>% arrange(desc(m)) %>% select(-m)

# df <- na.omit(df) # remove rows with NA
df_norm <- read.csv('/Users/hunglo/Documents/inscopix_csv/opto_csv/processed_csv/opto_data_norm.csv',header=TRUE)
# df_norm <- na.omit(df_norm) # remove rows with NA

## Set fig output folder
fig_output_folder <- '/Users/hunglo/Documents/inscopix_csv/opto_csv/fig_output/'

# ---- replace actual mouse ids with M1,M2...for simplicity  ----
# for(i in  seq_along(unique(df$mouse_id))){
#   df$mouse_id[df$mouse_id==unique(df$mouse_id)[i]] <- glue('M{i}')
# }
# for(i in  seq_along(unique(df_norm$mouse_id))){
#   df_norm$mouse_id[df_norm$mouse_id==unique(df_norm$mouse_id)[i]] <- glue('M{i}')
# }

#
my_color_map <- c('#56b4e9',
                  '#e69f00',
                  '#009e73',
                  '#f0e442',
                  '#0072b2',
                  '#d55e00',
                  '#cc79a7',
                  '#999999',
                  '#000000')
theme_set(theme_gray(base_size = 8))
# palette.colors(palette = "Okabe-Ito") ## where the colors from
# ---- mutate baseline appetite value to exact value (include 300 sec)  ----
## if baseline with the 300 sec baseline
#df$baseline.time <- df$baseline.time + 300


# ---- Quick plotting of the data  ----

## center title
theme_update(plot.title = element_text(hjust = 0.5))

## plotting for LED state vs milk consumption with individual mice as a column
ggplot(data = df, aes(x=LED.state,y=milk.consumption,color=virus))+
  # geom_boxplot()+
  geom_point(alpha=0.75)+
  # geom_beeswarm(alpha=0.75,aes(x=LED.state),size=1,dodge.width = 1)+
  # geom_beeswarm()+
  stat_summary(geom = "line", fun = mean, group = 1)+
  facet_wrap(~mouse_id,ncol=length(unique(df$mouse_id))/2)+
  scale_color_manual(values=my_color_map)


## same but plot all together
dodge <- position_dodge(0.75)
rect <- data.frame(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf)
virus.labs <- c("tdTomato","eOPN3")
names(virus.labs) <- c("0_tdTomato", "1_eOPN3")

# var <- "Total feeding duration after closed-loop trigger"


fig1 <- ggplot(data = df, aes(x= LED.state, y=milk.consumption))+
# fig1 <- ggplot(data = df, aes(x= reorder(LED.state,-milk.consumption,FUN=mean), y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "# of pump delivery")+
  scale_color_manual(values=c('#e69f00','#e69f00','#e69f00','#e69f00','#56b4e9','#56b4e9','#56b4e9','#56b4e9'))+ # 
  ggtitle("Total Ensure consumption")

fig1

ggsave(paste(fig_output_folder,'milk_consum.pdf',sep=''),fig1,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)



# ggsave(paste(fig_output_folder,'feeding_duration_after_trigger.png',sep=''),fig1,units="in",width=4,height=2.5,scale=1.25)


fig2 <- ggplot(data = df_norm, aes(x=LED.state,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  # stat_summary(aes(color=mouse_id),alpha=1,geom = "pointrange",fatten=-1,fun.data= mean_se,group=1,position=dodge)+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "# of pump delivery (normalized)")+
  scale_color_manual(values=my_color_map)+
  ggtitle("Total Ensure consumption (norm.)")
fig2
ggsave(paste(fig_output_folder,'milk_consum_norm.pdf',sep=''),fig2,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)


#### save other plots like licks & feeding duration
fig1 <- ggplot(data = df, aes(x=LED.state,y=feeding.duration))+
# fig1 <- ggplot(data = df, aes(x=LED.state,y=feeding.duration.after.closed.loop.trigger))+
# fig1 <- ggplot(data = df, aes(x=LED.state,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "feeding duration (second)")+
  scale_color_manual(values=my_color_map)+
  ggtitle("Total feeding duration")

fig1
ggsave(paste(fig_output_folder,'feeding_duration.pdf',sep=''),fig1,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)


fig1 <- ggplot(data = df_norm, aes(x=LED.state,y=feeding.duration))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=feeding.duration.after.closed.loop.trigger))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "feeding duration (normalized)")+
  scale_color_manual(values=my_color_map)+
  ggtitle("Total feeding duration (norm.)")

fig1
ggsave(paste(fig_output_folder,'feeding_duration_norm.pdf',sep=''),fig1,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)



fig1 <- ggplot(data = df, aes(x=LED.state,y=lick))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=feeding.duration.after.closed.loop.trigger))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "lick events")+
  scale_color_manual(values=my_color_map)+
  ggtitle("Total licks")

fig1
ggsave(paste(fig_output_folder,'lick.pdf',sep=''),fig1,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)


fig1 <- ggplot(data = df_norm, aes(x=LED.state,y=lick))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=feeding.duration.after.closed.loop.trigger))+
  # fig1 <- ggplot(data = df, aes(x=LED.state,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_beeswarm(alpha=0.8,aes(color=mouse_id),size=3,dodge.width = 0.75,cex=3,stroke=NA)+
  geom_smooth(method='lm',aes(group=mouse_id,color=mouse_id),linewidth=0.5, linetype = "dashed",position=dodge,se=FALSE)+
  stat_summary(geom = "line", fun = mean, group = 1,linewidth=1,colour='black')+
  facet_wrap(~virus, labeller = labeller(virus=virus.labs))+
  theme(panel.background = element_blank())+
  scale_x_discrete(name ="Session type")+
  theme(axis.title.x=element_blank(), legend.position = "none")+
  scale_y_continuous(name = "lick events (normalized)")+
  scale_color_manual(values=my_color_map)+
  ggtitle("Total licks (norm.)")

fig1
ggsave(paste(fig_output_folder,'lick_norm.pdf',sep=''),fig1,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)





# fig1/fig2

## combined plots
# p2/p + 
#   labs(subtitle = 'Individual mice',va='center')+ 
#   plot_annotation(title = 'Closed-loop opto-suppression of aPC neurons upon feeding')
  # plot_layout(width = c(2,1.5))

## correlation between exp session and milk consumption
LED_on_session <- seq(2,12,2)
rect <- data.frame(xmin=LED_on_session-0.5,xmax=LED_on_session+0.5,ymin=-Inf,ymax=Inf)
x_ticks <- seq(1,12,1)
fig3<-ggplot(data=df_norm,aes(x=session,y=milk.consumption,color=mouse_id))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_point(alpha=0.5,stroke=NA,size=2)+
  geom_line(aes(x=session,y=milk.consumption),alpha=0.3)+
  stat_summary(aes(x=session,y=milk.consumption,color=virus),color='black',alpha=1,geom = "line",fun=mean,group=1)+
  stat_summary(aes(x=session,y=milk.consumption,color=virus),color='black',alpha=1,geom = "pointrange",fun.data= mean_se,group=1,size=0.5,stroke=NA)+
  # geom_smooth(method='lm',alpha=0.1)+
  facet_wrap(~virus,ncol = 1,labeller = labeller(virus=virus.labs))+
  scale_x_continuous(breaks = x_ticks)+
  theme(panel.background = element_blank(), legend.position = "none")+
  scale_y_continuous(name ="Total Ensure consumption (normalized)")
fig3
ggsave(paste(fig_output_folder,'consumption_overview_norm.pdf',sep=''),fig3,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)



LED_on_session <- seq(2,12,2)
rect <- data.frame(xmin=LED_on_session-0.5,xmax=LED_on_session+0.5,ymin=-Inf,ymax=Inf)
x_ticks <- seq(1,12,1)
fig3<-ggplot(data=df,aes(x=session,y=milk.consumption,color=mouse_id))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  geom_point(alpha=0.5,stroke=NA,size=2)+
  geom_line(aes(x=session,y=milk.consumption),alpha=0.3)+
  stat_summary(aes(x=session,y=milk.consumption,color=virus),color='black',alpha=1,geom = "line",fun=mean,group=1)+
  stat_summary(aes(x=session,y=milk.consumption,color=virus),color='black',alpha=1,geom = "pointrange",fun.data= mean_se,group=1,size=0.5,stroke=NA)+
  # geom_smooth(method='lm',alpha=0.1)+
  facet_wrap(~virus,ncol = 1,labeller = labeller(virus=virus.labs))+
  scale_x_continuous(breaks = x_ticks)+
  theme(panel.background = element_blank(), legend.position = "none")+
  scale_y_continuous(name ="Total Ensure consumption")
fig3
ggsave(paste(fig_output_folder,'consumption_overview.pdf',sep=''),fig3,units="in",width=3.5,height=2.5,scale=1,useDingbats = FALSE)



ggplot(data=df_norm,aes(x=session,y=milk.consumption))+
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='lightgreen',alpha=0.25,inherit.aes = FALSE)+
  # geom_point(alpha=0.5)+
  # geom_line(aes(x=session,y=milk.consumption,color=virus),alpha=0.3)+
  aes(colour=factor(virus))+
  stat_summary(geom = "line",fun=mean)+
  stat_summary(geom = "pointrange",fun.data= mean_se)+
  # geom_smooth(method='lm',alpha=0.1)+
  # facet_wrap(~virus,ncol = 1)+
  scale_x_continuous(breaks = x_ticks)+
  theme(panel.background = element_blank())+
  scale_y_continuous(name ="Total milk consumption")

## correlation between baseline time (proxy for appetite) and milk consumption
ggplot(data=df,aes(x=baseline.time,y=milk.consumption,color=mouse_id))+
  geom_point()+
  geom_smooth(method='lm',alpha=0.1)+
  facet_wrap(~virus,ncol = 1)  


# ## plotting for sex and milk consumption
# ggplot(data=df,aes(x=LED.state,y=milk.consumption,color=sex))+
#   geom_boxplot()+
#   geom_point()+
#   geom_smooth(method='lm',alpha=0.1)+
#   facet_wrap(~virus*sex,ncol = 2)

## plotting sex and milk consumption (regardless of LED states)
ggplot(data=df,aes(x=sex,y=milk.consumption,color=sex))+
  geom_boxplot()+
  geom_point()+
  geom_smooth(method='lm',alpha=0.1)+
  facet_wrap(~virus,ncol = 2)

# ## inner state (bodyweight vs milk consumption, only on no LED sessions)
# ggplot(data=df[df$LED.state=='LED off',],aes(x=bodyweight,y=milk.consumption,color=mouse_id))+
#   # geom_boxplot()+
#   geom_point()+
#   geom_smooth(method='lm',alpha=0.1)+
#   facet_wrap(~virus,ncol = 1)
# 
# ## inner state (bodyweight vs milk consumption, only on no LED sessions)
# ggplot(data=df[df$LED.state=='LED off',],aes(x=bodyweight,y=baseline.time,color=mouse_id))+
#   # geom_boxplot()+
#   geom_point()+
#   geom_smooth(method='lm',se=F)+
#   facet_wrap(~virus,ncol = 1)

# ---- Mixed effect model with lme4 ----

## First choose which formula to put in to our model
# ## original
formula_inter <- milk.consumption ~ 1 +session+bodyweight+baseline.time+sex+ virus + LED.state+ virus * LED.state + (1|mouse_id)
formula_null  <- milk.consumption ~ 1 +session+bodyweight+baseline.time+sex+ virus + LED.state + (1|mouse_id)

# formula_inter <- lick ~ 1 +session+bodyweight+baseline.time+sex+ virus * LED.state + (1|mouse_id)
# formula_null  <- lick ~ 1 +session+bodyweight+baseline.time+sex+ virus + LED.state + (1|mouse_id)

# ## standardized
# formula_inter <- milk.consumption ~ 1 + sex + session + bwScaled +baselineTimeScaled+ virus * LED.state + (1|mouse_id)
# formula_null  <- milk.consumption ~ 1 + sex+ session + bwScaled +baselineTimeScaled+ virus + LED.state + (1|mouse_id)
# 
# formula_inter <- feeding.duration ~ 1 +session+bodyweight+baseline.time+sex+ virus * LED.state + (1|mouse_id)
# formula_null  <- feeding.duration ~ 1 +session+bodyweight+baseline.time+sex+ virus + LED.state + (1|mouse_id)

# formula_inter <- feeding.duration.after.closed.loop.trigger ~ 1 +session+bodyweight+baseline.time+sex+ virus * LED.state + (1|mouse_id)
# formula_null  <- feeding.duration.after.closed.loop.trigger ~ 1 +session+bodyweight+baseline.time+sex+ virus + LED.state + (1|mouse_id)
# 


### Use LMM to test the effects of virus & LED
## compare interaction between virus & LED
# dfS <- data.frame(df, bwScaled = scale(df$bodyweight), baselineTimeScaled = scale(df$baseline.time))

# df.model = glmer(formula_inter, data = dfS, family = Gamma(link="log"))
## ^^ other distribution family didn't converge, use simpler methods

df.model = lmer(formula_inter, data=df,REML=FALSE)
df.null  = lmer(formula_null,  data=df,REML=FALSE)

## test if some variables can be dropped, just it's also fine to just keep them (?)
# drop1(df.model,test="Chisq")
# drop1(df.null,test="Chisq")

summary(df.model)
summary(df.null)

# performance::icc(df.model)#,by_group = TRUE)
# performance::icc(df.null)#,by_group = TRUE)

## test the interaction of virus and LED state has effects (by comparing 2 models)
## where you can report the significance of the interaction of virus type * LED state
anova(df.null,df.model)

## plotting residuals with QQ plot
# qqnorm(resid(df.model))
# qqline(resid(df.model))
# 
# qqnorm(resid(df.null))
# qqline(resid(df.null))

## plot residuals 
res <- residuals(df.model, type="deviance")
plot(fitted(df.model), res)
abline(0,0)
plot(density(res))

# xyplot(res~milk.consumption,groups=mouse_id,data=dfS)
# xyplot(res~milk.consumption|mouse_id,data=dfS,aspect=1,
#        panel = function(x, y) {
#          panel.grid(h=-1, v= -2)
#          panel.xyplot(x, y,pch=16)
#          # panel.loess(x,y, span=1)
#        })
## some viz for the estimates (slopes)
## list of labels for plotting
label_list = c("Virus (eOPN3) x Session type (LED on)",
               "Session type (LED on)",
               "Virus (eOPN3)", 
               "Sex (Male)", 
               "Time of 20th delivery of milk",
               "Body Weight")
               #"Session #")

sjPlot::plot_model(df.model,
                   show.values=TRUE, show.p=TRUE,)#
                   # axis.labels = label_list)


##

milk.emmn <- emmeans(df.model, c("LED.state", "virus"), type="response")
milk.cont <- contrast(milk.emmn, method="trt.vs.ctrl", ref="LED off 1_eOPN3")
milk.cont
plot(milk.cont,xlab='Difference in milk consumption')

## plot marginal effects (after fixing other parameters other than virus & LED.state)
library(ggeffects)
plot(ggpredict(df.model, c("virus","LED.state")),dodge=0.75)


## print table as figure
# sjPlot::tab_model(df.model, 
#                   show.re.var= TRUE, 
#                   pred.labels =c("(Intercept)", "Body Weight", "Sex", "Virus", "Trial types","Virus (eOPN3) x Trial types (LED on) "),
#                   dv.labels= "Effects of opto-suppression of aPC neurons upon feeding")

# ## extract coefficient
# beta <- round(summary(df.model)$coefficients[, 1], 4)
# sigma_e <- round(attr(VarCorr(df.model), "sc"), 2)
# coef(summary(df.model))[ , "Estimate"]
# coef(summary(df.model))[ , "Pr(>|t|)"]

## print random intercept for each mouse
# coef(df.model)$mouse_id

# ## Use mixed() from afex package for testing each of the fixed effects, though probably not necessary
# df.model_afex <- mixed(formula_inter, data=df,REML=FALSE)
# df.model_afex


# ---- Mixed effect model with Bayesian brms ----

brmsFitNull  <- brms::brm(formula_null, data=df, family = "gamma", cores =4, control = list(adapt_delta=0.95), iter = 3000, save_pars = save_pars(all = T))
brmsFitInter <- brms::brm(formula_inter,data=df, family = "gamma", cores =4, control = list(adapt_delta=0.95), iter = 3000, save_pars = save_pars(all = T))
bayesplot::mcmc_pairs(brmsFitInter, pars=c("shape", "Intercept"))
brms::bayes_factor(brmsFitInter,brmsFitNull)

brms::conditional_effects(x = brmsFitInter)

bayesplot::pp_check(brmsFitInter)

y_rep <- posterior_predict(brmsFitInter)
## change the column for plotting here
bayesplot::ppc_dens_overlay_grouped(y = df$milk.consumption, yrep = y_rep[1:50,], group = interaction(df$LED.state, df$virus))
print(brmsFitInter)

plot(brmsFitInter)
# p_brms <-plot(brmsFitInter, variable = "^b_", regex = TRUE) 
## ^^ only plot population-level (fixed) effects, no sigma & sd mouse_id

## Estimation plots
sjPlot::plot_model(brmsFitInter,
                   show.values=TRUE, 
                   axis.labels = label_list)

