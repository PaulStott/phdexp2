#NB: an additional 60 participants were recruited on 24/08 as a power issue was identified

#load packages
library(tidyverse)      #for tidying and wrangling
library(lme4)           #for GLMM building
library(lmerTest)       #for better output from GLMMs
library(emmeans)        #for our planned comparisons
library(fitdistrplus)   #for Cullen and Frey plots

#control scientific notation - makes reading p-values easier
options("scipen"=100000, "digits"=10)

#load in pre-tidied and wrangled dataset
exp2_tidied <- read_csv("exp2_tidied_ALL-PARTICIPANTS_data.csv")

#check participant no. = 118 as expected as had to drop 2
exp2_tidied %>%
  distinct(Participant) %>%
  nrow()

#quick check
head(exp2_tidied)

#filter NONCE control trials
exp2_tidied <- exp2_tidied %>%
  filter(!ProbeType == "NONCE")

#histogram to check distribution - looks like gamma
exp2_tidied %>%
  ggplot(aes(x = RT)) + #specify what data to show on x axis
  geom_histogram(binwidth = 200) + #produce histogram with bins that are 200ms wide
  labs(x = "RT (ms)", y ="Frequency") + #x and y axes labels
  ggtitle("Hist. of RTs (ms)") + #add title
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centred
        text = element_text(size = 12)) #adjust rest of chart text - size 12

#Cullen and Frey Graph to check distribution more quantitatively - definitely gamma
descdist(exp2_tidied$RT,
           discrete = FALSE, #we are not dealing with discrete data - it's continuous RT data
           method = "unbiased", #don't use bootstrapping for sampling - use the exact data without replacement
           graph = TRUE, #produce a graph!
           obs.col = "red", #the colour of the "observation" marker - red in this case
           obs.pch = 15) #the shape of the "observation" marker - square in this case

#factorising
exp2_tidied <- exp2_tidied %>%
  mutate(Participant = factor(Participant),
         ProbeType = factor(ProbeType), 
         ProbePosition = factor(ProbePosition), 
         StimID = factor(StimID),
         Condition = factor(Condition))

#summary stats.
exp2_tidied %>%
  group_by(Condition) %>% #group by Condition
  summarise("MeanRT" = mean(RT), #mean for each condition
            "MedianRT" = median(RT), #median for each condition
            "StdDevRT" = sd(RT)) %>% #SD for each condition
  arrange(desc(MeanRT)) #arrange descending

#viz - some really extreme datapoints but GLMER should be fine with these
set.seed(5675)
exp2_tidied %>%
  mutate(Condition = fct_reorder(Condition, log(RT), .desc = FALSE)) %>% #reorder log of RT by Condition from lowest to highest from bottom up
  ggplot(aes(x = Condition, y = log(RT), colour = ProbeType)) + #specify what data to show on x and y axes, 
  geom_violin() + #produce a violin plot
  geom_jitter(width = .2, alpha = .08) + #add jitter to minimise overlap of points
  stat_summary(fun.data = mean_se, colour = "black") + 
  guides() + #legend
  labs(x = "ProbeType x ProbePosition x SubjectType", #x axis label
       y = "LogRT (ms)") + #y axis label
  ggtitle("Effect of Probe Type x Probe Position Interaction on RTs") + #add title'
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centered
        text = element_text(size = 12)) + #adjust rest of chart text - size 12
  coord_flip() #flip horizontal

#deviation coding factors
contrasts(exp2_tidied$ProbeType) <- contr.sum(2) #removed NONCE trials as these are controls
contrasts(exp2_tidied$ProbePosition) <- contr.sum(2)

#looking at covariates:
#checking for repetition effect/effect of trial number: weak -VE correlation
ggscatter(exp2_tidied, x = "TrialN", y = "RT", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Trial Number", ylab = "RT (ms)")

#looking at covariates:
#checking for correlation between reading times and RTs: weak +VE correlation
ggscatter(exp2_tidied, x = "RT", y = "ReadingTime", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Reading Time (ms) between Overt Sentential Subject Onset and LDT Probe Onset", ylab = "RT (ms)")

#looking at covariates:
#checking for correlation between ProbeFrequency and RTs: weak -VE correlation (NB: not specified in prereg!)
ggscatter(exp2_tidied, x = "ProbeFrequency", y = "RT", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Zipf Frequency of Probe Word", ylab = "RT (ms)")

#MODEL BUILDING
#gamma model w/ identity link (see Lo & Andrews 2015)
#maximal model: "keep it maximal" (see Bates et al. 2013)
QNP.exp2.gamma.identity1 <- glmer(RT ~ ProbeType * ProbePosition + #two factors w/ 2-way ixn.
                      (1+ProbeType * ProbePosition | StimID) + #ranefx.
                      (1+ProbeType * ProbePosition | Participant), #ranefx.
                    family = Gamma(link = "identity"), #gamma model w/ "identity" link fxn.
                    data = exp2_tidied,
                    control=glmerControl(optCtrl=list(maxfun=2e5),
                                         optimizer = "bobyqa")) #failed to converge w/ default evaluations & optimizer
summary(QNP.exp2.gamma.identity1) #singular fit

#reduce ranefx. (see Matuschek et al. 2017; Barr et al. 2013)
QNP.exp2.gamma.identity2 <- glmer(RT ~ ProbeType * ProbePosition +
                                    (1+ProbeType * ProbePosition | StimID) + 
                                    (1+ ProbePosition | Participant), #drop ProbeType
                                  family = Gamma(link = "identity"),
                                  data = exp2_tidied,
                                  control=glmerControl(optCtrl=list(maxfun=2e5),
                                                       optimizer = "bobyqa"))
summary(QNP.exp2.gamma.identity2) #singular fit

#reduce ranefx. further
QNP.exp2.gamma.identity3 <- glmer(RT ~ ProbeType * ProbePosition +
                                    (1+ProbeType | StimID) + #drop ProbePosition
                                    (1+ ProbePosition | Participant),
                                  family = Gamma(link = "identity"),
                                  data = exp2_tidied,
                                  control=glmerControl(optCtrl=list(maxfun=2e5),
                                                       optimizer = "bobyqa"))
summary(QNP.exp2.gamma.identity3) #convergence warnings

#reduce ranefx. further
QNP.exp2.gamma.identity4 <- glmer(RT ~ ProbeType * ProbePosition +
                                    (1 | StimID) + #drop ProbeType
                                    (1+ ProbePosition | Participant),
                                  family = Gamma(link = "identity"),
                                  data = exp2_tidied,
                                  control=glmerControl(optCtrl=list(maxfun=2e5),
                                                       optimizer = "bobyqa"))
summary(QNP.exp2.gamma.identity4) #convergence warnings

#reduce ranefx. further
QNP.exp2.gamma.identity5 <- glmer(RT ~ ProbeType * ProbePosition +
                                    (1 | StimID) + #drop ProbeType
                                    (1 | Participant), #drop ProbePosition
                                  family = Gamma(link = "identity"),
                                  data = exp2_tidied,
                                  control=glmerControl(optCtrl=list(maxfun=2e5),
                                                       optimizer = "bobyqa"))
summary(QNP.exp2.gamma.identity5) #converges

#plot it
emmip(QNP.exp2.gamma.identity5, ProbePosition ~ ProbeType) #x-over seems to be there though looks weak

#EMMs
emmeans(QNP.exp2.gamma.identity5, pairwise ~ ProbeType * ProbePosition, adjust = "none") #predictions unsupported
emmeans(QNP.exp2.gamma.identity5, pairwise ~ ProbeType, adjust = "none") #significant ProbeType effect
emmeans(QNP.exp2.gamma.identity5, pairwise ~ ProbePosition, adjust = "none") #no significant ProbePosition effect

#CONCLUSION: differences in RTs are driven by ProbeType, with RELATED probes being responded to
#faster than UNRELATED probes, regardless of the position at which the probe is
#presented (AFTER-V vs. AFTER-VP). But *WHY* do we not see the same effect we saw in Exp1 for QNPs?

#something to do with Exp1 DNP items being blindly converted to QNP for Exp2?

exp1_items <- read_csv("exp1_experimental_items.csv") %>%
  dplyr::select(stim_id, sub_type) %>%
  rename(exp1_sub_type = sub_type) %>%
  rename(StimID = stim_id) %>%
  mutate(StimID = factor(StimID))

exp2_tidied <- exp2_tidied %>%
  full_join(slice(exp1_items, 1:60), by = "StimID")

view(exp1_exp2_dets)

determiner_type.lm <- lm(RT ~ exp1_sub_type, data = exp2_tidied)
summary(determiner_type.lm)

determiner_type.mod <- glmer(RT ~ ProbeType * ProbePosition + exp1_sub_type +
                                (1 | StimID) +
                                (1 | Participant),
                              family = Gamma(link = "identity"),
                              data = exp2_tidied,
                              control=glmerControl(optCtrl=list(maxfun=2e5),
                                                   optimizer = "bobyqa"))

summary(determiner_type.mod)

#doesn't look like that's had any effect on Exp2 results...

