#quantile survival analysis 
os_quan <- os %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")

require("survival")
sfit <- survfit(Surv(OS.Time,OS)~quan_group,data = os_quan)
sfit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = os_quan)
ggsurvplot(sfit,pval = TRUE)
ggsurvplot(
  sfit,
  data = os_quan,
  size = 1,                 # change line size
  palette =
    c("#FF0000", "#0000FF"),# custom color palettes
  conf.int = F,  # Add confidence interval
  conf.int.style = "step",
  pval = TRUE, # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  legend.title = "quartile",
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
) +
  ggtitle("PFI-HDAC6")

ggsurvplot(
  stage4fit,
  data = stage4,
  size = 1,       # change line size
  conf.int = F,  # Add confidence interval
  pval = TRUE, # Add p-value
  font.x = c(10), 
  font.y = c(10),
  font.tickslab = c(10),
  break.time.by = 1,
  xlim = c(0,5),
  xlab = "Time in years",
  ylab = "Overall survival",
  surv.scale="percent",
  palette = c("#FF0000", "#0000FF"),
  legend  = c(0.15,0.4),
  legend.title = "quartile",
  legend.labs = c("High", "Low"),
  conf.int.style = "step",
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.20, # Useful to change when you have multiple groups
  risk.table.title="Numbers of risk",
  ggtheme = theme_bw(),# Change ggplot2 theme
  tables.theme = theme_cleantable()
) +
  ggtitle("CESC stageIV-HDAC6")


ggsurvplot(
  stage4fit,
  data = stage4,
  size = 1,       # change line size
  conf.int = F,  # Add confidence interval
  pval = TRUE, # Add p-value
  font.x = c(10), 
  font.y = c(10),
  font.tickslab = c(10),
  break.time.by = 1,
  xlim = c(0,5),
  xlab = "Time in years",
  ylab = "Progression-free survival",
  surv.scale="percent",
  palette = c("#FC4E07","#2E9FDF"),
  legend  = c(0.15,0.4),
  legend.title = "quartile",
  legend.labs = c("High", "Low"),
  conf.int.style = "step",
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.20, # Useful to change when you have multiple groups
  risk.table.title="Numbers of risk",
  ggtheme = theme_bw(),# Change ggplot2 theme
  tables.theme = theme_cleantable()
) +
  ggtitle("CESC stageIV-HDAC6")

                       
                      
#stage -OS
fit <- coxph(Surv(OS.Time,OS)~stage_group+HDAC6,data = os_quan)
stage <- survfit(Surv(OS.Time,OS)~quan_group+strata(os_quan$stage_group),data = os_quan)
ggsurvplot(stage,pval = TRUE)
stage1 <- os_quan %>% dplyr::filter(stage_group == 1)
stage2 <- os_quan %>% dplyr::filter(stage_group == 2)
stage3 <- os_quan %>% dplyr::filter(stage_group == 3)
stage4 <- os_quan %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage1)
stage2fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage2)
stage3fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage3)
stage4fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage4)
