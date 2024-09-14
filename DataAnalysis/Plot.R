library(ggplot2)
library(tidyverse)
data.decay <- read.csv("DataAnalysis/Data/Cleaned_Data/decay2.csv")
data.GLMM <- read.csv("DataAnalysis/Data/Cleaned_Data/GLMM2.csv")
data.decay <- data.decay %>% 
  mutate(censor = ifelse(log10 == min(log10), 1, 0)) %>% 
  filter(! ((PATIENT == 26 & time > 13) | (PATIENT == 71 & time > 22.3)))

data.10days <- data.decay %>% 
  filter(time <= 10)
p1 <- ggplot(data.10days, aes(x=time, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=2, stroke=0) +
  geom_line(aes(group=PATIENT)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(values=c(16, 17),labels = c("No", "Yes"))+
  labs(shape = "Censored") + 
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

PATID <- unique(data.decay$PATIENT) 
subdata <- data.decay %>% 
  filter(PATIENT %in% PATID[c(1, 2, 3, 4, 5)])
p2 <- ggplot(subdata, aes(x=time, y=log10)) + 
  geom_point(aes(shape=factor(censor)),size=2, stroke=0) +
  geom_line(aes(group=PATIENT)) +
  scale_x_continuous("Day") + 
  scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
  scale_shape_manual(values=c(16, 17),labels = c("No", "Yes"))+
  labs(shape = "Censored") + 
  theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))

ggarrange(p1, p2, common.legend = TRUE, legend = "right")
