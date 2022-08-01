# raw data are included in the Source Data file sheet "fig6"

library(ggpubr)
library(BuenColors)

FKBP5 <- read.csv("path/to/FKBP5.csv")
TEAD3 <- read.csv("path/to/TEAD3.csv")
SPRK1 <- read.csv("path/to/SPRK1.csv")
PPARD <- read.csv("path/to/PPARD.csv")
MAPK14 <- read.csv("path/to/MAPK14.csv")

FKBP5_fig <- ggbarplot(FKBP5, x = "Treatment", y = "Expression", 
                         add = c("mean_sd", "jitter"),add.params = list(shape = "Subtype1"), color = "Treatment", palette = c("#00AFBB","#bb2600","#bb8d00"),
                        order = c("NT", "sgEnh1","sgEnh2"),
                         position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0)) #350x450

TEAD3_fig <- ggbarplot(TEAD3, x = "Treatment", y = "Expression", 
                       add = c("mean_sd", "jitter"),add.params = list(shape = "Subtype1"), color = "Treatment", palette = c("#00AFBB","#bb2600","#bb8d00"),
                       order = c("NT", "sgEnh1","sgEnh2"),
                       position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0)) 

SPRK1_fig <- ggbarplot(SPRK1, x = "Treatment", y = "Expression", 
                       add = c("mean_sd", "jitter"),add.params = list(shape = "Subtype1"), color = "Treatment", palette = c("#00AFBB","#bb2600","#bb8d00"),
                       order = c("NT", "sgEnh1","sgEnh2"),
                       position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0)) 

PPARD_fig <- ggbarplot(PPARD, x = "Treatment", y = "Expression", 
                       add = c("mean_sd", "jitter"),add.params = list(shape = "Subtype1"), color = "Treatment", palette = c("#00AFBB","#bb2600","#bb8d00"),
                       order = c("NT", "sgEnh1","sgEnh2"),
                       position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0)) 

MAPK14_fig <- ggbarplot(MAPK14, x = "Treatment", y = "Expression", 
                       add = c("mean_sd", "jitter"),add.params = list(shape = "Subtype1"), color = "Treatment", palette = c("#00AFBB","#bb2600","#bb8d00"),
                       order = c("NT", "sgEnh1","sgEnh2"),
                       position = position_dodge(0.8),ylim = c(0, 1.5))+scale_y_continuous(expand = c(0, 0)) 


library(multcomp)

FKBP5$Treatment = as.factor(FKBP5$Treatment)

res.aov <- aov(Expression ~ Treatment, data =FKBP5)
summary(res.aov)

# Summary of the analysis
glht.res <- glht(res.aov, linfct = mcp(Treatment = "Dunnett"))
summary(glht.res)
