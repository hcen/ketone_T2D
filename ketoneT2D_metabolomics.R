#install.packages("tidyverse")
library(tidyverse)

#install.packages("readxl")
library(readxl)

#install.packages("devtools") # most recent version of complexheatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) 
#install.packages("gridtext")
library(gridtext)
library(scales)

#install.packages('matrixStats')
library(matrixStats) # for rowMedians

#install.packages("ComplexUpset")
library(ComplexUpset)

library(ggbeeswarm)

library(ggrepel) # for label text in ggplot2 # https://ggrepel.slowkow.com/articles/examples.html

library(ggfortify) # for PCA plot (Not used?)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ##Set working directory to where this file is.
getwd()

#dir.create("input")
#dir.create("output")
#dir.create("figures")

# load data ========
plasma1 <- read_excel("input/AcuteKetoneT2D_metabolomics_mastersheet.xlsx",sheet = 1) 
plasma2 <- read_excel("input/AcuteKetoneT2D_metabolomics_mastersheet.xlsx",sheet = 3) 
View(plasma1)

pbmc <- read_excel("input/AcuteKetoneT2D_metabolomics_mastersheet.xlsx",sheet = 2)
traits <- read_excel("input/Participants characteristics_Acute ketone T2D_metabolomics_final.xlsx",sheet = 1) 
View(traits)
traits <- traits[-c(15:16)]

View(pbmc)
View(plasma1)
View(plasma2)
plasma1$batch <- "one"
plasma2$batch <- "two"

# combine and format plasma data =====

plasma1 <- as.data.frame(t(plasma1))
plasma2 <- as.data.frame(t(plasma2))

# identify uncommon metabolites in batch one and two
batch.one <- rownames(plasma1)[(rownames(plasma1) %in% rownames(plasma2))==FALSE]
batch.two <- rownames(plasma2)[(rownames(plasma2) %in% rownames(plasma1))==FALSE]

max.len = max(length(batch.one), length(batch.two))
batch.one = c(batch.one, rep(NA, max.len - length(batch.one)))
batch.two = c(batch.two, rep(NA, max.len - length(batch.two)))
plasma.diff <- data.frame(batch.one,batch.two)
View(plasma.diff)

write.csv(plasma.diff,"output/plasma_batch_unique.csv", row.names = F)

# manually matched some metabolites, load and replace some names

plasma.match <- read.csv("output/plasma_batch_unique_match.csv",na.strings=c("","NA"))
View(plasma.match)

plasma.match <- na.omit(plasma.match)

View(plasma1)
df1 <- plasma1[plasma.match$batch.one,] 
df2 <- plasma1[c(plasma.match$batch.one),] 
df2 <- plasma1[(rownames(plasma1) %in% plasma.match$batch.one)==FALSE,] 
dim(df2)
dim(df1)
dim(plasma1)

rownames(df1) <- plasma.match$batch.two
rownames(plasma1[plasma.match$batch.one,])

plasma1 <- rbind(df2,df1)
dim(plasma1) # 48 columns

#
View(plasma2)
dim(plasma2) # 54 columns

plasma <- inner_join(rownames_to_column(plasma1),rownames_to_column(plasma2), by="rowname") %>%
  column_to_rownames()
View(plasma)
dim(plasma) #102 columns

plasma <- t(plasma) %>% as.data.frame()
plasma$Time <- gsub(" ", "", plasma$Time,fixed = TRUE)
plasma <- plasma %>% mutate(sample=paste(ID, Time, Condition,sep = "."))

rownames(plasma) <- plasma$sample
View(plasma)
dim(plasma) # After matching some names, it increased from 186 to 193 rows (with 5 rows of meta data, so 188 common metabolites)




# format traits data ============

plasma$ID <- as.character(plasma$ID)
traits$PARTICIPANT <- as.character(traits$PARTICIPANT)
traits <- plasma[,c("sample","ID","Time","Condition","batch")] %>%
  left_join(traits,by=c("ID"="PARTICIPANT"))

#traits$Condition <- gsub("A", "ketone",traits$Condition,fixed = TRUE)
#traits$Condition <- gsub("B", "placebo",traits$Condition,fixed = TRUE)
View(traits)

#traits.pair <- filter(traits, !ID %in% c(3,8))
#View(traits.pair)
#write.csv(traits.pair, file = "output/traits_pair.csv", row.names = F)


# (cont.) format plasma data ============

plasma.t <- t(select(plasma,!c("ID","Time","Condition","batch","sample"))) %>% as.data.frame()
View(plasma.t)
dim(plasma.t) # 188 metabolites
plasma.df <- as.data.frame(lapply(plasma.t, as.numeric),check.names=F)
rownames(plasma.df) <- rownames(plasma.t)
str(plasma.df)
View(plasma.df)

plasma.df <- plasma.df[rowSums(plasma.df==0) < (ncol(plasma.df)*0.5), ] # remove metabolites with 0 in >=50% samples
plasma.df <- plasma.df[rowSums(is.na(plasma.df)) < (ncol(plasma.df)*0.5), ] # remove metabolites with NA in >= 50% samples
plasma.df <- plasma.df[, colSums(is.na(plasma.df)) != (nrow(plasma.df))] # remove missing samples (all NA)
dim(plasma.df) # 184 metabolites, 98 samples
dim(plasma.t) # 188 metabolites, 102 samples

write.csv(plasma.df, file = "output/plasma.csv")

#plasma.pair <- filter(plasma, !ID %in% c(3,8))
#View(plasma.pair)
#plasma.pair.df <- t(select(plasma.pair,!c("ID","Time","Condition","batch","sample"))) %>% as.data.frame()
#View(plasma.pair.df)

#write.csv(plasma.pair.df, file = "output/plasma.pair.csv")



# (cont.) format traits data ============
colnames(plasma.df)
traits$sample %in% colnames(plasma.df)
traits <- traits[traits$sample %in% colnames(plasma.df),]

dim(traits) # 98 x 19
View(traits)
write.csv(traits, file = "output/traits.csv", row.names = F)

# PCA ===========
#df <- as.data.frame(t(xcms[,-c(1:23)]))
df <- as.data.frame(t(plasma.df)) # each row should be a sample
View(df)
str(df)
dim(df) # 98 x 184
#grep("Pre",row.names(df))
#grep("Post",row.names(df))
#df$group <- c(rep("Pre",19),rep("Post",19))
#dim(df)

df.pca <- prcomp(df, center = TRUE,scale. = TRUE) # each row should be a sample

summary(df.pca)
str(df.pca)
View(summary(df.pca))

pc <- as.data.frame(summary(df.pca)[["importance"]])
View(pc)
pc
pc1.v <- round(pc["Proportion of Variance","PC1"],digits=3)*100 # variance PC1
pc2.v <- round(pc["Proportion of Variance","PC2"],digits=3)*100 # variance PC2

cbPalette <- c(
               "#E69F00", #lightorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#999999", #grey
               "#F0E442" #yellow
               )
col <- cbPalette

PCi<-data.frame(df.pca$x,batch=traits$batch,time=traits$Time,condition=traits$Condition)
PCi$batch <- factor(PCi$batch, levels = c("one","two"))
PCi$time <- factor(PCi$time, levels = c("0","90","180"))
PCi$condition <- factor(PCi$condition, levels = c("ketone","placebo"))
View(PCi)

#autoplot(xcms.pca,data = df,colour = 'group')+

ggplot(PCi,aes(x=PC1,y=PC2,
               #col=batch
               col=condition,shape=time
               ))+
  geom_text_repel(aes(label=row.names(PCi)),size=2,
                  color="grey50",
                  box.padding   = 0.4,
                  point.padding = 0,
                  #force=1,
                  #force_pull=10,
                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
                  segment.color = 'darkgrey')+
  geom_point(size=2)+
  scale_color_manual(values=col) +
  labs(x = paste0("PC1 (",pc1.v,"%)"), y = paste0("PC1 (",pc2.v,"%)"))+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        #legend.title = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16)
        
  )+
  theme(aspect.ratio=1/1)

ggsave(filename="figures/PCA_batch.png",width=15,height=12,units="cm",dpi=400)

ggsave(filename="figures/PCA_time_condition.png",width=15,height=12,units="cm",dpi=400)


# cluster
sampleDists <- dist(df,  method = "euclidean") # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters=hclust(sampleDists)
plot(clusters)


# overall heatmap ===============

m <- as.matrix(as.data.frame(lapply(plasma.df, as.numeric),check.names=F))

m.z <- t(scale(t(m))) #%>% as.data.frame()
View(m.z)
colnames(m.z)

m.z[m.z>5] <- NA
m.z <- t(scale(t(m.z)))
m.z[is.na(m.z)] <- 5

ceiling(max(abs(m.z)))

colnames(m.z)
#summary(m.z)

number_of_batch1 <- length(grep("one", traits$batch))
number_of_batch2 <- length(grep("two", traits$batch))

end_index_batch1 <- grep("one", traits$batch)[number_of_batch1]
end_index_batch2 <- grep("two", traits$batch)[number_of_batch2]

number_of_batch1
end_index_batch1
number_of_batch2
end_index_batch2

###


heatmap.all <- Heatmap(m.z, #matrix_Z_score_total,
                       name = "Z score",
                       show_row_names = FALSE,
                       show_column_names = TRUE,
                       show_row_dend = TRUE,
                       #row_labels = gt_render(rownames(plasma.df)),
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       column_names_side = "bottom",
                       row_dend_side = "left",
                       column_dend_side = "bottom",
                  #     layer_fun = function(j, i, x, y, width, height, fill) {
                  #       mat = restore_matrix(j, i, x, y)
                  #       ind = unique(c(mat[, c(end_index_pre, 
                  #                              end_index_post
                  #       )]))
                  #       grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                  #                 y = y[ind], 
                  #                 width = unit(0.03, "inches"), 
                  #                 height = unit(1/nrow(m.z), "npc"),
                  #                 gp = gpar(col = "white")
                  #       )
                  #     },
                       #col = colorRamp2(c(-ceiling(max(abs(m.z))), 0, ceiling(max(abs(m.z)))), c("blue", "white", "red")),
                       col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
                       
                       top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                            , height = unit(12, "mm")
                       )),
                       #   left_annotation = 
                       #     rowAnnotation(`Significant CW vs CS` = m.anno$CS_CW,
                       
                       #                   `Significant KW vs KS` = m.anno$KS_KW,
                       
                       #                   `Significant CW vs KW` = m.anno$KW_CW,
                       
                       #                   `Significant CS vs KS` = m.anno$KS_CS,
                       
                       #                   col = list(`Significant CW vs CS` = c("TRUE" = hue_pal()(4)[1]),
                       #                              `Significant KW vs KS` = c("TRUE" = hue_pal()(4)[2]),
                       #                              `Significant CW vs KW` = c("TRUE" = hue_pal()(4)[3]),
                       #                              `Significant CS vs KS` = c("TRUE" = hue_pal()(4)[4])
                       #                   ),
                       #                   annotation_legend_param = list(
                       
                       #                     `Significant CW vs CS` = list(
                       #                       title = c(""),
                       #                       labels = expression(italic("Ins1")^italic("+/+")*" Sucrose vs Water") #expression() issue - https://github.com/jokergoo/ComplexHeatmap/issues/678
                       #                     ),
                       
                       #                     `Significant KW vs KS` = list(
                       #                       title = c(""),
                       #                       labels = expression(italic("Ins1")^italic("+/-")*" Sucrose vs Water")
                       #                     ),
                       
                       #                     `Significant CW vs KW` = list(
                       #                       title = c(""),
                       #                       labels = expression("Water "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
                       #                     ),
                       
                       #                     `Significant CS vs KS` = list(
                       #                       title = c(""),
                       #                       labels = expression("Sucrose "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
                       #                     )
                       #                   ),
                       #                   show_annotation_name = FALSE,
                       #                   show_legend=FALSE,
                       #                   na_col = "white"
                       
                       #     ),
                       column_order = 1:ncol(m.z),
                       height = 
                         
                         unit(150, "mm"), # multi AA
                       
                       width = ncol(m.z)*unit(2.5, "mm"),
                       border_gp = gpar(col = "black"),
                       show_heatmap_legend = TRUE,
                       heatmap_legend_param = list(
                         title = "Z-score",
                         title_position = "topleft",
                         legend_height = unit(4, "cm")))

draw(heatmap.all)
#

png(file = "figures/heatmap_all_batch.png",
     width = 11, 
     height = 7, 
     units = "in", res = 600)



draw(heatmap.all)
#

seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))


####Condition labels

#Condition label 1


grid.rect(x = (loc2$x - loc1$x)*(end_index_batch1)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_batch1)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = cbPalette[1],
                    col = cbPalette[1]
          )
)
grid.text("batch one", 
          x = (loc2$x - loc1$x)*(end_index_batch1)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2

grid.rect(x = (loc2$x - loc1$x)*(end_index_batch1 + 
                                   end_index_batch2)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_batch2)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = cbPalette[2],
                    col = cbPalette[2]
          )
)

grid.text("batch two", 
          x = (loc2$x - loc1$x)*(end_index_batch1 + 
                                   end_index_batch2)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


dev.off()

## end of overall heatmap =======



###Top label gaps


#Vertical lines
grid.rect(x = end_index_pre/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.03, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))


dev.off()

## End of heatmap=====


# format PBMC metabolomics data ===============
pbmc <- pbmc %>% mutate(sample=paste(ID, Time, Condition,sep = "."))
rownames(pbmc) <- pbmc$sample
View(pbmc)
length(pbmc$sample)
nrow(pbmc)

pbmc.t <- t(select(pbmc,!c("ID","Time","Condition","sample"))) %>% as.data.frame()
View(pbmc.t)
colnames(pbmc.t) <- rownames(pbmc)
dim(pbmc.t) # 67 metabolites
str(pbmc.t)
pbmc.df <- pbmc.t %>% mutate_at(1:ncol(pbmc.t),as.numeric)
#pbmc.df <- as.data.frame(lapply(pbmc.t, as.numeric),check.names=F)
str(pbmc.df)
View(pbmc.df)

#pbmc.df <- pbmc.df[rowSums(pbmc.df==0) < (ncol(pbmc.df)*0.5), ] # remove metabolites with 0 in >=50% samples
pbmc.df <- pbmc.df[rowSums(is.na(pbmc.df)) < (ncol(pbmc.df)*0.5), ] # remove metabolites with NA in >= 50% samples
pbmc.df <- pbmc.df[, colSums(is.na(pbmc.df)) != (nrow(pbmc.df))] # remove missing samples (all NA)
dim(pbmc.df) # 66 metabolites, 45 samples
dim(pbmc.t) # 67 metabolites, 48 samples

write.csv(pbmc.df, file = "output/pbmc.csv")

# format traits.pbmc

ncol(pbmc.df)==nrow(traits) # false
all(colnames(pbmc.df) %in% traits$sample) # true

traits.pbmc <- traits[traits$sample %in% colnames(pbmc.df),]

dim(traits.pbmc) # 45 x 19
View(traits.pbmc)

all(colnames(pbmc.df) == traits.pbmc$sample) #True

write.csv(traits.pbmc, file = "output/traits.pbmc.csv", row.names = F)



# save plasma spreadsheets for MetaboAnalyst.ca ===============

plasma.all <- read.csv("output/plasma_upload.csv", check.names = F)
View(plasma.all)
colnames(plasma.all)
plasma.all[plasma.all==0] <- NA
plasma.drink.time <- plasma.all[!rowSums(is.na(plasma.all))>ncol(plasma.all)*0.4,!colSums(is.na(plasma.all))>nrow(plasma.all)*0.4]
View(plasma.drink.time)
write.csv(plasma.drink.time,"output/plasma_drink_time.csv",row.names = F)

plasma.ketone <- plasma.all %>% select(all_of(grep("\\.A|sample",colnames(plasma.all))))
colnames(plasma.ketone)

plasma.0min <- plasma.all %>% 
  select(all_of(grep("\\.0|sample",colnames(plasma.all))))
colnames(plasma.0min)
plasma.90min <- plasma.all %>% 
  select(all_of(grep("\\.90|sample",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.90.A"))) # 2.90.B is a missing sample, so removed subject 2
colnames(plasma.90min)
plasma.180min <- plasma.all %>% select(all_of(grep("\\.180|sample",colnames(plasma.all)))) %>%
  select(-all_of(c("4.180.B","4.180.A","6.180.B","6.180.A"))) # 4.180.B, 6.180.A, 6.180.B are missing samples, so removed subject 4 & 6
colnames(plasma.180min)

plasma.drink.0.90min <- plasma.all %>% select(all_of(grep("\\.0|\\.90|sample",colnames(plasma.all))))
colnames(plasma.drink.0.90min)
plasma.drink.0.180min <- plasma.all %>% select(all_of(grep("\\.0|\\.180|sample",colnames(plasma.all))))
colnames(plasma.drink.0.180min)

plasma.ketone.0.90 <- plasma.all %>% 
  select(all_of(grep("sample|\\.90.A|\\.0.A",colnames(plasma.all)))) 
colnames(plasma.ketone.0.90)
plasma.ketone.0.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.0.A",colnames(plasma.all)))) %>%
  select(-all_of(c("6.180.A","6.0.A"))) # 6.180.A is a missing sample, so removed subject 6
colnames(plasma.ketone.0.180)
plasma.ketone.90.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.90.A",colnames(plasma.all)))) %>%
  select(-all_of(c("6.180.A","6.90.A"))) # 6.180.A is a missing sample, so removed subject 6
colnames(plasma.ketone.90.180)

plasma.ctrl.0.90 <- plasma.all %>% 
  select(all_of(grep("sample|\\.90.B|\\.0.B",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.0.B"))) # 2.90.B is a missing sample, so removed subject 2
colnames(plasma.ctrl.0.90)
plasma.ctrl.0.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.0.B",colnames(plasma.all)))) %>%
  select(-all_of(c("4.180.B","4.0.B","6.180.B","6.0.B"))) # 4.180.B & 6.180.B are missing samples, so removed subject 4 & 6
colnames(plasma.ctrl.0.180)
plasma.ctrl.90.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.90.B",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.180.B","4.180.B","4.90.B","6.180.B","6.90.B"))) # 2.90.B, 4.180.B & 6.180.B are missing samples, so removed subject 2, 4 & 6
colnames(plasma.ctrl.90.180)


meta.all <- read.csv("output/traits_upload.csv", check.names = F)
View(meta.all)
colnames(meta.all)
meta.all$Time <- paste0("T",meta.all$Time)
meta.all$Subject <- paste0("S",meta.all$Subject)
meta.drink.time <- meta.all %>% 
  filter(sample %in% colnames(plasma.drink.time))
write.csv(meta.all.save,"output/traits_drink_time.csv",row.names = F)

meta.ketone <- meta.all %>% filter(sample %in% colnames(plasma.ketone))
meta.ketone$sample

meta.0min <- meta.all %>% filter(sample %in% colnames(plasma.0min))
meta.0min$sample
meta.90min <- meta.all %>% filter(sample %in% colnames(plasma.90min))
meta.90min$sample
meta.180min <- meta.all %>% filter(sample %in% colnames(plasma.180min))
meta.180min$sample

meta.drink.0.90min <- meta.all %>% filter(sample %in% colnames(plasma.drink.0.90min))
meta.drink.0.90min$sample
meta.drink.0.180min <- meta.all %>% filter(sample %in% colnames(plasma.drink.0.180min))
meta.drink.0.180min$sample

meta.ketone.0.90 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.0.90))
meta.ketone.0.90$sample
meta.ketone.0.180 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.0.180))
meta.ketone.0.180$sample
meta.ketone.90.180 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.90.180))
meta.ketone.90.180$sample

meta.ctrl.0.90 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.0.90))
meta.ctrl.0.90$sample
meta.ctrl.0.180 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.0.180))
meta.ctrl.0.180$sample
meta.ctrl.90.180 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.90.180))
meta.ctrl.90.180$sample

write.csv(plasma.ketone,"output/plasma_ketone.csv",row.names = F)

write.csv(plasma.0min,"output/plasma_0min.csv",row.names = F)
write.csv(plasma.90min,"output/plasma_90min.csv",row.names = F)
write.csv(plasma.180min,"output/plasma_180min.csv",row.names = F)

write.csv(plasma.drink.0.90min,"output/plasma_drink_0min.90min.csv",row.names = F)
write.csv(plasma.drink.0.180min,"output/plasma_drink_0min.180min.csv",row.names = F)

write.csv(plasma.ketone.0.90,"output/plasma_ketone_0min.90min.csv",row.names = F)
write.csv(plasma.ketone.0.180,"output/plasma_ketone_0min.180min.csv",row.names = F)
write.csv(plasma.ketone.90.180,"output/plasma_ketone_90min.180min.csv",row.names = F)

write.csv(plasma.ctrl.0.90,"output/plasma_ctrl_0min.90min.csv",row.names = F)
write.csv(plasma.ctrl.0.180,"output/plasma_ctrl_0min.180min.csv",row.names = F)
write.csv(plasma.ctrl.90.180,"output/plasma_ctrl_90min.180min.csv",row.names = F)

#
write.csv(meta.ketone,"output/traits_ketone.csv",row.names = F)

write.csv(meta.0min,"output/traits_0min.csv",row.names = F)
write.csv(meta.90min,"output/traits_90min.csv",row.names = F)
write.csv(meta.180min,"output/traits_180min.csv",row.names = F)

write.csv(meta.drink.0.90min,"output/traits_drink_0min.90min.csv",row.names = F)
write.csv(meta.drink.0.180min,"output/traits_drink_0min.180min.csv",row.names = F)

write.csv(meta.ketone.0.90,"output/traits_ketone_0min.90min.csv",row.names = F)
write.csv(meta.ketone.0.180,"output/traits_ketone_0min.180min.csv",row.names = F)
write.csv(meta.ketone.90.180,"output/traits_ketone_90min.180min.csv",row.names = F)

write.csv(meta.ctrl.0.90,"output/traits_ctrl_0min.90min.csv",row.names = F)
write.csv(meta.ctrl.0.180,"output/traits_ctrl_0min.180min.csv",row.names = F)
write.csv(meta.ctrl.90.180,"output/traits_ctrl_90min.180min.csv",row.names = F)

#

# save pbmc spreadsheets for MetaboAnalyst.ca ===============

pbmc.all <- read.csv("output/pbmc.csv", check.names = F, row.names = 1)
View(pbmc.all)
colnames(pbmc.all)
pbmc.all[pbmc.all==0] <- NA
pbmc.all <- pbmc.all[!rowSums(is.na(pbmc.all))>ncol(pbmc.all)*0.4,!colSums(is.na(pbmc.all))>nrow(pbmc.all)*0.4] %>%
  rownames_to_column(var = "sample")
View(pbmc.all)
write.csv(pbmc.all,"output/pbmc_drink_time.csv",row.names = F)

pbmc.ketone <- pbmc.all %>% select(all_of(grep("\\.A|sample",colnames(pbmc.all))))
colnames(pbmc.ketone)

pbmc.0min <- pbmc.all %>% 
  select(all_of(grep("\\.0|sample",colnames(pbmc.all)))) %>%
  select(-all_of(c("3.0.B","8.0.A")))  # 3.0.A, 8.0.B is missing, so removed subject 3, 8
colnames(pbmc.0min)
pbmc.90min <- pbmc.all %>% 
  select(all_of(grep("\\.90|sample",colnames(pbmc.all)))) %>%
  select(-all_of(c("2.90.A","3.90.B","8.90.A"))) # 2.90.B, 3.90.A, 8.90.B is a missing sample, so removed subject 2,3,8
colnames(pbmc.90min)
pbmc.180min <- pbmc.all %>% select(all_of(grep("\\.180|sample",colnames(pbmc.all)))) %>%
  select(-all_of(c("3.180.B","5.180.B","6.180.B","8.180.A"))) # 3.180.A, 5.180.A, 6.180.A, 8.180.B are missing samples, so removed subject 3,5,6,8
colnames(pbmc.180min)

pbmc.drink.0.90min <- pbmc.all %>% select(all_of(grep("\\.0|\\.90|sample",colnames(pbmc.all))))
colnames(pbmc.drink.0.90min)
pbmc.drink.0.180min <- pbmc.all %>% select(all_of(grep("\\.0|\\.180|sample",colnames(pbmc.all))))
colnames(pbmc.drink.0.180min)

pbmc.ketone.0.90 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.90.A|\\.0.A",colnames(pbmc.all)))) 
colnames(pbmc.ketone.0.90)
pbmc.ketone.0.180 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.0.A",colnames(pbmc.all)))) %>%
  select(-all_of(c("5.0.A","6.0.A"))) # 5.180.A,6.180.A are missing samples, so removed subject 5,6
colnames(pbmc.ketone.0.180)
pbmc.ketone.90.180 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.90.A",colnames(pbmc.all)))) %>%
  select(-all_of(c("5.90.A","6.90.A"))) # 5.180.A,6.180.A are missing samples, so removed subject 5,6
colnames(pbmc.ketone.90.180)

pbmc.ctrl.0.90 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.90.B|\\.0.B",colnames(pbmc.all)))) %>%
  select(-all_of(c("2.0.B"))) # 2.90.B is a missing sample, so removed subject 2
colnames(pbmc.ctrl.0.90)
pbmc.ctrl.0.180 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.0.B",colnames(pbmc.all)))) #%>%
colnames(pbmc.ctrl.0.180)
pbmc.ctrl.90.180 <- pbmc.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.90.B",colnames(pbmc.all)))) %>%
  select(-all_of(c("2.180.B"))) # 2.90.B is a missing sample, so removed subject 2
colnames(pbmc.ctrl.90.180)


meta.all <- read.csv("output/traits.pbmc.csv", check.names = F)
View(meta.all)
colnames(meta.all)
meta.all$Time <- paste0("T",meta.all$Time)
meta.all$Subject <- paste0("S",meta.all$ID)
meta.all <- meta.all %>% 
  relocate(Subject,.before = (ID)) %>%
  select(-ID)
meta.drink.time <- meta.all %>%
  filter(sample %in% colnames(pbmc.drink.time))
View(meta.drink.time)
write.csv(meta.drink.time,"output/traits_pbmc_drink_time.csv",row.names = F)

meta.ketone <- meta.all %>% filter(sample %in% colnames(pbmc.ketone))
meta.ketone$sample

meta.0min <- meta.all %>% filter(sample %in% colnames(pbmc.0min))
meta.0min$sample
meta.90min <- meta.all %>% filter(sample %in% colnames(pbmc.90min))
meta.90min$sample
meta.180min <- meta.all %>% filter(sample %in% colnames(pbmc.180min))
meta.180min$sample

meta.drink.0.90min <- meta.all %>% filter(sample %in% colnames(pbmc.drink.0.90min))
meta.drink.0.90min$sample
meta.drink.0.180min <- meta.all %>% filter(sample %in% colnames(pbmc.drink.0.180min))
meta.drink.0.180min$sample

meta.ketone.0.90 <- meta.all %>% filter(sample %in% colnames(pbmc.ketone.0.90))
meta.ketone.0.90$sample
meta.ketone.0.180 <- meta.all %>% filter(sample %in% colnames(pbmc.ketone.0.180))
meta.ketone.0.180$sample
meta.ketone.90.180 <- meta.all %>% filter(sample %in% colnames(pbmc.ketone.90.180))
meta.ketone.90.180$sample

meta.ctrl.0.90 <- meta.all %>% filter(sample %in% colnames(pbmc.ctrl.0.90))
meta.ctrl.0.90$sample
meta.ctrl.0.180 <- meta.all %>% filter(sample %in% colnames(pbmc.ctrl.0.180))
meta.ctrl.0.180$sample
meta.ctrl.90.180 <- meta.all %>% filter(sample %in% colnames(pbmc.ctrl.90.180))
meta.ctrl.90.180$sample

write.csv(pbmc.ketone,"output/pbmc_ketone.csv",row.names = F)

write.csv(pbmc.0min,"output/pbmc_0min.csv",row.names = F)
write.csv(pbmc.90min,"output/pbmc_90min.csv",row.names = F)
write.csv(pbmc.180min,"output/pbmc_180min.csv",row.names = F)

write.csv(pbmc.drink.0.90min,"output/pbmc_drink_0min.90min.csv",row.names = F)
write.csv(pbmc.drink.0.180min,"output/pbmc_drink_0min.180min.csv",row.names = F)

write.csv(pbmc.ketone.0.90,"output/pbmc_ketone_0min.90min.csv",row.names = F)
write.csv(pbmc.ketone.0.180,"output/pbmc_ketone_0min.180min.csv",row.names = F)
write.csv(pbmc.ketone.90.180,"output/pbmc_ketone_90min.180min.csv",row.names = F)

write.csv(pbmc.ctrl.0.90,"output/pbmc_ctrl_0min.90min.csv",row.names = F)
write.csv(pbmc.ctrl.0.180,"output/pbmc_ctrl_0min.180min.csv",row.names = F)
write.csv(pbmc.ctrl.90.180,"output/pbmc_ctrl_90min.180min.csv",row.names = F)

#
write.csv(meta.ketone,"output/traits_pbmc_ketone.csv",row.names = F)

write.csv(meta.0min,"output/traits_pbmc_0min.csv",row.names = F)
write.csv(meta.90min,"output/traits_pbmc_90min.csv",row.names = F)
write.csv(meta.180min,"output/traits_pbmc_180min.csv",row.names = F)

write.csv(meta.drink.0.90min,"output/traits_pbmc_drink_0min.90min.csv",row.names = F)
write.csv(meta.drink.0.180min,"output/traits_pbmc_drink_0min.180min.csv",row.names = F)

write.csv(meta.ketone.0.90,"output/traits_pbmc_ketone_0min.90min.csv",row.names = F)
write.csv(meta.ketone.0.180,"output/traits_pbmc_ketone_0min.180min.csv",row.names = F)
write.csv(meta.ketone.90.180,"output/traits_pbmc_ketone_90min.180min.csv",row.names = F)

write.csv(meta.ctrl.0.90,"output/traits_pbmc_ctrl_0min.90min.csv",row.names = F)
write.csv(meta.ctrl.0.180,"output/traits_pbmc_ctrl_0min.180min.csv",row.names = F)
write.csv(meta.ctrl.90.180,"output/traits_pbmc_ctrl_90min.180min.csv",row.names = F)

# organize MetaboAnalyst results ===================

ketone.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.90min_subject.csv",
         row.names = 1)
ketone.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.180min_subject.csv",
                        row.names = 1)
ketone.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_90min.180min_subject.csv",
                        row.names = 1)
ctrl.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.90min_subject.csv",
                        row.names = 1)
ctrl.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.180min_subject.csv",
                         row.names = 1)
ctrl.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_90min.180min_subject.csv",
                          row.names = 1)
ketone.ctrl.0 <- read.csv("input/MetaboAnalyst_results/covariate_result_0min_ctrl.ketone_subject.csv",
                             row.names = 1)
ketone.ctrl.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_90min_ctrl.ketone_subject.csv",
                             row.names = 1)
ketone.ctrl.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_180min_ctrl.ketone_subject.csv",
                             row.names = 1)

# organize ketone.0.90
plasma.ketone.0.90 <- read.csv("output/plasma_ketone_0min.90min.csv", check.names = F)
plasma.ketone.0.90.proc <- plasma.ketone.0.90[!rowSums(is.na(plasma.ketone.0.90))>
                                                  ncol(plasma.ketone.0.90)*0.4,
                                                !colSums(is.na(plasma.ketone.0.90))>
                                                  nrow(plasma.ketone.0.90)*0.4]
View(plasma.ketone.0.90)
colnames(ketone.0.90)
ketone.0.90.anno <- rownames_to_column(ketone.0.90) %>% 
  left_join(plasma.ketone.0.90.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.0.90.anno)
View(ketone.0.90.anno)
colnames(ketone.0.90.anno)
ketone.0.90.anno <- ketone.0.90.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ketone.0.90.anno[,grep("\\.0\\.",colnames(ketone.0.90.anno))],na.rm = T),
         mean.90min=rowMeans(ketone.0.90.anno[,grep("\\.90\\.",colnames(ketone.0.90.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ketone.0.90.anno[,grep("\\.0\\.",colnames(ketone.0.90.anno))]),na.rm = T),
         median.90min=rowMedians(as.matrix(ketone.0.90.anno[,grep("\\.90\\.",colnames(ketone.0.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.90min/mean.0min),
         logFC.median=log2(median.90min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ketone.0.90.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.0.90.anno[is.na(ketone.0.90.anno)] <- 0

colnames(ketone.0.90.anno)
ketone.0.90.up <- ketone.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.0.90.up)

ketone.0.90.down <- ketone.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.0.90.down)

# organize ketone.0.180
plasma.ketone.0.180 <- read.csv("output/plasma_ketone_0min.180min.csv", check.names = F)
plasma.ketone.0.180.proc <- plasma.ketone.0.180[!rowSums(is.na(plasma.ketone.0.180))>
                                                  ncol(plasma.ketone.0.180)*0.4,
                                                !colSums(is.na(plasma.ketone.0.180))>
                                                  nrow(plasma.ketone.0.180)*0.4]
View(plasma.ketone.0.180)
colnames(ketone.0.180)
ketone.0.180.anno <- rownames_to_column(ketone.0.180) %>% 
  left_join(plasma.ketone.0.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.0.180.anno)
View(ketone.0.180.anno)
colnames(ketone.0.180.anno)
ketone.0.180.anno <- ketone.0.180.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ketone.0.180.anno[,grep("\\.0\\.",colnames(ketone.0.180.anno))],na.rm = T),
         mean.180min=rowMeans(ketone.0.180.anno[,grep("\\.180\\.",colnames(ketone.0.180.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ketone.0.180.anno[,grep("\\.0\\.",colnames(ketone.0.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ketone.0.180.anno[,grep("\\.180\\.",colnames(ketone.0.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.0min),
         logFC.median=log2(median.180min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ketone.0.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.0.180.anno[is.na(ketone.0.180.anno)] <- 0

colnames(ketone.0.180.anno)
ketone.0.180.up <- ketone.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.0.180.up)

ketone.0.180.down <- ketone.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.0.180.down)

# organize ketone.90.180
plasma.ketone.90.180 <- read.csv("output/plasma_ketone_90min.180min.csv", check.names = F)
View(plasma.ketone.90.180)
plasma.ketone.90.180.proc <- plasma.ketone.90.180[!rowSums(is.na(plasma.ketone.90.180))>
                                                    ncol(plasma.ketone.90.180)*0.4,
                                                  !colSums(is.na(plasma.ketone.90.180))>
                                                    nrow(plasma.ketone.90.180)*0.4]
View(plasma.ketone.90.180.proc)
colnames(ketone.90.180)
ketone.90.180.anno <- rownames_to_column(ketone.90.180) %>% 
  left_join(plasma.ketone.90.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.90.180.anno)
View(ketone.90.180.anno)
colnames(ketone.90.180.anno)
ketone.90.180.anno <- ketone.90.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.90min=rowMeans(ketone.90.180.anno[,grep("\\.90\\.",colnames(ketone.90.180.anno))],na.rm = T),
         mean.180min=rowMeans(ketone.90.180.anno[,grep("\\.180\\.",colnames(ketone.90.180.anno))],na.rm = T),
         median.90min=rowMedians(as.matrix(ketone.90.180.anno[,grep("\\.90\\.",colnames(ketone.90.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ketone.90.180.anno[,grep("\\.180\\.",colnames(ketone.90.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.90min),
         logFC.median=log2(median.180min/median.90min)) %>%
  relocate(mean.90min:logFC.median,.after = logFC)
View(ketone.90.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.90.180.anno[is.na(ketone.90.180.anno)] <- 0

colnames(ketone.90.180.anno)
ketone.90.180.up <- ketone.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.90.180.up) # no metabolite

ketone.90.180.down <- ketone.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.90.180.down) # 1 metabolite


# organize ctrl.0.90
plasma.ctrl.0.90 <- read.csv("output/plasma_ctrl_0min.90min.csv", check.names = F)
plasma.ctrl.0.90.proc <- plasma.ctrl.0.90[!rowSums(is.na(plasma.ctrl.0.90))>
                                                ncol(plasma.ctrl.0.90)*0.4,
                                              !colSums(is.na(plasma.ctrl.0.90))>
                                                nrow(plasma.ctrl.0.90)*0.4]
View(plasma.ctrl.0.90)
colnames(ctrl.0.90)
ctrl.0.90.anno <- rownames_to_column(ctrl.0.90) %>% 
  left_join(plasma.ctrl.0.90.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.0.90.anno)
View(ctrl.0.90.anno)
colnames(ctrl.0.90.anno)
ctrl.0.90.anno <- ctrl.0.90.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ctrl.0.90.anno[,grep("\\.0\\.",colnames(ctrl.0.90.anno))],na.rm = T),
         mean.90min=rowMeans(ctrl.0.90.anno[,grep("\\.90\\.",colnames(ctrl.0.90.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ctrl.0.90.anno[,grep("\\.0\\.",colnames(ctrl.0.90.anno))]),na.rm = T),
         median.90min=rowMedians(as.matrix(ctrl.0.90.anno[,grep("\\.90\\.",colnames(ctrl.0.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.90min/mean.0min),
         logFC.median=log2(median.90min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ctrl.0.90.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.0.90.anno[is.na(ctrl.0.90.anno)] <- 0

colnames(ctrl.0.90.anno)
ctrl.0.90.up <- ctrl.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.0.90.up) # no metabolite

ctrl.0.90.down <- ctrl.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.0.90.down) # 5 metabolites

# organize ctrl.0.180
plasma.ctrl.0.180 <- read.csv("output/plasma_ctrl_0min.180min.csv", check.names = F)
plasma.ctrl.0.180.proc <- plasma.ctrl.0.180[!rowSums(is.na(plasma.ctrl.0.180))>
                                                  ncol(plasma.ctrl.0.180)*0.4,
                                                !colSums(is.na(plasma.ctrl.0.180))>
                                                  nrow(plasma.ctrl.0.180)*0.4]
View(plasma.ctrl.0.180)
colnames(ctrl.0.180)
ctrl.0.180.anno <- rownames_to_column(ctrl.0.180) %>% 
  left_join(plasma.ctrl.0.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.0.180.anno)
View(ctrl.0.180.anno)
colnames(ctrl.0.180.anno)
ctrl.0.180.anno <- ctrl.0.180.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ctrl.0.180.anno[,grep("\\.0\\.",colnames(ctrl.0.180.anno))],na.rm = T),
         mean.180min=rowMeans(ctrl.0.180.anno[,grep("\\.180\\.",colnames(ctrl.0.180.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ctrl.0.180.anno[,grep("\\.0\\.",colnames(ctrl.0.180.anno))],na.rm = T)),
         median.180min=rowMedians(as.matrix(ctrl.0.180.anno[,grep("\\.180\\.",colnames(ctrl.0.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.0min),
         logFC.median=log2(median.180min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ctrl.0.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.0.180.anno[is.na(ctrl.0.180.anno)] <- 0

colnames(ctrl.0.180.anno)
ctrl.0.180.up <- ctrl.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.0.180.up) # no metabolite

ctrl.0.180.down <- ctrl.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.0.180.down) # 6 metabolite

# organize ctrl.90.180

plasma.ctrl.90.180 <- read.csv("output/plasma_ctrl_90min.180min.csv", check.names = F)
View(plasma.ctrl.90.180)
plasma.ctrl.90.180.proc <- plasma.ctrl.90.180[!rowSums(is.na(plasma.ctrl.90.180))>
                                                    ncol(plasma.ctrl.90.180)*0.4,
                                                  !colSums(is.na(plasma.ctrl.90.180))>
                                                    nrow(plasma.ctrl.90.180)*0.4]
View(plasma.ctrl.90.180.proc)
colnames(ctrl.90.180)
ctrl.90.180.anno <- rownames_to_column(ctrl.90.180) %>% 
  left_join(plasma.ctrl.90.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.90.180.anno)
View(ctrl.90.180.anno)
colnames(ctrl.90.180.anno)
ctrl.90.180.anno <- ctrl.90.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.90min=rowMeans(ctrl.90.180.anno[,grep("\\.90\\.",colnames(ctrl.90.180.anno))],na.rm = T),
         mean.180min=rowMeans(ctrl.90.180.anno[,grep("\\.180\\.",colnames(ctrl.90.180.anno))],na.rm = T),
         median.90min=rowMedians(as.matrix(ctrl.90.180.anno[,grep("\\.90\\.",colnames(ctrl.90.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ctrl.90.180.anno[,grep("\\.180\\.",colnames(ctrl.90.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.90min),
         logFC.median=log2(median.180min/median.90min)) %>%
  relocate(mean.90min:logFC.median,.after = logFC)
View(ctrl.90.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.90.180.anno[is.na(ctrl.90.180.anno)] <- 0

colnames(ctrl.90.180.anno)
ctrl.90.180.up <- ctrl.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.90.180.up) # 1 beta-Hydroxybutyric acid

ctrl.90.180.down <- ctrl.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.90.180.down) # no metabolite


# organize ketone.ctrl.0
plasma.0min <- read.csv("output/plasma_0min.csv", check.names = F)
View(plasma.0min)
plasma.0min.proc <- plasma.0min[!rowSums(is.na(plasma.0min))>
                                                ncol(plasma.0min)*0.4,
                                              !colSums(is.na(plasma.0min))>
                                                nrow(plasma.0min)*0.4]
View(plasma.0min.proc)
colnames(ketone.ctrl.0)
ketone.ctrl.0.anno <- rownames_to_column(ketone.ctrl.0) %>% 
  left_join(plasma.0min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.0.anno)
View(ketone.ctrl.0.anno)
colnames(ketone.ctrl.0.anno)
ketone.ctrl.0.anno <- ketone.ctrl.0.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.0.anno[,grep("\\.B",colnames(ketone.ctrl.0.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.0.anno[,grep("\\.A",colnames(ketone.ctrl.0.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.0.anno[,grep("\\.B",colnames(ketone.ctrl.0.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.0.anno[,grep("\\.A",colnames(ketone.ctrl.0.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.0.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.0.anno[is.na(ketone.ctrl.0.anno)] <- 0

colnames(ketone.ctrl.0.anno)
ketone.ctrl.0.up <- ketone.ctrl.0.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.0.up) # no metabolite

ketone.ctrl.0.down <- ketone.ctrl.0.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.0.down) # no metabolite


# organize ketone.ctrl.90
plasma.90min <- read.csv("output/plasma_90min.csv", check.names = F)
View(plasma.90min)
plasma.90min.proc <- plasma.90min[!rowSums(is.na(plasma.90min))>
                                  ncol(plasma.90min)*0.4,
                                !colSums(is.na(plasma.90min))>
                                  nrow(plasma.90min)*0.4]
View(plasma.90min.proc)
colnames(ketone.ctrl.90)
ketone.ctrl.90.anno <- rownames_to_column(ketone.ctrl.90) %>% 
  left_join(plasma.90min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.90.anno)
View(ketone.ctrl.90.anno)
colnames(ketone.ctrl.90.anno)
ketone.ctrl.90.anno <- ketone.ctrl.90.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.90.anno[,grep("\\.B",colnames(ketone.ctrl.90.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.90.anno[,grep("\\.A",colnames(ketone.ctrl.90.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.90.anno[,grep("\\.B",colnames(ketone.ctrl.90.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.90.anno[,grep("\\.A",colnames(ketone.ctrl.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.90.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.90.anno[is.na(ketone.ctrl.90.anno)] <- 0

colnames(ketone.ctrl.90.anno)
ketone.ctrl.90.up <- ketone.ctrl.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.90.up) 

ketone.ctrl.90.down <- ketone.ctrl.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.90.down) 

# organize ketone.ctrl.180
plasma.180min <- read.csv("output/plasma_180min.csv", check.names = F)
View(plasma.180min)
plasma.180min.proc <- plasma.180min[!rowSums(is.na(plasma.180min))>
                                    ncol(plasma.180min)*0.4,
                                  !colSums(is.na(plasma.180min))>
                                    nrow(plasma.180min)*0.4]
View(plasma.180min.proc)
colnames(ketone.ctrl.180)
ketone.ctrl.180.anno <- rownames_to_column(ketone.ctrl.180) %>% 
  left_join(plasma.180min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.180.anno)
View(ketone.ctrl.180.anno)
colnames(ketone.ctrl.180.anno)
ketone.ctrl.180.anno <- ketone.ctrl.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.180.anno[,grep("\\.B",colnames(ketone.ctrl.180.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.180.anno[,grep("\\.A",colnames(ketone.ctrl.180.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.180.anno[,grep("\\.B",colnames(ketone.ctrl.180.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.180.anno[,grep("\\.A",colnames(ketone.ctrl.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.180.anno[is.na(ketone.ctrl.180.anno)] <- 0

colnames(ketone.ctrl.180.anno)
ketone.ctrl.180.up <- ketone.ctrl.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.180.up) 

ketone.ctrl.180.down <- ketone.ctrl.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.180.down) 

# plot altered metebolits in different comparisons onto UpSet plot ======================
sets1 <- ketone.0.90.up %>% mutate(ketone_0min.90min_up=1) %>%
  select(all_of(c("metabolite","ketone_0min.90min_up")))
sets2 <- ketone.0.90.down %>% mutate(ketone_0min.90min_down=1) %>%
  select(all_of(c("metabolite","ketone_0min.90min_down")))
sets3 <- ketone.0.180.up %>% mutate(ketone_0min.180min_up=1) %>%
  select(all_of(c("metabolite","ketone_0min.180min_up")))
sets4 <- ketone.0.180.down %>% mutate(ketone_0min.180min_down=1) %>%
  select(all_of(c("metabolite","ketone_0min.180min_down")))
sets5 <- ketone.90.180.up %>% mutate(ketone_90min.180min_up=1) %>%
  select(all_of(c("metabolite","ketone_90min.180min_up")))
sets6 <- ketone.90.180.down %>% mutate(ketone_90min.180min_down=1) %>%
  select(all_of(c("metabolite","ketone_90min.180min_down")))

sets7 <- ctrl.0.90.up %>% mutate(ctrl_0min.90min_up=1) %>%
  select(all_of(c("metabolite","ctrl_0min.90min_up")))
sets8 <- ctrl.0.90.down %>% mutate(ctrl_0min.90min_down=1) %>%
  select(all_of(c("metabolite","ctrl_0min.90min_down")))
sets9 <- ctrl.0.180.up %>% mutate(ctrl_0min.180min_up=1) %>%
  select(all_of(c("metabolite","ctrl_0min.180min_up")))
sets10 <- ctrl.0.180.down %>% mutate(ctrl_0min.180min_down=1) %>%
  select(all_of(c("metabolite","ctrl_0min.180min_down")))
sets11 <- ctrl.90.180.up %>% mutate(ctrl_90min.180min_up=1) %>%
  select(all_of(c("metabolite","ctrl_90min.180min_up")))
sets12 <- ctrl.90.180.down %>% mutate(ctrl_90min.180min_down=1) %>%
  select(all_of(c("metabolite","ctrl_90min.180min_down")))

sets13 <- ketone.ctrl.0.up %>% mutate(ketone.ctrl_0min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_0min_up")))
sets14 <- ketone.ctrl.0.down %>% mutate(ketone.ctrl_0min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_0min_down")))
sets15 <- ketone.ctrl.90.up %>% mutate(ketone.ctrl_90min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_90min_up")))
sets16 <- ketone.ctrl.90.down %>% mutate(ketone.ctrl_90min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_90min_down")))
sets17 <- ketone.ctrl.180.up %>% mutate(ketone.ctrl_180min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_180min_up")))
sets18 <- ketone.ctrl.180.down %>% mutate(ketone.ctrl_180min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_180min_down")))

# combine all sets
sets <- sets1 %>%
  full_join(sets2) %>%
  full_join(sets3) %>%
  full_join(sets4) %>%
  full_join(sets5) %>%
  full_join(sets6) %>%
  full_join(sets7) %>%
  full_join(sets8) %>%
  full_join(sets9) %>%
  full_join(sets10) %>%
  full_join(sets11) %>%
  full_join(sets12) %>%
  full_join(sets13) %>%
  full_join(sets14) %>%
  full_join(sets15) %>%
  full_join(sets16) %>%
  full_join(sets17) %>%
  full_join(sets18)
View(sets)
sets[,-1] = sets[,-1] == 1 # convert 1 to TRUE
View(sets)

# combine ketone 0, 90, 180 min
sets.ketone.time <- sets1 %>%
  full_join(sets2) %>%
  full_join(sets3) %>%
  full_join(sets4) %>%
  full_join(sets5) %>%
  full_join(sets6) 
sets.ketone.time[,-1] = sets.ketone.time[,-1] == 1 # convert 1 to TRUE
View(sets.ketone.time)
# combine ctrl 0, 90, 180
sets.ctrl.time <- sets7 %>%
  full_join(sets8) %>%
  full_join(sets9) %>%
  full_join(sets10) %>%
  full_join(sets11) %>%
  full_join(sets12) 
View(sets.ctrl.time)
sets.ctrl.time[,-1] = sets.ctrl.time[,-1] == 1 # convert 1 to TRUE
View(sets.ctrl.time)
#https://krassowski.github.io/complex-upset/articles/Examples_R.html

## select the data columns!!!!!!
data.sets <- colnames(sets)[-1]
data.sets <- colnames(sets.ctrl.time)[-1]
data.sets <- colnames(sets.ketone.time)[-1]

ComplexUpset::upset(

  ## select the data !!!!!!!!!!!!
                    sets,
                    #sets.ctrl.time,
                    #sets.ketone.time,
                    
                    data.sets, name='', 
                    sort_intersections_by='degree',
                    width_ratio=0.2,
                    height_ratio=0.6, 
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(text=element_text(size=12)),
                        'overall_sizes'=theme(axis.text.x=element_text(angle=90))
                      )),
                    
                    # min_size=500,
                    
                    base_annotations = list(
                      'Intersection size'=(
                        intersection_size(
                          
                          bar_number_threshold=1000000,
                          #text=list(vjust=1.1)
                          #counts=FALSE,
                          #mapping=aes(fill=general.type)
                        )
                        
                        + expand_limits(y=50)
                    #    + scale_fill_manual(values=cbPalette)
                        + theme(legend.title = element_blank(),
                                axis.text.y =element_blank(),
                                axis.title.y = element_text(margin = margin(r = -250)),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank())
                        
                      )),
                    set_sizes=(
                      upset_set_size()
                      + ylab('Set size')
                      + geom_text(aes(label=after_stat(count)), hjust=1, stat='count')
                      + expand_limits(y=80)
                      + theme(axis.title.x = element_text(margin = margin(t = -10)),
                              axis.text.x=element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())
                    ),
                    
                    stripes = c("white")) + 
  theme(panel.grid = element_blank(),
        text = element_text(size=14))+
  
  patchwork::plot_layout(heights=c(0.6, 1))

ggsave(filename="figures/upset.png",width=8,height=5,units="in",dpi=600)
ggsave(filename="figures/upset.pdf",width=8,height=5,units="in")
ggsave(filename="figures/upset_ctrl.pdf",width=8,height=5,units="in")
ggsave(filename="figures/upset_ketone.pdf",width=8,height=5,units="in")

# plot lines ketone down ==============
colnames(ketone.0.180.down)
colnames(ketone.0.90.down)
#setdiff(colnames(ketone.0.90.down),colnames(ketone.0.180.down))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.down),colnames(ketone.0.180.down)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ketone.0.90.180.down <- ketone.0.90.down[,1:13] %>%
  full_join(ketone.0.180.down[,1:13], suffix = c("_0.90","_0.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ketone.0.90.180.down)
colnames(ketone.0.90.180.down)

# normalized values of each metabolite to mean ketone 0 min

ketone.down.norm <- ketone.0.90.180.down[,c(1,7,19,26:123)] 
colnames(ketone.down.norm)
ketone.down.norm$mean.ketone.0min <- rowMeans(ketone.down.norm[,grep("\\.0.A",colnames(ketone.down.norm))],na.rm = T)
ketone.down.norm <- ketone.down.norm %>%
  mutate_at(vars(4:101),list(norm=~./mean.ketone.0min))
ketone.down.norm <- ketone.down.norm[,-(4:102)]
colnames(ketone.down.norm) <- gsub("_norm","",colnames(ketone.down.norm))
View(ketone.down.norm)
colnames(ketone.down.norm)

ketone.down.norm.long <- ketone.down.norm[,c(1:3,4:101)] %>%
  pivot_longer(-(1:3),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ketone.down.norm.long$drink<-factor(ketone.down.norm.long$drink,
                                    levels=c("Placebo","Ketone"))
levels(ketone.down.norm.long$drink)
View(ketone.down.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ketone.down.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)
p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
    geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  #geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
  #                 method = "quasirandom",##"pseudorandom",
                   #width = 3, #nbins = 10, #bandwidth = 1,#varwidth = TRUE, # different width based on dot density
  #                 dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
                                  color="black", size=14),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
                                   color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
                                   color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
                                   color="black",size=14), 
        legend.title = element_blank(),
        legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
                                  color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("ketone.down_placebo_ribon_lines.pdf",path="./figures/",width = 18, height = 10, units = "cm")
ggsave("ketone.down_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.down_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")


# plot lines ketone up ==============

colnames(ketone.0.180.up)
colnames(ketone.0.90.up)
colnames(ketone.90.180.up)
#setdiff(colnames(ketone.0.90.up),colnames(ketone.0.180.up))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.up),colnames(ketone.0.180.up)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ketone.0.90.180.up <- ketone.0.90.up[,1:13] %>%
  full_join(ketone.0.180.up[,1:13], suffix = c("","_0.180"),by="metabolite") %>%
  full_join(ketone.90.180.down[,1:13], suffix = c("_0.90","_90.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ketone.0.90.180.up)
colnames(ketone.0.90.180.up)

# normalized values of each metabolite to mean ketone 0 min

ketone.up.norm <- ketone.0.90.180.up[,c(1,7,19,31,38:135)] 
colnames(ketone.up.norm)
ketone.up.norm$mean.ketone.0min <- rowMeans(ketone.up.norm[,grep("\\.0.A",colnames(ketone.up.norm))],na.rm = T)
ketone.up.norm <- ketone.up.norm %>%
  mutate_at(vars(5:102),list(norm=~./mean.ketone.0min))
colnames(ketone.up.norm)
ketone.up.norm <- ketone.up.norm[,-(5:103)]
colnames(ketone.up.norm) <- gsub("_norm","",colnames(ketone.up.norm))
View(ketone.up.norm)
colnames(ketone.up.norm)

ketone.up.norm.long <- ketone.up.norm %>%
  pivot_longer(-(1:4),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ketone.up.norm.long$drink<-factor(ketone.up.norm.long$drink,
                                    levels=c("Placebo","Ketone"))
levels(ketone.up.norm.long$drink)
View(ketone.up.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ketone.up.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)

p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  #geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
                   method = "quasirandom",##"pseudorandom",
                   #width = 3,#nbins = 10,#bandwidth = 1, #varwidth = TRUE, # different width based on dot density
                   dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(
  
# (edit here!!!!!) two metabolite in ketone group is very high, so set limits to make other dots and lines more visible.     
    limits = c(0.5,1.5),
    breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("ketone.up_placebo_ribon_lines.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_lines_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_dots_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")

# plot lines ctrl ==============

colnames(ctrl.0.180.down)
colnames(ctrl.0.90.down)
colnames(ctrl.90.180.up)
#setdiff(colnames(ketone.0.90.up),colnames(ketone.0.180.up))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.up),colnames(ketone.0.180.up)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ctrl.0.90.180 <- ctrl.0.90.down[,1:13] %>%
  full_join(ctrl.0.180.down[,1:13], suffix = c("","_0.180"),by="metabolite") %>%
  full_join(ctrl.90.180.up[,1:13], suffix = c("_0.90","_90.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ctrl.0.90.180)
colnames(ctrl.0.90.180)

# normalized values of each metabolite to mean ketone 0 min

ctrl.norm <- ctrl.0.90.180[,c(1,7,19,31,38:135)] 
colnames(ctrl.norm)
ctrl.norm$mean.ketone.0min <- rowMeans(ctrl.norm[,grep("\\.0.A",colnames(ctrl.norm))],na.rm = T)
ctrl.norm <- ctrl.norm %>%
  mutate_at(vars(5:102),list(norm=~./mean.ketone.0min))
colnames(ctrl.norm)
ctrl.norm <- ctrl.norm[,-(5:103)]
colnames(ctrl.norm) <- gsub("_norm","",colnames(ctrl.norm))
View(ctrl.norm)
colnames(ctrl.norm)

ctrl.norm.long <- ctrl.norm %>%
  pivot_longer(-(1:4),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ctrl.norm.long$drink<-factor(ctrl.norm.long$drink,
                                  levels=c("Placebo","Ketone"))
levels(ctrl.norm.long$drink)
View(ctrl.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ctrl.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)

p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  #geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
  #                 method = "quasirandom",##"pseudorandom",
  #                 #width = 3,#nbins = 10,#bandwidth = 1, #varwidth = TRUE, # different width based on dot density
  #                 dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(

# (edit here!!!!!) one metabolite in ketone group is very high, so set limits to make other dots and lines more visible.     
    limits = c(0.6,1.5),
    breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("placebo.updown_ketone_ribon_lines_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("placebo.updown_ketone_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("placebo.updown_ketone_ribon_dots_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")

# plot all metabolites trends ===========
colnames(ketone.down.norm.long)

####### select data to plot !!!!!!!!!!
dat <- ketone.down.norm.long %>%
  arrange(logFC.mean_0.180,logFC.mean_0.90)
dat <- ketone.up.norm.long %>%
  arrange(-logFC.mean_0.180,-logFC.mean_0.90)
dat <- ctrl.norm.long %>%
  arrange(-logFC.mean_0.180,-logFC.mean_0.90)

dat$metabolite <- factor(dat$metabolite,levels=unique(dat$metabolite))
dat$drink<-factor(dat$drink,levels=c("Placebo","Ketone"))
View(dat)
length(unique(dat$metabolite))

df <- dat[which(dat$Value==0),]
View(df)

datTab <- dat %>%
  group_by(metabolite,time,drink)  %>%
  summarise(UQ=quantile(Value, 0.75,na.rm=T),LQ=quantile(Value,0.25,na.rm=T),
            Median=median(Value,na.rm=T),Mean=mean(Value,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)
p<-ggplot(dat,aes(x=time,y=Value,group=drink)) +
  geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(Subject,drink),colour=drink),
            linetype="solid",linewidth=0.5,alpha=0.4 )+ #
#  geom_quasirandom(aes(group=interaction(Subject,drink),colour=drink),
#                   method = "quasirandom",##"pseudorandom", #width = 3, #nbins = 10, #bandwidth = 1,#varwidth = TRUE, # different width based on dot density
#                   dodge.width=10,cex=0.8, alpha=0.3)+
  
  ##geom_beeswarm(aes(group=interaction(Subject,drink),
  ##                  colour=drink),
  ##              data=dat[dat$drink=="Ketone",],
  ##              side=1L,
  ##              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,linewidth=1.2,alpha=1)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "none") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(
    facets=~metabolite,
    
####### change # of rows and columes here !!!!!!!!!!!!!
    #nrow = 8, ncol = 8, # for ketone down
    #nrow = 1, ncol = 6, # for ketone up
    nrow = 1, ncol = 8, # for placebo
    scales = "free",
    shrink = TRUE,
    labeller = "label_value",
    as.table = TRUE,
    #switch = deprecated(),
    drop = TRUE,
    dir = "h",
    strip.position = "top"
  )


pdf("figures/ketone.down_placebo_ribon_allMetabolites_lines.pdf",width = 48,height=40)
pdf("figures/ketone.down_placebo_ribon_allMetabolites.pdf",width = 48,height=40)
pdf("figures/ketone.down_placebo_ribon_allMetabolites_dots.pdf",width = 48,height=40)

pdf("figures/ketone.up_placebo_ribon_allMetabolites_lines.pdf",width = 36, height=6)
pdf("figures/ketone.up_placebo_ribon_allMetabolites.pdf",width = 36, height=6)
pdf("figures/ketone.up_placebo_ribon_allMetabolites_dots.pdf",width = 36, height=6)

pdf("figures/placebo_ribon_allMetabolites_lines.pdf",width = 48, height=6)
pdf("figures/placebo_ribon_allMetabolites.pdf",width = 48, height=6)
pdf("figures/placebo_ribon_allMetabolites_dots.pdf",width = 48, height=6)


p

dev.off()

rm(p)
#ggsave("ketone.down_placebo_ribon_lines.pdf",path="./figures/",width = 18, height = 10, units = "cm")
#ggsave("ketone.down_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
#ggsave("ketone.down_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")



# volcano plots ketone vs control at 90 min or 180 min ==============
write.csv(ketone.ctrl.0.anno,
          "output/covariate_result_0min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
write.csv(ketone.ctrl.90.anno,
          "output/covariate_result_90min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
write.csv(ketone.ctrl.180.anno,
          "output/covariate_result_180min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
ketone.ctrl.0.anno <- read.csv("output/covariate_result_0min_ctrl.ketone_subject_annotate.csv",
                          check.names = F)
ketone.ctrl.90.anno <- read.csv("output/covariate_result_90min_ctrl.ketone_subject_annotate.csv",
                               check.names = F)
ketone.ctrl.180.anno <- read.csv("output/covariate_result_180min_ctrl.ketone_subject_annotate.csv",
                               check.names = F)
View(ketone.ctrl.0.anno) # no significant
View(ketone.ctrl.90.anno)
View(ketone.ctrl.180.anno)

View(plasma)
colnames(ketone.ctrl.90.anno)

df <- ketone.ctrl.90.anno %>%
  arrange(adj.P.Val)
## or: 
df <- ketone.ctrl.180.anno %>%
  arrange(adj.P.Val)

View(df)
colnames(df)

up <- df %>% 
  filter(adj.P.Val<0.05 & logFC.mean>0) #%>%
#  head(20)
head(up)

down <- df %>% 
  filter(adj.P.Val<0.05 & logFC.mean<0) #%>%
  #head(20)
head(down)


df <- df %>%
  mutate(sig= ifelse(adj.P.Val<0.05 & logFC.mean>0, "Up", 
                     ifelse(adj.P.Val<0.05 & logFC.mean<0,"Down","NS") 
  ),
  label= ifelse(sig!="NS" & 
                  metabolite %in% c(up$metabolite,down$metabolite), 
                metabolite, ""))


df$sig <- factor(df$sig, levels = c("Up","Down","NS"))
View(df)
#expression(italic("Insr")^italic("-/-")~"vs"~italic("Insr")^italic("+/+"))

# Volcano 
ggplot(data=df, 
       aes(x=logFC.mean, y=-log10(adj.P.Val), 
           col= sig,
           label=label)) +
  scale_color_manual(values=c("firebrick", # nothing positively correlate with leucine or Glu+FA
                              "dodgerblue3",#"blue",
                              "grey"), #c("#F8766D", "#00BFC4", "grey"),
                     guide = "none"
  ) +
  ylab(expression(-log[10]~(Adj~P.value))) +
  xlab(expression(log[2](Fold~Change)))+
  geom_point(aes(),
             size= 1.5) +
  #  scale_color_gradient2(high = "red3",mid = "white",low = "blue3",
  #                        midpoint = 0,na.value = "grey80"
  #                        #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
  #  )+
  #  scale_size_continuous(range = c(0.1, 4))+
  #  theme_minimal() +
  theme_bw()+
  geom_text_repel(
    data=subset(df,logFC.mean>0), 
    color="firebrick",
    size=4.5,
    #box.padding = 2,
    #  force = 5,
    #  force_pull = 3,
    min.segment.length = 0,
    nudge_x = 4 - subset(df,logFC.mean>0)$logFC.mean,
    segment.size=0.3, segment.color="grey", 
    direction="y", 
    hjust=0,
    max.overlaps = 200)+
  
  geom_text_repel(
    data=subset(df,logFC.mean<0), 
    color="dodgerblue3",
    size=4.3, ######### used 4.3 for 180 min, used 4.5 for 90 min###
    #box.padding = 2,
    #  force = 5,
    #  force_pull = 10,
    min.segment.length = 0,
    nudge_x = -2 - subset(df,logFC.mean<0)$logFC.mean,
    segment.size=0.3, segment.color="grey", 
    direction="y", 
    hjust=1,
    max.overlaps = 300)+
  
  #labs(size="Abundance (log10)",color="CoV (log10)")+
  scale_x_continuous(#breaks = c(-1,-0.5,0,0.5,1), 
    limits = c(-4, 4)
  )+
  scale_y_continuous(#breaks = c(-1,-0.5,0,0.5,1), 
    limits = c(0, 5)
  )+
 # annotate(geom="text",x=4, y=0, 
#           label=expression("Up in "*italic("Insr")^"f/f"),
#           color="firebrick",size=5) +
#  annotate(geom="text",x=-4, y=0, 
#           label=expression("Down in "*italic("Insr")^"f/f"),
#           color="steelblue3",size=5) +
  
  #xlim(-1.5,1.5) +
  #annotate(geom="text", x=2.1,y=-0.1, label="Up in T2D", color="#F8766D", size=7) +
  #annotate(geom="text", x=-2.1,y=-0.1, label="Up in ND", color="#00BFC4", size=7) +
  ggtitle(
    
########### change title here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #"Ketone vs Placebo 90 min"
    "Ketone vs Placebo 180 min"
    ) +
  theme(
    plot.title = element_text(color="black", size=18,hjust = 0.5),
    legend.position = "right",
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text=element_text(size=17),
    axis.text=element_text(size=17),
    axis.title=element_text(size=17),
    legend.text=element_text(size=17),
    legend.title=element_text(size=17),
    aspect.ratio = 1/1, panel.grid.major = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave(filename="figures/volcano_ketone.placebo_90min.png",width=6.5,height=6.8,units="in",dpi=600)
ggsave(filename="figures/volcano_ketone.placebo_90min.pdf",width=6.5,height=6.8,units="in")

ggsave(filename="figures/volcano_ketone.placebo_180min.png",width=6.5,height=6.8,units="in",dpi=600)
ggsave(filename="figures/volcano_ketone.placebo_180min.pdf",width=6.5,height=6.8,units="in")






# save spreadsheets for MetaboAnalyst.ca ===============

plasma.all <- read.csv("output/plasma_upload.csv", check.names = F)
View(plasma.all)
colnames(plasma.all)
plasma.all[plasma.all==0] <- NA
plasma.drink.time <- plasma.all[!rowSums(is.na(plasma.all))>ncol(plasma.all)*0.4,!colSums(is.na(plasma.all))>nrow(plasma.all)*0.4]
View(plasma.drink.time)
write.csv(plasma.drink.time,"output/plasma_drink_time.csv",row.names = F)

plasma.ketone <- plasma.all %>% select(all_of(grep("\\.A|sample",colnames(plasma.all))))
colnames(plasma.ketone)

plasma.0min <- plasma.all %>% 
  select(all_of(grep("\\.0|sample",colnames(plasma.all))))
colnames(plasma.0min)
plasma.90min <- plasma.all %>% 
  select(all_of(grep("\\.90|sample",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.90.A"))) # 2.90.B is a missing sample, so removed subject 2
colnames(plasma.90min)
plasma.180min <- plasma.all %>% select(all_of(grep("\\.180|sample",colnames(plasma.all)))) %>%
  select(-all_of(c("4.180.B","4.180.A","6.180.B","6.180.A"))) # 4.180.B, 6.180.A, 6.180.B are missing samples, so removed subject 4 & 6
colnames(plasma.180min)

plasma.drink.0.90min <- plasma.all %>% select(all_of(grep("\\.0|\\.90|sample",colnames(plasma.all))))
colnames(plasma.drink.0.90min)
plasma.drink.0.180min <- plasma.all %>% select(all_of(grep("\\.0|\\.180|sample",colnames(plasma.all))))
colnames(plasma.drink.0.180min)

plasma.ketone.0.90 <- plasma.all %>% 
  select(all_of(grep("sample|\\.90.A|\\.0.A",colnames(plasma.all)))) 
colnames(plasma.ketone.0.90)
plasma.ketone.0.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.0.A",colnames(plasma.all)))) %>%
  select(-all_of(c("6.180.A","6.0.A"))) # 6.180.A is a missing sample, so removed subject 6
colnames(plasma.ketone.0.180)
plasma.ketone.90.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.A|\\.90.A",colnames(plasma.all)))) %>%
  select(-all_of(c("6.180.A","6.90.A"))) # 6.180.A is a missing sample, so removed subject 6
colnames(plasma.ketone.90.180)

plasma.ctrl.0.90 <- plasma.all %>% 
  select(all_of(grep("sample|\\.90.B|\\.0.B",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.0.B"))) # 2.90.B is a missing sample, so removed subject 2
colnames(plasma.ctrl.0.90)
plasma.ctrl.0.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.0.B",colnames(plasma.all)))) %>%
  select(-all_of(c("4.180.B","4.0.B","6.180.B","6.0.B"))) # 4.180.B & 6.180.B are missing samples, so removed subject 4 & 6
colnames(plasma.ctrl.0.180)
plasma.ctrl.90.180 <- plasma.all %>% 
  select(all_of(grep("sample|\\.180.B|\\.90.B",colnames(plasma.all)))) %>%
  select(-all_of(c("2.90.B","2.180.B","4.180.B","4.90.B","6.180.B","6.90.B"))) # 2.90.B, 4.180.B & 6.180.B are missing samples, so removed subject 2, 4 & 6
colnames(plasma.ctrl.90.180)


meta.all <- read.csv("output/traits_upload.csv", check.names = F)
View(meta.all)
colnames(meta.all)
meta.all$Time <- paste0("T",meta.all$Time)
meta.all$Subject <- paste0("S",meta.all$Subject)
meta.drink.time <- meta.all %>% 
  filter(sample %in% colnames(plasma.drink.time))
write.csv(meta.all.save,"output/traits_drink_time.csv",row.names = F)

meta.ketone <- meta.all %>% filter(sample %in% colnames(plasma.ketone))
meta.ketone$sample

meta.0min <- meta.all %>% filter(sample %in% colnames(plasma.0min))
meta.0min$sample
meta.90min <- meta.all %>% filter(sample %in% colnames(plasma.90min))
meta.90min$sample
meta.180min <- meta.all %>% filter(sample %in% colnames(plasma.180min))
meta.180min$sample

meta.drink.0.90min <- meta.all %>% filter(sample %in% colnames(plasma.drink.0.90min))
meta.drink.0.90min$sample
meta.drink.0.180min <- meta.all %>% filter(sample %in% colnames(plasma.drink.0.180min))
meta.drink.0.180min$sample

meta.ketone.0.90 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.0.90))
meta.ketone.0.90$sample
meta.ketone.0.180 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.0.180))
meta.ketone.0.180$sample
meta.ketone.90.180 <- meta.all %>% filter(sample %in% colnames(plasma.ketone.90.180))
meta.ketone.90.180$sample

meta.ctrl.0.90 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.0.90))
meta.ctrl.0.90$sample
meta.ctrl.0.180 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.0.180))
meta.ctrl.0.180$sample
meta.ctrl.90.180 <- meta.all %>% filter(sample %in% colnames(plasma.ctrl.90.180))
meta.ctrl.90.180$sample

write.csv(plasma.ketone,"output/plasma_ketone.csv",row.names = F)

write.csv(plasma.0min,"output/plasma_0min.csv",row.names = F)
write.csv(plasma.90min,"output/plasma_90min.csv",row.names = F)
write.csv(plasma.180min,"output/plasma_180min.csv",row.names = F)

write.csv(plasma.drink.0.90min,"output/plasma_drink_0min.90min.csv",row.names = F)
write.csv(plasma.drink.0.180min,"output/plasma_drink_0min.180min.csv",row.names = F)

write.csv(plasma.ketone.0.90,"output/plasma_ketone_0min.90min.csv",row.names = F)
write.csv(plasma.ketone.0.180,"output/plasma_ketone_0min.180min.csv",row.names = F)
write.csv(plasma.ketone.90.180,"output/plasma_ketone_90min.180min.csv",row.names = F)

write.csv(plasma.ctrl.0.90,"output/plasma_ctrl_0min.90min.csv",row.names = F)
write.csv(plasma.ctrl.0.180,"output/plasma_ctrl_0min.180min.csv",row.names = F)
write.csv(plasma.ctrl.90.180,"output/plasma_ctrl_90min.180min.csv",row.names = F)

#
write.csv(meta.ketone,"output/traits_ketone.csv",row.names = F)

write.csv(meta.0min,"output/traits_0min.csv",row.names = F)
write.csv(meta.90min,"output/traits_90min.csv",row.names = F)
write.csv(meta.180min,"output/traits_180min.csv",row.names = F)

write.csv(meta.drink.0.90min,"output/traits_drink_0min.90min.csv",row.names = F)
write.csv(meta.drink.0.180min,"output/traits_drink_0min.180min.csv",row.names = F)

write.csv(meta.ketone.0.90,"output/traits_ketone_0min.90min.csv",row.names = F)
write.csv(meta.ketone.0.180,"output/traits_ketone_0min.180min.csv",row.names = F)
write.csv(meta.ketone.90.180,"output/traits_ketone_90min.180min.csv",row.names = F)

write.csv(meta.ctrl.0.90,"output/traits_ctrl_0min.90min.csv",row.names = F)
write.csv(meta.ctrl.0.180,"output/traits_ctrl_0min.180min.csv",row.names = F)
write.csv(meta.ctrl.90.180,"output/traits_ctrl_90min.180min.csv",row.names = F)

# organize MetaboAnalyst results ===================

ketone.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.90min_subject.csv",
                        row.names = 1)
ketone.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.180min_subject.csv",
                         row.names = 1)
ketone.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_90min.180min_subject.csv",
                          row.names = 1)
ctrl.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.90min_subject.csv",
                      row.names = 1)
ctrl.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.180min_subject.csv",
                       row.names = 1)
ctrl.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_90min.180min_subject.csv",
                        row.names = 1)
ketone.ctrl.0 <- read.csv("input/MetaboAnalyst_results/covariate_result_0min_ctrl.ketone_subject.csv",
                          row.names = 1)
ketone.ctrl.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_90min_ctrl.ketone_subject.csv",
                           row.names = 1)
ketone.ctrl.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_180min_ctrl.ketone_subject.csv",
                            row.names = 1)

# organize ketone.0.90
plasma.ketone.0.90 <- read.csv("output/plasma_ketone_0min.90min.csv", check.names = F)
plasma.ketone.0.90.proc <- plasma.ketone.0.90[!rowSums(is.na(plasma.ketone.0.90))>
                                                ncol(plasma.ketone.0.90)*0.4,
                                              !colSums(is.na(plasma.ketone.0.90))>
                                                nrow(plasma.ketone.0.90)*0.4]
View(plasma.ketone.0.90)
colnames(ketone.0.90)
ketone.0.90.anno <- rownames_to_column(ketone.0.90) %>% 
  left_join(plasma.ketone.0.90.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.0.90.anno)
View(ketone.0.90.anno)
colnames(ketone.0.90.anno)
ketone.0.90.anno <- ketone.0.90.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ketone.0.90.anno[,grep("\\.0\\.",colnames(ketone.0.90.anno))],na.rm = T),
         mean.90min=rowMeans(ketone.0.90.anno[,grep("\\.90\\.",colnames(ketone.0.90.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ketone.0.90.anno[,grep("\\.0\\.",colnames(ketone.0.90.anno))]),na.rm = T),
         median.90min=rowMedians(as.matrix(ketone.0.90.anno[,grep("\\.90\\.",colnames(ketone.0.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.90min/mean.0min),
         logFC.median=log2(median.90min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ketone.0.90.anno)

write.csv(ketone.0.90.anno,
          "output/covariate_result_ketone_0.90min_subject_annotate.csv",
          row.names = F)

# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.0.90.anno[is.na(ketone.0.90.anno)] <- 0

colnames(ketone.0.90.anno)
ketone.0.90.up <- ketone.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.0.90.up)

ketone.0.90.down <- ketone.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.0.90.down)

# organize ketone.0.180
plasma.ketone.0.180 <- read.csv("output/plasma_ketone_0min.180min.csv", check.names = F)
plasma.ketone.0.180.proc <- plasma.ketone.0.180[!rowSums(is.na(plasma.ketone.0.180))>
                                                  ncol(plasma.ketone.0.180)*0.4,
                                                !colSums(is.na(plasma.ketone.0.180))>
                                                  nrow(plasma.ketone.0.180)*0.4]
View(plasma.ketone.0.180)
colnames(ketone.0.180)
ketone.0.180.anno <- rownames_to_column(ketone.0.180) %>% 
  left_join(plasma.ketone.0.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.0.180.anno)
View(ketone.0.180.anno)
colnames(ketone.0.180.anno)
ketone.0.180.anno <- ketone.0.180.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ketone.0.180.anno[,grep("\\.0\\.",colnames(ketone.0.180.anno))],na.rm = T),
         mean.180min=rowMeans(ketone.0.180.anno[,grep("\\.180\\.",colnames(ketone.0.180.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ketone.0.180.anno[,grep("\\.0\\.",colnames(ketone.0.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ketone.0.180.anno[,grep("\\.180\\.",colnames(ketone.0.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.0min),
         logFC.median=log2(median.180min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ketone.0.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.0.180.anno[is.na(ketone.0.180.anno)] <- 0

colnames(ketone.0.180.anno)
ketone.0.180.up <- ketone.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.0.180.up)

ketone.0.180.down <- ketone.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.0.180.down)

write.csv(ketone.0.180.anno,
          "output/covariate_result_ketone_0.180min_subject_annotate.csv",
          row.names = F)


# organize ketone.90.180
plasma.ketone.90.180 <- read.csv("output/plasma_ketone_90min.180min.csv", check.names = F)
View(plasma.ketone.90.180)
plasma.ketone.90.180.proc <- plasma.ketone.90.180[!rowSums(is.na(plasma.ketone.90.180))>
                                                    ncol(plasma.ketone.90.180)*0.4,
                                                  !colSums(is.na(plasma.ketone.90.180))>
                                                    nrow(plasma.ketone.90.180)*0.4]
View(plasma.ketone.90.180.proc)
colnames(ketone.90.180)
ketone.90.180.anno <- rownames_to_column(ketone.90.180) %>% 
  left_join(plasma.ketone.90.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.90.180.anno)
View(ketone.90.180.anno)
colnames(ketone.90.180.anno)
ketone.90.180.anno <- ketone.90.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.90min=rowMeans(ketone.90.180.anno[,grep("\\.90\\.",colnames(ketone.90.180.anno))],na.rm = T),
         mean.180min=rowMeans(ketone.90.180.anno[,grep("\\.180\\.",colnames(ketone.90.180.anno))],na.rm = T),
         median.90min=rowMedians(as.matrix(ketone.90.180.anno[,grep("\\.90\\.",colnames(ketone.90.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ketone.90.180.anno[,grep("\\.180\\.",colnames(ketone.90.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.90min),
         logFC.median=log2(median.180min/median.90min)) %>%
  relocate(mean.90min:logFC.median,.after = logFC)
View(ketone.90.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.90.180.anno[is.na(ketone.90.180.anno)] <- 0

colnames(ketone.90.180.anno)
ketone.90.180.up <- ketone.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.90.180.up) # no metabolite

ketone.90.180.down <- ketone.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.90.180.down) # 1 metabolite

write.csv(ketone.90.180.anno,
          "output/covariate_result_ketone_90.180min_subject_annotate.csv",
          row.names = F)

# organize ctrl.0.90
plasma.ctrl.0.90 <- read.csv("output/plasma_ctrl_0min.90min.csv", check.names = F)
plasma.ctrl.0.90.proc <- plasma.ctrl.0.90[!rowSums(is.na(plasma.ctrl.0.90))>
                                            ncol(plasma.ctrl.0.90)*0.4,
                                          !colSums(is.na(plasma.ctrl.0.90))>
                                            nrow(plasma.ctrl.0.90)*0.4]
View(plasma.ctrl.0.90)
colnames(ctrl.0.90)
ctrl.0.90.anno <- rownames_to_column(ctrl.0.90) %>% 
  left_join(plasma.ctrl.0.90.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.0.90.anno)
View(ctrl.0.90.anno)
colnames(ctrl.0.90.anno)
ctrl.0.90.anno <- ctrl.0.90.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ctrl.0.90.anno[,grep("\\.0\\.",colnames(ctrl.0.90.anno))],na.rm = T),
         mean.90min=rowMeans(ctrl.0.90.anno[,grep("\\.90\\.",colnames(ctrl.0.90.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ctrl.0.90.anno[,grep("\\.0\\.",colnames(ctrl.0.90.anno))]),na.rm = T),
         median.90min=rowMedians(as.matrix(ctrl.0.90.anno[,grep("\\.90\\.",colnames(ctrl.0.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.90min/mean.0min),
         logFC.median=log2(median.90min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ctrl.0.90.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.0.90.anno[is.na(ctrl.0.90.anno)] <- 0

colnames(ctrl.0.90.anno)
ctrl.0.90.up <- ctrl.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.0.90.up) # no metabolite

ctrl.0.90.down <- ctrl.0.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.0.90.down) # 5 metabolites

write.csv(ctrl.0.90.anno,
          "output/covariate_result_ctrl_0.90min_subject_annotate.csv",
          row.names = F)

# organize ctrl.0.180
plasma.ctrl.0.180 <- read.csv("output/plasma_ctrl_0min.180min.csv", check.names = F)
plasma.ctrl.0.180.proc <- plasma.ctrl.0.180[!rowSums(is.na(plasma.ctrl.0.180))>
                                              ncol(plasma.ctrl.0.180)*0.4,
                                            !colSums(is.na(plasma.ctrl.0.180))>
                                              nrow(plasma.ctrl.0.180)*0.4]
View(plasma.ctrl.0.180)
colnames(ctrl.0.180)
ctrl.0.180.anno <- rownames_to_column(ctrl.0.180) %>% 
  left_join(plasma.ctrl.0.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.0.180.anno)
View(ctrl.0.180.anno)
colnames(ctrl.0.180.anno)
ctrl.0.180.anno <- ctrl.0.180.anno%>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.0min=rowMeans(ctrl.0.180.anno[,grep("\\.0\\.",colnames(ctrl.0.180.anno))],na.rm = T),
         mean.180min=rowMeans(ctrl.0.180.anno[,grep("\\.180\\.",colnames(ctrl.0.180.anno))],na.rm = T),
         median.0min=rowMedians(as.matrix(ctrl.0.180.anno[,grep("\\.0\\.",colnames(ctrl.0.180.anno))],na.rm = T)),
         median.180min=rowMedians(as.matrix(ctrl.0.180.anno[,grep("\\.180\\.",colnames(ctrl.0.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.0min),
         logFC.median=log2(median.180min/median.0min)) %>%
  relocate(mean.0min:logFC.median,.after = logFC)
View(ctrl.0.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.0.180.anno[is.na(ctrl.0.180.anno)] <- 0

colnames(ctrl.0.180.anno)
ctrl.0.180.up <- ctrl.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.0.180.up) # no metabolite

ctrl.0.180.down <- ctrl.0.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.0.180.down) # 6 metabolite

write.csv(ctrl.0.180.anno,
          "output/covariate_result_ctrl_0.180min_subject_annotate.csv",
          row.names = F)

# organize ctrl.90.180

plasma.ctrl.90.180 <- read.csv("output/plasma_ctrl_90min.180min.csv", check.names = F)
View(plasma.ctrl.90.180)
plasma.ctrl.90.180.proc <- plasma.ctrl.90.180[!rowSums(is.na(plasma.ctrl.90.180))>
                                                ncol(plasma.ctrl.90.180)*0.4,
                                              !colSums(is.na(plasma.ctrl.90.180))>
                                                nrow(plasma.ctrl.90.180)*0.4]
View(plasma.ctrl.90.180.proc)
colnames(ctrl.90.180)
ctrl.90.180.anno <- rownames_to_column(ctrl.90.180) %>% 
  left_join(plasma.ctrl.90.180.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ctrl.90.180.anno)
View(ctrl.90.180.anno)
colnames(ctrl.90.180.anno)
ctrl.90.180.anno <- ctrl.90.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.90min=rowMeans(ctrl.90.180.anno[,grep("\\.90\\.",colnames(ctrl.90.180.anno))],na.rm = T),
         mean.180min=rowMeans(ctrl.90.180.anno[,grep("\\.180\\.",colnames(ctrl.90.180.anno))],na.rm = T),
         median.90min=rowMedians(as.matrix(ctrl.90.180.anno[,grep("\\.90\\.",colnames(ctrl.90.180.anno))]),na.rm = T),
         median.180min=rowMedians(as.matrix(ctrl.90.180.anno[,grep("\\.180\\.",colnames(ctrl.90.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.180min/mean.90min),
         logFC.median=log2(median.180min/median.90min)) %>%
  relocate(mean.90min:logFC.median,.after = logFC)
View(ctrl.90.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ctrl.90.180.anno[is.na(ctrl.90.180.anno)] <- 0

colnames(ctrl.90.180.anno)
ctrl.90.180.up <- ctrl.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ctrl.90.180.up) # 1 beta-Hydroxybutyric acid

ctrl.90.180.down <- ctrl.90.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ctrl.90.180.down) # no metabolite

write.csv(ctrl.90.180.anno,
          "output/covariate_result_ctrl_90.180min_subject_annotate.csv",
          row.names = F)

# organize ketone.ctrl.0
plasma.0min <- read.csv("output/plasma_0min.csv", check.names = F)
View(plasma.0min)
plasma.0min.proc <- plasma.0min[!rowSums(is.na(plasma.0min))>
                                  ncol(plasma.0min)*0.4,
                                !colSums(is.na(plasma.0min))>
                                  nrow(plasma.0min)*0.4]
View(plasma.0min.proc)
colnames(ketone.ctrl.0)
ketone.ctrl.0.anno <- rownames_to_column(ketone.ctrl.0) %>% 
  left_join(plasma.0min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.0.anno)
View(ketone.ctrl.0.anno)
colnames(ketone.ctrl.0.anno)
ketone.ctrl.0.anno <- ketone.ctrl.0.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.0.anno[,grep("\\.B",colnames(ketone.ctrl.0.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.0.anno[,grep("\\.A",colnames(ketone.ctrl.0.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.0.anno[,grep("\\.B",colnames(ketone.ctrl.0.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.0.anno[,grep("\\.A",colnames(ketone.ctrl.0.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.0.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.0.anno[is.na(ketone.ctrl.0.anno)] <- 0

colnames(ketone.ctrl.0.anno)
ketone.ctrl.0.up <- ketone.ctrl.0.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.0.up) # no metabolite

ketone.ctrl.0.down <- ketone.ctrl.0.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.0.down) # no metabolite


# organize ketone.ctrl.90
plasma.90min <- read.csv("output/plasma_90min.csv", check.names = F)
View(plasma.90min)
plasma.90min.proc <- plasma.90min[!rowSums(is.na(plasma.90min))>
                                    ncol(plasma.90min)*0.4,
                                  !colSums(is.na(plasma.90min))>
                                    nrow(plasma.90min)*0.4]
View(plasma.90min.proc)
colnames(ketone.ctrl.90)
ketone.ctrl.90.anno <- rownames_to_column(ketone.ctrl.90) %>% 
  left_join(plasma.90min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.90.anno)
View(ketone.ctrl.90.anno)
colnames(ketone.ctrl.90.anno)
ketone.ctrl.90.anno <- ketone.ctrl.90.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.90.anno[,grep("\\.B",colnames(ketone.ctrl.90.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.90.anno[,grep("\\.A",colnames(ketone.ctrl.90.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.90.anno[,grep("\\.B",colnames(ketone.ctrl.90.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.90.anno[,grep("\\.A",colnames(ketone.ctrl.90.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.90.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.90.anno[is.na(ketone.ctrl.90.anno)] <- 0

colnames(ketone.ctrl.90.anno)
ketone.ctrl.90.up <- ketone.ctrl.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.90.up) 

ketone.ctrl.90.down <- ketone.ctrl.90.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.90.down) 

# organize ketone.ctrl.180
plasma.180min <- read.csv("output/plasma_180min.csv", check.names = F)
View(plasma.180min)
plasma.180min.proc <- plasma.180min[!rowSums(is.na(plasma.180min))>
                                      ncol(plasma.180min)*0.4,
                                    !colSums(is.na(plasma.180min))>
                                      nrow(plasma.180min)*0.4]
View(plasma.180min.proc)
colnames(ketone.ctrl.180)
ketone.ctrl.180.anno <- rownames_to_column(ketone.ctrl.180) %>% 
  left_join(plasma.180min.proc, by=c("rowname"="sample")) %>%
  rename(metabolite=rowname) %>%
  mutate_at(-1,as.numeric)
str(ketone.ctrl.180.anno)
View(ketone.ctrl.180.anno)
colnames(ketone.ctrl.180.anno)
ketone.ctrl.180.anno <- ketone.ctrl.180.anno %>%
  # mean and median were calculated without NA. Later, NA will be given 0 for plotting.
  mutate(mean.ctrl=rowMeans(ketone.ctrl.180.anno[,grep("\\.B",colnames(ketone.ctrl.180.anno))],na.rm = T),
         mean.ketone=rowMeans(ketone.ctrl.180.anno[,grep("\\.A",colnames(ketone.ctrl.180.anno))],na.rm = T),
         median.ctrl=rowMedians(as.matrix(ketone.ctrl.180.anno[,grep("\\.B",colnames(ketone.ctrl.180.anno))]),na.rm = T),
         median.ketone=rowMedians(as.matrix(ketone.ctrl.180.anno[,grep("\\.A",colnames(ketone.ctrl.180.anno))]),na.rm = T)
  ) %>%
  mutate(logFC.mean=log2(mean.ketone/mean.ctrl),
         logFC.median=log2(median.ketone/median.ctrl)) %>%
  relocate(mean.ctrl:logFC.median,.after = logFC)
View(ketone.ctrl.180.anno)
# missing values NA were assigned 0 for plotting purpose. (In Metabolomics, they were imputed as 1/5 lowest values.)
ketone.ctrl.180.anno[is.na(ketone.ctrl.180.anno)] <- 0

colnames(ketone.ctrl.180.anno)
ketone.ctrl.180.up <- ketone.ctrl.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean>0)
View(ketone.ctrl.180.up) 

ketone.ctrl.180.down <- ketone.ctrl.180.anno %>%
  filter(adj.P.Val<0.05 & logFC.mean<0)
View(ketone.ctrl.180.down) 

# plot altered metebolits in different comparisons onto UpSet plot ======================
sets1 <- ketone.0.90.up %>% mutate(ketone_0min.90min_up=1) %>%
  select(all_of(c("metabolite","ketone_0min.90min_up")))
sets2 <- ketone.0.90.down %>% mutate(ketone_0min.90min_down=1) %>%
  select(all_of(c("metabolite","ketone_0min.90min_down")))
sets3 <- ketone.0.180.up %>% mutate(ketone_0min.180min_up=1) %>%
  select(all_of(c("metabolite","ketone_0min.180min_up")))
sets4 <- ketone.0.180.down %>% mutate(ketone_0min.180min_down=1) %>%
  select(all_of(c("metabolite","ketone_0min.180min_down")))
sets5 <- ketone.90.180.up %>% mutate(ketone_90min.180min_up=1) %>%
  select(all_of(c("metabolite","ketone_90min.180min_up")))
sets6 <- ketone.90.180.down %>% mutate(ketone_90min.180min_down=1) %>%
  select(all_of(c("metabolite","ketone_90min.180min_down")))

sets7 <- ctrl.0.90.up %>% mutate(ctrl_0min.90min_up=1) %>%
  select(all_of(c("metabolite","ctrl_0min.90min_up")))
sets8 <- ctrl.0.90.down %>% mutate(ctrl_0min.90min_down=1) %>%
  select(all_of(c("metabolite","ctrl_0min.90min_down")))
sets9 <- ctrl.0.180.up %>% mutate(ctrl_0min.180min_up=1) %>%
  select(all_of(c("metabolite","ctrl_0min.180min_up")))
sets10 <- ctrl.0.180.down %>% mutate(ctrl_0min.180min_down=1) %>%
  select(all_of(c("metabolite","ctrl_0min.180min_down")))
sets11 <- ctrl.90.180.up %>% mutate(ctrl_90min.180min_up=1) %>%
  select(all_of(c("metabolite","ctrl_90min.180min_up")))
sets12 <- ctrl.90.180.down %>% mutate(ctrl_90min.180min_down=1) %>%
  select(all_of(c("metabolite","ctrl_90min.180min_down")))

sets13 <- ketone.ctrl.0.up %>% mutate(ketone.ctrl_0min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_0min_up")))
sets14 <- ketone.ctrl.0.down %>% mutate(ketone.ctrl_0min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_0min_down")))
sets15 <- ketone.ctrl.90.up %>% mutate(ketone.ctrl_90min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_90min_up")))
sets16 <- ketone.ctrl.90.down %>% mutate(ketone.ctrl_90min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_90min_down")))
sets17 <- ketone.ctrl.180.up %>% mutate(ketone.ctrl_180min_up=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_180min_up")))
sets18 <- ketone.ctrl.180.down %>% mutate(ketone.ctrl_180min_down=1) %>%
  select(all_of(c("metabolite","ketone.ctrl_180min_down")))

# combine all sets
sets <- sets1 %>%
  full_join(sets2) %>%
  full_join(sets3) %>%
  full_join(sets4) %>%
  full_join(sets5) %>%
  full_join(sets6) %>%
  full_join(sets7) %>%
  full_join(sets8) %>%
  full_join(sets9) %>%
  full_join(sets10) %>%
  full_join(sets11) %>%
  full_join(sets12) %>%
  full_join(sets13) %>%
  full_join(sets14) %>%
  full_join(sets15) %>%
  full_join(sets16) %>%
  full_join(sets17) %>%
  full_join(sets18)
View(sets)
sets[,-1] = sets[,-1] == 1 # convert 1 to TRUE
View(sets)

# combine ketone 0, 90, 180 min
sets.ketone.time <- sets1 %>%
  full_join(sets2) %>%
  full_join(sets3) %>%
  full_join(sets4) %>%
  full_join(sets5) %>%
  full_join(sets6) 
sets.ketone.time[,-1] = sets.ketone.time[,-1] == 1 # convert 1 to TRUE
View(sets.ketone.time)
# combine ctrl 0, 90, 180
sets.ctrl.time <- sets7 %>%
  full_join(sets8) %>%
  full_join(sets9) %>%
  full_join(sets10) %>%
  full_join(sets11) %>%
  full_join(sets12) 
View(sets.ctrl.time)
sets.ctrl.time[,-1] = sets.ctrl.time[,-1] == 1 # convert 1 to TRUE
View(sets.ctrl.time)
#https://krassowski.github.io/complex-upset/articles/Examples_R.html

## select the data columns!!!!!!
data.sets <- colnames(sets)[-1]
data.sets <- colnames(sets.ctrl.time)[-1]
data.sets <- colnames(sets.ketone.time)[-1]

ComplexUpset::upset(
  
  ## select the data !!!!!!!!!!!!
  sets,
  #sets.ctrl.time,
  #sets.ketone.time,
  
  data.sets, name='', 
  sort_intersections_by='degree',
  width_ratio=0.2,
  height_ratio=0.6, 
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=12)),
      'overall_sizes'=theme(axis.text.x=element_text(angle=90))
    )),
  
  # min_size=500,
  
  base_annotations = list(
    'Intersection size'=(
      intersection_size(
        
        bar_number_threshold=1000000,
        #text=list(vjust=1.1)
        #counts=FALSE,
        #mapping=aes(fill=general.type)
      )
      
      + expand_limits(y=50)
      #    + scale_fill_manual(values=cbPalette)
      + theme(legend.title = element_blank(),
              axis.text.y =element_blank(),
              axis.title.y = element_text(margin = margin(r = -250)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
    )),
  set_sizes=(
    upset_set_size()
    + ylab('Set size')
    + geom_text(aes(label=after_stat(count)), hjust=1, stat='count')
    + expand_limits(y=80)
    + theme(axis.title.x = element_text(margin = margin(t = -10)),
            axis.text.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  ),
  
  stripes = c("white")) + 
  theme(panel.grid = element_blank(),
        text = element_text(size=14))+
  
  patchwork::plot_layout(heights=c(0.6, 1))

ggsave(filename="figures/upset.png",width=8,height=5,units="in",dpi=600)
ggsave(filename="figures/upset.pdf",width=8,height=5,units="in")
ggsave(filename="figures/upset_ctrl.pdf",width=8,height=5,units="in")
ggsave(filename="figures/upset_ketone.pdf",width=8,height=5,units="in")

# plot lines ketone down ==============
colnames(ketone.0.180.down)
colnames(ketone.0.90.down)
#setdiff(colnames(ketone.0.90.down),colnames(ketone.0.180.down))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.down),colnames(ketone.0.180.down)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ketone.0.90.180.down <- ketone.0.90.down[,1:13] %>%
  full_join(ketone.0.180.down[,1:13], suffix = c("_0.90","_0.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ketone.0.90.180.down)
colnames(ketone.0.90.180.down)

# normalized values of each metabolite to mean ketone 0 min

ketone.down.norm <- ketone.0.90.180.down[,c(1,7,19,26:123)] 
colnames(ketone.down.norm)
ketone.down.norm$mean.ketone.0min <- rowMeans(ketone.down.norm[,grep("\\.0.A",colnames(ketone.down.norm))],na.rm = T)
ketone.down.norm <- ketone.down.norm %>%
  mutate_at(vars(4:101),list(norm=~./mean.ketone.0min))
ketone.down.norm <- ketone.down.norm[,-(4:102)]
colnames(ketone.down.norm) <- gsub("_norm","",colnames(ketone.down.norm))
View(ketone.down.norm)
colnames(ketone.down.norm)

ketone.down.norm.long <- ketone.down.norm[,c(1:3,4:101)] %>%
  pivot_longer(-(1:3),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ketone.down.norm.long$drink<-factor(ketone.down.norm.long$drink,
                                    levels=c("Placebo","Ketone"))
levels(ketone.down.norm.long$drink)
View(ketone.down.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ketone.down.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)
p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  #geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
  #                 method = "quasirandom",##"pseudorandom",
  #width = 3, #nbins = 10, #bandwidth = 1,#varwidth = TRUE, # different width based on dot density
  #                 dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("ketone.down_placebo_ribon_lines.pdf",path="./figures/",width = 18, height = 10, units = "cm")
ggsave("ketone.down_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.down_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")


# plot lines ketone up ==============

colnames(ketone.0.180.up)
colnames(ketone.0.90.up)
colnames(ketone.90.180.up)
#setdiff(colnames(ketone.0.90.up),colnames(ketone.0.180.up))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.up),colnames(ketone.0.180.up)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ketone.0.90.180.up <- ketone.0.90.up[,1:13] %>%
  full_join(ketone.0.180.up[,1:13], suffix = c("","_0.180"),by="metabolite") %>%
  full_join(ketone.90.180.down[,1:13], suffix = c("_0.90","_90.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ketone.0.90.180.up)
colnames(ketone.0.90.180.up)

# normalized values of each metabolite to mean ketone 0 min

ketone.up.norm <- ketone.0.90.180.up[,c(1,7,19,31,38:135)] 
colnames(ketone.up.norm)
ketone.up.norm$mean.ketone.0min <- rowMeans(ketone.up.norm[,grep("\\.0.A",colnames(ketone.up.norm))],na.rm = T)
ketone.up.norm <- ketone.up.norm %>%
  mutate_at(vars(5:102),list(norm=~./mean.ketone.0min))
colnames(ketone.up.norm)
ketone.up.norm <- ketone.up.norm[,-(5:103)]
colnames(ketone.up.norm) <- gsub("_norm","",colnames(ketone.up.norm))
View(ketone.up.norm)
colnames(ketone.up.norm)

ketone.up.norm.long <- ketone.up.norm %>%
  pivot_longer(-(1:4),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ketone.up.norm.long$drink<-factor(ketone.up.norm.long$drink,
                                  levels=c("Placebo","Ketone"))
levels(ketone.up.norm.long$drink)
View(ketone.up.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ketone.up.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)

p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  #geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
                   method = "quasirandom",##"pseudorandom",
                   #width = 3,#nbins = 10,#bandwidth = 1, #varwidth = TRUE, # different width based on dot density
                   dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(
    
    # (edit here!!!!!) two metabolite in ketone group is very high, so set limits to make other dots and lines more visible.     
    limits = c(0.5,1.5),
    breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("ketone.up_placebo_ribon_lines.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_lines_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("ketone.up_placebo_ribon_dots_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")

# plot lines ctrl ==============

colnames(ctrl.0.180.down)
colnames(ctrl.0.90.down)
colnames(ctrl.90.180.up)
#setdiff(colnames(ketone.0.90.up),colnames(ketone.0.180.up))

colnames(plasma.drink.time)
#plasma.plot <- plasma.drink.time %>%
#  select(any_of(c("sample",colnames(ketone.0.90.up),colnames(ketone.0.180.up)))) 
#View(plasma.ketone.plot)
#colnames(plasma.ketone.plot)

ctrl.0.90.180 <- ctrl.0.90.down[,1:13] %>%
  full_join(ctrl.0.180.down[,1:13], suffix = c("","_0.180"),by="metabolite") %>%
  full_join(ctrl.90.180.up[,1:13], suffix = c("_0.90","_90.180"),by="metabolite") %>%
  left_join(plasma.drink.time,by=c("metabolite"="sample")) %>%
  mutate_at(-1,as.numeric)
View(ctrl.0.90.180)
colnames(ctrl.0.90.180)

# normalized values of each metabolite to mean ketone 0 min

ctrl.norm <- ctrl.0.90.180[,c(1,7,19,31,38:135)] 
colnames(ctrl.norm)
ctrl.norm$mean.ketone.0min <- rowMeans(ctrl.norm[,grep("\\.0.A",colnames(ctrl.norm))],na.rm = T)
ctrl.norm <- ctrl.norm %>%
  mutate_at(vars(5:102),list(norm=~./mean.ketone.0min))
colnames(ctrl.norm)
ctrl.norm <- ctrl.norm[,-(5:103)]
colnames(ctrl.norm) <- gsub("_norm","",colnames(ctrl.norm))
View(ctrl.norm)
colnames(ctrl.norm)

ctrl.norm.long <- ctrl.norm %>%
  pivot_longer(-(1:4),
               names_to = "sample",
               values_to = "Value"
  ) %>%
  left_join(meta.drink.time, by="sample") %>%
  mutate(time=gsub("T","",Time),
         drink=ifelse(Phenotype=="A", "Ketone", "Placebo")) %>%
  relocate(time:drink,.after = Time) %>%
  mutate_at("time", as.numeric) 
ctrl.norm.long$drink<-factor(ctrl.norm.long$drink,
                             levels=c("Placebo","Ketone"))
levels(ctrl.norm.long$drink)
View(ctrl.norm.long)

colnames(meta.drink.time)


## plot general pattern lines ---------
dat <- ctrl.norm.long %>%
  group_by(metabolite,time,drink) %>%
  summarise(Value.UQ=quantile(Value, 0.75,na.rm=T),Value.LQ=quantile(Value,0.25,na.rm=T),
            Value.Median=median(Value,na.rm=T),Value.Mean=mean(Value,na.rm=T))
View(dat)
datTab <- dat %>%
  group_by(time,drink) %>%
  summarise(UQ=quantile(Value.Mean, 0.75,na.rm=T),LQ=quantile(Value.Mean,0.25,na.rm=T),
            Median=median(Value.Mean,na.rm=T),Mean=mean(Value.Mean,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)

p<-ggplot(dat,aes(x=time,y=Value.Mean,group=drink))
p+geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(metabolite,drink),colour=drink),linetype="solid",size=0.5,alpha=0.4 )+ #
  #geom_beeswarm(aes(group=interaction(metabolite,drink),
  #                  colour=drink),
  #              data=dat[dat$drink=="Ketone",],
  #              side=1L,
  #              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,size=1.2,alpha=1)+
  #geom_quasirandom(aes(group=interaction(metabolite,drink),colour=drink),
  #                 method = "quasirandom",##"pseudorandom",
  #                 #width = 3,#nbins = 10,#bandwidth = 1, #varwidth = TRUE, # different width based on dot density
  #                 dodge.width=10,cex=0.8, alpha=0.3)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(
    
    # (edit here!!!!!) one metabolite in ketone group is very high, so set limits to make other dots and lines more visible.     
    limits = c(0.6,1.5),
    breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "right") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("placebo.updown_ketone_ribon_lines_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("placebo.updown_ketone_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
ggsave("placebo.updown_ketone_ribon_dots_limits.pdf",path="./figures/",width = 15, height = 10, units = "cm")

# plot all metabolites trends ===========
colnames(ketone.down.norm.long)

####### select data to plot !!!!!!!!!!
dat <- ketone.down.norm.long %>%
  arrange(logFC.mean_0.180,logFC.mean_0.90)
dat <- ketone.up.norm.long %>%
  arrange(-logFC.mean_0.180,-logFC.mean_0.90)
dat <- ctrl.norm.long %>%
  arrange(-logFC.mean_0.180,-logFC.mean_0.90)

dat$metabolite <- factor(dat$metabolite,levels=unique(dat$metabolite))
dat$drink<-factor(dat$drink,levels=c("Placebo","Ketone"))
View(dat)
length(unique(dat$metabolite))

df <- dat[which(dat$Value==0),]
View(df)

datTab <- dat %>%
  group_by(metabolite,time,drink)  %>%
  summarise(UQ=quantile(Value, 0.75,na.rm=T),LQ=quantile(Value,0.25,na.rm=T),
            Median=median(Value,na.rm=T),Mean=mean(Value,na.rm=T))
datTab$drink<-factor(datTab$drink,levels=c("Placebo","Ketone"))

View(datTab)

##### Plot and save figure #####
cbPalette <- c(
  "#999999", #grey
  "#E69F00", #lightorange
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#56B4E9", #blue
  "#CC79A7", #magenta
  "grey90",
  
  "#009E73", #green
  "#F0E442" #yellow
)
p<-ggplot(dat,aes(x=time,y=Value,group=drink)) +
  geom_ribbon(aes(x=time,y=Median,ymin=LQ,ymax=UQ,fill=drink), 
              data=datTab,alpha=0.2)+
  geom_line(aes(group=interaction(Subject,drink),colour=drink),
            linetype="solid",linewidth=0.5,alpha=0.4 )+ #
  #  geom_quasirandom(aes(group=interaction(Subject,drink),colour=drink),
  #                   method = "quasirandom",##"pseudorandom", #width = 3, #nbins = 10, #bandwidth = 1,#varwidth = TRUE, # different width based on dot density
  #                   dodge.width=10,cex=0.8, alpha=0.3)+
  
  ##geom_beeswarm(aes(group=interaction(Subject,drink),
  ##                  colour=drink),
  ##              data=dat[dat$drink=="Ketone",],
  ##              side=1L,
  ##              size=0.5,alpha=0.4 ) +
  geom_line(aes(y=Median,colour = drink,group=drink), #
            data=datTab,linewidth=1.2,alpha=1)+
  geom_point(aes(y=Median,colour = drink,group=drink,#shape=drink,
                 fill=drink), 
             data=datTab,size=2,alpha=1)+
  #coord_cartesian(ylim=c(0,12))+
  scale_x_continuous(breaks=unique(dat$time))+
  scale_y_continuous(breaks=pretty_breaks(n=6))+
  labs(x="Time (minutes)",y="Relative concentration")+
  scale_colour_manual(name="Drink", values=cbPalette)+
  #scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Drink", values=cbPalette)+
  #geom_hline(yintercept = 1,colour="black",linetype=3)+
  theme(axis.title = element_text(#family = "Roboto Light", 
    color="black", size=14),
    axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(#family = "Roboto Light",
    color="black",size=14))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(#family = "Roboto Light",
    color="black",size=14))+
  theme(legend.text = element_text(#family = "Roboto Light",
    color="black",size=14), 
    legend.title = element_blank(),
    legend.position = "none") + #c(0.15,0.9))+
  theme(strip.text = element_text(#family="Roboto Light",
    color="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(
    facets=~metabolite,
    
    ####### change # of rows and columes here !!!!!!!!!!!!!
    #nrow = 8, ncol = 8, # for ketone down
    #nrow = 1, ncol = 6, # for ketone up
    nrow = 1, ncol = 8, # for placebo
    scales = "free",
    shrink = TRUE,
    labeller = "label_value",
    as.table = TRUE,
    #switch = deprecated(),
    drop = TRUE,
    dir = "h",
    strip.position = "top"
  )


pdf("figures/ketone.down_placebo_ribon_allMetabolites_lines.pdf",width = 48,height=40)
pdf("figures/ketone.down_placebo_ribon_allMetabolites.pdf",width = 48,height=40)
pdf("figures/ketone.down_placebo_ribon_allMetabolites_dots.pdf",width = 48,height=40)

pdf("figures/ketone.up_placebo_ribon_allMetabolites_lines.pdf",width = 36, height=6)
pdf("figures/ketone.up_placebo_ribon_allMetabolites.pdf",width = 36, height=6)
pdf("figures/ketone.up_placebo_ribon_allMetabolites_dots.pdf",width = 36, height=6)

pdf("figures/placebo_ribon_allMetabolites_lines.pdf",width = 48, height=6)
pdf("figures/placebo_ribon_allMetabolites.pdf",width = 48, height=6)
pdf("figures/placebo_ribon_allMetabolites_dots.pdf",width = 48, height=6)


p

dev.off()

rm(p)
#ggsave("ketone.down_placebo_ribon_lines.pdf",path="./figures/",width = 18, height = 10, units = "cm")
#ggsave("ketone.down_placebo_ribon.pdf",path="./figures/",width = 15, height = 10, units = "cm")
#ggsave("ketone.down_placebo_ribon_dots.pdf",path="./figures/",width = 15, height = 10, units = "cm")



# volcano plots ketone vs control at 90 min or 180 min ==============
write.csv(ketone.ctrl.0.anno,
          "output/covariate_result_0min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
write.csv(ketone.ctrl.90.anno,
          "output/covariate_result_90min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
write.csv(ketone.ctrl.180.anno,
          "output/covariate_result_180min_ctrl.ketone_subject_annotate.csv",
          row.names = F)
ketone.ctrl.0.anno <- read.csv("output/covariate_result_0min_ctrl.ketone_subject_annotate.csv",
                               check.names = F)
ketone.ctrl.90.anno <- read.csv("output/covariate_result_90min_ctrl.ketone_subject_annotate.csv",
                                check.names = F)
ketone.ctrl.180.anno <- read.csv("output/covariate_result_180min_ctrl.ketone_subject_annotate.csv",
                                 check.names = F)
View(ketone.ctrl.0.anno)
View(ketone.ctrl.90.anno)
View(ketone.ctrl.180.anno)

View(plasma)
colnames(ketone.ctrl.90.anno)

df <- ketone.ctrl.90.anno %>%
  arrange(adj.P.Val)
## or: 
df <- ketone.ctrl.180.anno %>%
  arrange(adj.P.Val)

View(df)
colnames(df)

up <- df %>% 
  filter(adj.P.Val<0.05 & logFC.mean>0) #%>%
#  head(20)
head(up)

down <- df %>% 
  filter(adj.P.Val<0.05 & logFC.mean<0) #%>%
#head(20)
head(down)


df <- df %>%
  mutate(sig= ifelse(adj.P.Val<0.05 & logFC.mean>0, "Up", 
                     ifelse(adj.P.Val<0.05 & logFC.mean<0,"Down","NS") 
  ),
  label= ifelse(sig!="NS" & 
                  metabolite %in% c(up$metabolite,down$metabolite), 
                metabolite, ""))


df$sig <- factor(df$sig, levels = c("Up","Down","NS"))
View(df)
#expression(italic("Insr")^italic("-/-")~"vs"~italic("Insr")^italic("+/+"))

# Volcano 
ggplot(data=df, 
       aes(x=logFC.mean, y=-log10(adj.P.Val), 
           col= sig,
           label=label)) +
  scale_color_manual(values=c("firebrick", # nothing positively correlate with leucine or Glu+FA
                              "dodgerblue3",#"blue",
                              "grey"), #c("#F8766D", "#00BFC4", "grey"),
                     guide = "none"
  ) +
  ylab(expression(-log[10]~(Adj~P.value))) +
  xlab(expression(log[2](Fold~Change)))+
  geom_point(aes(),
             size= 1.5) +
  #  scale_color_gradient2(high = "red3",mid = "white",low = "blue3",
  #                        midpoint = 0,na.value = "grey80"
  #                        #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
  #  )+
  #  scale_size_continuous(range = c(0.1, 4))+
  #  theme_minimal() +
  theme_bw()+
  geom_text_repel(
    data=subset(df,logFC.mean>0), 
    color="firebrick",
    size=4.5,
    #box.padding = 2,
    #  force = 5,
    #  force_pull = 3,
    min.segment.length = 0,
    nudge_x = 4 - subset(df,logFC.mean>0)$logFC.mean,
    segment.size=0.3, segment.color="grey", 
    direction="y", 
    hjust=0,
    max.overlaps = 200)+
  
  geom_text_repel(
    data=subset(df,logFC.mean<0), 
    color="dodgerblue3",
    size=4.3, ######### used 4.3 for 180 min, used 4.5 for 90 min###
    #box.padding = 2,
    #  force = 5,
    #  force_pull = 10,
    min.segment.length = 0,
    nudge_x = -2 - subset(df,logFC.mean<0)$logFC.mean,
    segment.size=0.3, segment.color="grey", 
    direction="y", 
    hjust=1,
    max.overlaps = 300)+
  
  #labs(size="Abundance (log10)",color="CoV (log10)")+
  scale_x_continuous(#breaks = c(-1,-0.5,0,0.5,1), 
    limits = c(-4, 4)
  )+
  scale_y_continuous(#breaks = c(-1,-0.5,0,0.5,1), 
    limits = c(0, 5)
  )+
  # annotate(geom="text",x=4, y=0, 
  #           label=expression("Up in "*italic("Insr")^"f/f"),
  #           color="firebrick",size=5) +
  #  annotate(geom="text",x=-4, y=0, 
  #           label=expression("Down in "*italic("Insr")^"f/f"),
  #           color="steelblue3",size=5) +
  
  #xlim(-1.5,1.5) +
  #annotate(geom="text", x=2.1,y=-0.1, label="Up in T2D", color="#F8766D", size=7) +
  #annotate(geom="text", x=-2.1,y=-0.1, label="Up in ND", color="#00BFC4", size=7) +
  ggtitle(
    
    ########### change title here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #"Ketone vs Placebo 90 min"
    "Ketone vs Placebo 180 min"
  ) +
  theme(
    plot.title = element_text(color="black", size=18,hjust = 0.5),
    legend.position = "right",
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text=element_text(size=17),
    axis.text=element_text(size=17),
    axis.title=element_text(size=17),
    legend.text=element_text(size=17),
    legend.title=element_text(size=17),
    aspect.ratio = 1/1, panel.grid.major = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave(filename="figures/volcano_ketone.placebo_90min.png",width=6.5,height=6.8,units="in",dpi=600)
ggsave(filename="figures/volcano_ketone.placebo_90min.pdf",width=6.5,height=6.8,units="in")

ggsave(filename="figures/volcano_ketone.placebo_180min.png",width=6.5,height=6.8,units="in",dpi=600)
ggsave(filename="figures/volcano_ketone.placebo_180min.pdf",width=6.5,height=6.8,units="in")

# plot individual fold changes from baseline on PCA and heatmap ============
plasma.drink.time <- read.csv("output/plasma_drink_time.csv", check.names = F) %>%
  mutate_at(-1, as.numeric)

str(plasma.drink.time)
View(plasma.drink.time)
colnames(plasma.drink.time) # missing 2.90.B, 4.180.B, 6.180.A, 6.180.B

plasma.A0 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.0.A", colnames(plasma.drink.time)))))
View(plasma.A0)

plasma.A90 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.90.A", colnames(plasma.drink.time)))))

plasma.A180 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.180.A", colnames(plasma.drink.time)))))

plasma.B0 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.0.B", colnames(plasma.drink.time)))))
View(plasma.B0)

plasma.B90 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.90.B", colnames(plasma.drink.time)))))

plasma.B180 <- plasma.drink.time %>%
  select(all_of(c(grep("\\.180.B", colnames(plasma.drink.time)))))

plasma.fc.A90.0 <- plasma.A90/plasma.A0
View(plasma.fc.A90.0)
colnames(plasma.A90)
colnames(plasma.A0)

plasma.fc.A180.0 <- plasma.A180/(plasma.A0 %>% select(-'6.0.A'))
colnames(plasma.A180)
colnames((plasma.A0 %>% select(-'6.0.A')))

plasma.fc.B90.0 <- plasma.B90/(plasma.B0 %>% select(-'2.0.B'))
View(plasma.fc.B90.0)

plasma.fc.B180.0 <- plasma.B180/(plasma.B0 %>% select(-c('4.0.B','6.0.B')))
View(plasma.fc.B180.0)

plasma.fc <- cbind(plasma.fc.B90.0, plasma.fc.A90.0, plasma.fc.B180.0, plasma.fc.A180.0)
rownames(plasma.fc) <- plasma.drink.time$sample

## baseline fold change PCA -----------
#install.packages("ggfortify")

#df <- as.data.frame(t(xcms[,-c(1:23)]))
df <- as.data.frame(t(plasma.fc)) #%>% na.omit() # each row should be a sample
df[is.na(df)] <- 0
View(df)
str(df)
dim(df) # 98 x 184
#grep("Pre",row.names(df))
#grep("Post",row.names(df))
#df$group <- c(rep("Pre",19),rep("Post",19))
#dim(df)

df.pca <- prcomp(df, center = TRUE,scale. = TRUE) # each row should be a sample

summary(df.pca)
str(df.pca)
View(summary(df.pca))

pc <- as.data.frame(summary(df.pca)[["importance"]])
View(pc)
pc
pc1.v <- round(pc["Proportion of Variance","PC1"],digits=3)*100 # variance PC1
pc2.v <- round(pc["Proportion of Variance","PC2"],digits=3)*100 # variance PC2

cbPalette <- c(
  "#E69F00", #lightorange
  "#009E73", #green
  "#CC79A7", #magenta
  "#56B4E9", #blue
  "#0072B2", #darkblue
  "#D55E00", #darkorange
  "#999999", #grey
  "#F0E442" #yellow
)
rownames(df)
#df$group <- c(rep("Pre",19),rep("Post",19))
PCi<-data.frame(df.pca$x,group=c(rep("Placebo 90/0 min", 
                                     length(grep("\\.90.B", rownames(df)))),
                                 rep("Ketone 90/0 min", 
                                     length(grep("\\.90.A", rownames(df)))),
                                 rep("Placebo 180/0 min", 
                                     length(grep("\\.180.B", rownames(df)))),
                                 rep("Ketone 180/0 min", 
                                     length(grep("\\.180.A", rownames(df))))
                                 ))
PCi$group <- factor(PCi$group, levels = c("Placebo 90/0 min",
                                          "Placebo 180/0 min",
                                          "Ketone 90/0 min",
                                          "Ketone 180/0 min"))

View(PCi)

#autoplot(xcms.pca,data = df,colour = 'group')+

ggplot(PCi,aes(x=PC1,y=PC2,
               col=group #shape=time
))+
#  geom_text_repel(aes(label=row.names(PCi)),size=2,
#                  color="grey50",
#                  box.padding   = 0.4,
#                  point.padding = 0,
#                  #force=1,
#                  #force_pull=10,
#                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
#                  segment.color = 'darkgrey')+
  geom_point(size=2)+
  scale_color_manual(values=cbPalette) +
  labs(x = paste0("PC1 (",pc1.v,"%)"), 
       y = paste0("PC2 (",pc2.v,"%)"),
       color='Individual fold change')+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        #legend.title = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16)
        
  )+
  theme(aspect.ratio=1/1)

ggsave(filename="figures/PCA_plasma_foldchange.png",width=15,height=12,units="cm",dpi=400)
ggsave(filename="figures/PCA_plasma_foldchange.pdf",width=15,height=12,units="cm")

## heatmap individual fold changes ---------------
m <- plasma.fc

ketone.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.90min_subject.csv",
                        row.names = 1) %>%
  filter(adj.P.Val<0.05)
ketone.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_0min.180min_subject.csv",
                         row.names = 1) %>%
  filter(adj.P.Val<0.05)
#ketone.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ketone_90min.180min_subject.csv",
#                          row.names = 1) %>%
#  filter(adj.P.Val<0.05)
ctrl.0.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.90min_subject.csv",
                      row.names = 1) %>%
  filter(adj.P.Val<0.05)
ctrl.0.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_0min.180min_subject.csv",
                       row.names = 1) %>%
  filter(adj.P.Val<0.05)
#ctrl.90.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_ctrl_90min.180min_subject.csv",
#                        row.names = 1) %>%
#  filter(adj.P.Val<0.05)
ketone.ctrl.0 <- read.csv("input/MetaboAnalyst_results/covariate_result_0min_ctrl.ketone_subject.csv",
                          row.names = 1) %>%
  filter(adj.P.Val<0.05)
View(ketone.ctrl.0)
ketone.ctrl.90 <- read.csv("input/MetaboAnalyst_results/covariate_result_90min_ctrl.ketone_subject.csv",
                           row.names = 1)%>%
  filter(adj.P.Val<0.05)
ketone.ctrl.180 <- read.csv("input/MetaboAnalyst_results/covariate_result_180min_ctrl.ketone_subject.csv",
                            row.names = 1)%>%
  filter(adj.P.Val<0.05)

View(ketone.0.180)                      
View(plasma.fc)                           
plasma.fc.anno <- plasma.fc %>% mutate(metabolite=rownames(plasma.fc)) %>%
  filter(metabolite %in% c(rownames(ketone.0.90), rownames(ketone.0.180), rownames(ctrl.0.90), rownames(ctrl.0.180)))
plasma.fc.anno <- plasma.fc.anno %>% 
  mutate(`ketone_90.0`= ifelse(metabolite %in% rownames(ketone.0.90), "TRUE",""), #"Ketone 90 vs 0 min", ""),
         `ketone_180.0`= ifelse(metabolite %in% rownames(ketone.0.180), "TRUE",""), # "Ketone 180 vs 0 min", ""),
         `placebo_90.0`= ifelse(metabolite %in% rownames(ctrl.0.90), "TRUE",""), # "Placebo 90 vs 0 min", ""),
         `placebo_180.0`= ifelse(metabolite %in% rownames(ctrl.0.180), "TRUE",""), # "Placebo 180 vs 0 min", ""),
         )
View(plasma.fc.anno)

plasma.fc.anno.de <- plasma.fc.anno %>%
  filter(metabolite %in% c(rownames(ketone.ctrl.90), rownames(ketone.ctrl.180)))
View(plasma.fc.anno.de)

# all significant metabolites in comparison to baseline 0 min
m <- plasma.fc.anno[,1:64] 
m.anno <- plasma.fc.anno

# all significant metabolites in comparison to baseline 0 min, then filtered for significant ones at ketone vs placebo 90 & 180 min (some metabolites are similarly changes in ketone and placebo group)
m <- plasma.fc.anno.de[,1:64] 
m.anno <- plasma.fc.anno.de

#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)

m.z <- log2(m)

#m.z <- t(scale(t(m))) 
View(m.z)
colnames(m.z)
ceiling(max(abs(m.z)))
#m.z[m.z>5] <- NA
#m.z <- t(scale(t(m.z)))
#m.z[is.na(m.z)] <- 5
dim(m.z)
max(na.omit(m.z))
min(na.omit(m.z))

number_of_p90 <- length(grep("\\.90.B", colnames(m.z)))
number_of_p180 <- length(grep("\\.180.B", colnames(m.z)))
number_of_k90 <- length(grep("\\.90.A", colnames(m.z)))
number_of_k180 <- length(grep("\\.180.A", colnames(m.z)))

end_index_p90 <- grep("\\.90.B", colnames(m.z))[number_of_p90]
end_index_p180 <- grep("\\.180.B", colnames(m.z))[number_of_p180]
end_index_k90 <- grep("\\.90.A", colnames(m.z))[number_of_k90]
end_index_k180 <- grep("\\.180.A", colnames(m.z))[number_of_k180]

heatmap.all <- Heatmap(m.z, #matrix_Z_score_total,
                      name = "Z score",
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      show_row_dend = TRUE,
                      
                      row_labels = gt_render(m.anno$metabolite),
                      
                      row_names_gp = gpar(fontsize = 11),
                      #row_names_gp = gpar(fontsize = 9),
                      
                      column_names_gp = gpar(fontsize = 8),
                      
                      column_names_side = "bottom",
                      column_dend_side = "bottom",
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "ward.D2",
                      row_dend_side = "left",
                      row_dend_width = unit(5, "mm"),
                      layer_fun = function(j, i, x, y, width, height, fill) {
                        mat = restore_matrix(j, i, x, y)
                        ind = unique(c(mat[, c(
                          end_index_p90,
                          end_index_k90,
                          end_index_p180
                        )]))
                        grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                  y = y[ind], 
                                  width = unit(0.03, "inches"), 
                                  height = unit(1/nrow(m.z), "npc"),
                                  gp = gpar(col = "white", fill = "white")
                        )
                      },
                      col = colorRamp2(#c(-3,0,3), 
                        c(-5,0,5),
                        c("blue", "white", "red")),
                      top_annotation = columnAnnotation(
                        empty = anno_empty(
                          border = FALSE, 
                          height = unit(6, "mm")
                        )),
                      
                      column_order = 1:ncol(m.z),
                      height = #nrow(m.z)*unit(6, "mm"),
                        unit(180, "mm"), 
                      
                      width = ncol(m.z)*unit(2.5, "mm"),
                      border_gp = gpar(col = "black"),
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(
                        title = "log2 fold change",
                        title_position = "leftcenter-rot",
                        legend_height = unit(4, "cm")))


draw(heatmap.all)

heatmap.de <- Heatmap(m.z, #matrix_Z_score_total,
                      name = "Z score",
                      show_row_names = TRUE,
                      show_column_names = TRUE,
                      show_row_dend = TRUE,
                      
                      row_labels = gt_render(m.anno$metabolite),
                      
                      row_names_gp = gpar(fontsize = 11),
                      #row_names_gp = gpar(fontsize = 9),
                      
                      column_names_gp = gpar(fontsize = 8),
                      
                      column_names_side = "bottom",
                      column_dend_side = "bottom",
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "ward.D2",
                      row_dend_side = "left",
                      row_dend_width = unit(5, "mm"),
                      layer_fun = function(j, i, x, y, width, height, fill) {
                        mat = restore_matrix(j, i, x, y)
                        ind = unique(c(mat[, c(
                          end_index_p90,
                          end_index_k90,
                          end_index_p180
                        )]))
                        grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                  y = y[ind], 
                                  width = unit(0.03, "inches"), 
                                  height = unit(1/nrow(m.z), "npc"),
                                  gp = gpar(col = "white", fill = "white")
                        )
                      },
                      col = colorRamp2(#c(-3,0,3), 
                        c(-5,0,5),
                                       c("blue", "white", "red")),
                          top_annotation = columnAnnotation(
                            empty = anno_empty(
                              border = FALSE, 
                              height = unit(6, "mm")
                              )),
                      
                         left_annotation = 
                           rowAnnotation(`Placebo 90 vs 0 min` = m.anno$placebo_90.0,
                                         `Ketone 90 vs 0 min` = m.anno$ketone_90.0,
                                         `Placebo 180 vs 0 min` = m.anno$placebo_180.0,
                                         `Ketone 180 vs 0 min` = m.anno$ketone_180.0,
                                         
                                         col = list(`Placebo 90 vs 0 min` = c("TRUE" = alpha(cbPalette[1],0.5)),#hue_pal()(4)[3]),
                                                    `Ketone 90 vs 0 min` = c("TRUE" = alpha(cbPalette[3],0.5)),#hue_pal()(4)[1]),
                                                    `Placebo 180 vs 0 min` = c("TRUE" = alpha(cbPalette[2],0.5)),#hue_pal()(4)[4])
                                                    `Ketone 180 vs 0 min` = c("TRUE" = alpha(cbPalette[4],0.5))#hue_pal()(4)[2]),
                                         ),
#                                         annotation_legend_param = list(
#                      
#                                           `Significant CW vs CS` = list(
#                                             title = c(""),
#                                             labels = expression(italic("Ins1")^italic("+/+")*" Sucrose vs Water") #expression() issue - https://github.com/jokergoo/ComplexHeatmap/issues/678
#                                           ),
                      
#                                           `Significant KW vs KS` = list(
#                                             title = c(""),
#                                             labels = expression(italic("Ins1")^italic("+/-")*" Sucrose vs Water")
#                                           ),
                      
#                                           `Significant CW vs KW` = list(
#                                             title = c(""),
#                                             labels = expression("Water "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
#                                           ),
                      
#                                           `Significant CS vs KS` = list(
#                                             title = c(""),
#                                             labels = expression("Sucrose "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
#                                           )
#                                         ),
                                         show_annotation_name = FALSE,
                                         show_legend=FALSE,
                                         na_col = "white",
                                         width = unit(8, "mm")
                      
                           ),
                      
                      column_order = 1:ncol(m.z),
                      height = #nrow(m.z)*unit(6, "mm"),
                       unit(180, "mm"), 
                      
                      width = ncol(m.z)*unit(2.5, "mm"),
                      border_gp = gpar(col = "black"),
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(
                        title = "log2 fold change",
                        title_position = "topcenter",
                        legend_height = unit(4, "cm")))


draw(heatmap.de)



pdf(file = "figures/heatmap_foldchange_plasma_log2.pdf",
     width = 7.5, 
     height = 8)
draw(heatmap.all)

pdf(file = "figures/heatmap_foldchange_plasma_de_log2_anno.pdf",
    width = 11, 
    height = 8)

png(file = "figures/heatmap_foldchange_plasma_de_log2_anno.png",
    width = 11, 
    height = 8,
    units="in",
    res=600)

draw(heatmap.de)

pdf(file = "figures/heatmap_foldchange_plasma_de_log2_anno_baseline_ketone.placebo.pdf",
    width = 11, 
    height = 8)

png(file = "figures/heatmap_foldchange_plasma_de_log2_anno_baseline_ketone.placebo.png",
    width = 11, 
    height = 8,
    units="in",
    res=600)

draw(heatmap.de)

#
seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))


####Condition labels

#Condition label 1


grid.rect(x = (loc2$x - loc1$x)*(end_index_p90)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_p90)/ncol(m.z), 
          height = (loc2$y - loc1$y),
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[1],0.5),
                    col = alpha(cbPalette[1],0.5)
          )
)
grid.text("Placebo 90/0 min", 
          x = (loc2$x - loc1$x)*(end_index_p90)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2

grid.rect(x = (loc2$x - loc1$x)*(end_index_p90 + 
                                   end_index_k90)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_k90)/ncol(m.z), 
          height = (loc2$y - loc1$y),
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[3],0.5),
                    col = alpha(cbPalette[3],0.5)
          )
)

grid.text("Ketone 90/0 min", 
          x = (loc2$x - loc1$x)*(end_index_p90 + 
                                   end_index_k90)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

#Condition label 3

grid.rect(x = (loc2$x - loc1$x)*(end_index_k90 + 
                                   end_index_p180)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_p180)/ncol(m.z), 
          height = (loc2$y - loc1$y),
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[2],0.5),
                    col = alpha(cbPalette[2],0.5)
          )
)

grid.text("Placebo 180/0 min", 
          x = (loc2$x - loc1$x)*(end_index_k90 + 
                                   end_index_p180)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Condition label 4

grid.rect(x = (loc2$x - loc1$x)*(end_index_k180 + 
                                   end_index_p180)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_k180)/ncol(m.z), 
          height = (loc2$y - loc1$y),
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[4],0.5),
                    col = alpha(cbPalette[4],0.5)
          )
)

grid.text("Ketone 180/0 min", 
          x = (loc2$x - loc1$x)*(end_index_k180 + 
                                   end_index_p180)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Vertical lines
grid.rect(x = end_index_p90/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.5, "mm"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_k90/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.5, "mm"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_p180/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.5, "mm"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

dev.off()


# plot PCA and heatmap for PBMC ===============
pbmc <- read.csv("output/pbmc_drink_time.csv", check.names = F) %>%
  column_to_rownames("sample")
View(pbmc)

meta <- read.csv("output/traits_pbmc_drink_time.csv", check.names = F)
meta <- meta %>%
  mutate(Condition=ifelse(Condition=="A", "Ketone", "Placebo"),
         Time=ifelse(Time=="T0", "0min", ifelse(Time=="T90", "90min", "180min"))) %>%
  mutate(group=paste0(Condition," ",Time))
View(meta)

all(colnames(pbmc)==meta$sample)

df <- as.data.frame(t(pbmc)) # each row should be a sample
View(df)
str(df)

#grep("Pre",row.names(df))
#grep("Post",row.names(df))
#df$group <- c(rep("Pre",19),rep("Post",19))
#dim(df)

df[is.na(df)] <- 0
df.pca <- prcomp(df, center = TRUE,scale. = TRUE) # each row should be a sample

summary(df.pca)
str(df.pca)
View(summary(df.pca))

pc <- as.data.frame(summary(df.pca)[["importance"]])
View(pc)
pc
pc1.v <- round(pc["Proportion of Variance","PC1"],digits=3)*100 # variance PC1
pc2.v <- round(pc["Proportion of Variance","PC2"],digits=3)*100 # variance PC2

cbPalette <- c(
  "#E69F00", #lightorange
  "#009E73", #green
  "#CC79A7", #magenta
  "#0072B2", #darkblue
  "#56B4E9", #blue
  "#D55E00", #darkorange
  "#999999", #grey
  "#F0E442" #yellow
)
View(df.pca$x)
PCi<-data.frame(df.pca$x,group=meta$group,time=meta$Time,condition=meta$Condition, sex=meta$gender)
PCi$group <- factor(PCi$group, levels = unique(PCi$group)[c(4,5,6,1,2,3)])
PCi$time <- factor(PCi$time, levels = c("0min","90min","180min"))
PCi$condition <- factor(PCi$condition, levels = c("Ketone","Placebo"))
View(PCi)

#autoplot(xcms.pca,data = df,colour = 'group')+

ggplot(PCi,aes(x=PC1,y=PC2,
               col=sex
               #col=condition,shape=time
))+
#  geom_text_repel(aes(label=row.names(PCi)),size=2,
#                  color="grey50",
#                  box.padding   = 0.4,
#                  point.padding = 0,
                  #force=1,
                  #force_pull=10,
#                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
#                  segment.color = 'darkgrey')+
  geom_point(size=2)+
  
##############  
#  scale_color_manual(values=hue_pal()(9)[c(1,2,3,8,6,7)]) + #cbPalette
  labs(x = paste0("PC1 (",pc1.v,"%)"), y = paste0("PC1 (",pc2.v,"%)"),
       color="")+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        #legend.title = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 16)
        
  )+
  theme(aspect.ratio=1/1)
ggsave(filename="figures/PCA_PBMC_metabolomics_label.pdf",width=15,height=12,units="cm")
ggsave(filename="figures/PCA_PBMC_metabolomics_group.pdf",width=15,height=12,units="cm")
ggsave(filename="figures/PCA_PBMC_metabolomics_sex.pdf",width=15,height=12,units="cm")

ggsave(filename="figures/PCA_time_condition.png",width=15,height=12,units="cm",dpi=400)

# heatmap PBMC --------------
m <- pbmc
m <- m[,c(grep("\\.0.B",colnames(m)),
          grep("\\.0.A",colnames(m)),
          grep("\\.90.B",colnames(m)),
          grep("\\.90.A",colnames(m)),
          grep("\\.180.B",colnames(m)),
          grep("\\.180.A",colnames(m)))]

#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)

m.z <- t(scale(t(m))) 
View(m.z)
colnames(m.z)
ceiling(max(abs(m.z)))
#m.z[m.z>5] <- NA
#m.z <- t(scale(t(m.z)))
#m.z[is.na(m.z)] <- 5
dim(m.z)
max(na.omit(m.z))
min(na.omit(m.z))

number_of_p0 <- length(grep("\\.0.B", colnames(m.z)))
number_of_p90 <- length(grep("\\.90.B", colnames(m.z)))
number_of_p180 <- length(grep("\\.180.B", colnames(m.z)))
number_of_k0 <- length(grep("\\.0.A", colnames(m.z)))
number_of_k90 <- length(grep("\\.90.A", colnames(m.z)))
number_of_k180 <- length(grep("\\.180.A", colnames(m.z)))

end_index_p0 <- grep("\\.0.B", colnames(m.z))[number_of_p0]
end_index_k0 <- grep("\\.0.A", colnames(m.z))[number_of_k0]
end_index_p90 <- grep("\\.90.B", colnames(m.z))[number_of_p90]
end_index_k90 <- grep("\\.90.A", colnames(m.z))[number_of_k90]
end_index_p180 <- grep("\\.180.B", colnames(m.z))[number_of_p180]
end_index_k180 <- grep("\\.180.A", colnames(m.z))[number_of_k180]

heatmap.all <- Heatmap(m.z, #matrix_Z_score_total,
                      name = "Z score",
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      show_row_dend = TRUE,
                      
                      #row_labels = gt_render(m.anno$metabolite),
                      
                      row_names_gp = gpar(fontsize = 11),
                      #row_names_gp = gpar(fontsize = 9),
                      
                      column_names_gp = gpar(fontsize = 8),
                      
                      column_names_side = "bottom",
                      column_dend_side = "bottom",
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "ward.D2",
                      row_dend_side = "left",
                      row_dend_width = unit(5, "mm"),
                      layer_fun = function(j, i, x, y, width, height, fill) {
                        mat = restore_matrix(j, i, x, y)
                        ind = unique(c(mat[, c(
                          end_index_p0,
                          end_index_k0,
                          end_index_p90,
                          end_index_k90,
                          end_index_p180
                        )]))
                        grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                  y = y[ind], 
                                  width = unit(0.03, "inches"), 
                                  height = unit(1/nrow(m.z), "npc"),
                                  gp = gpar(col = "white", fill = "white")
                        )
                      },
                      col = colorRamp2(#c(-3,0,3), 
                        c(-5,0,5),
                        c("blue", "white", "red")),
                      top_annotation = columnAnnotation(
                        empty = anno_empty(
                          border = FALSE, 
                          height = unit(6, "mm")
                        )),
                      
                     
                      column_order = 1:ncol(m.z),
                      height = #nrow(m.z)*unit(6, "mm"),
                        unit(180, "mm"), 
                      
                      width = ncol(m.z)*unit(2.5, "mm"),
                      border_gp = gpar(col = "black"),
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(
                        title = "z-score",
                        title_position = "topleft",
                        legend_height = unit(4, "cm")))


pdf(file = "figures/heatmap_all_PBMC.pdf",
    width = 6, 
    height = 8)


draw(heatmap.all)

#
seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))


####Condition labels
grid.text("Placebo 0", 
          x = (loc2$x - loc1$x)*(end_index_p0)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))


#Condition label 2

grid.text("Ketone 0", 
          x = (loc2$x - loc1$x)*(end_index_p0 + end_index_k0)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))

#Condition label 3

grid.text("Placebo 90", 
          x = (loc2$x - loc1$x)*(end_index_p90 + end_index_k0)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))

#Condition label 4

grid.text("Ketone 90", 
          x = (loc2$x - loc1$x)*(end_index_p90 + end_index_k90)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))

#Condition label 5

grid.text("Placebo 180", 
          x = (loc2$x - loc1$x)*(end_index_p180 + end_index_k90)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))

#Condition label 6

grid.text("Ketone 180", 
          x = (loc2$x - loc1$x)*(end_index_p180 + end_index_k180)/2/ncol(m.z),
          y = 0.5,
          just = c("center", "center"),
          gp = gpar(fontsize = 9,
                    col = "black"))


dev.off()

# MetaboAnalyst R (not fully successful) =========================
#Function to download packages
metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# install RTools: https://cran.r-project.org/bin/windows/Rtools/

# had issues with rlang, installed new version
#install.packages("rlang")
#packageVersion("rlang")

# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

# Load the MetaboAnalystR package
library("MetaboAnalystR")


# upload data as concentration, multiple factors/covariate, (samples in columns), paired is false (will set subject as covariate later)

#rm(mSet)
mSet<-InitDataObjects("conc", "mf", FALSE)
?InitDataObjects
#mSet<-SetDesignType(mSet, "time")
mSet<-SetDesignType(mSet, "multi")
?SetDesignType

mSet<-Read.TextData(mSet, "output/plasma_drink_180min.csv", "colmf" #"colts", "cont"
                  )
mSet<-ReadMetaData(mSet, "output/traits_drink_180min.csv")
??ReadMetaData
View(mSet)

# data integrity check

mSet<-SanityCheckData(mSet)

## missing values are replace by LoDs (1/5 of the minimum positive value of each variable)
## remove features with >50 missing values
#mSet<-ReplaceMin(mSet) #replace zero/missing values by half of the smallest positive value in the original dataset
??ReplaceMin
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeMissingVar(mSet, method="min") #Replace missing variables by min/mean/median/KNN/BPCA/PPCA/svdImpute.
??ImputeMissingVar

mSet<-SanityCheckMeta(mSet, 1)

## changed subject form continuous to categorical
## default - meta data with missing value will be removed

mSet<-SetDataTypeOfMeta(mSet)

# data filtering 
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, 
                     filter="none", 
                     -1,  # what's this?
                     qcFilter="F", 
                     rsd=25, # won't be used in function unless qcFilter=T
                     F # what's this
                     )
??FilterVariable

# Normalization

##
mSet<-PreparePrenormData(mSet)
## Sample normlaization - None
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
## sample normalization - median
mSet<-Normalization(mSet, 
                    "MedianNorm", # sample normalization
                    "LogNorm", # data transformation
                    "NULL",  # data scaling
                    ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)


adj.vec <<- [Subject]
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_0_dpi72.png", "png", "Time", "0", "batch" , 0.05,"anova")
??CovariateScatter.Anal
adj.vec <<- [Subject]
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_4_dpi72.png", "png", "Time", "0", "batch" , 0.05,"90")
Xia Lab @ McGill   (last updated 2023-03-24)
# PID of current job: 2187579
mSet<-InitDataObjects("conc", "mf", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, "output/plasma_upload.csv", "colmf");
mSet<-ReadMetaData(mSet, "output/traits_upload.csv");
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeMissingVar(mSet, method="min")
mSet<-SanityCheckMeta(mSet, 1)
mSet<-SetDataTypeOfMeta(mSet)

mSet<-FilterVariable(mSet, "none", "F", 25)
?FilterVariable
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_0_", "png", 72, width=NA, 5, "ID", "Time")
mSet<-iPCA.Anal(mSet, "ipca_3d_0_.json")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_1_", "png", 72, width=NA, 5, "Time", "Condition")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_2_", "png", 72, width=NA, 2, "null", "null")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_3_", "png", 72, width=NA, 5, "Time", "Condition")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_4_", "png", 72, width=NA, 5, "Condition", "Time")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_4_", "png", 300, width=NA, 5, "Condition", "Time")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_4_", "pdf", 72, width=NA, 5, "Condition", "Time")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_5_", "png", 72, width=NA, 5, "Time", "Condition")
mSet<-PlotPCAPairSummaryMeta(mSet, "pca_pair_meta_5_", "pdf", 72, width=NA, 5, "Time", "Condition")
meta.vec.hm2 <- [ID, Time, Condition, batch]

mSet<-PlotHeatMap2(mSet, "heatmap2_0_", "norm", "row", "png", 72, width=NA, "euclidean","ward.D","bwm", 8, "overview","mean",2000, F, F,  F, T, T, F)
meta.vec.hm2 <- [ID, Time, Condition, batch]
mSet<-PlotHeatMap2(mSet, "heatmap2_1_", "raw", "row", "png", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, T)
meta.vec.hm2 <- [ID, Time, Condition, batch] 
mSet<-PlotHeatMap2(mSet, "heatmap2_2_", "raw", "row", "png", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, F)
meta.vec.hm2 <- [ID, Time, Condition, batch]
mSet<-PlotHeatMap2(mSet, "heatmap2_3_", "raw", "row", "png", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, F)
mSet<-PlotHeatMap2(mSet, "heatmap2_3_", "raw", "row", "pdf", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, F)
meta.vec.hm2 <- [ID, Time, Condition, batch]
mSet<-PlotHeatMap2(mSet, "heatmap2_4_", "raw", "row", "png", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, F)

mSet<-PlotHeatMap2(mSet, "heatmap2_4_", "raw", "row", "pdf", 72, width=NA, "euclidean","ward.D","bwm", 10, "overview","mean",25, F, F,  F, T, T, F)

adj.vec <<- [ID, Time, gender, age, weight (kg), BMI (kg/m2), A1C, ketones (g), medications]

mSet<-CovariateScatter.Anal(mSet, "covariate_plot_0_dpi72.png", "png", 72, "default", "Condition", "B", "batch" , 0.05, FALSE)

adj.vec <<- [ID, Time, gender, age, weight (kg), BMI (kg/m2), A1C]

mSet<-CovariateScatter.Anal(mSet, "covariate_plot_1_dpi72.png", "png", 72, "default", "Condition", "B", "batch" , 0.05, FALSE)
adj.vec <<- [Time, gender, age, weight (kg), BMI (kg/m2), A1C]
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_2_dpi72.png", "png", 72, "default", "Condition", "B", "batch" , 0.05, FALSE)
adj.vec <<- [Time, gender, age]
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_3_dpi72.png", "png", 72, "default", "Condition", "B", "batch" , 0.05, FALSE)
mSet<-PlotCmpdSummary(mSet, "beta-Hydroxybutyric acid","NA", 0, "png", 72, width=NA)
mSet<-PlotCmpdSummary(mSet, "beta-Hydroxybutyric acid","Condition", 0, "png", 72, width=NA)
adj.vec <<- [Time]
mSet<-CovariateScatter.Anal(mSet, "covariate_plot_4_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_5_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_6_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age, weight (kg)]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_7_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age, BMI (kg/m2)]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_8_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_9_dpi72.png", "png", 72, "default", "Condition", "A", "batch" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, gender, age]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_10_dpi72.png", "png", 72, "default", "Condition", "A", "NA" , 0.05, FALSE)
                                                                                     adj.vec <<- [Time, batch, gender, age]
                                                                                     mSet<-CovariateScatter.Anal(mSet, "covariate_plot_11_dpi72.png", "png", 72, "default", "Condition", "A", "NA" , 0.05, FALSE)
                                                                                     mSet<-PlotCovariateMap(mSet, "ggplot", "covariate_plot_11_dpi72.pdf", "pdf", 72)
                                                                                     mSet<-PlotCovariateMap(mSet, "default", "covariate_plot_11_dpi72.pdf", "pdf", 72)
                                                                                     meta.vec.aov <- [Time, Condition]
                                                                                     mSet<-ANOVA2.Anal(mSet, 0.05, "fdr", "multi", "between")
                                                                                     mSet<-PlotANOVA2(mSet, "aov2_1_", "png", 72, width=NA)
                                                                                     mSet<-PlotANOVA2(mSet, "aov2_1_", "pdf", 72, width=NA)
                                                                                     meta.vec.asca <- [Time, Condition]
                                                                                     mSet<-Perform.ASCA(mSet, 1, 1, 2, 2)
                                                                                     mSet<-PlotModelScree(mSet, "asca_scree_0_", "png", 72, width=NA)
                                                                                     mSet<-CalculateImpVarCutoff(mSet, 0.05, 0.9)
                                                                                     mSet<-PlotAscaImpVar(mSet, "asca_impa_0_", "png", 72, width=NA, "a")
                                                                                     mSet<-PlotAscaImpVar(mSet, "asca_impb_0_", "png", 72, width=NA, "b")
                                                                                     mSet<-PlotAscaImpVar(mSet, "asca_impab_0_", "png", 72, width=NA, "ab")
                                                                                     mSet<-PlotASCAModel(mSet, "asca_fa_0_", "png", 72, width=NA, "a",FALSE)
                                                                                     mSet<-PlotASCAModel(mSet, "asca_fb_0_", "png", 72, width=NA, "b",FALSE)
                                                                                     mSet<-PlotInteraction(mSet, "asca_fab_0_", "png", 72,FALSE, width=NA)
                                                                                     meta.vec.rf <- [Time]
                                                                                     mSet<-RF.AnalMeta(mSet, 500,7,1, "Condition")
                                                                                     mSet<-PlotRF.ClassifyMeta(mSet, "rf_cls_1_", "png", 72, width=NA)
                                                                                     mSet<-PlotRF.VIPMeta(mSet, "rf_imp_1_", "png", 72, width=NA)
                                                                                     mSet<-PlotRF.Outlier(mSet, "rf_outlier_1_", "png", 72, width=NA)
                                                                                     mSet<-SaveTransformedData(mSet)
                                                                                     