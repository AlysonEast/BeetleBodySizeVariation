library(ggplot2)
library(dplyr)
library(ggplot2)
library(sf)
library(maps)
library(dplyr)
library(ggrepel)

getwd()
setwd("/home/aly/Beetles/BeetleBodySizeVariation")
df<-read.csv("./Data/BeetleMeasurements.csv")
meta<-read.csv("./Data/NEON_Field_Site_Metadata_20260130.csv")
puum<-read.csv("./Data/trait_annotations.csv")

df<-merge(df, meta, by.x="siteID", by.y="site_id", all.x=TRUE)

head(df)

df_sub<-subset(df, user_name=="IsaFluck")

head(df_sub)
str(df_sub)

ElytraLength<-subset(df_sub, structure=="ElytraLength")

hist(ElytraLength$dist_cm)

site_spp<-table(ElytraLength$species, ElytraLength$siteID)
site_spp<-as.data.frame(site_spp)

ElytraLength<-subset(ElytraLength, siteID!="DELA" & siteID!="DSNY")

#### Explore the data ####
# Get unique siteID levels ordered by latitude
site_order <- unique(ElytraLength[order(ElytraLength$latitude, decreasing = TRUE), "siteID"])

# Convert siteID into a factor with the correct order
ElytraLength$siteID <- factor(ElytraLength$siteID, levels = site_order)
ElytraLength$log_dist_cm<-log10(ElytraLength$dist_cm)

library(ggpubr)
#### 
ElytraLength$Spp_plot<-paste0(ElytraLength$scientificName, ElytraLength$plotID)
forPlot<-subset(ElytraLength, domain_id=="D01")
forPlot<-subset(forPlot, scientificName=="Carabus goryi"|
                  scientificName=="Cymindis neglecta"|
                  scientificName=="Pterostichus pensylvanicus"|
                  scientificName=="Pterostichus tristis"|
                  scientificName=="Sphaeroderus stenostomus"|
                  scientificName=="Synuchus impunctatus")
forPlot$scientificName <- factor(
  forPlot$scientificName,
  levels = c("Cymindis neglecta",
             "Pterostichus pensylvanicus",
             "Pterostichus tristis",
             "Synuchus impunctatus",
             "Sphaeroderus stenostomus",
             "Carabus goryi"), ordered = TRUE)
table(forPlot$Spp_plot)
ggarrange(nrow=1,
          ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
            geom_density(alpha=0.5) +
            scale_fill_manual(values=c(rep("#028A9B",3), 
                                       rep("#BFA194",5), 
                                       rep("#F8CF9B",4),
                                       rep("#757943",1),
                                       rep("#DE5A5A",1),
                                       rep("#BF3F76",16))) +
            #xlim(0.4,1.6) +
            facet_wrap(.~scientificName, ncol = 1) +
            theme_pubr() +
            theme(legend.position="none"),
          ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot, colour = Spp_plot)) +
            geom_density(alpha=.5, position="stack") +
            scale_fill_manual(values=c(rep("#028A9B",3), 
                                       rep("#BFA194",5), 
                                       rep("#F8CF9B",4),
                                       rep("#757943",1),
                                       rep("#DE5A5A",1),
                                       rep("#BF3F76",16))) +
            scale_color_manual(values=c("black",rep("#028A9B",2), 
                                        "black",rep("#BFA194",4), 
                                        "black",rep("#F8CF9B",3),
                                        rep("black",1),
                                        rep("black",1),
                                        "black",rep("#BF3F76",15))) +
            #xlim(0.4,1.6) +
            facet_wrap(.~scientificName, scales = "free_y", ncol = 1) +
            theme_pubr() +
            theme(legend.position="none"),
          ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
            geom_density(alpha=.5, position="stack", col="#00000000") +
            scale_fill_manual(values=c(rep("#028A9B",3), 
                                       rep("#BFA194",5), 
                                       rep("#F8CF9B",4),
                                       rep("#757943",1),
                                       rep("#DE5A5A",1),
                                       rep("#BF3F76",16))) +
            scale_color_manual(values=c(rep("#028A9B",3), 
                                        rep("#BFA194",5), 
                                        rep("#F8CF9B",4),
                                        rep("#757943",1),
                                        rep("#DE5A5A",1),
                                        rep("#BF3F76",16))) +
            #xlim(0.4,1.6) +
            facet_wrap(.~siteID, ncol = 1) +
            theme_pubr() +
            theme(legend.position="none")
)

#D00###
forPlot<-subset(ElytraLength, domain_id=="D10")
table(forPlot$scientificName, forPlot$siteID)[,12:13]
table(forPlot$scientificName, forPlot$plotID)
forPlot$scientificName <- factor(
  forPlot$scientificName,
  levels = c("Discoderus parallelus",
             "Euryderus grossus",
             "Axinopalpus biplagiatus",
             "Selenophorus planipennis",
             "Pasimachus elongatus",
             "Cicindela punctulata",
             "Cratacanthus dubius",
             "Harpalus paratus",
             "Harpalus caliginosus",
             "Anisodactylus rusticus",
             "Harpalus desertus",
             "Amara carinata",
             "Cyclotrachelus torvus",
             "Harpalus pensylvanicus",
             "Poecilus scitulus",
             "Pterostichus protractus",
             "Calathus advena",
             "Pterostichus restrictus"), ordered = TRUE)
table(forPlot$Spp_plot)

TPD_plot_Spp<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
  geom_density(alpha=0.5) +
  geom_rug(alpha= 0.9, linewidth = 0.2, length = unit(.4, "npc"), aes(colour = Spp_plot))+
  # geom_vline(aes(xintercept = dist_cm)) +
  scale_fill_manual(values=c(rep("#1B9E77",3),
                             rep("#D95F02",2),
                             rep("#7570B3",1),
                             rep("#E7298A",4),
                             rep("#66A61E",1),
                             rep("#E6AB02",9),
                             rep("#A6761D",1),
                             rep("#666666",3),
                             rep("#1F78B4",2),
                             rep("#B2DF8A",2),
                             rep("#FB9A99",1),
                             rep("#CAB2D6",3),
                             rep("#FDBF6F",3),
                             rep("#6A3D9A",4),
                             rep("#B15928",1),
                             rep("#17BECF",4),
                             rep("#9B2F5D",2),
                             rep("#2E8B57",2))) +
  scale_color_manual(values=c(rep("#1B9E77",3),
                              rep("#D95F02",2),
                              rep("#7570B3",1),
                              rep("#E7298A",4),
                              rep("#66A61E",1),
                              rep("#E6AB02",9),
                              rep("#A6761D",1),
                              rep("#666666",3),
                              rep("#1F78B4",2),
                              rep("#66A61E",2),
                              rep("#FB9A99",1),
                              rep("#CAB2D6",3),
                              rep("#FDBF6F",3),
                              rep("#6A3D9A",4),
                              rep("#B15928",1),
                              rep("#17BECF",4),
                              rep("#9B2F5D",2),
                              rep("#2E8B57",2))) +
  #xlim(0.4,1.6) +
  facet_wrap(.~scientificName, ncol = 1, scales = "free_y") +
  theme_pubr() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

TPD_Spp<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot, colour = Spp_plot)) +
  geom_density(alpha=.5, position="stack") +
  scale_fill_manual(values=c(rep("#1B9E77",3),
                             rep("#D95F02",2),
                             rep("#7570B3",1),
                             rep("#E7298A",4),
                             rep("#66A61E",1),
                             rep("#E6AB02",9),
                             rep("#A6761D",1),
                             rep("#666666",3),
                             rep("#1F78B4",2),
                             rep("#B2DF8A",2),
                             rep("#FB9A99",1),
                             rep("#CAB2D6",3),
                             rep("#FDBF6F",3),
                             rep("#6A3D9A",4),
                             rep("#B15928",1),
                             rep("#17BECF",4),
                             rep("#9B2F5D",2),
                             rep("#2E8B57",2))) +
  scale_color_manual(values=c("black",rep("#1B9E77",2),
                              "black",rep("#D95F02",1),
                              rep("black",1),
                              "black",rep("#E7298A",3),
                              rep("black",1),
                              "black",rep("#E6AB02",8),
                              rep("black",1),
                              "black",rep("#666666",2),
                              "black",rep("#1F78B4",1),
                              "black",rep("#B2DF8A",1),
                              rep("black",1),
                              "black",rep("#CAB2D6",2),
                              "black",rep("#FDBF6F",2),
                              "black",rep("#6A3D9A",3),
                              rep("black",1),
                              "black",rep("#17BECF",3),
                              "black",rep("#9B2F5D",1),
                              "black",rep("#2E8B57",1))) +
  #xlim(0.4,1.6) +
  facet_wrap(.~scientificName, scales = "free_y", ncol = 1) +
  theme_pubr() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

forPlot$plotID
forPlot$plotID <- factor(
  forPlot$plotID,
  levels = c("CPER_008", "CPER_003", "CPER_004", "CPER_006", "CPER_002", "CPER_009",
             "STER_033", "STER_026", "STER_028", "STER_034", "STER_035", "STER_027", 
             "STER_031", "STER_006", "STER_029", "STER_032", 
             "RMNP_016", "RMNP_002", "RMNP_003", "RMNP_007", "RMNP_008", "RMNP_012", "RMNP_014"), 
  ordered = TRUE)

TPD_plot<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
  geom_density(position="stack", linewidth = 1.5, aes(x=dist_cm, colour = plotID))  +
  geom_density(position="stack", col="#00000000")  +
  scale_fill_manual(values=c(rep("#8DCEBB",3),
                             rep("#ECAF80",2),
                             rep("#BAB8D9",1),
                             rep("#F394C4",4),
                             rep("#B2D28E",1),
                             rep("#F2D580",9),
                             rep("#D2BA8E",1),
                             rep("#B2B2B2",3),
                             rep("#8FBCDA",2),
                             rep("#D8EFC4",2),
                             rep("#FDCCCC",1),
                             rep("#E4D8EA",3),
                             rep("#FEDFB7",3),
                             rep("#B49ECC",4),
                             rep("#D8AC94",1),
                             rep("#8BDEE7",4),
                             rep("#CD97AE",2),
                             rep("#96C5AB",2))) +
  scale_color_manual(values=c("#4C78A8",
                             "#F58518",    
                             "#54A24B",
                             "#E45756",
                             "#72B7B2",
                             "#EECA3B",
                             "#B279A2",
                             "#FF9DA6",
                             "#9D755D",
                             "#BAB0AC",
                             "#2F4B7C",
                             "#A05195",
                             "#D45087",
                             "#F95D6A",
                             "#FF7C43",
                             "#FFA600",
                             "#665191",
                             "#00876C",
                             "#7A5195",
                             "#EF5675",
                             "#003F5C",
                             "#BC5090",
                             "#FFA600")) +
  
  #xlim(0.4,1.6) +
  facet_wrap(.~plotID, scale  = "free_y", ncol = 1) +
  theme_pubr() +
  xlab("Elytra Length (cm)") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.8)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

forPlot$plotSpp<-paste0(forPlot$plotID,forPlot$scientificName)
table(forPlot$plotSpp)
dim(table(forPlot$plotSpp))
forPlot$plotSpp<-sort(forPlot$plotSpp, decreasing = FALSE)

TPD_site<-ggplot(data = forPlot, aes(x=dist_cm, fill = plotID)) +
  geom_density(position="stack", linewidth = 3, alpha=0, aes(x=dist_cm, fill = Spp_plot, colour = siteID))  +
  geom_density(position="stack", linewidth = 0, aes(x=dist_cm, colour = Spp_plot), fill="#FFFFFF")  +
  geom_density(position="stack", linewidth = 0, alpha=0.7, aes(x=dist_cm, colour = Spp_plot))  +
  # geom_density(position="stack", col="#00000000")  +
  scale_fill_manual(values=c("#4C78A8",
                              "#F58518",    
                              "#54A24B",
                              "#E45756",
                              "#72B7B2",
                              "#EECA3B",
                              "#B279A2",
                              "#FF9DA6",
                              "#9D755D",
                              "#BAB0AC",
                              "#2F4B7C",
                              "#A05195",
                              "#D45087",
                              "#F95D6A",
                              "#FF7C43",
                              "#FFA600",
                              "#665191",
                              "#00876C",
                              "#7A5195",
                              "#EF5675",
                              "#003F5C",
                              "#BC5090",
                              "#FFA600", rep("#00000000",48))) +
  # scale_color_manual(values=c("black",rep("#00000000",25),"black", rep("#00000000",10+25))) +
  # scale_fill_manual(values=c(rep("#8DCEBB",3),
  #                            rep("#ECAF80",2),
  #                            rep("#BAB8D9",1),
  #                            rep("#F394C4",4),
  #                            rep("#B2D28E",1),
  #                            rep("#F2D580",9),
  #                            
  #                            rep("#B2B2B2",3),
  #                            rep("#8FBCDA",2),
  #                            rep("#D8EFC4",2),
  #                            rep("#FDCCCC",1),
  #                            rep("#E4D8EA",3),
  #                            rep("#FEDFB7",3),
  #                            rep("#B49ECC",4),
  #                            rep("#D8AC94",1),
  #                            rep("#8BDEE7",4),
  #                            rep("#CD97AE",2),
  #                            rep("#96C5AB",2))) +
  # scale_fill_manual(values=c(rep("#8DCEBB",3),
  #                            rep("#ECAF80",2),
  #                            "#E45756",
  #                            rep("#F394C4",4),
  #                            "#B279A2",
  #                            rep("#F2D580",9),
  #                            "#F95D6A", #rep("#D2BA8E",1),
  #                            c("#4C78A8","#F58518","#54A24B"),
  #                            "#F58518","#54A24B",
  #                            "#F95D6A","#B279A2",
  #                            "#A05195",#rep("#FDCCCC",1),
  #                            "#B279A2","#FF9DA6","#F95D6A",
  #                            "#F95D6A","#FF7C43","#FFA600",#rep("#FEDFB7",3),
  #                            "#E45756","#72B7B2","#EECA3B","#B279A2",
  #                            "#A05195",#rep("#D8AC94",1),
  #                            rep("#8BDEE7",4),
  #                            rep("#CD97AE",2),
  #                            "#F58518","#EECA3B")) +
  scale_color_manual(values=c(rep("#AA3377",1),
                             rep("#CCBB44",1),
                             rep("#4477AA",1),
                             rep("#000000",52))) +
  # #xlim(0.4,1.6) +
  facet_wrap(.~siteID, ncol = 1) +
  theme_pubr() +
  xlab("Elytra Length (cm)") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

png("./Figures/TPDexampleV2.png", units = "in", res = 300, height = 9, width = 15)
annotate_figure(
  ggarrange(nrow=1,
            TPD_plot_Spp, TPD_Spp, TPD_plot, TPD_site),
  bottom = text_grob("Elytra Length (cm)", color = "black", size = 14))
dev.off()

forPlot$plotSpp<-paste0(forPlot$plotID,forPlot$scientificName)
table(forPlot$plotSpp)
forPlot$plotSpp<-sort(forPlot$plotSpp, decreasing = FALSE)

png("./Figures/TPDexampleDomain.png", units = "in", res = 300, height = 4, width = 4)
ggplot(data = forPlot, aes(x=dist_cm, fill = plotSpp, color=plotSpp)) +
  geom_density(position="stack", alpha = 0.5, aes(x=dist_cm))  +
  scale_fill_manual(values=c(rep("#AA3377",12),
                             rep("#CCBB44",10),
                             rep("#4477AA",26))) +
  scale_color_manual(values=c("black",rep("#00000000",(12+10+25)))) +
  #xlim(0.4,1.6) +
  theme_pubr() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )
dev.off()

#D10, 2 sites
forPlot<-subset(forPlot, siteID!="RMNP")
forPlot<-subset(forPlot, plotID!="STER_028" & plotID!="STER_034")

table(forPlot$scientificName, forPlot$siteID)[,12:13]
table(forPlot$scientificName, forPlot$plotID)
forPlot$scientificName <- factor(
  forPlot$scientificName,
  levels = c("Discoderus parallelus",
             "Euryderus grossus",
             "Axinopalpus biplagiatus",
             "Selenophorus planipennis",
             "Pasimachus elongatus",
             "Cicindela punctulata",
             "Harpalus paratus",
             "Harpalus caliginosus",
             "Cratacanthus dubius",
             "Anisodactylus rusticus",
             "Amara carinata",
             "Cyclotrachelus torvus",
             "Harpalus pensylvanicus",
             "Poecilus scitulus",
             "Pterostichus protractus",
             "Calathus advena",
             "Pterostichus restrictus",
             "Harpalus desertus"), ordered = TRUE)
table(forPlot$Spp_plot)

TPD_plot_Spp<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
  geom_density(alpha=0, linewidth=0) +
  geom_rug(alpha= 0.9, linewidth = 0.2, length = unit(.7, "npc"), aes(colour = Spp_plot))+
  # geom_vline(aes(xintercept = dist_cm)) +
  scale_fill_manual(values=c(rep("#1B9E77",3),
                             rep("#D95F02",3),
                             rep("#7570B3",1),
                             rep("#1F78B4",7),
                             rep("#66A61E",1),
                             rep("#E6AB02",3),
                             rep("#A6761D",2),
                             rep("#666666",2),
                             rep("#E7298A",1),
                             rep("#9B2F5D",3),
                             rep("#FB9A99",3),
                             rep("#CAB2D6",4),
                             rep("#17BECF",3))) +
  scale_color_manual(values=c(rep("#1B9E77",3),
                              rep("#D95F02",3),
                              rep("#7570B3",1),
                              rep("#1F78B4",7),
                              rep("#66A61E",1),
                              rep("#E6AB02",3),
                              rep("#A6761D",2),
                              rep("#666666",2),
                              rep("#E7298A",1),
                              rep("#9B2F5D",3),
                              rep("#FB9A99",3),
                              rep("#CAB2D6",4),
                              rep("#17BECF",3))) +
  #xlim(0.4,1.6) +
  facet_wrap(.~scientificName, ncol = 1, scales = "free_y") +
  theme_pubr() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

TPD_Spp<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot, colour = Spp_plot)) +
  geom_density(alpha=.5, position="stack") +
  scale_fill_manual(values=c(rep("#1B9E77",3),
                           rep("#D95F02",3),
                           rep("#7570B3",1),
                           rep("#1F78B4",7),
                           rep("#66A61E",1),
                           rep("#E6AB02",3),
                           rep("#A6761D",2),
                           rep("#666666",2),
                           rep("#E7298A",1),
                           rep("#9B2F5D",3),
                           rep("#FB9A99",3),
                           rep("#CAB2D6",4),
                           rep("#17BECF",3))) +
  scale_color_manual(values=c("black",rep("#1B9E77",2),
                              "black",rep("#D95F02",2),
                              "black",
                              "black",rep("#1F78B4",6),
                              "black",
                              "black",rep("#E6AB02",2),
                              "black",rep("#A6761D",1),
                              "black",rep("#666666",1),
                              "black",
                              "black",rep("#9B2F5D",2),
                              "black",rep("#FB9A99",2),
                              "black",rep("#CAB2D6",3),
                              "black",rep("#17BECF",2))) +
  facet_wrap(.~scientificName, scales = "free_y", ncol = 1) +
  theme_pubr() +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

forPlot$plotID
forPlot$plotID <- factor(
  forPlot$plotID,
  levels = c("CPER_008", "CPER_003", "CPER_004", "CPER_006", "CPER_002", "CPER_009",
             "STER_033", "STER_026", "STER_028", "STER_034", "STER_035", 
             "STER_006", "STER_029", "STER_031", "STER_032", "STER_027"), 
  ordered = TRUE)

TPD_plot<-ggplot(data = forPlot, aes(x=dist_cm, fill = Spp_plot)) +
  geom_density(position="stack", linewidth = 1.5, alpha=0.5, aes(x=dist_cm, colour = plotID))  +
  # geom_density(position="stack", col="#00000000", alpha=0.5)  +
  scale_fill_manual(values=c(rep("#1B9E77",3),
                             rep("#D95F02",3),
                             rep("#7570B3",1),
                             rep("#1F78B4",7),
                             rep("#66A61E",1),
                             rep("#E6AB02",3),
                             rep("#A6761D",2),
                             rep("#666666",2),
                             rep("#E7298A",1),
                             rep("#9B2F5D",3),
                             rep("#FB9A99",3),
                             rep("#CAB2D6",4),
                             rep("#17BECF",3))) +
  scale_color_manual(values=c("#4C78A8",
                              "#F58518",    
                              "#54A24B",
                              "#E45756",
                              "#72B7B2",
                              "#EECA3B",
                              "#B279A2",
                              "#00876C",
                              "#9D755D",
                              "#BAB0AC",
                              "#2F4B7C",
                              "#A05195",
                              "#FFA600",
                              "#F95D6A")) +
  
  #xlim(0.4,1.6) +
  facet_wrap(.~plotID, scale  = "free_y", ncol = 1) +
  theme_pubr() +
  xlab("Elytra Length (cm)") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.8)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

forPlot$plotSpp<-paste0(forPlot$plotID,forPlot$scientificName)
table(forPlot$plotSpp)
dim(table(forPlot$plotSpp))
forPlot$plotSpp<-sort(forPlot$plotSpp, decreasing = FALSE)

TPD_site<-ggplot(data = forPlot, aes(x=dist_cm, fill = plotID)) +
  geom_density(position="stack", linewidth = 3, alpha=0, aes(x=dist_cm, fill = Spp_plot, colour = siteID))  +
  geom_density(position="stack", linewidth = 0, aes(x=dist_cm, colour = Spp_plot), fill="#FFFFFF")  +
  geom_density(position="stack", linewidth = 0, alpha=0.7, aes(x=dist_cm, colour = Spp_plot))  +
  # geom_density(position="stack", col="#00000000")  +
  scale_fill_manual(values=c("#4C78A8",
                             "#F58518",    
                             "#54A24B",
                             "#E45756",
                             "#72B7B2",
                             "#EECA3B",
                             "#B279A2",
                             "#00876C",
                             "#9D755D",
                             "#BAB0AC",
                             "#2F4B7C",
                             "#A05195",
                             "#FFA600",
                             "#F95D6A", rep("#00000000",48))) +
  scale_color_manual(values=c(rep("#AA3377",1),
                              rep("#CCBB44",1),
                              rep("#4477AA",1),
                              rep("#000000",52))) +
  # #xlim(0.4,1.6) +
  facet_wrap(.~siteID, ncol = 1) +
  theme_pubr() +
  xlab("Elytra Length (cm)") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 0.5)))+
  theme(legend.position="none") +
  theme(strip.background = element_blank(), # Removes the background box
        strip.text = element_blank(),        # Removes the text
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
  )

png("./Figures/TPDexampleV3.png", units = "in", res = 300, height = 9, width = 15)
annotate_figure(
  ggarrange(nrow=1,
            TPD_plot_Spp, TPD_Spp, TPD_plot, TPD_site),
  bottom = text_grob("Elytra Length (cm)", color = "black", size = 14))
dev.off()

png("./Figures/TPDexampleV3Site.png", units = "in", res = 300, height = 6, width = 5)
TPD_site +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(add = c(0, 10)))
dev.off()
