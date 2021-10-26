###########################################
########      CITP diversity         ######
########       Fig 1 and 2          #######
###########################################



library(ggmap)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("/Users/anastasia/Documents/CITP")

##################################
######   Fig1. Strain map      ###
##################################



#register_google(key="mygoogleAPIkey")
species<-c(rep("C. elegans",8),rep("C. briggsae",8),rep("C. tropicalis",6))

# coordinates for the CITP strains
names<-c("N2-PD1073",
         "CB4856",
         "ED3040",
         "JU1088",
         "JU1652",
         "JU775",
         "MY16",
         "QX1211",
         "AF16",
         "HK104",
         "JU1264", 
         "JU1348",
         "JU726",
         "NIC20",
         "QR25",
         "ED3092",
         "JU1373",
         "JU1630",
         "NIC122",
         "NIC58",
         "QG131",
         "QG834")



	
myLocation <-
  data.frame(    lon = c(
    -2.587900,
    -158.000100,
    28.047300,
    138.015000,
    -56.164500,
    -9.148500,
    7.570400,
    -122.4331,
    72.571400,
    133.919800,
    49.1260,
    77.236800,
    109.640500,
    120.2269,
    -73.567300,
    36.798500,
    55.6885,
    -25.1222,
    -61.6563,
    -52.6787,
    -159.476,
    -79.8308
  ),
    lat = c(
      51.4545,
      21.4389,
      -26.2041,
      34.7612,
      -34.9011,
      38.7183,
      51.9247,
      37.7502,
      23.0225,
      34.6555,
      49.126,
      9.4622,
      25.8952,
      22.999,
      45.5017,
      -1.3245,
      -21.0473,
      17.1064,
      16.1779,
      4.08731,
      22.22955,
      9.15385
    )
  )


samples<-data.frame(myLocation,Species=species,Names=names)
samples$Species <- factor(samples$Species , levels = c("C. elegans","C. briggsae","C. tropicalis"))



center = paste(min(myLocation$lat)+(max(myLocation$lat)-min(myLocation$lat))/2,
               min(myLocation$lon)+(max(myLocation$lon)-min(myLocation$lon))/2, sep=" ")

citpmap<-get_map(location=center,zoom=1,maptype ="terrain")
mapl <-    ggmap(citpmap,extent = "device",fullpage=T) +
  geom_text_repel(
    data = samples,
    aes(label = Names),
    box.padding = unit(0.29, "lines"),
    colour = "#000000",
    size = 2
  ) + geom_point(data = samples,
                 aes(
                   x = lon,
                   y = lat,
                   shape = Species,
                   color = Species
                 ),
                 size = 2)  + theme(panel.background = element_blank())+
  theme(axis.line = element_line(color = NA))  + labs(title="A", x = "longitude", y = "latitude") + theme_bw() + theme_bw(base_size = 12) +
  theme(legend.position ="bottom")  + theme(legend.text = element_text(face = "italic", size = rel(1))) + scale_color_manual(values = c("#000000", "#f2b230", "#078a8a")) + theme(legend.title = element_blank()) + scale_y_continuous(limits = c(-74, 74)) + 
  scale_x_continuous(limits = c(-180, 180))  + theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  theme(legend.position = c(0.5,1.07)) + guides(col=guide_legend(ncol=3,byrow=FALSE)) + theme(plot.title = element_text(face = "bold")) 



##################################
####### Fig2. popVAE results  ####
##################################


CE<-read.csv("CE_popvae_FIN_latent_coords.txt", header=TRUE, sep="\t")
CB<-read.csv("CB_popvae_FIN_latent_coords.txt", header=TRUE, sep="\t")
CT<-read.csv("CT_popvae_FIN_latent_coords.txt", header=TRUE, sep="\t")

CE$species <-"C. elegans"
CB$species<-"C. briggsae"
CT$species<-"C. tropicalis"



COMBO<-rbind(CE,CB)
COMBO<-rbind(COMBO,CT)

COMBO$sampleID<-gsub("^C[EBT]_","",perl=T,COMBO$sampleID)
COMBO$sampleID<-gsub("PD1073","N2-PD1073",COMBO$sampleID)
COMBO$species <- factor(COMBO$species , levels = c("C. elegans","C. briggsae","C. tropicalis"))
COMBO$core<-"no"
COMBO[COMBO$sampleID %in% c("JU775",
                            "MY16",
                            "N2-PD1073", "AF16",
                            "ED3092",
                            "HK104","JU1373",
                            "JU1630","QG834"),]$core<-"core"

###LD shows the latent dimentions
pl<-ggplot(COMBO, aes(x = mean1, y = mean2, color=species, alpha=core, ordered = FALSE))  +
  labs(title="B", x = "LD1", y = "LD2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_wrap(.~species,scales = "free") +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) +
  geom_text_repel(aes(label = sampleID), force = 1, min.segment.length = 0, seed = 42, box.padding = unit(0.8, "lines"),colour = "#000000",size=3) +
  geom_point(size=2) + scale_color_manual(values = c("#000000","#f2b230","#078a8a")) + scale_alpha_manual(values = c(1,0.4))


png(filename = "CITP_popvae.png",width = 5,height = 3.25, res=500, units = "in")
plot(pl)
dev.off()


############################
###### Nucleotide diversity
############################

my_files<-list.files(pattern="*_pi_100kb_FIN.bed")
my_names<-gsub("_pi_100kb_FIN.bed","", perl=F, my_files)

STATS<-c();for (i in 1:(length(my_files)) ) { B=read.csv(file=my_files[i], header=F, sep="\t"); STATS<-rbind(STATS,data.frame(B,sample=my_names[i]));}

STATS$Species<-gsub("^(C[BET])_.*", "\\1", perl=T, STATS$sample)
STATS$Core<-gsub("_CITP", "", perl=T, STATS$sample)
STATS$Core<-gsub(".*.3-300", "All", perl=F, STATS$Core)
STATS$Core<-gsub(".*.core", "Core", perl=F, STATS$Core)

STATS$Species<-gsub("CB", "C. briggsae",STATS$Species)
STATS$Species<-gsub("CE", "C. elegans",STATS$Species)
STATS$Species<-gsub("CT", "C. tropicalis",STATS$Species)
STATS$V4<-as.numeric(as.character(STATS$V4))

STATS$Species <- factor(STATS$Species , levels = c("C. elegans","C. briggsae","C. tropicalis"))




pipl<-ggplot(STATS, aes(x=Core, y = V4, color=Species, fill=Species, alpha=Core, ordered = FALSE))   + geom_jitter(size=0.25) + geom_violin() + labs(title="C", x = "", y = expression(pi)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  facet_wrap(.~Species,scales = "free") +
  theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt"))   + scale_color_manual(values = c("#000000","#f2b230","#078a8a")) + scale_fill_manual(values = c("#000000","#f2b230","#078a8a")) + scale_alpha_manual(values = c(0.3,0.7)) + geom_boxplot(width=0.1, color="#4d4c4c",outlier.shape = NA)




#combine all plots
complot<-  mapl | (pl / pipl)

png(filename = "CITP_Fig1.png",width = 12,height = 5.4, res=500, units = "in")
plot(complot)
dev.off()


