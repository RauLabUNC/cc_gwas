# Plots to make for IMGS 
# 1. Manhattan plot for delta Wall thickness day 28
# 2. Manhattan plots for Fibs + CMs at control
# 3. Locus zoom plots for each 1 and 2 at chr 15
# 4. Dot plot for Arid2 vs Wall.Thicknessd.28/LVDWs.29

library(tidyverse)
library(qqman)

#### Manhattan for delta wall thickness ####
# Load in GWAS results for delta
traits_delta <- read.csv("data/processed/christophGWAS/PerChange_pvals.csv")

WallThickd28Delta <- traits_delta |> 
  select(Name:MAF, Wall.Thicknessd.28) |> 
  filter(Chr %in% 1:20) |> 
  mutate(Chr = as.numeric(Chr),
         Pos_mm39 = as.numeric(Pos_mm39)) 


don <- WallThickd28Delta %>% 
  
  # Compute chromosome size
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos_mm39)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(WallThickd28Delta, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos_mm39) %>%
  mutate( BPcum=Pos_mm39+tot)

axisdf = don %>%
  group_by(Chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


wallPlot <- ggplot(don, aes(x=BPcum, y= -log10(Wall.Thicknessd.28))) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  geom_hline(yintercept = -log10(10^-5), linetype = "dashed", color = "#feb078") +
  scale_color_manual(values = rep(c("#b73779", "#f1605d"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +
  scale_y_continuous(breaks = c(3,6,9),limits = c(0, 9)) +  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 8, color = "#fcfdbf"),
    axis.text = element_text(size = 8, color = "#fcfdbf"),
    axis.ticks = element_line(color = "#fcfdbf")
  ) + 
  labs(x = "Chromosome", y = "-log10(p-value)") 

ggsave("results/imgs/deltaWallThickManhattan.png", wallPlot, 
       height = 6, width = 14, units = "cm", dpi = 600)

#### Plot manhattan for FBs and CMs ####

# Load in GWAS results for delta
cells <- read.csv("data/processed/christophGWAS/RoughPercentCellTypes_pvals.csv")

cells_clean <- cells |> 
  select(Name:MAF, Cardiomyocytes_Control, Fibroblast_Control, Endothelial.Cells_Iso) |> 
  filter(Chr %in% 1:20 & !is.na(Pos_mm39)) |> 
  mutate(Chr = as.numeric(Chr),
         Pos_mm39 = as.numeric(Pos_mm39)) 



FBs <- cells_clean %>% 
  select(-Cardiomyocytes_Control, -Endothelial.Cells_Iso) |> 
  # Compute chromosome size
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos_mm39)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(cells_clean, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos_mm39) %>%
  mutate( BPcum=Pos_mm39+tot)

axisdf = FBs %>%
  group_by(Chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


fbPlot <- ggplot(FBs, aes(x=BPcum, y= -log10(Fibroblast_Control))) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  geom_hline(yintercept = -log10(10^-5), linetype = "dashed", color = "#feb078") +
  scale_color_manual(values = rep(c("#b73779", "#f1605d"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +   # remove space between plot area and x axis
  scale_y_continuous(breaks = c(2,4,6),limits = c(0, 6)) +  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 8, color = "#fcfdbf"),
    axis.text = element_text(size = 8, color = "#fcfdbf"),
    axis.ticks = element_line(color = "#fcfdbf")
  ) + 
  labs(x = "Chromosome", y = "-log10(p-value)") 

ggsave("results/imgs/controlFBManhattan.png", fbPlot, 
       height = 6, width = 14, units = "cm", dpi = 600)



#### Plot CM abundance #####

CMs <- cells_clean %>% 
  select(-Endothelial.Cells_Iso, -Fibroblast_Control) |> 
  # Compute chromosome size
  group_by(Chr) %>% 
  summarise(chr_len=max(Pos_mm39)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(cells_clean, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, Pos_mm39) %>%
  mutate( BPcum=Pos_mm39+tot)

axisdf = CMs %>%
  group_by(Chr) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cmPlot <- ggplot(CMs, aes(x=BPcum, y= -log10(Cardiomyocytes_Control))) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  geom_hline(yintercept = -log10(10^-5), linetype = "dashed", color = "#feb078") +
  scale_color_manual(values = rep(c("#b73779", "#f1605d"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +   # remove space between plot area and x axis
  scale_y_continuous(breaks = c(2,4,6),limits = c(0, 6)) +  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 8, color = "#fcfdbf"),
    axis.text = element_text(size = 8, color = "#fcfdbf"),
    axis.ticks = element_line(color = "#fcfdbf")
  ) + 
  labs(x = "Chromosome", y = "-log10(p-value)") 

ggsave("results/imgs/controlCMManhattan.png", cmPlot, 
       height = 6, width = 14, units = "cm", dpi = 600)



#### Make locus zoom plots ####
library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
# Load genome assembly
mm39 <- assembly(Genome = "mm39_GRCm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene", OrgDb = "org.Mm.eg.db")

wallPG <- WallThickd28Delta |> 
  mutate(pos = Pos_mm39,
         chrom = paste0("chr", Chr),
         p = Wall.Thicknessd.28) |> 
  dplyr::select(pos, chrom, p) |> 
  filter(!is.na(p))

fbPG <- FBs |> 
  mutate(pos = Pos_mm39,
         chrom = paste0("chr", Chr),
         p = Fibroblast_Control) |> 
  dplyr::select(pos, chrom, p) |> 
  filter(!is.na(p))
# Set plot boundaries

upper.lim <- 96497549 - 5*10^5
lower.lim <- 96497549 + 5*10^5
# Create a plotgardener page
pageCreate(
  width = 6, height = 5, default.units = "inches",
  showGuides = FALSE, xgrid = 0, ygrid = 0
)

params <- pgParams(assembly = mm39, quiet = T, 
                   just = c("center", "center"), 
                   default.units = "inches",
                   chromstart = upper.lim, chromend = lower.lim)
# Plot LOD score data using plotSignal
manhattanPlot <- plotManhattan(
  data = wallPG,
  chrom = "chr15", 
  trans = "-log10", sigVal = 10^(-5),
  x = 3, y = 0, width = 6, height = 1.2,
  just = c("center", "top"), 
  xfield = "pos", yfield = "p", fill = "black", sigCol = "darkorange",
  sigLine = T, baseline = T,
  yscale_reverse = F, params = params
)

  ## Annotate y-axis
annoYaxis(
  plot = manhattanPlot, at = seq(0, 6),
  axisLine = TRUE, fontsize = 8, params = params)

manhattanPlot <- plotManhattan(
  data = fbPG,
  chrom = "chr15", 
  trans = "-log10", sigVal = 10^(-5),
  x = 3, y = 1.4, width = 6, height = 1.2,
  just = c("center", "top"), 
  xfield = "pos", yfield = "p", fill = "black", sigCol = "darkorange",
  sigLine = T, baseline = T,
  yscale_reverse = F, params = params
)

## Annotate y-axis
annoYaxis(
  plot = manhattanPlot, at = seq(0, 6),
  axisLine = TRUE, fontsize = 8, params = params)


## Plot text label

plotGenes(
  chrom = "chr15", chromstart = upper.lim, chromend = lower.lim,
  assembly = mm39,
  x = 3, y = 2.5, width = 6, height = 2,
  just = c("center", "top"), default.units = "inches",
  params = params
)

plotGenomeLabel(
  chrom = "chr15", chromstart = upper.lim, chromend = lower.lim,
  assembly = mm39,
  x = 3, y = 4.5, length = 6, scale = "Mb",
  just = c("center", "top"), default.units = "inches",
  params = params, quiet = T
)




