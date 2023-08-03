#!/bin/bash

# Data processing using AWK and SED
for file in *.MR.txt;
do
  awk -F"\t" '{ print $6,"\t",$4,"\t",$3,"\t",$5,"\t",$7,"\t",$7-1.96*$8,"\t", $7+1.96*$8,"\t",$9,"\t",$7,"\t",1.96*$8 }' "$file" > "tmp";
  sed -e '1s/0/lower.confidence.interval/' -e '1s/0/upper.confidence.interval/' -e '1s/0/ci/' "tmp" > "$file.Results.txt";
done


# Read the Previous results files and filter rows with "Inverse variance weighted" in the methods column
results_files <- list.files(pattern = "*Results.txt")
all_rows <- do.call(rbind, lapply(results_files, function(file) 
{data <- read.delim(file, header = TRUE, sep = "\t")
  subset(data, methods == "Inverse variance weighted")}
                                 ))

# Save the filtered rows into a new file
write.csv(all_rows, file = "all.traits.IVW.results.csv", sep = "\t", row.names = FALSE, col.names = TRUE)



# R script - Scatterplot generation
Rscript -e 'library(ggplot2);
library(RColorBrewer);
library(viridis);

ds 
# Read the file
ds <- 
read.csv("/mnt/isilon/thom_lab/thomlab_MR/shivakums1/Bitarello.Multi.Blood/all.traits.IVW.results.csv")

# Extract the trait name
ds$id.exposure <- sub("^[^.]+\\.(.*?)\\..*$", "\\1", ds$id.exposure)

# Setting as factors
ds$Exposure <- factor(ds$id.exposure, levels = unique(ds$id.exposure))

# Plot
bp <- ggplot(ds, aes(x = b, y = Exposure, color = Exposure)) + 
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_errorbarh(aes(xmin = lower.confidence.interval, xmax = upper.confidence.interval), height = 0.3) +
  labs(title = "DrnksWk EUR -> EUR",y = NULL, x = "Effect (SD units)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
#, panel.grid = element_blank()) + 
  scale_y_discrete(limits = rev(levels(ds$Exposure))) +
  xlim(-0.4, 0.4) +  
  geom_text(aes(label = ifelse(pval < 0.05, "*", "")),vjust = -0.5, color 
= "black") +
 geom_vline(xintercept = 0, linetype = "solid", color = "grey")  


# Color code y-axis labels
my_colors <- viridis_pal(option = "viridis")(16)
bp <- bp + scale_color_manual(values = my_colors)

# Save the plot
ggsave("scatterplot.pdf", plot = bp)

