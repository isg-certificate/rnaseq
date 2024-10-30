library("tidyverse")
library("DESeq2")

# get a list of count files
directory <- "../../results/04_counts/counts"
sampleFiles <- list.files(directory, pattern = ".*counts$")

# read in the metadata table linking SRA accessions to experimental conditions
meta <- read.csv("../../metadata/SraRunTable.txt") %>%
  mutate(population=str_replace(population, "Eliza.*", "ER")) %>% # change population names to be simpler
  mutate(population=str_replace(population, "King.*", "KC")) %>% # change population names to be simpler
  mutate(pcb_dosage = case_when(pcb_dosage == 0 ~ "control", pcb_dosage > 0 ~ "exposed")) %>% # change pcb dosage to control vs exposed
  filter(population %in% c("ER","KC")) %>%
  select(Run, population, pcb_dosage)

# ensure that sampleFiles and metadata table are in the same order
all( str_remove(sampleFiles, ".counts") == meta[,1] )

# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
  sampleName = meta$Run,
  fileName = sampleFiles,
  population = meta$population,
  dose = meta$pcb_dosage
)

# create the DESeq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable, 
  directory = directory, 
  design = ~ population*dose
)

# Reset population factor levels so KC is base level
treatments <- c("KC","ER")
ddsHTSeq$population <- factor(ddsHTSeq$population, levels = treatments)

# Run statistical analysis
dds <- DESeq(ddsHTSeq)

# Extract significance tests for the dose factor for control vs exposed for King's Creek population (KC)
res <- results(dds, contrast=c("dose","exposed","control"))

# get a quick summary of results
summary(res)

