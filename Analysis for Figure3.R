#--------------------Figure 3 Impact of contaminants----------------------------

#1.---------------- Load required libraries-----------------------------------
library(tidyverse)
library(openxlsx)
library(patchwork)

#2.-----------------Figure 3A: Theoretical peptides-----------------------------

# Set working directory for read Theoretical Peptides
setwd()

# Read files into a list
TheoreticalPepList <- lapply(1:12, function(i) {
  read.csv(paste0('TheoreticalPep_Mix-', i, ".txt"), sep = "\t")
})

# Combine list and process
Theoretical_all <- bind_rows(TheoreticalPepList) %>%
  group_by(Sequence) %>%
  summarise(across(starts_with("Mix"), ~ ifelse(sum(., na.rm = TRUE) == 0, NA, sum(., na.rm = TRUE))),
            Protein.Group = first(Protein.Group)) %>%
  mutate(Count = rowSums(!is.na(select(., -Sequence, -Protein.Group)))) %>%
  select(Protein.Group, Sequence, starts_with("Mix"), Count) %>%  # Reorder columns
  mutate(PepLength = nchar(Sequence)) %>%  # Calculate sequence length
  filter(PepLength <= 30 & PepLength >= 7)  # DIA-NN default 7-30

# Calculate frequency/count
SharePep_Theoretical <- Theoretical_all %>%
  count(Count) %>%
  mutate(Ratio = round(n / sum(n) * 100, 2),
         AccPercentage = cumsum(Ratio),
         log10 = log10(n))

# Create figure
figure_1a <- ggplot(SharePep_Theoretical) +
  geom_bar(aes(x = Count, y = log10), stat = "identity", width = 0.85, fill = "#CCCCCC") +
  labs(x = "Peptides presence frequency", y = "log10, peptides") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 6)) +
  geom_text(aes(x = Count, y = log10, label = n), vjust = 0.5, hjust = 1.2, size = 3.5, color = "black", angle = 90)

#3.-----------------Figure 3B: DIA-NN 1.8.1 contaminants-------------------

# Set working directory
setwd() # for output of diann 1.8.1

# Read pr.matrix report
prList181 <- lapply(1:12, function(i) {
  read.csv(paste0("Mix-", i, ".pr_matrix.tsv"), sep = "\t") %>%
    filter(Proteotypic == 1) %>%
    select(Stripped.Sequence, Precursor.Id) %>%
    mutate(Mix = paste0("Mix", i))
})

# Combine and process data
prall_181 <- bind_rows(prList181) %>%
  pivot_wider(names_from = Mix, values_from = Mix) %>%
  group_by(Precursor.Id) %>%
  summarise(across(starts_with("Mix"), ~ ifelse(sum(!is.na(.), na.rm = TRUE) == 0, NA, sum(!is.na(.), na.rm = TRUE))),
            Stripped.Sequence = first(Stripped.Sequence)) %>%
  mutate(Count = rowSums(!is.na(select(., -Stripped.Sequence, -Precursor.Id)))) %>%
  rename(Sequence = Stripped.Sequence) %>%
  select(Precursor.Id, Sequence, starts_with("Mix"), Count)

# Calculate frequency/counts
SharedPep_181 <- prall_181 %>%
  distinct(Sequence, .keep_all = TRUE) %>%
  count(Count) %>%
  mutate(Ratio = n / sum(n) * 100,
         AccPercentage = cumsum(Ratio),
         log10 = log10(n))

# Create figure
figure_1b <- ggplot(SharedPep_181) +
  geom_bar(aes(x = Count, y = log10), stat = "identity", width = 0.85, fill = "#FFCD92") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 6)) +
  geom_text(aes(x = Count, y = log10, label = n), vjust = 0.5, hjust = 1.2, size = 3.5, color = "black", angle = 90)

#4.-----------------Figure 3C: DIA-NN 1.9 contaminants--------------------------

# Set working directory
setwd() # for output of diann 1.9

# Read pr.matrix report
prList_19 <- lapply(1:12, function(i) {
  read.csv(paste0("Mix-", i, ".pr_matrix.tsv"), sep = "\t") %>%
    filter(Proteotypic == 1) %>%
    select(Stripped.Sequence, Precursor.Id) %>%
    mutate(Mix = paste0("Mix", i))
})

# Combine and process data
prall_19 <- bind_rows(prList_19) %>%
  pivot_wider(names_from = Mix, values_from = Mix) %>%
  group_by(Precursor.Id) %>%
  summarise(across(starts_with("Mix"), ~ ifelse(sum(!is.na(.), na.rm = TRUE) == 0, NA, sum(!is.na(.), na.rm = TRUE))),
            Stripped.Sequence = first(Stripped.Sequence)) %>%
  mutate(Count = rowSums(!is.na(select(., -Stripped.Sequence, -Precursor.Id)))) %>%
  rename(Sequence = Stripped.Sequence) %>%
  select(Precursor.Id, Sequence, starts_with("Mix"), Count)

# Calculate frequency/counts
SharedPep_19 <- prall_19 %>%
  distinct(Sequence, .keep_all = TRUE) %>%
  count(Count) %>%
  mutate(Ratio = n / sum(n) * 100,
         AccPercentage = cumsum(Ratio),
         log10 = log10(n))

# Create figure
figure_1c <- ggplot(SharedPep_19) +
  geom_bar(aes(x = Count, y = log10), stat = "identity", width = 0.85, fill = "#91DA91") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 6)) +
  geom_text(aes(x = Count, y = log10, label = n), vjust = 0.5, hjust = 1.2, size = 3.5, color = "black", angle = 90)

#5.-----------------Figure 3D: FDR calculations---------------------------------

# 5.1 Function to create peptide lists for each mix (with carryovers)
createMixList <- function(prall) {
  lapply(1:12, function(i) {
    Mixname <- paste0("Mix", i)
    prall %>%
      filter(!is.na(!!sym(Mixname))) %>%
      select(Sequence, all_of(Mixname), Count) %>%
      arrange(desc(Count)) %>%
      distinct(Sequence, .keep_all = TRUE)
  })
}

# Function to calculate FDR
calculateFDR <- function(mixList, theoreticalPepList) {
  FDR_list <- lapply(1:12, function(countLimit) {
    temp_list <- lapply(1:length(mixList), function(i) {
      df <- mixList[[i]]
      fasta <- theoreticalPepList[[i]]
      df %>%
        filter(Count <= countLimit) %>%
        mutate(Exist = ifelse(Sequence %in% fasta$Sequence, 1, 0)) %>%
        summarise(Alldiscovery = n(), Fdiscovery = sum(Exist == 0)) %>%
        mutate(SharedNumber = countLimit)
    })
    bind_rows(temp_list)
  })
  list(FDR_list = FDR_list,
       FDR_summary = bind_rows(FDR_list) %>%
         mutate(FDR = round(Fdiscovery / Alldiscovery * 100, 2)) %>%
         group_by(SharedNumber) %>%
         summarise(FDR = round(mean(FDR), 2)))
}

# Generate MixLists and calculate FDR for both DIA-NN 1.8.1 and 1.9
MixList_181_withcarry <- createMixList(prall_181)
MixList_19_withcarry <- createMixList(prall_19)

FDR_result_181_withcarry <- calculateFDR(MixList_181_withcarry, TheoreticalPepList)
FDR_result_19_withcarry <- calculateFDR(MixList_19_withcarry, TheoreticalPepList)

# Extract FDR summary and FDR list
FDR_summary_181_withcarry <- FDR_result_181_withcarry$FDR_summary
FDR_list_181_withcarry <- FDR_result_181_withcarry$FDR_list

FDR_summary_19_withcarry <- FDR_result_19_withcarry$FDR_summary
FDR_list_19_withcarry <- FDR_result_19_withcarry$FDR_list

# Combine FDR results (with carryovers)
FDR_summary_withcarry <- bind_cols(FDR_summary_181_withcarry, FDR_summary_19_withcarry) %>%
  select(1, 2, 4) %>%
  setNames(c("Count", "diann181", "diann19"))

# 5.2 Function to create peptide lists without carryovers
createMixListWithoutCarryovers <- function(prall) {
  Mix1 <- prall %>%
    filter(!is.na(Mix1)) %>%
    select(Sequence, Mix1, Count) %>%
    arrange(desc(Count)) %>%
    distinct(Sequence, .keep_all = TRUE)
  
  Mix2to12 <- lapply(2:12, function(i) {
    Mixname <- paste0("Mix", i)
    Mixname_1 <- paste0("Mix", i - 1)
    prall %>%
      filter(!is.na(!!sym(Mixname))) %>%
      select(Sequence, !!sym(Mixname_1), !!sym(Mixname), Count) %>%
      arrange(desc(Count)) %>%
      distinct(Sequence, .keep_all = TRUE) %>%
      filter(is.na(!!sym(Mixname_1))) %>%
      select(Sequence, !!sym(Mixname), Count)
  })
  c(list(Mix1), Mix2to12)
}

# Calculate FDR without carryovers for both DIA-NN 1.8.1 and 1.9
MixList_181_withoutcarry <- createMixListWithoutCarryovers(prall_181)
MixList_19_withoutcarry <- createMixListWithoutCarryovers(prall_19)

FDR_result_181_withoutcarry <- calculateFDR(MixList_181_withoutcarry, TheoreticalPepList)
FDR_result_19_withoutcarry <- calculateFDR(MixList_19_withoutcarry, TheoreticalPepList)

# Extract FDR summary and FDR list
FDR_summary_181_withoutcarry <- FDR_result_181_withoutcarry$FDR_summary
FDR_list_181_withoutcarry <- FDR_result_181_withoutcarry$FDR_list

FDR_summary_19_withoutcarry <- FDR_result_19_withoutcarry$FDR_summary
FDR_list_19_withoutcarry <- FDR_result_19_withoutcarry$FDR_list

# Combine FDR results (without carryovers)
FDR_summary_withoutcarry <- bind_cols(FDR_summary_181_withoutcarry, FDR_summary_19_withoutcarry) %>%
  select(1, 2, 4) %>%
  setNames(c("Count", "diann181", "diann19"))

# 5.3 Combine FDR results for plotting
FDR_figure <- bind_cols(FDR_summary_withcarry, FDR_summary_withoutcarry) %>%
  select(1, 2, 3, 5, 6) %>%
  setNames(c("Count",
             "DIA-NN 1.8.1 with carryovers",
             "DIA-NN 1.9 with carryovers",
             "DIA-NN 1.8.1 without carryovers",
             "DIA-NN 1.9 without carryovers"))

# Convert data to long format for plotting
FDR_figure_long <- pivot_longer(FDR_figure, cols = -Count,
                                names_to = "version",
                                values_to = "FDR")

# Set factor levels for Count
FDR_figure_long$Count <- factor(as.character(FDR_figure_long$Count), levels = as.character(12:1))

# Create figure
figure_1d <- ggplot(FDR_figure_long, aes(x = Count, y = FDR, group = version, color = version)) +
  geom_line(aes(linetype = ifelse(grepl("with carryovers", version), "solid", "dashed")),
            linewidth = 1) +  
  geom_point(size = 2, shape = 21, aes(fill = version), color = "black", stroke = 0.8) +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c(
    "DIA-NN 1.8.1 without carryovers" = "#FFCD92",
    "DIA-NN 1.8.1 with carryovers" = "#FFCD92",
    "DIA-NN 1.9 without carryovers" = "#91DA91",
    "DIA-NN 1.9 with carryovers" = "#91DA91"
  )) +
  scale_fill_manual(values = c(
    "DIA-NN 1.8.1 without carryovers" = "#FFCD92",
    "DIA-NN 1.8.1 with carryovers" = "#FFCD92",
    "DIA-NN 1.9 without carryovers" = "#91DA91",
    "DIA-NN 1.9 with carryovers" = "#91DA91"
  )) +  
  theme(
    legend.position = c(0.7, 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10.1),
    axis.text.y = element_text(size = 10.1),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 5, color = "black"),
    axis.ticks = element_line(linewidth = 0.25)
  ) +
  labs(x = "Peptides detection frequency", y = "Average FDR (%)")

#6. ----------------Patch the Figure 3------------------------------------------
(figure_1a + figure_1b) / (figure_1c + figure_1d)