# script generates multiple individual participant reports for gut (or any other) microbiome
# this is intended to be for fun or general interest only NOT diagnostic or clinical
# text and data are synthetic for illustration purposes
# family-level on page 1; named species on page 2

# load packages
library(dplyr)
library(tibble)
library(phyloseq)
library(ggplot2)
library(grid)

# folder for the generated participant reports
output_dir <- "participant_reports"
dir.create(output_dir)

# import microbiome data - based off of phyloseq format
# gut_asv: ASV abundance table (rows = ASVs (or OTUs), columns = sample IDs)
# gut_tax: taxonomic classifications for each ASV (e.g., Kingdom, Phylum, Family)
# metadata: participant information (e.g., sample ID and study visit)

asv_df  <- read.delim("gut_asv.txt", header = TRUE, check.names = FALSE)
tax_df  <- read.delim("gut_tax.txt", header = TRUE)
meta_df <- read.delim("metadata.txt", header = TRUE)

# prepare data for phyloseq object
asv_df <- asv_df %>% tibble::column_to_rownames('asv')
asv_df[] <- lapply(asv_df, function(x) as.numeric(as.character(x)))
asv_mat <- as.matrix(asv_df)
if (anyNA(asv_mat) || any(!is.finite(asv_mat)) || any(asv_mat < 0)) stop('ASV table has missing or negative values.')
tax_mat <- as.matrix(tibble::column_to_rownames(tax_df, 'asv'))
meta_df <- tibble::column_to_rownames(meta_df, 'sample_id')

# keep only IDs shared across tables
shared_taxa <- intersect(rownames(asv_mat), rownames(tax_mat))
shared_samples <- intersect(colnames(asv_mat), rownames(meta_df))
if (!length(shared_taxa) || !length(shared_samples)) stop('there are no matching taxa or sample IDs') # stop running if no matching ASVs or samples.
asv_mat <- asv_mat[shared_taxa, shared_samples, drop = FALSE]
tax_mat <- tax_mat[shared_taxa, , drop = FALSE]
meta_df <- meta_df[shared_samples, , drop = FALSE]
if (!all(c('Family', 'Genus', 'Species') %in% colnames(tax_mat))) stop('taxonomy needs Family, Genus and Species columns.') # check that taxonomy columns exist

# create phyloseq object
Stool <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_df)
)

# remove samples with no bacterial abundance
Stool <- prune_samples(sample_sums(Stool) > 0, Stool)

if (nsamples(Stool) == 0) {
  stop('No samples remain, 0 abundance.')
}

# Calculate relative abundances for each sample (e.g. counts or proportions)
Stool_rel <- transform_sample_counts(Stool, function(x) x / sum(x))
long <- psmelt(Stool_rel) %>%
  mutate(Abundance = Abundance * 100,
         Family = trimws(sub("^f__", "", as.character(Family))),
         Genus = trimws(sub("^g__", "", as.character(Genus))),
         Species = trimws(sub("^s__", "", as.character(Species))))

# summarise families - unnamed classifications are not displayed
valid_name <- function(x) {
  !is.na(x) & nzchar(trimws(x)) &
    !tolower(trimws(x)) %in% c("na", "unknown", "unclassified", "unidentified", "uncultured", "f__", "g__", "s__")
}
family_data <- long %>%
  filter(valid_name(Family)) %>%
  group_by(Sample, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance)) %>%
  ungroup()

# combine ASVs with the same named species for the species chart
make_species_label <- function(genus, species) {
  if (!valid_name(genus) || !valid_name(species)) return(NA_character_)
  species <- gsub("_", " ", trimws(species), fixed = TRUE)
  if (grepl("unclassified|unidentified|metagenome|^sp\\.?$", species, ignore.case = TRUE)) return(NA_character_)
  if (startsWith(tolower(species), paste0(tolower(genus), " "))) return(species)
  if (grepl(" ", species, fixed = TRUE)) return(NA_character_)
  paste(genus, species)
}
long$Species_label <- mapply(make_species_label, long$Genus, long$Species, USE.NAMES = FALSE)
species_data <- long %>%
  filter(!is.na(Species_label)) %>%
  group_by(Sample, Species_label) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance)) %>%
  ungroup()

# you can add in participant-friendly facts/interpretations. This is to give some context to the families, but make sure to not diagnose!
# the "facts" below are synthetic
# these might be examples of what SOME members can do & not claims about every species.
family_facts <- c(
  Prevotellaceae = 'This family includes Prevotella. Some Prevotellaceae members help break down fibre and other complex carbohydrates from plant foods. The abundance can be very different from person to person.',
  Bacteroidaceae = 'This family includes Bacteroides and Phocaeicola. Some members can tackle food components that our own digestive enzymes cannot break down.',
  Lachnospiraceae = 'Some members turn carbohydrates into butyrate, a substance that helps fuel the cells lining your gut.',
  Oscillospiraceae = 'Some members use nutrients or by-products made by other microbes. Gut bacteria often depend on one another!',
  Rikenellaceae = 'This family includes Alistipes. Researchers are exploring how its members respond to the food and conditions inside the gut.',
  Ruminococcaceae = 'Some members can break apart tough plant fibres that our own digestive enzymes struggle to digest.')

safe_id <- function(x) gsub('[^A-Za-z0-9_-]', '_', x)

# make report
make_report <- function(sample_id) {
  df <- family_data %>% filter(Sample == sample_id, Abundance > 0) %>% arrange(desc(Abundance))
  sp <- species_data %>% filter(Sample == sample_id, Abundance > 0) %>% arrange(desc(Abundance)) %>% slice_head(n = 10)
  if (nrow(df) == 0) { warning('No named families for ', sample_id); return(invisible(NULL)) }
  top20 <- head(df, 20)
  top5 <- head(df, 5)
  other <- max(0, 100 - sum(top5$Abundance))
  comp <- bind_rows(top5, data.frame(Sample = sample_id, Family = 'Other identified families', Abundance = other))
  palette <- c('#4776BA', '#53A9A1', '#E8AF54', '#987AC4', '#D87D88')
  fills <- c(setNames(palette[seq_len(nrow(top5))], top5$Family), 'Other identified families' = '#B7BEC6')
  
  p1 <- ggplot(top20, aes(x = reorder(Family, Abundance), y = Abundance)) +
    geom_col(fill = '#4776BA', width = .7) +
    geom_text(aes(label = sprintf('%.1f%%', Abundance)), hjust = -.12, size = 2.6) +
    coord_flip(clip = 'off') +
    scale_y_continuous(limits = c(0, max(top20$Abundance) * 1.3), expand = c(0, 0)) +
    labs(title = 'Your top 20 bacterial families', x = NULL, y = 'Share of identified family-level bacteria (%)') +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_text(size = 7.5), panel.grid.major.y = element_blank(), plot.title = element_text(face = 'bold'),
          plot.margin = margin(8, 12, 8, 12))
  
  p2 <- ggplot(comp, aes(x = 'Your sample', y = Abundance, fill = Family)) +
    geom_col(width = .92) + coord_flip() +
    scale_fill_manual(values = fills, breaks = comp$Family) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.05)) +
    labs(title = 'Your whole bacterial community', x = NULL, y = 'Share of identified family-level bacteria (%)', fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = 'bold'), legend.position = 'bottom',
          legend.text = element_text(size = 8), panel.grid.major.y = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  
  facts <- vapply(top5$Family, function(g) {
    fact <- if (g %in% names(family_facts)) unname(family_facts[g]) else
      'Scientists are still learning about this family.'
    paste0(g, ' (', sprintf('%.1f%%', top5$Abundance[match(g, top5$Family)]), '): ', fact)
  }, character(1))
  wrapped <- paste(vapply(facts, function(x) paste(strwrap(x, width = 72), collapse = '\n'), character(1)), collapse = '\n\n')
  
  filename <- file.path(output_dir, paste0('Gut_Microbiome_Families_', safe_id(sample_id), '.pdf'))
  grDevices::pdf(filename, width = 8.27, height = 11.69, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  # PAGE ONE: establish an inset content area with 1.5 cm margins.
  grid.newpage()
  pushViewport(viewport(x = unit(1.5, 'cm'), y = unit(1.5, 'cm'),
                        width = unit(1, 'npc') - unit(3, 'cm'),
                        height = unit(1, 'npc') - unit(3, 'cm'),
                        just = c('left', 'bottom')))
  pushViewport(viewport(layout = grid.layout(5, 1, heights = unit(c(1.8, 4.5, 2.8, 4.2, 0.5), 'null'))))
  # Title and parent-friendly introduction above the family chart
  pushViewport(viewport(layout.pos.row = 1))
  grid.text('Meet Your Gut Microbiome!',
            x = unit(.5, 'npc'), y = unit(.83, 'npc'),
            gp = gpar(fontsize = 19, fontface = 'bold'))
  grid.text(paste0('Sample: ', sample_id),
            x = unit(.5, 'npc'), y = unit(.64, 'npc'),
            gp = gpar(fontsize = 10, col = 'grey40'))
  intro_text <- paste(
    'Your gut is home to many different types of bacteria, which naturally live together as a community.',
    'Some help break down parts of our food, while others use substances made by their neighbours.',
    'Diet is one of several things that can influence this mixture, which can change over time.'
  )
  grid.text(paste(strwrap(intro_text, width = 85), collapse = '\n'),
            x = unit(.5, 'npc'), y = unit(.30, 'npc'),
            gp = gpar(fontsize = 10, col = 'grey30'),
            just = 'center')
  popViewport()
  print(p1, vp = viewport(layout.pos.row = 2), newpage = FALSE)
  print(p2, vp = viewport(layout.pos.row = 3), newpage = FALSE)
  pushViewport(viewport(layout.pos.row = 4))
  grid.text(paste0('GET TO KNOW YOUR BACTERIAL FAMILIES\n\n', wrapped), x = unit(.05, 'npc'),
            y = unit(.95, 'npc'), just = c('left', 'top'), gp = gpar(fontsize = 8.5))
  popViewport()
  # each popViewport() closes a different layer as this pdf has different nested layout sections
  popViewport()
  popViewport() # leave the page-one margin viewport
  
  # PAGE TWO
  grid.newpage()
  pushViewport(viewport(x = unit(1.5, 'cm'), y = unit(1.5, 'cm'),
                        width = unit(1, 'npc') - unit(3, 'cm'),
                        height = unit(1, 'npc') - unit(3, 'cm'),
                        just = c('left', 'bottom')))
  if (nrow(sp) > 0) {
    p3 <- ggplot(sp, aes(x = reorder(Species_label, Abundance), y = Abundance)) +
      geom_col(fill = '#53A9A1', width = .72) +
      geom_text(aes(label = sprintf('%.2f%%', Abundance)), hjust = -.12, size = 2.8) +
      coord_flip(clip = 'off') +
      scale_y_continuous(limits = c(0, max(sp$Abundance) * 1.35), expand = c(0, 0)) +
      labs(title = 'Your top 10 named bacterial species',
           subtitle = paste('Sample:', sample_id),
           x = NULL, y = 'Share of identified species-level bacteria (%)',
           caption = 'The chart shows the 10 most common named species; percentages include all identified species.') +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = 'bold', size = 16),
            axis.text.y = element_text(face = 'italic', size = 8.5),
            panel.grid.major.y = element_blank(),
            axis.title.x = element_text(hjust = 0),
             plot.caption = element_text(hjust = 0),
             plot.margin = margin(20, 15, 20, 15))
    print(p3, newpage = FALSE)
  } else {
    grid.text('A closer look at your species\n\nNo confidently named species were available for this sample.',
              gp = gpar(fontsize = 15))
  }
  popViewport() # leave the page-two margin viewport
  # footer text page 2
  grid.text('A snapshot of bacterial DNA detected in your sample. This is not a clinical assessment.',
            x = unit(0.5, 'npc'), y = unit(0.65, 'cm'),
            gp = gpar(fontsize = 8, col = 'grey40'), just = 'center')
  message('Created: ', filename)
  invisible(filename)
}

# create a PDF for every sample
all_sample_ids <- as.character(unique(family_data$Sample))
reports_created <- 0
for (id in all_sample_ids) {
  result <- make_report(id)
  if (!is.null(result)) reports_created <- reports_created + 1
}
cat('Reports generated:', reports_created, '\n')

cat('\nSamples included:', nsamples(Stool), '\n')
cat('PDF folder:', output_dir, '\n')
