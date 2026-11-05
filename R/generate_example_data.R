
# generate synthetic gut microbiome data for demonstration
# all sample IDs and abundance values are fictional

set.seed(123)

# example bacterial species
species <- c(
  "Bacteroides uniformis",
  "Bacteroides vulgatus",
  "Bacteroides thetaiotaomicron",
  "Phocaeicola dorei",
  "Prevotella copri",
  "Faecalibacterium prausnitzii",
  "Agathobacter rectalis",
  "Roseburia intestinalis",
  "Blautia obeum",
  "Blautia wexlerae",
  "Anaerostipes hadrus",
  "Coprococcus comes",
  "Dorea longicatena",
  "Ruminococcus bromii",
  "Ruminococcus bicirculans",
  "Subdoligranulum variabile",
  "Alistipes putredinis",
  "Alistipes finegoldii",
  "Parabacteroides distasonis",
  "Parabacteroides merdae",
  "Akkermansia muciniphila",
  "Bifidobacterium longum",
  "Bifidobacterium adolescentis",
  "Collinsella aerofaciens",
  "Escherichia coli",
  "Eubacterium siraeum",
  "Odoribacter splanchnicus",
  "Barnesiella intestinihominis",
  "Dialister invisus",
  "Veillonella parvula",
  "Streptococcus salivarius",
  "Clostridium butyricum",
  "Turicibacter sanguinis",
  "Sutterella wadsworthensis",
  "Christensenella minuta"
)

# family classifications matching the species above
family <- c(
  "Bacteroidaceae", "Bacteroidaceae", "Bacteroidaceae",
  "Bacteroidaceae", "Prevotellaceae",
  "Oscillospiraceae", "Lachnospiraceae",
  "Lachnospiraceae", "Lachnospiraceae",
  "Lachnospiraceae", "Lachnospiraceae",
  "Lachnospiraceae", "Lachnospiraceae",
  "Ruminococcaceae", "Oscillospiraceae",
  "Oscillospiraceae", "Rikenellaceae",
  "Rikenellaceae", "Tannerellaceae",
  "Tannerellaceae", "Akkermansiaceae",
  "Bifidobacteriaceae", "Bifidobacteriaceae",
  "Coriobacteriaceae", "Enterobacteriaceae",
  "Eubacteriaceae", "Odoribacteraceae",
  "Barnesiellaceae", "Veillonellaceae",
  "Veillonellaceae", "Streptococcaceae",
  "Clostridiaceae", "Turicibacteraceae",
  "Sutterellaceae", "Christensenellaceae"
)

# synthetic taxonomy (taxa) table
tax_df <- data.frame(
  asv = paste0("ASV", seq_along(species)),
  Family = family,
  Genus = sub(" .*", "", species),
  Species = species
)

# create 5 sample IDs
sample_ids <- paste0("EXAMPLE_", sprintf("%03d", 1:5))

# give some bacterial species higher abundances
weights <- c(
  12, 10, 8, 6, 7, 12, 8, 6, 6, 5,
  4, 3, 3, 5, 3, 4, 3, 2, 4, 3,
  4, 3, 2, 2, 1, 2, 2, 2, 1, 1,
  1, 1, 1, 1, 1
)

# generate random relative abundances
abundance <- replicate(length(sample_ids), {
  x <- rgamma(length(weights), shape = weights, rate = 1)
  x / sum(x) # Each sample sums to 1 (100%)
})

# ASV abundance table
asv_df <- data.frame(
  asv = tax_df$asv,
  abundance,
  check.names = FALSE
)

names(asv_df)[-1] <- sample_ids

# fake participant metadata
meta_df <- data.frame(
  sample_id = sample_ids,
  Subject_ID = paste0("Subject_", 1:5),
  Visit = c(1, 1, 2, 1, 2)
)

# save the tab-separated text files
write.table(asv_df, "gut_asv.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(tax_df, "gut_tax.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(meta_df, "metadata.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("created synthetic microbiome data\n")
