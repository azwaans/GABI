#!/usr/bin/env Rscript

#function to write nexus
write_gestalt_nexus = function(cell_by_site_dat, output_file){
  #based on S.Seidel's write_nexus in SciPhy
  
  # Basic checks on input data
  if (!is.matrix(cell_by_site_dat) && !is.data.frame(cell_by_site_dat)) {
    stop("Error: 'cell_by_site_dat' must be a matrix or a data frame.")
  }
  
  # Test if we can write to the output file or if it exists
  if (file.exists(output_file) && file.access(output_file, mode = 2) != 0) {
    stop("Error: The output file is not writable.")
  }
  
  n_taxa = nrow(cell_by_site_dat)
  tax_labels = cell_by_site_dat$cell_id
  n_sites = length(unlist(strsplit(cell_by_site_dat$barcode[1],"-")))
  
  #Write nexus header
  write(x = "#NEXUS", output_file)
  write(x = " ", output_file, append = T)
  

  
  # write data block
  write(x="begin characters;", output_file, append = T)
  write(x=paste0("dimensions nchar=", n_sites, ";"), output_file, append=T)
  write(x="format datatype=gstaltData;", output_file, append = T)
  write(x="matrix", output_file, append = T) 
  
  for (taxon_label in tax_labels){
    sequence_concatenated = gsub("-",",",cell_by_site_dat[cell_by_site_dat$cell_id==taxon_label,'barcode'])
    write(x = paste(taxon_label, sequence_concatenated), output_file, append = T)
  }
  write(x=";", output_file, append = T)
  write(x="end;", output_file, append = T)
}

# Get command-line arguments
args = commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript convert_cell_id_barcode_to_beast_input.R <input_data_file> <output_file>")
}

# Read arguments
input_data_file = args[1]
output_file = args[2]

# Load the data (assuming it's in CSV format, but adjust as needed)
cell_id_barcode = read.csv(input_data_file)

write_gestalt_nexus(cell_id_barcode, output_file)
