# Assign HUMAnN/MetaCyc pathway names to broad functional categories for custom GMT / g:Profiler.

pathway_gmt_id <- function(pathway_name) {
  gsub("[^A-Za-z0-9]+", "_", pathway_name)
}

classify_pathway_category <- function(pathway_name) {
  n <- tolower(pathway_name)
  rules <- list(
    SCFA_and_fermentation = "butanoate|butyr|propanoate|propion|acetate|ferment|succinate.*butano",
    Amino_acid_metabolism = "arginine|ornithine|lysine|glutam|aspart|proline|histidine|tryptoph|phenylalan|tyrosine|branched|isoleucine|leucine|valine|serine|glycine|alanine|citrulline|stickland",
    Nucleotide_metabolism = "purine|pyrimidine|nucleotide|guanosine|adenosine|thymidine|uridine|cytidine",
    Carbohydrate_degradation = "glucose|fructose|galactose|mannose|starch|glycogen|cellulose|chitin|fucose|rhamnose|xylose|arabinose|sucrose|maltose|trehalose|chondroitin|heparin",
    Energy_and_respiration = "respiration|oxidation|tca|glycolysis|gluconeogenesis|nadh|electron|cytochrome|atp synth",
    Cell_wall_and_peptidoglycan = "peptidoglycan|lipopolysaccharide|o.antigen|teichoic|murein",
    Lipid_and_fatty_acid = "fatty acid|lipid|cholesterol|bile|phospholipid|sphingo",
    Cofactor_and_vitamin = "thiamine|riboflavin|folate|biotin|cobalamin|menaquinone|ubiquinone|nad |nadp|coenzyme|porphyrin|heme"
  )
  for (lab in names(rules)) {
    if (grepl(rules[[lab]], n)) return(lab)
  }
  "Other_metabolism"
}
