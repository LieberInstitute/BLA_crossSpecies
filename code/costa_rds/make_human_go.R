library(org.Hs.eg.db)
library(GO.db)

make_go_human = function() {
  go_sets = as.list(org.Hs.egGO2ALLEGS)
  all_genes = unique(unlist(go_sets))
  entrez_to_symbols = unlist(mget(all_genes, org.Hs.egSYMBOL))
  go_sets = lapply(go_sets, function(gs) { unique(entrez_to_symbols[gs]) })
  names(go_sets) = sapply(names(go_sets), full_go_name)
  saveRDS(go_sets, "go_human.rds")
}

full_go_name = function(go_id) {
  go_term = GOTERM[[go_id]]
  return(paste(go_id, Term(go_term), Ontology(go_term), sep = "|"))
}

make_go_human()
