# gene_netwoRks.R
# contact: pdeweird@broadinstitute.org
# Use biogrid, STRING and depmap to create gene interaction networks
library(tidygraph)
library(ggraph)
library(here)
library(tidyverse)
library(STRINGdb)
library(KEGGREST)
library(igraph)

## Biogrid

get_biogrid_ppi <- function(gene_tibble, throughput = 'any') {
  # Query biogrid for interactions between genes
  #
  # Args:
  #   genes: list of query genes
  #   throughput: experiment type. Can be high, low, or any
  piped_genes <- paste(gene_tibble$`Gene ID`, collapse = '|')
  ppi_publications <- read_tsv(url(
    paste('https://webservice.thebiogrid.org/interactions/?searchIDs=true&geneList=',
          piped_genes,'&taxId=9606&includeInteractors=false&includeHeader=true&throughputTag=', 
          throughput,
          'format=tab2&,selfcitationsExcluded=false&accesskey=ef396dc62ddfc51b05547bf9dd6feb32',
          sep = ''))) 
  publication_count <- ppi_publications %>%
    mutate(A = ifelse(`Entrez Gene Interactor A` < `Entrez Gene Interactor B`, 
                      `Entrez Gene Interactor A`, `Entrez Gene Interactor B`),
           B = ifelse(`Entrez Gene Interactor A` > `Entrez Gene Interactor B`, 
                      `Entrez Gene Interactor A`, `Entrez Gene Interactor B`)) %>%
    group_by(A, B) %>%
    summarise(publications = n())
  graph_tbl <- as_tbl_graph(publication_count) %>%
    activate(nodes) %>%
    inner_join(gene_tibble, by = c('name' = 'Gene ID')) %>%
    mutate(`Gene ID` = name, name = `Gene Symbol`) %>%
    select(-`Gene Symbol`) %>%
    activate(edges) %>%
    mutate(weight = publications)
  return(graph_tbl)
}

plot_biogrid_ppi <- function(network) {
  graph = ggraph(network, layout = 'fr') +
    geom_edge_link(color = 'grey') +
    geom_node_point(color = 'darkgrey', fill = 'white', size = 3, pch = 21) +
    geom_node_text(aes(label = name), size = 3) +
    theme_graph() 
  return(graph)
}

## Depmap

load_depmap_data <- function() {
  depmap_data <- read_csv(here('networks', 'depmap', 
                               'avana-public-tentative-19q4_v2-gene-effect.csv')) %>%
    rename(cell = X1)
  return(depmap_data)
}

get_gene_overlap <- function(reference, genes) {
  # return reference genes that are in a gene set of interest
  overlap <- reference[reference %in% genes]
}

get_hugo_annotion <- function(gene, hugo_annotations) {
  mapped_symbol = hugo_annotations$`Approved symbol`[
    hugo_annotations$Input == gene]
  if (length(mapped_symbol) > 1) {
    mapped_symbol = first(mapped_symbol)
  } else if (length(mapped_symbol) == 0) {
    mapped_symbol = gene
  }
  return(mapped_symbol)
}

extract_depmap_data <- function(genes) {
  # extract the relvant columns from the depmap data
  depmap_data <- load_depmap_data()
  depmap_ids <- str_extract(colnames(depmap_data), '(?<=\\()[0-9]*')
  relevant_columns <- depmap_ids %in% genes$`Gene ID`
  relevant_depmap_ids <- depmap_ids[relevant_columns]
  depmap_id_df <- tibble(`Gene ID` = relevant_depmap_ids) %>%
    mutate(row = row_number())
  joined_id_genes_df <- inner_join(depmap_id_df, genes) %>%
    arrange(row)
  relevant_depmap_gene_names <- joined_id_genes_df$`Gene Symbol`
  relevant_data <- depmap_data[relevant_columns]
  colnames(relevant_data) <- relevant_depmap_gene_names
  diff = nrow(genes) - ncol(relevant_data)
  if (diff > 0) {
    warning(paste('Could not find', as.character(diff), 'genes:', 
                  paste(setdiff(genes$`Gene Symbol`, relevant_depmap_gene_names), collapse = ', ')))
  }
  return(relevant_data)
}

infer_depmap_cor_net <- function(genes, cutoff = 0.2) {
  # Use a defined correlation cutoff to infer a 
  # network from genes in the depmap data.
  # 
  # Args: 
  #   genes: list of genes
  #   cutoff: correlation cutoff
  #   abs: Use absolute value to infer edges
  # 
  # Returns: 
  #   network object from tidygraph
  relevant_data <- extract_depmap_data(genes)
  correlation_mat <- cor(relevant_data, use ='pairwise.complete.obs')
  edges <- abs(correlation_mat) > cutoff
  graph_matrix <- correlation_mat*edges
  diag(graph_matrix) <- 0
  graph_tbl <- as_tbl_graph(graph_matrix, directed = FALSE) %>%
    activate(edges) %>%
    mutate(cor = weight, abs_cor = abs(weight), weight = abs_cor) %>%
    activate(nodes) %>%
    inner_join(genes, by = c('name' = 'Gene Symbol'))
  return(graph_tbl)
}

plot_stringdb_net <- function(network, node_weight, names = TRUE, 
                              node_range = c(4,7), edge_range = c(0.3,2)) {
  
  lay = 'fr'
  graph = ggraph(network, layout = lay) +
    geom_edge_link(aes(width = weight), color = 'grey') +
    scale_edge_width(range = edge_range, breaks = c(0.4,0.65,0.9)) +
    theme_void() +
    theme(text = element_text(size = 10, family = 'Arial')) +
    labs(edge_width = 'Combined Score') +
    geom_node_point(aes(size = abs(!!as.name(node_weight)), 
                        color = as.factor(sign(!!as.name(node_weight)))), 
                    fill = 'white', pch = 21) +
    scale_color_manual(values = c(`-1` = '#ef8a62', `1` = '#67a9cf')) +
    scale_size(range = node_range, breaks = c(4,6)) + 
    labs(color = 'Node Sign', size = 'Score') +
    guides(edge_width = guide_legend(order = 1, ncol = 1), 
           size = guide_legend(order = 2, ncol = 2), 
           color = guide_legend(order = 3, ncol = 2)) 
  if (names) {
    graph = graph + geom_node_text(aes(label = name), size = 3, 
                   repel = TRUE, point.padding = 0, min.segment.length = 0.1) 
  } else {
    graph = graph + geom_node_text(aes(label = name), size = 0.5,
                                   repel = TRUE, point.padding = 0, 
                                   box.padding = 0.01, min.segment.length = 0.1) 
  }
  return(graph)
}

plot_depmap_cor_net <- function(network, node_weight = FALSE, names = TRUE, 
                                node_range = c(4,7), edge_range = c(0.3,2), 
                                all = FALSE) {
  lay = 'fr'
  graph <- ggraph(network, layout = lay) +
    geom_edge_link(aes(color = as.factor(sign(cor)), 
                       width = abs_cor)) +
    scale_edge_width(range = edge_range, breaks = c(0.2,0.4,0.6)) +
    scale_edge_color_manual(values = c(`-1` = 'pink', `1` = 'skyblue')) +
    theme_void() +
    theme(text = element_text(size = 10, family = 'Arial')) +
    labs(edge_width = 'Correlation', edge_color = 'Edge Sign') 
  if (node_weight != FALSE) {
    graph <- graph + geom_node_point(aes(size = abs(!!as.name(node_weight)), 
                                         color = as.factor(sign(!!as.name(node_weight)))), 
                                     fill = 'white', pch = 21) +
      scale_color_manual(values = c(`-1` = '#ef8a62', `1` = '#67a9cf')) +
      scale_size(range = node_range, breaks = c(2,4,6)) + 
      labs(color = 'Node Sign', size = 'Score') +
      guides(edge_width = guide_legend(order = 1, ncol = 2), 
             size = guide_legend(order = 2, ncol = 2), 
             edge_color = guide_legend(order = 3, ncol = 2), 
             color = guide_legend(order = 4, ncol = 2))
      
  } else {
    graph <- graph + geom_node_point(color = 'darkgrey', fill = 'white', 
                                     pch = 21) 
  }
  if (names) {
    graph <- graph + geom_node_text(aes(label = name), size = 3, 
                                    repel = TRUE, point.padding = 0, 
                                    min.segment.length = 0.1)
  } else {
    graph = graph + geom_node_text(aes(label = name), size = 0.5,
                                   repel = TRUE, point.padding = 0, 
                                   box.padding = 0.01, 
                                   min.segment.length = 0.1) 
  }
  return(graph)
}

plot_combined_network <- function(network, node_weight) {
  graph <- ggraph(network, layout = 'kk') +
    geom_edge_fan(aes(color = Source), width = 1) +
    scale_edge_color_manual(values = c(`string` = 'grey', `depmap` = '#8da0cb')) +
    #scale_edge_linetype_manual(values = c(`string` = 'solid', `depmap` = 'dotdash')) +
    theme_graph(base_family = 'Helvetica') +
    geom_node_point(aes(size = abs(!!as.name(node_weight)), 
                        color = as.factor(sign(!!as.name(node_weight)))), 
                    fill = 'white', pch = 21) +
    scale_color_manual(values = c(`-1` = 'red', `1` = 'blue')) +
    scale_size(range = c(4,7)) + 
    labs(color = 'Node Sign', size = 'Score') +
    guides(size = guide_legend(order = 1),
           color = guide_legend(order = 2)) +
    geom_node_text(aes(label = name), size = 3, 
                   repel = TRUE, point.padding = 0) 
  return(graph)
}

plot_combined_weighted_network <- function(network, node_weight) {
  graph = ggraph(network, layout = 'gem') +
    geom_edge_link(aes(color = depmap_delta, width = weight)) +
    labs(edge_color = 'Depmap Delta', edge_width = 'Combined Score') +
    scale_edge_width(range = c(0.3,1.5), limits = c(0,1000)) +
    scale_edge_color_gradient(low = 'grey', high = 'royalblue', limits = c(0,500)) +
    theme_graph(base_family = 'Helvetica') +
    geom_node_point(aes(size = abs(!!as.name(node_weight)), 
                        color = as.factor(sign(!!as.name(node_weight)))), 
                    fill = 'white', pch = 21) +
    scale_color_manual(values = c(`-1` = 'red', `1` = 'blue')) +
    scale_size(range = c(4,7)) + 
    labs(color = 'Node Sign', size = 'Score') +
    guides(size = guide_legend(order = 1, ncol = 2),
           color = guide_legend(order = 2), 
           edge_color = guide_legend(order = 3, ncol = 2), 
           edge_width = guide_legend(order = 4, ncol = 2)) +
    geom_node_text(aes(label = name), size = 3, 
                   repel = TRUE, point.padding = 0)
  return(graph)
}

get_depmap_null_degree <- function(resamples, n_genes, cor_cutoff) {
  depmap_data <- load_depmap_data()
  degrees <- vector(mode = 'numeric', length = resamples*n_genes)
  for(i in 1:resamples) {
    random_columns <- sample(2:ncol(depmap_data), size = n_genes)
    random_data <- depmap_data[,random_columns]
    correlations <- cor(random_data, use = 'pairwise.complete.obs')
    edges <- (correlations > cor_cutoff)*upper.tri(correlations)
    degree_dist <- colSums(edges)
    degrees[((i-1)*n_genes + 1):(i*n_genes)] <- degree_dist
  }
  degree_df <- tibble(degree = degrees) %>%
    group_by(degree) %>%
    summarise(n = n()/resamples) 
  return(degree_df)
}

infer_depmap_huge_net <- function(genes) {
  # Use coessentiality correlations to create a network. 
  # Infer sparse network using the package 'huge'
  #
  # Args:
  #   genes - list of genes
  #   
  # Returns:
  #   tidygraph object
  relevant_data <- extract_depmap_data(genes)
  huge_estimations <- huge(scale(as.matrix(relevant_data)))
  huge_graph <- huge.select(huge_estimations, criterion = 'ric')
  graph_tbl <- as_tbl_graph(as.matrix(huge_graph$refit), directed = FALSE) %>%
    activate(nodes) %>%
    mutate(name = colnames(relevant_data))
  return(graph_tbl)
}


plot_depmap_huge_net <- function(network) {
  graph <- ggraph(network, layout = 'fr') +
    geom_edge_link(color = 'grey') +
    geom_node_point(color = 'darkgrey', fill = 'white', size = 3, pch = 21) +
    geom_node_text(aes(label = name), size = 3) +
    theme_graph()
  return(graph)
}

## STRING

get_stringdb_net <- function(genes, cutoff = 150) {
  # Use the stringdb R package to get a network of string interactions from a 
  # gene list
  #
  # Args:
  #   genes: list of genes
  #   cutoff: stringdb score to cut the network off at
  #
  # Returns:
  #   network object from tidygraph
  string_db <- STRINGdb$new(version = '10', species = 9606, 
                            score_threshold = cutoff)
  gene_df <- data.frame(gene = genes$`Gene Symbol`)
  string_id_df <- string_db$map(gene_df, 'gene') %>%
    drop_na()
  interactions <- string_db$get_interactions(string_id_df$STRING_id)
  graph_tbl <- as_tbl_graph(interactions, directed = FALSE) %>%
    activate(nodes) %>%
    right_join(string_id_df, by = c('name' = 'STRING_id')) %>%
    mutate(string_id = name, name = gene) %>%
    inner_join(genes, by = c('name' = 'Gene Symbol')) %>%
    select(-gene) %>%
    activate(edges) %>%
    mutate(weight = combined_score)
  return(graph_tbl)
}

get_string_null_degree <- function(resamples, n_genes, cutoff) {  
  print('Setting Up')
  string_db <- STRINGdb$new(version = '10', species = 9606, 
                            score_threshold = cutoff)
  all_aliases <- string_db$get_aliases()
  degrees <- vector(mode = 'numeric', length = resamples*n_genes)
  print('Generating Null')
  for (i in 1:resamples) {
    random_columns <- sample(1:nrow(all_aliases), size = n_genes)
    random_string_tibble <- data.frame(all_aliases[random_columns,]) %>%
      drop_na()
    rownames(random_string_tibble) <- NULL
    random_string_ids <- random_string_tibble$STRING_id
    interactions <- string_db$get_interactions(random_string_ids)
    
    degree_tibble = as_tbl_graph(interactions, directed = FALSE) %>%
      activate(nodes) %>%
      mutate(degree = centrality_degree()) %>%
      as_tibble()
    full_degree_tibble = right_join(degree_tibble, random_string_tibble, 
                                    by = c('name' = 'STRING_id')) %>%
      mutate(degree = replace_na(degree, 0))
    degrees[((i-1)*n_genes + 1):(i*n_genes)] <- full_degree_tibble$degree
  }
  degree_df <- tibble(degree = degrees) %>%
    group_by(degree) %>%
    summarise(n = n()/resamples) 
  return(degree_df)
}

# Network Utitilites

get_edgelist <- function(graph) {
  # Extract an edge list from tidygraph object 
  edge_tib <- graph %>%
    activate(edges) %>%
    as_tibble()
  node_tib <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(index = row_number())
  edgelist <- inner_join(edge_tib, node_tib, by = c('from' = 'index')) %>%
    inner_join(node_tib, by = c('to' = 'index'), suffix = c('.from', '.to')) %>%
    select(name.from, name.to) %>%
    mutate(A = ifelse(name.from < name.to, name.from, name.to), 
           B = ifelse(name.from > name.to, name.from, name.to),
           combo = paste(A, B, sep = '_'))
  return(edgelist)
}

get_joined_edgelists <- function(network_list, network_names) {
  network_edgelists = list()
  for (i in 1:length(network_list)) {
    name = network_names[[i]]
    network <- network_list[[i]] 
    edgelist <- get_edgelist(network)
    network_edgelist <- edgelist %>%
      mutate(!!name := 1)
    network_edgelists[[i]] = network_edgelist
  }
  joined_edgelists <- reduce(network_edgelists, full_join,
                             by = c('name.from', 'name.to')) %>%
    mutate(combined_name = paste(name.from, name.to, sep = '_')) %>%
    dplyr::select(-c(name.from, name.to)) %>%
    arrange_at(vars(-combined_name)) %>%
    mutate(row = row_number(), 
           combined_name = fct_reorder(as.factor(combined_name), row)) %>%
    dplyr::select(-row)
  
  return(joined_edgelists)
}

## KEGG

get_kegg_network <- function(gene_tibble, cutoff = 200) {
  kegg_pathways = keggLink('pathway', 'hsa')
  kegg_pathways_tibble <- tibble(`Gene ID` = word(names(kegg_pathways), 2, sep = ':'),
                                 pathway = kegg_pathways) %>%
    group_by(pathway) %>%
    filter(n() < cutoff) %>%
    ungroup()
  genes_of_interest_tibble <- kegg_pathways_tibble %>%
    filter(`Gene ID` %in% gene_tibble$`Gene ID`)
  kegg_gene_connections_tibble <- inner_join(genes_of_interest_tibble, 
                                             genes_of_interest_tibble, 
                                             by = 'pathway') %>%
    select(`Gene ID.x`, `Gene ID.y`) %>%
    mutate(A = ifelse(`Gene ID.x` < `Gene ID.y`, `Gene ID.x`, `Gene ID.y`),
           B = ifelse(`Gene ID.x` > `Gene ID.y`, `Gene ID.x`, `Gene ID.y`)) %>%
    select(A, B) %>%
    distinct() %>%
    filter(A != B)
  kegg_network <- as_tbl_graph(kegg_gene_connections_tibble, directed = FALSE) %>%
    activate(nodes) %>%
    inner_join(gene_tibble, by = c('name' = 'Gene ID')) %>%
    mutate(`Gene ID` = name, name = `Gene Symbol`) %>%
    select(-`Gene Symbol`)
  return(kegg_network)
}

get_edge_meta <- function(net) {
  edge_tib <- net %>%
    activate(edges) %>%
    as_tibble()
  node_tib <- net %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(index = row_number())
  edgelist <- inner_join(edge_tib, node_tib, by = c('from' = 'index')) %>% 
    inner_join(node_tib, by = c('to' = 'index'), suffix = c('.from', '.to')) %>%
    mutate(A = ifelse(name.from < name.to, name.from, name.to), 
           B = ifelse(name.from > name.to, name.from, name.to))
  return(edgelist)
}

rank_edges <- function(network, ranked_hits) {
  ranked_edgelist <- get_edge_meta(network) %>%
    inner_join(ranked_hits, by = c('A' = 'Gene Symbol')) %>%
    inner_join(ranked_hits, by = c('B' = 'Gene Symbol'), suffix = c('.A', '.B')) %>%
    mutate(max_rank = pmax(rank.A, rank.B)) 
}

get_ranked_edge_density <- function(network, gene_hits) {
  ranked_hits <- gene_hits %>%
    mutate(rank = rank(-abs(all), ties.method = 'min'))
  edges_ranked <- rank_edges(network, ranked_hits)
  n_hits <- nrow(ranked_hits)
  densities <- vector(mode = 'numeric', n_hits - 1)
  for (i in 2:n_hits) {
    subnetwork <- edges_ranked %>%
      filter(max_rank <= i)
    subnetwork_density <- nrow(subnetwork)/choose(i, 2)
    densities[i - 1] = subnetwork_density
  }
  null_density <- nrow(edges_ranked)/choose(n_hits, 2)
  ranked_edge_density <- tibble(rank = 2:n_hits, density = densities) %>%
    mutate(normed_density = density/null_density)
  return(list(density = ranked_edge_density, ranked_edges = edges_ranked))
}

plot_clusters <- function(clustered_net, gene, z_cutoff, edge_cutoff, database) {
  clusters = clustered_net %>% 
    activate(nodes) %>%
    as_tibble() %>%
    group_by(group) %>% 
    filter(n() > 1) %>% 
    select(group)  %>% 
    distinct() %>%
    arrange(group)
  curr_dir = here('manuscript_figures', gene, 'clusters')
  dir.create(curr_dir)
  for (interesting_cluster in clusters$group) {
    print(interesting_cluster)
    subnetwork = clustered_net %>% filter(group == interesting_cluster)
    if (database == 'depmap') {
      p = plot_depmap_cor_net(subnetwork,'all') 
    } else if (database == 'string') {
      p = plot_stringdb_net(subnetwork, 'all') 
    } else if (database == 'depmap_string_simple') {
      p = plot_combined_network(subnetwork, 'all') 
    } else if (database == 'depmap_string_weighted') {
      p = plot_combined_weighted_network(subnetwork, 'all')
    }
    ggsave(file.path(curr_dir, paste(gene, '_node_cutoff_',
                                           as.character(z_cutoff), 
                                           '_edge_cutoff_', as.character(edge_cutoff),
                                           database, '_cluster_', 
                                     as.character(interesting_cluster), 
                                     '.svg', sep = '')), p, 
           width = 8.7, height = 7.7, units = 'cm')
    
  }
}

get_clusters <- function(network, hits, by = c('name' = 'Gene Symbol', 'Gene ID')) {
  clusters <- network %>%
    activate(nodes) %>%
    inner_join(hits, by = by) %>%
    mutate(group = group_louvain(weights = weight)) 
  barcodes <- clusters %>%
    activate(nodes) %>%
    as_tibble() %>%
    select(`group`,n_guides:all) %>%
    select(-n_guides) %>%
    gather(assay, score, -group)
  p = ggplot(barcodes %>% group_by(group, assay) %>% filter(n() > 2)) +
    aes(x = score, y = assay, color = assay) +
    geom_point(pch = 124, size = 2) +
    facet_wrap('group') +
    theme_classic() + 
    guides(color = FALSE) +
    geom_vline(xintercept = 0, linetype = 'dashed')
  return(list(clusters = clusters, barcode_p = p))
}

get_cluster_enrichment <- function(network) {
  string_db <- STRINGdb$new(version = '10', species = 9606)
  network_groups <- network %>%
    group_by(group) %>%
    filter(n() > 1)
  nodes <- network_groups %>% activate(nodes) %>% as_tibble() %>% data.frame()
  string_mapped_nodes <- nodes %>%
    string_db$map('name')
  enriched_clusters <- string_mapped_nodes %>% 
    group_by(group) %>%
    nest() %>%
    mutate(bound_enrichment = map(data,
                                  function(df) {
                                    comp = string_db$get_enrichment(df$STRING_id, category = 'Component')
                                    proc = string_db$get_enrichment(df$STRING_id, category = 'Process')
                                    func = string_db$get_enrichment(df$STRING_id, category = 'Function')
                                    bound_enrichment = bind_rows(list(comp, proc, func))
                                    return(bound_enrichment)
                                  })
    ) %>%
    unnest(bound_enrichment)
  most_enriched_terms <- enriched_clusters %>%
    filter(proteins < 100) %>%
    group_by(group) %>%
    top_n(3, -pvalue) %>%
    arrange(group, pvalue)
  return(most_enriched_terms)
}

plot_cluster_barcodes <- function(clusters) {
  network_nodes <- clusters %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(community = as.factor(group),
           community = fct_reorder(community, -abs(all), mean))
  
  p = ggplot(network_nodes %>% group_by(group) %>% filter(n() > 1)) +
    aes(x = community, y = all) +
    geom_point(pch = 95, size = 8, alpha = 0.5) +
    geom_hline(yintercept = 2, linetype = 'dashed') +
    geom_hline(yintercept = -2, linetype = 'dashed') +
    theme_classic() +
    theme(text = element_text(size = 12)) +
    ylab('Avg. Z-score') +
    xlab('Community')
  return(p)
}

get_long_degree <- function(observed, null, resamples) {
  observed_long = map_df(observed, rep, observed$n)
  null_long = map_df(null, rep, null$n*resamples)
  long_degree = bind_rows(observed_long, null_long) %>%
    select(-n)
  return(long_degree)
}

get_long_degree_ks <- function(long_degree) {
  ks_test = tidy(ks.test(long_degree$degree[long_degree$type == 'observed'], 
                         long_degree$degree[long_degree$type == 'null'], 
                         alternative = 'less'))
  return(ks_test)
}

plot_degree_dist <- function(network, type, resamples, cutoff, gene) {
  degree_distribution <- network %>%
    activate(nodes) %>%
    mutate(degree = centrality_degree(loops = FALSE)) %>%
    as_tibble() %>%
    group_by(degree) %>%
    summarise(n = n()) %>%
    mutate(type = 'observed')
  if (type == 'STRING') {
    null_degree <- get_string_null_degree(resamples, with_graph(network, graph_order()), cutoff) %>%
      mutate(type = 'null')
  } else if (type == 'co-essentiality') {
    null_degree <- get_depmap_null_degree(resamples, with_graph(network, graph_order()), 
                                          cutoff) %>%
      mutate(type = 'null')
  }
  
  bound_degree_dists <- bind_rows(degree_distribution, null_degree) %>%
    group_by(type) %>%
    mutate(fraction = n/sum(n)) %>%
    ungroup()
  
  long_degree = get_long_degree(degree_distribution, null_degree, resamples)
  ks = get_long_degree_ks(long_degree)
  ks_chars = paste('ks: ', as.character(format(ks$statistic, digits = 3)), '\n',
                   'p: ',  as.character(format(ks$p.value, scientific = TRUE, digits = 2)), 
                   sep = '')
  label_info = tibble(label = ks_chars, degree = max(bound_degree_dists$degree),
                      fraction = max(bound_degree_dists$fraction))
  
  degree_distribution <- ggplot(bound_degree_dists) +
    aes(x = degree, y = fraction, fill = type) +
    geom_bar(stat = 'identity', position = 'identity', alpha = 0.5) +
    scale_fill_manual(values = c('#d95f02', '#7570b3')) +
    theme_grey() +
    theme(text = element_text(size = 10, family = 'Arial')) +
    ggtitle(paste(gene, type)) +
    xlab('Number of Interactors') +
    ylab('Fraction of Genes') +
    theme(legend.key.size = unit(0.4, 'cm')) +
    facet_zoom(xlim = c(0, 15), zoom.size = 3) +
    geom_text(aes(label = label, fill = NULL), data = label_info, hjust = 1, vjust = 1, 
              size = 3) 
  return(degree_distribution)
}

