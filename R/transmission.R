library(tidyverse)
library(tidytree)
library(ggraph)
library(tidygraph)

cs <- 10

# code for the transmission analysis
# takes two inputs - metadata and distances

meta <- tibble(id=character(), S_SC_STc=character(), gp_specific=character(), gp=character(), gp_abbv=character())
dis <- tibble(id_1=character(), id_2=character(), Core=integer())

dis <- dis %>%
  filter(id_1 != "Reference" & id_2 != "Reference") %>%
  filter(id_1 != id_2) %>%
  rowwise() %>%
  mutate(id=paste(sort(c(id_1, id_2)), collapse="___")) %>%
  ungroup() %>%
  arrange(id) %>%
  separate(id, c("id_f", "id_s"), sep="___") %>%
  filter(id_1 == id_f & id_2 == id_s) %>%
  select(id_1, id_2, Core) %>%
  filter(Core <= cs)

dis_meta <- dis %>%
  left_join(., meta, by=c("id_1"="id")) %>%
  rename(gp_1=gp) %>%
  rename(gp_abbv_1=gp_abbv) %>%
  rename(gp_specific_1=gp_specific) %>%
  rename(S_SC_STc_1=S_SC_STc) %>%
  left_join(., meta, by=c("id_2"="id")) %>%
  rename(gp_2=gp) %>%
  rename(gp_abbv_2=gp_abbv) %>%
  rename(gp_specific_2=gp_specific) %>%
  rename(S_SC_STc_2=S_SC_STc) %>%
  filter(!is.na(gp_1) & !is.na(gp_2)) %>%
  filter(! is.na(S_SC_STc_1) & ! is.na(S_SC_STc_2)) %>%
  filter(S_SC_STc_1 == S_SC_STc_2) %>%
  select(-S_SC_STc_2) %>%
  rename(S_SC_STc=S_SC_STc_1) %>%
  rowwise() %>%
  mutate(gp=paste(sort(c(gp_1, gp_2)), collapse="~~~")) %>%
  ungroup()

sc_count <- meta %>%
  filter(S_SC_STc %in% dis_meta$S_SC_STc) %>%
  group_by(S_SC_STc) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  arrange(desc(count))

all_events <- NULL
all_events_detailed <- NULL
all_events_super_detailed <- NULL
for(sc in sc_count$S_SC_STc){
  
  sc_ids_clus <- NULL
  
  n_list <- meta %>%
    filter(sc == S_SC_STc)
  
  e_list <- dis_meta %>%
    filter(S_SC_STc == sc)
  
  ne <- tbl_graph(nodes=n_list, edges=e_list, directed=FALSE)
  
  ne_com <- ne %>%
    morph(to_components)
  
  sc_events <- 0
  e_list_events <- NULL
  n_list_ids_clus <- NULL
  for(com_i in 1:length(ne_com)){
    
    n_list_com <- ne_com[[com_i]] %>%
      activate(nodes) %>%
      as_tibble()
    
    n_list_com_ids_clus <- n_list_com %>%
      select(id) %>%
      mutate(cluster=com_i,
             sketch_size=ss,
             core_snps=cs,
             cat=paste(ss, cs, sep="_"))
    
    n_list_ids_clus <- bind_rows(n_list_ids_clus, n_list_com_ids_clus)
    
    e_list_com_events <- e_list %>%
      filter(id_1 %in% n_list_com$id & id_2 %in% n_list_com$id) %>%
      mutate(spark_id_1=str_replace(id_1, "_C\\d+", ""),
             spark_id_2=str_replace(id_2, "_C\\d+", "")) %>%
      filter(spark_id_1 != spark_id_2) %>%
      filter(gp_specific_1 != gp_specific_2) %>%
      select(-spark_id_1, -spark_id_2) %>%
      mutate(events="no")
    
    if(nrow(e_list_com_events) > 0){
      e_list_com_events <- e_list_com_events %>%
        group_by(gp) %>%
        sample_n(1) %>%
        ungroup() %>%
        mutate(events="yes")
    }
    
    sc_events <- sc_events + nrow(e_list_com_events)
    
    e_list_events <- bind_rows(e_list_events, e_list_com_events)
  }
  
  sc_ids_clus <- bind_rows(sc_ids_clus, n_list_ids_clus)
  
  all_events_tmp <- e_list_events %>%
    select(S_SC_STc, gp)
  
  all_events <- bind_rows(all_events, all_events_tmp)
  all_events_detailed <- bind_rows(all_events_detailed, e_list_events)
  
  e_list <- left_join(e_list, e_list_events) %>%
    mutate(events=ifelse(is.na(events), "no", events)) %>%
    arrange(events) %>%
    mutate(events=factor(events, levels=unique(events)))
  
  all_events_super_detailed <- bind_rows(all_events_super_detailed, e_list)
  
}

# dis_ran <- tibble(id_1=sample(meta$id, 500, replace=TRUE),
#        id_2=sample(meta$id, 500, replace=TRUE),
#        Core=sample(1:20, 500, replace=TRUE))

# meta_ran <- tibble(id=paste("ID_", sample(1:20, 20, replace=TRUE), "_C1", sep=""),
#           S_SC_STc=paste("SC_", sample(1:4, 20, replace=TRUE), sep=""),
#           gp_specific=paste(sample(c("cow", "pig", "sheep"), 20, replace=TRUE), "_", sample(1:4, 20, replace=TRUE), sep=""),
#           gp=str_replace(gp_specific, "_\\d+", ""),
#           gp_abbv=str_extract(gp, "^\\S"))
