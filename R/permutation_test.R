library(tidyverse)

# code for the permutation test

in_file <- "Table_S2.csv"

sg <- read_csv(in_file, guess_max=10000) %>%
  filter(plate_type != "diagnostic" & !is.na(group_summary) & !is.na(S_SC_STc)) %>%
  select(species_abbv, group_summary)

s_all <- sg$species_abbv
g_all <- sg$group_summary

sg_r <- sg %>%
  group_by(species_abbv, group_summary) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  spread(group_summary, count) %>%
  gather("group_summary", "count", -species_abbv) %>%
  mutate(count=replace_na(count, 0))

l <- nrow(sg)
reps <- 10000

sg_base_r <- tibble(species_abbv=rep(sg_r$species_abbv, reps),
                    group_summary=rep(sg_r$group_summary, reps),
                    rep=rep(1:reps, each=nrow(sg_r)))

r_ran <- rep(1:reps, each=l)
s_ran <- rep(s_all, reps)
g_ran <- as.vector(replicate(reps, sample(g_all, replace=FALSE)))

sg_ran_r <- tibble(rep=r_ran, species_abbv=s_ran, group_summary=g_ran) %>%
  group_by(rep, species_abbv, group_summary) %>%
  summarise(count=n()) %>%
  left_join(sg_base_r, .) %>%
  mutate(count=replace_na(count, 0))

ggplot() +
  geom_histogram(data=sg_ran_r, aes(x=count), bins=20) +
  geom_vline(data=sg_r, aes(xintercept=count), linetype="dashed", colour="red") +
  facet_grid(species_abbv~group_summary, scales="fixed") +
  theme(strip.text.x=element_text(angle=90, hjust=0),
        strip.text.y=element_text(angle=0, hjust=0))

sg_ran_r_c <- sg_ran_r %>%
  rename(ran_count=count) %>%
  left_join(., sg_r) %>%
  rename(obs_count=count)

p_vals <- sg_ran_r_c %>%
  group_by(species_abbv, group_summary) %>%
  summarise(n_gr=sum(obs_count <= ran_count),
            p_gr=n_gr/10000,
            n_le=sum(obs_count >= ran_count),
            p_le=n_le/10000) %>%
  ungroup() %>%
  mutate(n=ifelse(n_gr < n_le, n_gr, n_le),
         p=ifelse(p_gr < p_le, p_gr, p_le),
         cat=ifelse(p_gr < p_le, "over", "under"),
         p_gr=ifelse(p_gr == 0, 1/10000, p_gr),
         p_le=ifelse(p_le == 0, 1/10000, p_le),
         p=ifelse(p == 0, 1/10000, p),
         p_adj=p.adjust(p, method="BH"),
         result=ifelse(p_adj < 0.025, "sig", "non-sig"),
         p_adj_lab=ifelse(result == "sig", round(p_adj, digits=4), NA)) %>%
  arrange(p_adj)
