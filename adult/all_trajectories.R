# Code to identify trajectory profiles using multiple imputation and the EM algorithm (flexmix)
# Wednesday July 30, 2025

# load packages and set preferences via conflicts package

curr_fxns <- as.vector(lsf.str())

pckgs <- c(
  "glue", "conflicted", "tidyverse", "tibble"
)

all_pckgs <- .packages(all.available=T)

pckgs_yes <- pckgs[pckgs %in% all_pckgs]
pckgs_no <- pckgs[!pckgs %in% all_pckgs]

silent_load <- function(package){
  suppressPackageStartupMessages(do.call("library", list(package)))
  return()
}

no_ret <- sapply(pckgs_yes, silent_load)

conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::select())
conflicts_prefer(dplyr::lag())
conflicts_prefer(dplyr::first())
conflicts_prefer(dplyr::last())
conflicts_prefer(dplyr::summarize())
conflicts_prefer(lubridate::month)
conflicts_prefer(lubridate::year)
conflicts_prefer(lubridate::day)
conflicts_prefer(base::cbind())
conflicts_prefer(rlang::`:=`)
conflicts_prefer(tidyr::replace_na)

library(mice)
library(miceadds)
library(ggmice)
library(flexmix)
library(ggrepel)
library(gtools)
library(irr)

# The input dataset for this analysis is a "long" dataset, with one row per visit, and 5 visits/rows per participant.
# Overall there are 3659x5 = 18295 rows.
# The variables are as follows:
# # record_id           Unique patient identifier
# # visit_month_touse   3, 6, 9, 12, or 15
# # pasc_score          Long COVID Research Index, continuous
# # reinf               Indicator variable of whether the visit was between 7 days before and 30 days after an on-study reinfection
# # age_inf             Age at the time of first SARS-CoV-2 infection (years), continuous, not rounded to nearest year
# # biosex_an           Sex assigned at birth, binary
# # race_unique_an      Race/ethnicity, categorical
# # ref_type_an2        Referral source, categorical
# # gen_1H_inf          Hospitalization status during acute phase of first SARS-CoV-2 infection, binary
# # fvacc_index_inf_an  Vaccination status at first SARS-CoV-2 infection, binary
# # marital             Marital status, categorical
# # cc_imm_inf          Comorbidity: Immunocompromised condition, such as a transplant, HIV, or an immune deficiency
# # cc_autoimm_inf      Comorbidity: Autoimmune disease, binary
# # cc_cancer_inf       Comorbidity: Cancer, binary
# # cc_liver_inf        Comorbidity: Chronic liver disease, binary
# # cc_obesity_inf      Comorbidity: Obesity, binary
# # cc_diabetes_inf     Comorbidity: Diabetes, binary
# # cc_renal_inf        Comorbidity: Kidney disease, binary
# # cc_cvd_inf          Comorbidity: Cardiovascular disease, binary
# # cc_stroke_inf       Comorbidity: Stroke or TIA, binary
# # cc_asthma_inf       Comorbidity: Asthma, binary
# # cc_lung_inf         Comorbidity: Lung disease, binary
# # cc_dementia_inf     Comorbidity: Dementia, binary
# # cc_mh_inf           Comorbidity: Mental health disorder, binary
# # cc_fibromyalgia_inf Comorbidity: Chronic pain syndrome or fibromyalgia, binary
# # cc_cfs_inf          Comorbidity: Myalgic encephalomyelitis/chronic fatigue syndrome, binary
# # cc_pots_inf         Comorbidity: Postural orthostatic tachycardia syndrome, binary
# # cc_neurO_inf        Comorbidity: Neurological condition, binary
# # sd_homeless         Social determinant of health: Experiencing homelessness or transient housing, binary
# # sd_disability       Social determinant of health: Disability, binary
# # sd_unemploy         Social determinant of health: Unemployment, binary
# # sd_medicaid         Social determinant of health: Medicaid, binary
# # sd_uninsured        Social determinant of health: Uninsured, binary
# # sd_lostinsur        Social determinant of health: Lost insurance due to pandemic, binary
# # sd_income           Social determinant of health: 2019 household income, categorical
# # sd_moneyshort       Social determinant of health: Difficulty covering expenses in the previous month, binary
# # sd_docvisit         Social determinant of health: Most recent doctor's visit before index was within the last 5 years, binary
# # sd_skipcare         Social determinant of health: Skipping care in the 12 months before enrollment, binary
# # sd_food             Social determinant of health: Food insecurity in the 12 months before enrollment, binary
# # sd_english          Social determinant of health: English as primary language, binary
# # educ_coll           Social determinant of health: College degree or higher, binary

sdl <- readRDS("sdl_traj.rds") # named SDL for "scaled data long"

sdw <- sdl %>%
  pivot_wider(id_cols = c("record_id"), names_from = "visit_month_touse", 
              values_from = c("pasc_score"), names_glue = "pasc_score_{visit_month_touse}")

# Multiple imputation using 'mice'
meths <- make.method(sdl)

meths["pasc_score"] <- "2l.pmm"
meths[names(meths) %!in% c("record_id", "visit_month_touse", "pasc_score", "race_unique_an", "reinf")] <- "2lonly.function"

pm <- make.predictorMatrix(sdl)

pm[, "record_id"] <- -2
pm[c("record_id", "visit_month_touse", "race_unique_an", "reinf", "ref_type_an2"),] <- 0
pm[rownames(pm) %!in% c("pasc_score"), "visit_month_touse"] <- 0
pm[rownames(pm) %!in% c("pasc_score"), "reinf"] <- 0
pm


imps <- 25
set.seed(904)

data_mice <- mice(sdl, m = imps, method = meths, pred = pm, 
                  imputationFunction = list("age_inf" = "pmm",
                                            "biosex_an" = "logreg",
                                            "ref_type_an2" = "logreg",
                                            "gen_1H_inf" = "logreg",
                                            "fvacc_index_inf_an" = "polyreg",
                                            "marital" = "logreg",
                                            "cc_imm_inf" = "logreg",
                                            "cc_autoimm_inf" = "logreg",
                                            "cc_cancer_inf" = "logreg",
                                            "cc_liver_inf" = "logreg",
                                            "cc_obesity_inf" = "logreg",
                                            "cc_diabetes_inf" = "logreg",
                                            "cc_renal_inf" = "logreg",
                                            "cc_cvd_inf" = "logreg",
                                            "cc_stroke_inf" = "logreg",
                                            "cc_asthma_inf" = "logreg",
                                            "cc_lung_inf" = "logreg",
                                            "cc_dementia_inf" = "logreg",
                                            "cc_mh_inf" = "logreg",
                                            "cc_fibromyalgia_inf" = "logreg",
                                            "cc_cfs_inf" = "logreg",
                                            "cc_pots_inf" = "logreg",
                                            "cc_neuro_inf" = "logreg",
                                            "sd_homeless" = "logreg",
                                            "sd_disability" = "logreg",
                                            "sd_unemploy" = "logreg",
                                            "sd_medicaid" = "logreg",
                                            "sd_uninsured" = "logreg",
                                            "sd_lostinsur" = "logreg",
                                            "sd_income" = "polyreg",
                                            "sd_moneyshort" = "polyreg",
                                            "sd_docvisit" = "polyreg",
                                            "sd_skipcare" = "logreg",
                                            "sd_food" = "logreg",
                                            "sd_english" = "logreg",
                                            "educ_coll" = "logreg"
                  ),
                  cluster_var = "record_id"
)

nc_list <- vector()
sd_list <- list()
bics_list <- list()
aics_list <- list()

set.seed(827)

# Calculating AIC and BIC using different numbers of clusters and choosing the 
# number of clusters with the smallest average BIC across imputations
for (i in 1:imps) {
  
  dat_run <- mice::complete(data_mice, action = i)
  
  sd_list[[i]] = dat_run
  
  all_fm <- stepFlexmix(pasc_score ~ visit_month_touse + I(visit_month_touse^2) | record_id, 
                        model = FLXMRglm(family = "poisson"), verbose = F,
                        data = dat_run, k = 2:10, nrep = 1)
  
  bics <- BIC(all_fm)
  aics <- AIC(all_fm)
  
  bics_list[[i]] <- bics
  aics_list[[i]] <- aics
  
  nc <- as.numeric(names(bics[which(bics == min(bics))]))
  
  nc_list[i] <- nc
  
}

a <- rowMeans(sapply(aics_list, unlist))
b <- rowMeans(sapply(bics_list, unlist))

k <- as.numeric(names(b[which(b==min(b))]))

print("K:")
print(k)

# Assigning each participant to a cluster in each imputed dataset
clusters_across <- sdw %>%
  select(record_id)

set.seed(828)

for (i in 1:imps) {
  
  print(i)
  
  varname <- paste0("cluster_", i)
  
  scaled_data <- sd_list[[i]] 
  
  k_used = 0
  
  while(k_used != k) {
    
    tes <- flexmix(pasc_score ~ visit_month_touse + I(visit_month_touse^2) | record_id, 
                   model = FLXMRglm(family = "poisson"),
                   data = scaled_data, k = k)
    
    k_used <- tes@k
    
  }
  
  kmeans_clust <- scaled_data %>% mutate(cluster = clusters(tes)) %>%
    select(record_id, cluster) %>%
    group_by(record_id) %>%
    filter(row_number() == 1)
  
  avg_scores <- scaled_data %>% mutate(cluster = clusters(tes)) %>%
    group_by(cluster) %>%
    summarize(avg = mean(pasc_score))
  
  lowest4 <- avg_scores %>%
    mutate(rank = rank(avg)) %>%
    filter(rank <= 4) %>%
    pull(cluster)
  
  kmeans_clust2 <- kmeans_clust
  
  scaled_data$cluster <- clusters(tes)
  
  dat_clusts <- scaled_data %>% 
    group_by(cluster) %>%
    summarize(avg = mean(pasc_score)) %>%
    ungroup() %>%
    mutate(cluster2 = rank(avg),
           clus = toupper(letters[cluster2]))
  
  dat_out2 <- scaled_data %>%
    left_join(dat_clusts) %>%
    filter(visit_month_touse >= 3, visit_month_touse <= 15) %>%
    group_by(cluster, visit_month_touse) %>%
    summarize(med = median(pasc_score, na.rm = T),
              up = quantile(pasc_score, 0.75, na.rm = T),
              down = quantile(pasc_score, 0.25, na.rm = T)
    )
  
  clusters_across <- clusters_across %>%
    left_join(kmeans_clust2 %>% 
                rename(!! quo_name(varname) := cluster), by = join_by(record_id))
  
}

# Consensus-based approach to assign participants to clusters
imp_list <- list()
maxk_list <- vector()

clusters_across_ind <- clusters_across 

for (i in 1:imps) {
  
  print(i)
  varname <- paste0("cluster_", i)
  
  perms <- as.data.frame(t(permutations(k, k, 1:k))) %>%
    mutate(!! quo_name(varname) := V1)
  
  clusters_across2 <- clusters_across_ind %>%
    left_join(perms, by = join_by(!!varname))
  
  kappa_list <- vector()  
  
  for (j in 1:(ncol(perms) - 1)) {
    varj <- paste0("V", j)
    kappa_input <- clusters_across2 %>% select("cluster_1", !!varj)
    kappa_list[j] <- kappa2(kappa_input)$value
  }
  
  maxk <- which(kappa_list == max(kappa_list, na.rm = T))
  maxk_list[i] = maxk
  imp_list[i] <- clusters_across2[,paste0("V", maxk)]
  ids <- clusters_across2[,"record_id"]
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

t <- as.data.frame(do.call(cbind, imp_list)) %>%
  cbind(ids) %>%
  rowwise() %>%
  mutate(mode = getmode(c_across(starts_with("V"))))

clus_ds <- sdl %>%
  left_join(t %>% select(record_id, cluster = mode), by = "record_id")

# Note that clusters were re-ordered A-H for presentation purposes to be in order of decreasing average LCRI across time




