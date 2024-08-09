# aggregate all raw pi/dxy calculations from various software packages
# KMS 2020-04-02
# KLK edited 2020-05-15
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#library("lavaan")
#library("polycor")
#library("emmeans")
#library("dplyr")
library("tidyverse")
library("officer")
library("aod")
library("rvg")
library("ggrastr")
library("patchwork")
library("reshape2")
library("gridExtra")
library("cocor")
#library("Cairo")

######################################## 
# read in and format data for plots
######################################## 
sim_dat <- read_rds("data/sim_dat_all.rds") %>%
  mutate(vcf_source = gsub("_invar.*", "", vcf_source))

# expected pi
Ne <- 1e6
Mu <- 1e-8
#exp_pi <- 4*Ne*Mu
exp_pi <- (12*Ne*Mu)/(3+(16*Ne*Mu))
#exp_pi <- 0

# this is the maxmium value of pi for a VCF
# i.e. the "true" per site estimate of pi for that sample
# with zero missing data
# used for scaling
max_stats <- sim_dat  %>%
  filter(missing_type == "sites") %>%
  filter(missing_data == 0.00) %>% 
  mutate(max_wt = avg_watterson_theta) %>% 
##  mutate(max_dxy = avg_dxy) %>%
  select(vcf_source, method, max_wt)

sim_dat <-  sim_dat %>%
  left_join(max_stats) %>%
  filter(missing_type != "accuracy") %>%
  mutate(wt_scaled = avg_watterson_theta/max_wt) #%>%
##  mutate(dxy_scaled = avg_dxy/max_dxy)

######################################## 
# statistical tests 
######################################## 

#stats_tests_pi <- scikit_dat %>%
#  filter(missing_data < 1) %>%
#  group_by(method,missing_type) %>%
##  do(model = lm(.$wt_scaled ~ .$missing_data)) %>%
##  broom::tidy(model) %>%
##  filter(term != "(Intercept)") %>%
#  arrange(missing_type) %>%
##  select(-term) %>%
#  mutate(stat = "pi")

#stats_tests_dxy <- scikit_dat %>%
#  filter(missing_data < 1) %>%
#  filter(!is.na(dxy_scaled)) %>%
#  group_by(method,missing_type) %>%
#  filter(method != "VCFtools") %>%
#  do(model = lm(.$dxy_scaled ~ .$missing_data)) %>%
#  broom::tidy(model) %>%
#  filter(term != "(Intercept)") %>%
#  arrange(missing_type) %>%
#  select(-term) %>%
#  mutate(stat = "dxy")

#bind_rows(stats_tests_pi, stats_tests_dxy) %>%
#  write.csv(file = "figures/TableS2.csv", row.names = FALSE, quote = FALSE)

######################################## 
# Figure 2, simulated data
######################################## 


# unscaled pi for all methods
# maybe S1?

wt_dat <- sim_dat[sim_dat$method != "vcftools", ]
td_dat <- sim_dat[sim_dat$method != "pixy.sites", ]

wt_sites <- wt_dat[wt_dat$missing_type == "sites", ]
wt_genos <- wt_dat[wt_dat$missing_type == "genotypes", ]

td_sites <- td_dat[td_dat$missing_type == "sites", ]
td_genos <- td_dat[td_dat$missing_type == "genotypes", ]
td_genos_subset <- td_dat[td_dat$missing_type == "genotypes" & td_dat$missing_data < 0.80, ]

#wt_sites %>% group_by(method) %>%
#  summarize(cor=cor(wt_scaled, missing_data))

for (data_set in c(wt_sites, wt_genos)) {
print(head(data_set$method))
}

head(wt_genos$method)

for (data_set in c(wt_sites, wt_genos)) {
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "pixy.sites", ], 
                    data_set[data_set$method == "popgenome", ]))
  
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "pixy", ], 
                    data_set[data_set$method == "popgenome", ]))
  
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "scikitallel", ], 
                    data_set[data_set$method == "popgenome", ]))
  
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "pixy.sites", ], 
                    data_set[data_set$method == "pixy", ]))
  
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "pixy.sites", ], 
                    data_set[data_set$method == "scikitallel", ]))
  
  cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
        data = list(data_set[data_set$method == "pixy", ], 
                    data_set[data_set$method == "scikitallel", ]))
}

cocor(~ missing_data + avg_watterson_theta | missing_data + avg_watterson_theta,
      data = list(wt_genos[wt_genos$method == "pixy", ], 
                  wt_genos[wt_genos$method == "popgenome", ]))

cocor(~ missing_data + tajima_d | missing_data + tajima_d,
     data = list(td_genos_subset[td_genos_subset$method == "pixy", ], 
                 td_genos_subset[td_genos_subset$method == "popgenome", ]))

wt_sites_model <- lm(missing_data ~ wt_scaled:method, data = wt_sites)
wt_sites_emmeans <- emmeans(wt_sites_model, ~ wt_scaled:method)
pairs(wt_sites_emmeans, by = "wt_scaled")
summary(wt_sites_model)

wt_genos_model <- lm(missing_data ~ wt_scaled:method, data = wt_genos)
wt_genos_emmeans <- emmeans(wt_genos_model, ~ wt_scaled:method)
pairs(wt_genos_emmeans, by = "wt_scaled")
summary(wt_genos_model)

td_genos_model <- lm(missing_data ~ tajima_d:method, data = td_genos_subset)
td_genos_emmeans <- emmeans(td_genos_model, ~ tajima_d:method)
pairs(td_genos_emmeans, by = "tajima_d")
summary(td_genos_model)

# VCFtools Tajima's D ----

vcftools_dat <- sim_dat[sim_dat$method == "vcftools", ]

vcftools_genos = subset(vcftools_dat, missing_type == "genotypes")
vcftools_sites = subset(vcftools_dat, missing_type == "sites")

vcftools_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, vcftools_genos))
vcftools_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, vcftools_sites))

vcftools_genos_tajimad_regression$coefficients[8]
vcftools_sites_tajimad_regression$coefficients[8]

vcftools_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
             missing_type = factor("genotypes", levels = c("genotypes","sites")), method = "vcftools")
vcftools_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
             missing_type = factor("sites", levels = c("genotypes","sites")), method = "vcftools")

# PopGenome Tajima's D ----

popgenome_dat <- sim_dat[sim_dat$method == "popgenome", ]

popgenome_genos = subset(popgenome_dat, missing_type == "genotypes")
popgenome_sites = subset(popgenome_dat, missing_type == "sites")

popgenome_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, popgenome_genos))
popgenome_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, popgenome_sites))

popgenome_genos_tajimad_regression$coefficients[8]
popgenome_sites_tajimad_regression$coefficients[8]

popgenome_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
              missing_type = factor("genotypes", levels = c("genotypes","sites")), method = "popgenome")
popgenome_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
              missing_type = factor("sites", levels = c("genotypes","sites")), method = "popgenome")

# Scikit-allel Tajima's D ----

scikit_dat <- sim_dat[sim_dat$method == "scikitallel", ]

scikit_genos = subset(scikit_dat, missing_type == "genotypes")
scikit_sites = subset(scikit_dat, missing_type == "sites")

scikit_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, scikit_genos))
scikit_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, scikit_sites))

scikit_genos_tajimad_regression$coefficients[8]
scikit_sites_tajimad_regression$coefficients[8]

scikit_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
                   missing_type = factor("genotypes",levels = c("genotypes","sites")), method = "scikitallel")
scikit_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
                   missing_type = factor("sites",levels = c("genotypes","sites")) , method = "scikitallel")

# Pixy Tajima's D ----

pixy_dat <- sim_dat[sim_dat$method == "pixy", ]

pixy_genos = subset(pixy_dat, missing_type == "genotypes")
pixy_sites = subset(pixy_dat, missing_type == "sites")

pixy_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, pixy_genos))
pixy_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, pixy_sites))

pixy_sites_tajimad_regression$coefficients[8]

pixy_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
                  missing_type = factor("genotypes",levels = c("genotypes","sites")), method = "pixy")
pixy_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -3,
                  missing_type = factor("sites",levels = c("genotypes","sites")), method = "pixy")

# PopGenome Watterson's Theta ----

popgenome_genos_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, popgenome_genos))
popgenome_sites_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, popgenome_sites))

popgenome_genos_theta_regression$coefficients[8]
popgenome_sites_theta_regression$coefficients[8]
.Machine$double.xmin

popgenome_genos_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 0.4,
                missing_type = factor("genotypes",levels = c("genotypes","sites")), method = "popgenome")
popgenome_sites_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 1.3,
                missing_type = factor("sites",levels = c("genotypes","sites")), method = "popgenome")

# Scikit-allel Watterson's Theta ----

scikit_genos_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, scikit_genos))
scikit_sites_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, scikit_sites))

scikit_genos_theta_regression$coefficients[8]
scikit_sites_theta_regression$coefficients[8]

scikit_genos_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 0.4,
                                     missing_type = factor("genotypes",levels = c("genotypes","sites")), method = "scikitallel")
scikit_sites_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 1.3,
                                     missing_type = factor("sites",levels = c("genotypes","sites")), method = "scikitallel")

# Pixy Watterson's Theta ----

pixy_genos_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy_genos))
pixy_sites_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy_sites))

pixy_genos_theta_regression$adj.r.squared
pixy_genos_theta_regression$coefficients[8]
pixy_sites_theta_regression$adj.r.squared
pixy_sites_theta_regression$coefficients[8]

pixy_genos_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 0.4,
              missing_type = factor("genotypes", levels = c("genotypes","sites")), method = "pixy")
pixy_sites_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 1.3,
              missing_type = factor("sites", levels = c("genotypes","sites")), method = "pixy")

# Pixy.sites Watterson's Theta ----

pixy.sites_dat <- sim_dat[sim_dat$method == "pixy.sites", ]

pixy.sites_genos = subset(pixy.sites_dat, missing_type == "genotypes")
pixy.sites_sites = subset(pixy.sites_dat, missing_type == "sites")

pixy.sites_genos_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy.sites_genos))
pixy.sites_sites_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy.sites_sites))

pixy.sites_genos_theta_regression$coefficients[8]


pixy.sites_genos_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 0.4,
                                   missing_type = factor("genotypes", levels = c("genotypes","sites")), method = "pixy.sites")
pixy.sites_sites_theta_ann <- data.frame(missing_data = 0.5, wt_scaled = 1.3,
                                   missing_type = factor("sites", levels = c("genotypes","sites")), method = "pixy.sites")

# Plot Watterson's Theta ----

wt <- wt_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = wt_scaled))+
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, color = "grey50")+
  #geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = 1, color = "black", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method)+
  xlab("Proportion of Data Missing")+
  ylab(expression("Scaled " * theta[w] * " Estimate")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  ylim(0, 1.6) + 
  geom_text(data = pixy_genos_theta_ann %>% filter(method == "pixy"), label = paste("R^2 ==", 
            round(pixy_genos_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = pixy_sites_theta_ann %>% filter(method == "pixy"), label = paste("R^2 ==", 
            round(pixy_sites_theta_regression$r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = pixy.sites_genos_theta_ann %>% filter(method == "pixy.sites"), label = paste("R^2 ==", 
             round(pixy.sites_genos_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = pixy.sites_sites_theta_ann %>% filter(method == "pixy.sites"), label = paste("R^2 ==", 
            round(pixy.sites_sites_theta_regression$r.squared, digits = 6)), parse = TRUE,
            size = 3) +
  geom_text(data = popgenome_genos_theta_ann %>% filter(method == "popgenome"), label = paste("R^2 ==", 
          round(popgenome_genos_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
          size = 3) +
  geom_text(data = popgenome_sites_theta_ann %>% filter(method == "popgenome"), label = paste("R^2 ==", 
            round(popgenome_sites_theta_regression$r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = scikit_genos_theta_ann %>% filter(method == "scikitallel"), label = paste("R^2 ==", 
          round(scikit_genos_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
          size = 3) +
  geom_text(data = scikit_sites_theta_ann %>% filter(method == "scikitallel"), label = paste("R^2 ==", 
            round(scikit_sites_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 3)
wt

# Plot Tajima's D ----

td <- td_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d)) +
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, color = "grey50") +
#  geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
  geom_smooth(color = "red", se = FALSE) +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method)+
  xlab("Proportion of Data Missing") +
  ylab("Tajima's D Estimate") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  ylim(-4, 4) + 
  geom_text(data = pixy_genos_tajimad_ann %>% filter(method == "pixy"), label = paste("R^2 ==", 
            round(pixy_genos_tajimad_regression$r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = pixy_sites_tajimad_ann %>% filter(method == "pixy"), label = paste("R^2 ==", 
            round(pixy_sites_tajimad_regression$r.squared, digits = 5)), parse = TRUE,
            size = 3) +
  geom_text(data = popgenome_genos_tajimad_ann %>% filter(method == "popgenome"), label = paste("R^2 ==", 
            round(popgenome_genos_tajimad_regression$r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = popgenome_sites_tajimad_ann %>% filter(method == "popgenome"), label = paste("R^2 ==", 
            round(popgenome_sites_tajimad_regression$r.squared, digits = 5)), parse = TRUE,
            size = 3) +
  geom_text(data = scikit_genos_tajimad_ann %>% filter(method == "scikitallel"), label = paste("R^2 ==", 
            round(scikit_genos_tajimad_regression$r.squared, digits = 3)), parse = TRUE,
            size = 3) +
  geom_text(data = scikit_sites_tajimad_ann %>% filter(method == "scikitallel"), label = paste("R^2 ==", 
            round(scikit_sites_tajimad_regression$r.squared, digits = 5)), parse = TRUE,
            size = 3) +
  geom_text(data = vcftools_genos_tajimad_ann %>% filter(method == "vcftools"), label = paste("R^2 ==", 
          round(scikit_genos_tajimad_regression$r.squared, digits = 3)), parse = TRUE,
          size = 3) +
  geom_text(data = vcftools_sites_tajimad_ann %>% filter(method == "vcftools"), label = paste("R^2 ==", 
            round(scikit_sites_tajimad_regression$r.squared, digits = 5)), parse = TRUE,
            size = 3)
td

# same for dxy
#dxy <- sim_dat %>%
#  filter(missing_data < 1) %>%
#  ggplot(aes(x = missing_data, y = dxy_scaled))+
#  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, color = "grey50", raster.height = 1, raster.width = 1)+
  #geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
#  geom_smooth(color = "red", se = FALSE)+
#  geom_hline(yintercept = 1, color = "black", size = 0.5, linetype = 2) +
#  facet_grid(missing_type~method)+
#  xlab("Proportion of Data Missing")+
#  ylab("Scaled Dxy Estimate") +
#  theme_bw()+
#  theme(panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(), 
#        strip.background = element_blank(), 
#        axis.line = element_line(colour = "black"),
#        axis.text.x = element_text(angle = 45, hjust = 1))+
#  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
#  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) 

# Final Plot ----

compound <- fig3 <-pi + #/ dxy + 
   plot_annotation(tag_levels = "A")
compound
ggsave("figures/Figure3_raw.pdf", plot = compound, device = "pdf", scale = 1 ,width = 6, height = 6)
