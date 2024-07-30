# aggregate all raw pi/dxy calculations from various software packages
# KMS 2020-04-02
# KLK edited 2020-05-15
setwd('pixy_analysis/')
library("tidyverse")
library("officer")
library("aod")
library("rvg")
library("ggrastr")
library("patchwork")
library("reshape2")
library("gridExtra")
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

vcftools_dat <- sim_dat[sim_dat$method == "vcftools", ]

vcftools_genos = subset(vcftools_dat, missing_type == "genotypes")
vcftools_sites = subset(vcftools_dat, missing_type == "sites")

vcftools_genos_regression = summary(lm(tajima_d ~ missing_data, vcftools_genos))
vcftools_sites_regression = summary(lm(tajima_d ~ missing_data, vcftools_sites))

vcftools_genos_ann <- data.frame(missing_data = 0.5, tajima_d = 4,
             missing_type = factor("genotypes", levels = c("genotypes","sites")))
vcftools_sites_ann <- data.frame(missing_data = 0.5, tajima_d = 4,
             missing_type = factor("sites", levels = c("genotypes","sites")))

vcftools_dat %>%
#  sample_frac(0.05) %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d))+
  geom_point_rast(size = 0.5, alpha = 1, shape = 16, color = "black")+
#  geom_point(size = 0.5, alpha = 1, shape = 16, color = "black")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = exp_pi, color = "blue", size = 0.75, linetype = 2) +
  facet_grid(method~missing_type) +
  xlab("Proportion of Data Missing") +
  ylab("VCFtools Tajima's D Estimate") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = vcftools_genos_ann, label = paste("R^2 ==", 
            round(vcftools_genos_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = vcftools_sites_ann, label = paste("R^2 ==", 
            round(vcftools_sites_regression$r.squared, digits = 5)), parse = TRUE,
            size = 5)

popgenome_dat <- sim_dat[sim_dat$method == "popgenome", ]

popgenome_genos = subset(popgenome_dat, missing_type == "genotypes")
popgenome_sites = subset(popgenome_dat, missing_type == "sites")

popgenome_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, popgenome_genos))
popgenome_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, popgenome_sites))

popgenome_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = 4,
              missing_type = factor("genotypes", levels = c("genotypes","sites")))
popgenome_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = 4,
              missing_type = factor("sites", levels = c("genotypes","sites")))

popgenome_dat %>%
  #  sample_frac(0.05) %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d))+
  geom_point_rast(size = 0.5, alpha = 1, shape = 16, color = "black")+
  #  geom_point(size = 0.5, alpha = 1, shape = 16, color = "black")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = exp_pi, color = "blue", size = 0.75, linetype = 2) +
  facet_grid(method~missing_type)+
  xlab("Proportion of Data Missing")+
  ylab("PopGenome Tajima's D Estimate") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = popgenome_genos_tajimad_ann, label = paste("R^2 ==", 
            round(popgenome_genos_tajimad_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = popgenome_sites_tajimad_ann, label = paste("R^2 ==", 
            round(popgenome_sites_tajimad_regression$r.squared, digits = 5)), parse = TRUE,
            size = 5)

scikit_dat <- sim_dat[sim_dat$method == "scikitallel", ]

scikit_genos = subset(scikit_dat, missing_type == "genotypes")
scikit_sites = subset(scikit_dat, missing_type == "sites")

scikit_genos_regression = summary(lm(tajima_d ~ missing_data, scikit_genos))
scikit_sites_regression = summary(lm(tajima_d ~ missing_data, scikit_sites))

scikit_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = 7,
                   missing_type = factor("genotypes",levels = c("genotypes","sites")))
scikit_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = 7,
                   missing_type = factor("sites",levels = c("genotypes","sites")))

scikit_dat %>%
  #  sample_frac(0.05) %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d)) +
  geom_point_rast(size = 0.5, alpha = 1, shape = 16, color = "black")+
  #  geom_point(size = 0.5, alpha = 1, shape = 16, color = "black")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = exp_pi, color = "blue", size = 0.75, linetype = 2) +
  facet_grid(method~missing_type)+
  xlab("Proportion of Data Missing")+
  ylab("scikit-allel Tajima's D Estimate") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = scikit_genos_tajimad_ann, label = paste("R^2 ==", 
            round(scikit_genos_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = scikit_sites_tajimad_ann, label = paste("R^2 ==", 
            round(scikit_sites_regression$r.squared, digits = 5)), parse = TRUE,
            size = 5)

pixy_dat <- sim_dat[sim_dat$method == "pixy", ]

pixy_genos = subset(pixy_dat, missing_type == "genotypes") #& missing_data < 0.9)
pixy_sites = subset(pixy_dat, missing_type == "sites")

#pixy_genos_subset <- pixy_genos[pixy_genos$missing_data < 0.9, ]

pixy_genos_tajimad_regression = summary(lm(tajima_d ~ missing_data, pixy_genos))
pixy_sites_tajimad_regression = summary(lm(tajima_d ~ missing_data, pixy_sites))

summary(lm(tajima_d ~ missing_data, pixy_genos))
summary(lm(tajima_d ~ missing_data, pixy_sites))

pixy_genos_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -50,
                  missing_type = factor("genotypes",levels = c("genotypes","sites")))
pixy_sites_tajimad_ann <- data.frame(missing_data = 0.5, tajima_d = -50,
                  missing_type = factor("sites",levels = c("genotypes","sites")))

pixy_dat %>%
  #  sample_frac(0.05) %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d))+
  geom_point_rast(size = 0.5, alpha = 1, shape = 16, color = "black")+
  #  geom_point(size = 0.5, alpha = 1, shape = 16, color = "black")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = 0, color = "blue", size = 0.75, linetype = 2) +
  facet_wrap(~missing_type)+
  xlab("Proportion of Data Missing")+
  ylab("Our Tajima's D Estimate") +
  theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = pixy_genos_tajimad_ann, label = paste("R^2 ==", 
            round(pixy_genos_tajimad_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = pixy_sites_tajimad_ann, label = paste("R^2 ==", 
            round(pixy_sites_tajimad_regression$r.squared, digits = 6)), parse = TRUE,
            size = 5)

popgenome_genos_wt_regression = summary(lm(avg_watterson_theta ~ missing_data, popgenome_genos))
popgenome_sites_wt_regression = summary(lm(avg_watterson_theta ~ missing_data, popgenome_sites))

summary(lm(avg_watterson_theta ~ missing_data, popgenome_genos))
summary(lm(avg_watterson_theta ~ missing_data, popgenome_sites))

popgenome_genos_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.2,
                missing_type = factor("genotypes",levels = c("genotypes","sites")))
popgenome_sites_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.2,
                missing_type = factor("sites",levels = c("genotypes","sites")))

popgenome_theta <- popgenome_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = wt_scaled)) +
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, 
                  color = "grey50") +
  #geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = 1, color = "blue", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method) +
  xlab("Proportion of Data Missing") +
  ylab("popgenome Watterson's Theta Estimate") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = popgenome_genos_theta_ann, label = paste("R^2 ==", 
              round(popgenome_genos_wt_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = popgenome_sites_theta_ann, label = paste("R^2 ==", 
              round(popgenome_sites_wt_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5)

popgenome_theta

scikit_genos_regression = summary(lm(avg_watterson_theta ~ missing_data, scikit_genos))
scikit_sites_regression = summary(lm(avg_watterson_theta ~ missing_data, scikit_sites))

scikit_genos_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.2,
                                     missing_type = factor("genotypes",levels = c("genotypes","sites")))
scikit_sites_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.2,
                                     missing_type = factor("sites",levels = c("genotypes","sites")))

# scaled pi for all methods
# figure 2B
scikit_theta <- scikit_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = wt_scaled))+
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, color = "grey50")+
  #geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = 1, color = "blue", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method)+
  xlab("Proportion of Data Missing")+
  ylab("scikit-allel Watterson's Theta Estimate") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = scikit_genos_theta_ann, label = paste("R^2 ==", 
                                                         round(scikit_genos_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = scikit_sites_theta_ann, label = paste("R^2 ==", 
                                                         round(scikit_sites_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5)
scikit_theta

pixy_genos_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy_genos))
pixy_sites_theta_regression = summary(lm(avg_watterson_theta ~ missing_data, pixy_sites))

pixy_genos_theta_regression$adj.r.squared
pixy_genos_theta_regression$coefficients

pixy_genos_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.7,
              missing_type = factor("genotypes",levels = c("genotypes","sites")))
pixy_sites_theta_ann <- data.frame(missing_data = 0.2, wt_scaled = 0.7,
              missing_type = factor("sites",levels = c("genotypes","sites")))

pixy_theta <- pixy_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = wt_scaled)) +
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, 
                  color = "grey50") +
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = 1, color = "blue", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method)+
  xlab("Proportion of Data Missing")+
  ylab("Our Watterson's Theta Estimate") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  geom_text(data = pixy_genos_theta_ann, label = paste("R^2 ==", 
            round(pixy_genos_theta_regression$adj.r.squared, digits = 3)), parse = TRUE,
            size = 5) +
  geom_text(data = pixy_sites_theta_ann, label = paste("R^2 ==", 
            round(pixy_sites_theta_regression$r.squared, digits = 6)), parse = TRUE,
            size = 5)
pixy_theta

wt_dat <- sim_dat[sim_dat$method != "vcftools", ]

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
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))
wt

td <- sim_dat %>%
  filter(missing_data < 1) %>%
  ggplot(aes(x = missing_data, y = tajima_d))+
  geom_point_rast(size = 0.25, alpha = 0.4, shape = 16, color = "grey50")+
  #geom_point(size = 0.5, alpha = 0.4, shape = 16, color = "grey50")+
  geom_smooth(color = "red", se = FALSE)+
  geom_hline(yintercept = , color = "black", size = 0.5, linetype = 2) +
  facet_grid(missing_type~method)+
  xlab("Proportion of Data Missing")+
  ylab("Tajima's D Estimate") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) #+
#  ylim(-3, 3)
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

compound <- fig3 <-pi + #/ dxy + 
   plot_annotation(tag_levels = "A")
compound
ggsave("figures/Figure3_raw.pdf", plot = compound, device = "pdf", scale = 1 ,width = 6, height = 6)
