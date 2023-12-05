#!/usr/bin/env Rscript

library(tidyverse)
library(Hmisc)

BQ_bin <- 5
inputfile <- commandArgs(TRUE)[1] 
#inputfile <- 'data/validation_N500_LR_with_MQ_ME5_MM3_MQ40_BQ0_NMD1_filterBQ.tsv'
output_prefix <- commandArgs(TRUE)[2]
d <- read_tsv(inputfile)

p <- d |>
  filter(NA37<10) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(BQ, muttype) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen > 0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=muttype, fill=muttype, ymin=lower, ymax=upper)) + 
  geom_ribbon(alpha=0.5) + geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 
ggsave(paste(output_prefix, "no_mismatch_BQ_muttype.png", sep="_"), p, width=10, height=6)

p <- d |>
  filter(NA37<10) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(BQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen > 0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.5) +
  geom_point() + geom_line() +
  scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 

ggsave(paste(output_prefix, "no_mismatch_BQ.png", sep="_"), p, width=10, height=6)


p <- d |>
  filter(NA37<10) |>
  #filter(max_oldBQ==37) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(BQ, muttype, oldBQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen>0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=muttype, fill=muttype, ymin=lower, ymax=upper)) + 
  #geom_ribbon(alpha=0.5) +
  geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() + facet_wrap(~oldBQ)

ggsave(paste(output_prefix, "no_mismatch_BQ_muttype_oldBQ.png", sep="_"), p, width=10, height=6)


p <- d |>
  filter(NA37<10) |>
  #filter(max_oldBQ==37) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(BQ, oldBQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen>0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=oldBQ, fill=oldBQ, ymin=lower, ymax=upper)) + 
  #geom_ribbon(alpha=0.5) +
  geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 

ggsave(paste(output_prefix, "no_mismatch_BQ_oldBQ.png", sep="_"), p, width=10, height=6)

p <- d |>
  filter(NA37<10) |>
  #filter(max_oldBQ==37) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin,
         max_oldBQ = factor(max_oldBQ)) |>
  group_by(BQ, max_oldBQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen>0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=max_oldBQ, fill=max_oldBQ, ymin=lower, ymax=upper)) + 
  #geom_ribbon(alpha=0.5) +
  geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 

ggsave(paste(output_prefix, "no_mismatch_BQ_maxoldBQ.png", sep="_"), p, width=10, height=6)


p <- d |>
  filter(NA37<10) |>
  #filter(max_oldBQ==37) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(muttype, oldBQ, BQ, validation) |>
  summarise(count = sum(count)) |>
  group_by(muttype, oldBQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count),
            BQ = sum(BQ*count)/sum(count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen > 0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=oldBQ, shape=muttype, fill=muttype, ymin=lower, ymax=upper)) + 
  geom_errorbar(alpha=0.5, width=0.2) +    
  geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw()

ggsave(paste(output_prefix, "no_mismatch_muttype_oldBQ.png", sep="_"), p, width=10, height=6)


p <- d |>
  filter(NA37<10) |>
  filter(n_mismatch==0) |>
  #filter(max_oldBQ==37) |>
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(oldBQ) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count),
            BQ = sum(BQ*count)/sum(count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen >0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=oldBQ, fill=oldBQ, ymin=lower, ymax=upper)) +
  geom_errorbar(alpha=0.5, width=0.2) + geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 

ggsave(paste(output_prefix, "no_mismatch_oldBQ.png", sep="_"), p, width=7, height=6)


p <- d |>
  filter(NA37<10) |>
  #filter(max_oldBQ==37) |>
  filter(n_mismatch==0) |> 
  mutate(BQ = as.integer(newBQ/BQ_bin)*BQ_bin) |>
  group_by(muttype) |>
  summarise(n_seen = sum(as.integer(validation==1)*count),
            n_notseen = sum(as.integer(validation==0)*count),
            BQ = sum(BQ*count)/sum(count)) |>
  filter(n_seen+n_notseen > 50,
         n_seen>0) |>
  mutate(rate = n_seen/(n_seen+n_notseen),
         lower = binconf(n_seen, n_seen+n_notseen)[,2],
         upper = binconf(n_seen, n_seen+n_notseen)[,3]) |>
  ggplot(aes(x=BQ, y=rate, color=muttype, fill=muttype, ymin=lower, ymax=upper)) + 
  geom_errorbar(alpha=0.5, width=0.2) + geom_point() + geom_line() + scale_y_log10("Probability allele seen at same site in HiFi data") +
  theme_bw() 

ggsave(paste(output_prefix, "no_mismatch_muttype.png", sep="_"), p, width=7, height=6)
