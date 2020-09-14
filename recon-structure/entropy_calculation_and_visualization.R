######
## Calculation, processing and visualization of the entropy 
## induced by the enzymatic reactions
######
library(stringr)
library(ggplot2)
library(tibble)

recon <- read.table('recon22.sfba', sep="\t")
enzymatic.reactions <- recon[recon$rules != "",]

no.reactions <- sort(table(enzymatic.reactions$rules), decreasing = T)
top_rules <- names(no.reactions)
top_entropies <- sapply(top_rules, function(rule){
  rows <- as.numeric(rownames(enzymatic.reactions[enzymatic.reactions$rules == rule,]))
  met.names <- select.unique.metab(enzymatic.reactions[enzymatic.reactions$rules == rule,])
  cols <- which(colnames(stoch.matrix) %in% met.names)
  entr.graph(stoch.matrix, list(cols = cols, rows = rows))
})

hist((top_entropies), breaks = 200, xlab="Entropy", main="Histogram of entropy in enzymatic reactions")
plot(as.numeric(no.reactions), top_entropies, pch = "*", xlab="Number of reactions", ylab="Entropy")

#########
## Generation of random entropies
#########
to_sample <- 1:nrow(enzymatic.reactions)

random.entropy.4 <- sapply(rep(1:300, 20), function(no.rec){
  print(no.rec)
  rows <- sample(to_sample, size = no.rec)
  met.names <- select.unique.metab(enzymatic.reactions[rows,])
  cols <- which(colnames(stoch.matrix) %in% met.names)
  entr.graph(stoch.matrix, list(cols = cols, rows = rows))
})

hist((random.entropy), breaks = 200, xlab="Entropy", main="Histogram of entropy in random reactions")

#########
## Summarizing the RECON entropy data
#########
entropy_data <- rbind(
  data.frame(no_reactions = rep(1:300, 20), 
             entropy = random.entropy.4, 
             type="Sampled"),
  data.frame(no_reactions = as.numeric(no.reactions), entropy = top_entropies, type="RECON"),
  data.frame(no_reactions = as.numeric(no.reactions[important_random_data]), entropy = top_entropies[important_random_data], type="Random Data"),
  data.frame(no_reactions = as.numeric(no.reactions[important_before_reduction]), entropy = top_entropies[important_before_reduction], type="Before adjustment"),
  data.frame(no_reactions = as.numeric(no.reactions[important_after_reduction]), entropy = top_entropies[important_after_reduction], type="After adjustment"),
  data.frame(no_reactions = as.numeric(no.reactions[poor_prognosis]), entropy = top_entropies[poor_prognosis], type="Poor prognosis"))

entropy_data <- entropy_data %>% 
  group_by(no_reactions, type) %>% 
  mutate(lower = quantile(entropy, 0.05), 
         upper = quantile(entropy, 0.95), 
         mean_ent = mean(entropy))

#########
## Visualization
#########
cols <- c("RECON" = "darkblue", "Sampled" = "lightgrey", "Random Data" = "darkgreen", 
          "Before adjustment" = "purple",  "After adjustment" = "orange", "Poor prognosis" = "red3")
my.shapes <-  c("RECON" = 18, "Sampled" = 18, "Random Data" = 18, 
                "Before adjustment" = 17 ,  "After adjustment" = 16 , "Poor prognosis" = 15)

entropy_subset <- entropy_data[entropy_data$mean_ent < 30 & entropy_data$no_reactions < 20,]
entropy_conf_int <- entropy_subset[entropy_subset$type %in% c("Sampled", "RECON") & 
                                     entropy_subset$lower != entropy_subset$upper,]
entropy_conf_outliers <- entropy_subset[entropy_subset$type %in% c("Sampled", "RECON") & 
                                          entropy_subset$lower == entropy_subset$upper,]
entropy_RCC <- entropy_subset[entropy_subset$type %in% c("Before adjustment", "Poor prognosis", "After adjustment") ,]
override.shape <- c(16, 17, 18, NA, NA)
override.fill <- c("orange", "purple", "red3", "white", "white")

g <-  ggplot(entropy_subset, 
             aes(x = round(no_reactions, 0), y = entropy, color = type), shape=16) + 
  geom_line(data = entropy_conf_int, 
            aes(x = no_reactions, y = mean_ent, group = type, colour = type)) +
  geom_ribbon(data = entropy_conf_int, 
              aes(x = no_reactions, ymin = lower, ymax = upper, fill = type, group = type), alpha = 0.3) +
  geom_point(data = entropy_conf_outliers, aes(x = no_reactions, y = entropy, fill = type), pch = 18) + 
  geom_point(data = entropy_RCC, aes(shape = type, color = type), size = 2.25, alpha = 0.5) + 
  scale_x_continuous(breaks =  seq(1, 25, by=2))  + # + c(0, 50, 100, 150, 297))
  ylab("Entropy") + xlab("Hub size") + 
  guides(shape = guide_legend( override.aes = list(colour = c("orange","red3"), 
                                                   shape = c(16,15),
                                                   alpha = 1)), 
         fill = guide_legend( override.aes = list(shape = c(18, NA), 
                                                  colour = c(alpha("darkblue", 0.8), alpha("lightgrey", 0.8)), 
                                                  fill = c(alpha("darkblue", 0.8), alpha("lightgrey", 0.3))))) +
  scale_color_manual(values = cols, guide = FALSE) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = my.shapes) +
  labs(shape="RCC discriminating", colour="Petal width label", fill = "Reaction Hubs") +
  theme_minimal() + 
  theme(legend.position = "Bottom", legend.box="None", legend.margin = margin(-5, 0, 0, 0),
                           strip.text.x = element_blank()) + 
  facet_grid( . ~ factor(no_reactions < 250, levels=c(T, F)), scales = "free_x", space = "free_x" )

g