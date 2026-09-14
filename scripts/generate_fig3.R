# This script plots Figure 3 -  an example of how the modelling framework can be
# used to simulate spillover, onward transmission, and surveillance for an emerging
# virus.

# It depends on:
#### results/simulations/fig3.rds
#### results/simulations/fig3_linelist.rds

# It produces:
### figures/fig3.png

# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)
library(purrr)
library(here)

# load simulations for plotting
det_clinical <- readRDS(here("results/simulations/fig3.rds"))
linelist <- readRDS(here("results/simulations/fig3_linelist.rds"))

# 4. Take a look at results ----------------------------------------------------
# ------------------------------------------------------------------------------

##### A. infections (by type) and viral load shed to wastewater

# Approximate start and mid points for each month (non-leap year)
month_start <- c(0, 31, 59, 90, 120, 151,
                 181, 212, 243, 273, 304, 334, 364)
p1 <- ggplot(linelist) +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(hjust = -0.66),) +
  geom_ribbon(data = det$sampled_data, aes(x = day, ymax = shedding_value*1.5, ymin = 0), fill = "#D3D3D3") +
  geom_histogram(aes(
    x = time_infection,
    fill = factor(type, levels = c("Human-to-human", "Animal-to-human"))
  ),
  binwidth = 1,
  boundary = 0
  ) +
  scale_fill_manual(values = c("#826699", "#9CAF88"), name = "Source") +
  xlab(NULL) +
  ylab("Incident infections") +
  scale_x_continuous(
    breaks = month_start, # tick marks at start of each month
    labels = c(month.abb, " "),
    limits = c(0, 365),
    expand = c(0.01,0)
  ) +
  scale_y_continuous(
    # secondary y-axis - rescale to account for shedding rel to sc2
    sec.axis = sec_axis(~ .*(1/1)/1.5, name = "Effective n individuals\nshedding into WW")
  )

##### B. probability of detection in wastewater & stochastic detection events
p2 <- ggplot(det$sampled_data)+
  geom_vline(data = det$sampled_data %>% filter(sampled == "Yes"),
             aes(xintercept = day), linewidth = 0.5, colour = "grey90")+
  geom_ribbon(data = det$sampled_data,
              aes(x = day, ymax = prob_detect, ymin = 0), fill = "salmon4", alpha = 0.6)+
  geom_line(aes(y = prob_detect, x = day), colour = "salmon4", alpha = 0.3)+
  geom_vline(data = det$sampled_data %>% filter(sampled == "Yes" & detect_draw>0),
             aes(xintercept = day), linewidth = 0.5, colour = "coral3")+

  scale_x_continuous(
    breaks = month_start, # tick marks at start of each month
    labels = c(NULL), # labels for months
    limits = c(0, 365),
    expand = c(0.01,0)
  ) +
  ylim(0,0.35)+
  xlab(NULL)+
  ylab("Probability detected\nin wastewater")+
  theme_classic()+
  theme(text = element_text(size = 15))

##### C. clinical stochastic detection events

p3 <- ggplot(det_clinical)+
  geom_ribbon(aes(x = day, ymax = prob_detection_clinical, ymin = 0), fill = "#5B7C99", alpha = 0.5)+
  geom_line(aes(y = prob_detection_clinical, x = day), colour = "#5B7C99", alpha = 0.7)+
  geom_vline(data = linelist%>%filter(!is.na(time_diagnosis)),
             aes(xintercept = time_diagnosis), linewidth = 0.5, colour = "coral3")+
  scale_x_continuous(
    breaks = month_start,        # tick marks at start of each month
    labels = c(month.abb, " "),           # labels for months
    limits = c(0, 365),
    expand = c(0.01,0)
  ) +
  scale_x_continuous(
    breaks = month_start,                 # tick marks at start of each month
    labels = c(month.abb, " "),
    limits = c(0, 365),
    expand = c(0.01,0)
  )+
  ylim(0,0.35)+
  ylab("Probability detected\nclinically ")+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = -0.66),
        text = element_text(size = 15))+
  labs(x = NULL)

##### D. plot transmission trees

# Create edges and vertices
linelistplot <- linelist %>%
  filter(time_infection<366)%>%
  mutate(
    id_outbreak = paste0(id, "-", spillover_ID),
    id_infector = paste0(infector, "-", spillover_ID),
  ) %>%
  select(id_outbreak, id_infector, symptomatic, diagnosis, type,
         time_infection, infection_generation) %>%
  rename(to = id_outbreak, from = id_infector)
edges <- linelistplot %>%
  select(from, to)

# Make sure vertices$"name" exists and covers ALL from/to
vertices <- linelistplot %>%
  transmute(name = to,
            symptomatic, diagnosis, type, time_infection,
            infection_generation)
# Add missing infector names (animal reservoir)
animal_names <- unique(linelistplot$from[grepl(linelistplot$from,
                                               pattern = "animal")])
animal_res <- linelistplot %>% filter(from %in% animal_names)%>%
  select(from, time_infection)%>%
  mutate(name = from,
         symptomatic = NA,
         diagnosis = NA,
         type = "index",
         infection_generation = 0)
vertices <- bind_rows(vertices, animal_res) %>%
  distinct(name, .keep_all = TRUE)  # remove duplicates if any

# Build edges with attributes of the 'to' node
edges_attr <- edges %>%
  left_join(vertices %>% select(name, type), by = c("to" = "name"))

# Build graph with edge attributes included
g <- graph_from_data_frame(d = edges_attr, vertices = vertices, directed = TRUE)

# Layout to plot by time and generation and add jitter to humans but not animals
vertices <- vertices %>%
  mutate(infection_generation = if_else(grepl(name, pattern = "animal"),
                                        infection_generation,
                                        infection_generation + rnorm(n(), 0, 0.1)))
# Flip the layout so gneration 0 is at the top of the y axis
flipped_layout <- as.matrix(vertices[, c("time_infection", "infection_generation")])
flipped_layout[, 2] <- -flipped_layout[, 2]  # Multiply the y-axis by -1

# Plot
p0 <- ggraph(g, layout = flipped_layout) +
  geom_edge_diagonal(arrow = arrow(length = unit(2, 'mm')),
                     end_cap = circle(1.2, 'mm'),
                     width = 0.5, colour = "grey40") +
  geom_node_point(aes(color = type), size = 4) +
  theme_classic() +
  ylab("Infection generation")+
  xlab(NULL)+
  scale_x_continuous(
    breaks = month_start,        # tick marks at start of each month
    labels = NULL,           # labels for months
    limits = c(0, 365),
    expand = c(0.01,0)
  ) +
  scale_y_continuous(labels = c(5,4,3,2,1,0), breaks = c(-5,-4,-3,-2,-1,0))+
  scale_color_manual(values = c("#9CAF88","#826699","grey40"),
                     name = "Transmission source")+
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    edge_color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(1.2, "cm"),
        axis.text.x = element_text(hjust = -0.75),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))

legend <- get_legend(p1)

p0123<- plot_grid(p0+theme(legend.position = "none"), p1+theme(legend.position = "none"),p2,p3, nrow = 4, align = "v")

ptot <- plot_grid(legend,p0123, nrow = 2, rel_heights = c(1,10))

ggsave(dpi = 300, ptot, file = here("figures/fig3.png", width = 10, height = 12))
