ward_rates %>% 
  left_join(newID_COVIDstat) %>% 
  group_by(newID) %>% 
  summarise(earliest = min(time)) %>% 
  mutate(weekday = weekdays(earliest)) %>% 
pull(weekday) %>% table

hist(ward_rates$contactIntensity)

ward_rates %>% 
  as_tibble %>% 
  ggplot(aes(x = contactIntensity)) + geom_histogram() + 
  facet_wrap(~ward_id)

ward_expRates <- ward_intervals %>%
  group_by(ward_id) %>%
  summarize(expRate = 1 / mean(contactIntensity, na.rm = T))

ward_expRates$expRate
typeID
unique(ward_intervals$ward_id)

distribution_contactIntensity <- ward_intervals %>% 
  left_join(typeID %>% select(ward_id, newID)) %>% 
  select(ward_id, newID, contactIntensity) %>% 
  filter(!is.na(contactIntensity), !is.na(newID))

ward_expRates <- distribution_contactIntensity %>% 
  group_by(ward_id, newID) %>% 
  summarize(expRate = 1 / mean(contactIntensity, na.rm = T))


line_tib <- expand.grid(list(newID = distribution_contactIntensity$newID %>% unique %>% sort, 
     contactIntensity = seq(0, distribution_contactIntensity$contactIntensity %>% max, by = 0.01))) %>% 
  as_tibble %>% 
  left_join(ward_expRates) %>% 
  mutate(density = dexp(x = contactIntensity, rate = expRate))


# Plot the histograms with exponential lines
distribution_contactIntensity %>% 
  ggplot(aes(x = contactIntensity)) + 
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue") + 
  geom_line(data = line_tib, aes(x = contactIntensity, y = density, colour = "Exponential"), linetype = "dashed", size = 0.5) + 
  facet_wrap(~newID, nrow = 2) +
  scale_colour_manual(values = c(Exponential = "black")) + 
  theme_bw() +
  labs(x = "Contact rate", y = "Density", colour = "") + 
  coord_cartesian(x = c(0, 1.0), y = c(0, 4.5)) + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, ""))

ggsave(filename = paste0("output/SuppFig1 - exponential rates.png")
       , width = 35, height = 12, units = "cm", dpi = 1000)



