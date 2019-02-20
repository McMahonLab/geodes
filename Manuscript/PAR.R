
enviro_data <- read.csv(paste(path2, "Desktop/geodes/environmental_data/for_humans/compiled_field_data_for_R.csv", sep = ""), header = T)


enviro_data$Time <- factor(enviro_data$Time, levels = c("9", "13", "17", "21", "1", "5"))
p <- ggplot(enviro_data[which(enviro_data$Depth == 0),], aes(x = Time, y = Light)) + geom_boxplot() + geom_point(aes(color = Lake))

enviro_data$condition <- "Day"
enviro_data$condition[which(enviro_data$Time == "21" | enviro_data$Time == "5" | enviro_data$Time == "1")] <- "Night"
p <- ggplot(enviro_data[which(enviro_data$Depth == 0),], aes(x = Time, y = log10(Light), fill = condition)) + geom_boxplot() + facet_wrap(~Lake, nrow = 1) + scale_fill_manual(values = c("yellow", "deepskyblue")) + labs(y = "Log of PAR (umol/s/m)", x = "Time") + background_grid(major = "xy")

save_plot("/Users/Alex/Desktop/geodes/Manuscript/supplemental_materials/PAR.pdf", p, base_height = 3, base_aspect_ratio = 2)
