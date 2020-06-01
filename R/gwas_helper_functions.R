

plot_qq <- function(p) {

    Observed <- -log10(sort(p, decreasing = FALSE, na.last = TRUE))
    Expected <- -log10(ppoints(length(p)))

    qplot(x = Expected,
          y = Observed) +
        geom_abline(intercept = 0, slope = 1) +
        theme_bw() +
        theme(panel.grid = element_blank())
}

plot_qq_raster <- function(p) {
    
    Observed <- -log10(sort(p, decreasing = FALSE, na.last = TRUE))
    Expected <- -log10(ppoints(length(p)))
    
    ggplot() +
        geom_point_rast(aes(x = Expected,
                            y = Observed)) +
        geom_abline(intercept = 0, slope = 1) +
        theme_bw() +
        theme(panel.grid = element_blank())
}


plot_manhattan <- function(data) {
    chr_numbers <- as.numeric(factor(data$chr,
                                     levels = unique(data$chr)))
    data$chr_indicator <- factor(chr_numbers %% 2)

    qplot(x = global_pos,
          y = -log10(p_wald),
          colour = chr_indicator,
          size = I(0.5),
          data = data) +
        scale_x_continuous(breaks = chr_breaks$global_pos_break,
                           labels = chr_breaks$chr_masked)
}

plot_manhattan_raster <- function(data) {
    chr_numbers <- as.numeric(factor(data$chr,
                                     levels = unique(data$chr)))
    data$chr_indicator <- factor(chr_numbers %% 2)
    
    ggplot() +
        geom_point_rast(aes(x = global_pos,
                            y = -log10(p_wald),
                            colour = chr_indicator),
                        data = data) +
        scale_x_continuous(breaks = chr_breaks$global_pos_break,
                           labels = chr_breaks$chr_masked)
}


flatten_coordinates <- function(chr,
                                pos,
                                chr_lengths) {
    pos_flat <- pos
    offset <- 0
 
    for (chr_ix in 1:nrow(chr_lengths)) {
        on_chr <- chr == chr_lengths$chr[chr_ix]
        pos_flat[on_chr] <- pos[on_chr] + offset
        offset <- offset + chr_lengths$length[chr_ix]
    }
 
    pos_flat
}
