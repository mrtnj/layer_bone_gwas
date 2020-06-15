
## Summarise results of conditional body weight GWAS

library(dplyr)
library(ggplot2)
library(readr)


main <- read_tsv("gwas/all_weight/output/all_weight.assoc.txt",
                 col_types = "ccnnccnnnnnnnnn")

conditional <- read_tsv("gwas/all_residual_weight/output/all_residual_weight.assoc.txt",
                        col_types = "ccnnccnnnnnnnnn")


combined <- rbind(transform(filter(main, chr == 4),
                            scan = "Default"),
                  transform(filter(conditional, chr == 4),
                            scan = "Conditional"))


formatting <- list(geom_hline(yintercept = -log10(5e-8),
                              colour = "red",
                              linetype = 2),
                   geom_hline(yintercept = -log10(1e-4),
                              colour = "blue",
                              linetype = 2),
                   theme_bw(),
                   theme(panel.grid = element_blank(),
                         legend.position = "right"),
                   scale_colour_manual(values = c("grey", "black"),
                                       name = ""),
                   xlab("Position (Mbp)"),
                   ylab(""))


plot_comparison <- qplot(x = ps/1e6, y = -log10(p_wald),
                         colour = scan,
                         data = combined) +
    formatting +
    ggtitle("Conditional GWAS of body weight on chromosome 4")


pdf("figures/plot_chr4_conditional.pdf")
print(plot_comparison)
dev.off()
