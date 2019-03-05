

require(stargazer)
#rownames(dd) <- NULL
# dd.tex <- stargazer(tex.dd, summary = F, digits = 2, perl = T, 
#                     title = "Conditional independence of determinants of protein production rate fold-change P",
#                     out = "~/Dropbox/PhD/thesis/TeX/items/tables/chapter5/table_results_conditional_independence_draft.tex"  )

gene.master.table       <- read.table("~/Documents/MRC16/results/master tables/gene.master.table.txt", header=1) # especially s-tAI scores and translation dynamics data from stochastic simulation


table.data <- with( subset(gene.master.table, valid == T  & time == 0 & total.translation_time < 1500 & n.events > 1, 
             select=c(name, n.codons, mRNA_abundance, av.initiation_time, av.elongation_time, total.translation_time)), 
             data.frame(check.names = F,
                        name = name,
                        'Length of CDS' = n.codons,
                        'Absolute mRNA abundance'= mRNA_abundance,
                        'Initiation frequency (min-1)' = 1/av.initiation_time*60,
                        'Elongation speed (s-1)' = n.codons / av.elongation_time,
                        'Total translation time of 1 protein (min)' = total.translation_time/60
                        ) 
      )

head(table.data)

table.tex.data <- do.call(rbind,lapply(table.data[,-1], function(x){ data.frame(mean = round(mean(x,na.rm=T),1), 
                                                               median = round(median(x,na.rm=T),1), 
                                                               sd = round(sd(x,na.rm=T),1), 
                                                               min = paste( round(min(x,na.rm=T),1), " (", as.character(table.data$name[which.min(x)])   ,")", sep  ="" ),
                                                               max = paste( round(max(x,na.rm=T),1), " (", as.character(table.data$name[which.max(x)])  ,")", sep = "" )
                                                               ) } ) )

stargazer(table.tex.data, summary = F, digits = 1, perl = T, title = "Model parameters in normal conditions",
          out = "~/Dropbox/PhD/thesis/TeX/items/tables/chapter5/table_results_parameters_normal_0.draft.tex")


subset( gene.master.table, total.translation_time < 5)

table( gene.master.table$total.translation_time > 1500 ) 
