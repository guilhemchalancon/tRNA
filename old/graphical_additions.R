data(XPR)


g <- ggplot(XPR$mRNA.abundance, aes(x= sqrt(mRNA.ab_YPD.rpkm)) ) + 
  geom_histogram( binwidth = 0.5) + theme_minimal()
ggsave(g , filename = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part2/mRNAab_distribution.pdf", useDingbats = F, width = 3, height = 2)

require(scales)

g <- ggplot( data=anticodon.master.table, aes(x=factor(1), y=total.tRNAab) ) + 
  geom_boxplot(fill="#99CC33", width=0.8, notch=T) + 
  geom_segment( aes(y = median(total.tRNAab), yend = median(total.tRNAab), xend=-Inf, x=1) ) +
  geom_text( data = subset(anticodon.master.table, total.tRNAab>4e5), 
             aes(label= paste( paste(anticodon, aa, sep="-"), paste( experiment, time,sep=" "), sep="\n")), size = 2.5 ) +
  geom_text(  data = data.frame(total.tRNAab = median(anticodon.master.table$total.tRNAab)), aes(label = round(total.tRNAab)),x=+Inf, vjust = 1.2, size = 3 ) +
  labs(y="",x="") +
  scale_y_continuous( labels = prettyNum ) +
#   scale_y_continuous( breaks = trans_breaks("log10", function(x) 10^x),
#                  labels = trans_format("log10", math_format(10^.x))
#   )  +
   coord_flip() + theme_minimal() + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())

ggsave(g , filename = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part2/anticodon.supply_distribution.pdf",
       useDingbats = F, width = 8, height = 1.5, dpi=300)



g <- ggplot( data=anticodon.master.table, aes(x=factor(1), y=demand.mRNAab) ) + 
  geom_boxplot(fill="#33CCFF", width=0.8, notch=T) + 
  geom_segment( aes(y = median(demand.mRNAab), yend = median(demand.mRNAab), xend=-Inf, x=1) ) +
  geom_point( data = subset(anticodon.master.table, demand.mRNAab>1.38e6), 
              size = 2 ) +
  geom_text( data = subset(anticodon.master.table, demand.mRNAab>1.38e6), 
             aes(label= paste( paste(anticodon, aa, sep="-"), paste( experiment, time,sep=" "), sep="\n")), size = 2.5 ) +
  geom_text(  data = data.frame(demand.mRNAab = median(anticodon.master.table$demand.mRNAab)), aes(label = round(demand.mRNAab)),x=+Inf, vjust = 1.2, size = 3 ) +
  labs(y="",x="") +
  scale_y_continuous( labels = prettyNum ) +
  coord_flip() + theme_minimal() + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())

ggsave(g , filename = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part2/anticodon.demand_distribution.pdf",
       useDingbats = F, width = 8, height = 1.5, dpi=300)
