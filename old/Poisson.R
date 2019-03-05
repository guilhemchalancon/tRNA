g1 <- ggplot(data.frame(x=c(0.5, 100)), aes(x)) + 
  geom_abline(slope=1, intercept=0, lty = 3) +
  stat_function(fun = function(x){sqrt(x)}) +
  coord_fixed(ratio = 1) + ylim(0,100) +
   labs(x=expression("mean:  "*lambda), y = expression("SD:  "*sqrt(lambda)), 
       title = "How SD scales with mean\n in a Poisson process" )

g2 <- ggplot(data.frame(x=c(0.5, 100)), aes(x)) + 
  stat_function(fun= function(x){sqrt(x)/x}) + 
  coord_fixed(ratio = 100) +
  labs(x=expression("mean:  "*lambda), y = expression("CV:  "*sqrt(lambda)/lambda), 
       title = "How CV scales with mean\n in a Poisson process" )


require(gridExtra)
grid.arrange(g1, g2, nrow=1)