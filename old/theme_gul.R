require(ggpplot2)
theme_gul <- theme_set( theme_gray() )

theme_gul <- theme_update( 
                # axis
                axis.text =  element_text( color = "gray25", size = rel(0.75)),
                axis.ticks = element_line( color = "gray70"),
                # background
                panel.background = element_rect(fill="gray94", color="gray25"),
                panel.grid.major = element_line(color="gray85", size = rel(0.6)),#,
                panel.grid.minor = element_line(color="gray88", size = rel(0.4)),
                # strips
                strip.background = element_rect(fill="gray55", color="gray25",size=0.2),
                strip.text       = element_text(color="white", size = rel(0.85) )
               )

