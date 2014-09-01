
library(ggplot2, quietly=TRUE)
library(RColorBrewer)

######################################################################
# Theme
######################################################################
immigration_theme <- theme(
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    plot.title = element_text(face="bold", size=22,colour="black"),
    axis.text.x = element_text(colour="black",size=13,angle=45),
    axis.text.y = element_text(colour="black",size=13,angle=0),
    axis.title.x = element_text(colour="black",size=20,angle=0),
    axis.title.y = element_text(colour="black",size=20,angle=90),
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white")
    #legend.title = element_blank()
    #axis.ticks = element_blank()
    #legend.position = "none"
    )

ng1 = theme(panel.background = element_rect(fill = "white",colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    axis.line = element_line(size = 1.2, colour="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black",size=15),
    axis.title = element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    strip.text.y = element_text(color="black",face="bold",size=15),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white")
    #legend.title = element_blank()
    )

ng2 = theme(panel.background = element_rect(fill = "grey90",colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    axis.line = element_line(size = 1.2, colour="black"),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(color="black",size=15),
    axis.text.x = element_text(colour="black",angle=45,hjust=1,vjust=1),
    axis.title=element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.text.y = element_text(size=13,face="bold"),
    strip.text.x = element_text(size=13,face="bold")
    )

hclusttheme <- theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(color="black",size=20),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.text.y = element_text(size=13,face="bold"),
    strip.text.x = element_text(size=13,face="bold")
    )
    #axis.text.x = element_text(colour="black",angle=45,hjust=1,vjust=1),
    #legend.position = "top", 
    #legend.direction="horizontal", 


polar_theme = theme(panel.background = element_rect(fill = "white",colour = "white"),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face="bold", size=20,colour="black"),
    panel.grid.minor = element_line(colour = NA),
    strip.text.y = element_text(color="black",face="bold",size=15),
    #legend.position = "top", 
    #legend.direction="horizontal", 
    legend.text = element_text(size=13),
    legend.key = element_rect(fill = "white")
    #legend.title = element_blank()
    )

#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  print(numPlots)
  print(cols*ceiling(numPlots/cols))
  print(ceiling(numPlots/cols))
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

