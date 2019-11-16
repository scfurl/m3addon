#' Segment a Range
#'
#' @description This function takes a ordered range of 2 values, fractionates the range based 
#' on input fraction and returns the range of the nth fraction
#' @param cds Input cell_data_set object.
#' @export

segment<-function(range, fraction=0.2, n=1){
  #this function takes a ordered range of 2 values, fractionates the range based on input fraction and returns the range of the nth fraction
  step_length<-(range[2]-range[1])*fraction
  n_steps<-1/fraction
  step_vec<-1:n_steps
  vals<-c(range[1],range[1]+step_vec*step_length)
  return(c(vals[n], vals[n+1]))
}

#' Plot a fixed size for plot and legend
#'
#' @description Creates a grid object comprised of both a plot (ggplot only) and a legend.  Function takes
#' a ggplot object as input and returns a grid arrange object with plot on left, legend on right and sizes as dictated
#' by width argument

#' @param a.gplot A ggplot2 object
#' @param widths A vector of length 2 specifying the size of the plot grid and legend grid, respectively.
#' @return grid arranged plot
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#' @references Puram, S. V. et al. Single-Cell Transcriptomic Analysis of Primary and 
#' Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 171, 1611.e1–1611.e24 (2017).


split_plot<-function(a.gplot, widths=c(4.5/6,1.5/6)){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  grid.arrange( a.gplot+ theme(legend.position = 'none'), legend, 
                ncol=2, nrow=1, widths=widths)
}


#' Reduced Dimensionality Plot
#'
#' @description Adds an icon to the lower left hand corner of a plot indicating reducted dimensionality method used
#' 

#' @param plot A ggplot2 object
#' @param method method indicating which icon to use.  Currently 'tSNE' and 'UMAP' supported
#' @param widths A vector of length 2 specifying the size of the plot grid and legend grid, respectively.
#' @param correction_x, integer to offset icon in x axis
#' @param correction_y, integer to offset icon in y axis
#' @return grid arranged plot using split_plot function
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom png readPNG
#' @importFrom grid rasterGrob
#' @export
#' @references Puram, S. V. et al. Single-Cell Transcriptomic Analysis of Primary and 
#' Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 171, 1611.e1–1611.e24 (2017).


red_dim_plot<-function(plot, method="UMAP", size=20, ICON_DIR=system.file("icons", "m3addon"), correction_x=0, correction_y=0){
  xs<-segment(layer_scales(plot)$x$range$range)
  ys<-segment(layer_scales(plot)$y$range$range)
  if(!is.null(plot$layers$title)){
    title=plot$layers$title
    plot_title=TRUE
  }else{
    plot_title=FALSE
  }
  plot+theme_void()+geom_segment(aes(x=xs[1], xend = xs[2] , y=ys[1], yend = ys[2]), size=1.5,
                              arrow = arrow(length = unit(0.2,"cm"))) 
  img_name<-paste0(method, "_", size, "x.png")
  img <- readPNG(file.path(ICON_DIR, img_name))
  g <- rasterGrob(img, interpolate=TRUE)
  if(!plot_title) {
    return(split_plot(plot+theme_void()+annotation_custom(g, xmin=xs[1]-correction_x, ymin=ys[1]-correction_y, xmax = xs[2], ymax=ys[2])))
  }else{
    return(split_plot(plot+theme_void()+ggtitle(label=title)+annotation_custom(g, xmin=xs[1]+correction_x, ymin=ys[1]+correction_y, xmax = xs[2], ymax=ys[2]))+ theme(plot.title = element_text(size=18)))
  }
}
