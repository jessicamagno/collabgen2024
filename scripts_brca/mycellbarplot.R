my_cell_bar_plot<- function(input, id = "Samples", title = "Cell Fraction", features = NULL, pattern = NULL, legend.position = "bottom",
                         coord_filp = TRUE, palette = 3, show_col = F, cols = NULL){
  
  input<-as.data.frame(input)
  colnames(input)[which(colnames(input)==id)]<-"Samples"
  
  if(is.null(features)){
    
    if(is.null(pattern)) stop(">>>=== The 'pattern' parameter must be defined...")
    feas <- colnames(input)[str_detect(colnames(input), pattern = pattern)]
  }else{
    
    feas <- features
  }
  
  input <- input[, c("Samples", feas)]
  
  input<-remove_names(input_df = input, variable = "colnames", patterns_to_na = patterns_to_na, patterns_space = "_")
  ##################
  if(legend.position == "top"|legend.position=="bottom") {
    legend.direction<-"horizontal"
  }else{
    legend.direction<-"vertical"
  }
  
  
  if(!is.null(cols)){
    cols<-cols
  }else{
    if(is.null(palette)){
      cols<-IOBR::palettes(category = "random", palette = 4, show_col = show_col, show_message = T)
    }else{
      cols<-IOBR::palettes(category = "random", palette = palette, show_col = show_col, show_message = T)
    }
  }
  
  
  if(coord_filp){
    pp<-input %>%
      tidyr::gather(cell_type,fraction, -Samples) %>%
      # plot as stacked bar chart
      ggplot(aes(x=Samples, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      theme_light()+
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = rev(levels(input)))+
      ggtitle(paste0(title))+
      labs(x = "Samples", y = "Fraction") +
      theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.text.x= element_text(face="plain",angle=0,hjust = 1,color="black"),
            axis.text.y= element_text(face="plain",angle= 30,hjust = 1,color="black"))+
      theme(legend.title = element_blank(),
            legend.position= legend.position,
            legend.direction= legend.direction,
            legend.justification=c(.5,.5),
            legend.box="horizontal",
            legend.box.just="top")
  }else{
    pp<- input %>%
      tidyr::gather(cell_type,fraction, -Samples) %>%
      # plot as stacked bar chart
      ggplot(aes(x=Samples, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      # coord_flip() +
      theme_minimal()+
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = rev(levels(input)))+
      ggtitle(paste0(title))+
      labs(x = "BRCA LumB Samples", y = "Fraction") +
      theme(plot.title=element_text(size=rel(2),hjust=0.5),
            axis.text.x = element_blank(),
            axis.ticks.x=element_blank(),
            # axis.text.x= element_text(face="plain",angle=0,hjust = 1,
            #                           color="white"),
            axis.text.y= element_text(face="plain",angle=30,hjust = 1,
                                      color="black"))+
      theme(legend.title = element_blank(),
            legend.position= legend.position,
            legend.direction= legend.direction,
            legend.justification=c(.5,.5),
            legend.box="horizontal",
            legend.box.just="top")
  }
  
  print(pp)
  return(pp)
  
}
