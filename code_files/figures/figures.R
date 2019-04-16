figure_ranges = function(file, plotfile = NULL) {
  load(file)
  data = get("pbdb")
  
  df = data.frame(id=1:length(data[,1]), min = data[,1], max = data[,2])
  library(ggplot2)
  pl = ggplot(df, aes(x=id))+ xlab("") + ylab("Time before the present (My)") +
    geom_linerange(aes(ymin=min,ymax=max),linetype=1,color="blue")+
    geom_point(aes(y=min),size=1,color="red")+
    geom_point(aes(y=max),size=1,color="red")+
    theme_bw()
  
  if(!is.null(plotfile)) ggsave(paste0(plotfile,"_ranges.pdf"))
  else print(pl)
}

plot_comparison = function(file, names = c("divtime", "diversificationRateFBD","turnoverFBD", "samplingProportionFBD"), outname = NULL) {
  library(ggplot2)
  load(file)
  if(is.null(outname)) outname = tools::file_path_sans_ext(basename(file))
  fullfile = get(tools::file_path_sans_ext(basename(file)))
  
  plotnames = c("Divergence time", "Diversification rate","Turnover", "Sampling proportion")
  names(plotnames) = c("extmrca.divtime", "diversificationRateFBD","turnoverFBD", "samplingProportionFBD")
  
  for(name in names) {
    full = fullfile[[name]]
    cats = c("correct_ages", "median_ages", "interval_ages", "random_ages", "norm_intvl_ages")
    if(name == "divtime") name = "extmrca.divtime"
    
    dfrel = data.frame(Age = c(rep("Correct age",100), rep("Median age",100), rep("*Interval ages",100), rep("Random age",100), rep("*Symmetric ages",100)), 
                       Relative_error = unlist(lapply(cats, function(x) {full[[paste0(x,".",name,".relative_error")]]})))
    
    dfcov = data.frame(Age = c("Correct age","Median age","*Interval ages","Random age","*Symmetric ages"), 
                       Coverage = sapply(cats, function(x) {mean(full[[paste0(x,".",name,".coverage")]])}))
    
    p2 <- ggplot(dfrel, aes(x = Age, y = Relative_error)) + geom_boxplot() + xlab(plotnames[name]) + #theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 15)) + 
      ylim(0,0.8) + ylab("Relative error") + scale_x_discrete(limits = c("Correct age","*Interval ages","Median age","Random age","*Symmetric ages"))
    
    ggsave(paste0("~/Work/FBD_figures/",outname, "_",name,"_rel_error.pdf"),width = 5,height = 5)
    
    p2 <- ggplot(dfcov, aes(x = Age, y = Coverage)) + geom_point() + xlab(plotnames[name]) + #theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), text = element_text(size = 15)) + 
      ylim(0,1) + scale_x_discrete(limits = c("Correct age","*Interval ages","Median age","Random age","*Symmetric ages"))
    
    ggsave(paste0("~/Work/FBD_figures/",outname, "_",name,"_coverage.pdf"),width = 5,height = 5)
  }
}

cetaceans_plot = function(outfolder, plotfile = NULL) {
  library(ggplot2)
  estnames = c("diversificationRateFBD","turnoverFBD", "samplingProportionFBD")
  pnames = c("Diversification rate","Turnover", "Sampling proportion")
  names = c("median_ages", "random_ages", "interval_ages")
  
  logs = list()
  for(x in names) {
    logs[[x]] = read.table(file.path(outfolder, paste0("Cetaceans_", x,".log")),header = T,comment.char = "#")
    n = nrow(logs[[x]])
    logs[[x]] = logs[[x]][(round(n/10)+1):n,]
  }
  
  cbPalette <- c("#56B4E9", "#009E73", "#E69F00", "#CC79A7")
  p = list()
  for(i in 1:length(estnames)) {
    e = estnames[i]
    df = data.frame(values = c(), ds = c())
    for(x in names) {
      df = rbind(df, data.frame(values = logs[[x]][[e]], ds = x))
    }
    px = ggplot(df, aes(x = ds, y = values, fill = ds)) + geom_boxplot(notch = TRUE) + theme_bw() +
      scale_x_discrete(labels = c("Median ages", "Random ages", "Interval ages")) + xlab("") + ylab("") + 
      ggtitle(pnames[i]) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values=cbPalette)
    p[[i]] = px
  }
  
  p = gridExtra::grid.arrange(grobs = p, ncol = 3)
  if(!is.null(plotfile)) ggsave(plotfile, plot = p,width = 10, height = 5)
  else print(p)
}

plot_node_comparison = function(outfolder, string, plotfile = NULL) {
  library(ggplot2)
  
  names = c("median_ages", "random_ages", "interval_ages")
  
  logs = list()
  for(x in names) {
    logs[[x]] = read.table(file.path(outfolder, paste0("Cetaceans_", x,".log")),header = T,comment.char = "#")
    n = nrow(logs[[x]])
    logs[[x]] = logs[[x]][(round(n/10)+1):n,]
    logs[[x]] = logs[[x]][[string]]
  }
  
  df = data.frame(values = logs[[names[3]]], Age = "Interval ages")
  df = rbind(df, data.frame(values = logs[[names[1]]], Age = "Median ages"))
  df = rbind(df, data.frame(values = logs[[names[2]]], Age = "Random ages"))
  
  # change default colors to colorblind-friendly
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
  
  p = ggplot(df, aes(x = values, y=..scaled.., col = Age)) + geom_density(n = 2^8) + theme_bw() + xlab("") + ggtitle("Divergence time of Phocoenidae") +
    geom_hline(color = "black", yintercept = 0, size = 0.5) + xlim(12,22) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Density (scaled to [0;1])") +
    scale_colour_manual(values=cbPalette)
  
  if(is.null(plotfile)) print(p)
  else ggsave(plotfile, height = 5, width = 10)
}