library('ape')
library('cowplot')
library('ggplot2')
library('ggtree')
library('ggnewscale')
library('phytools')
library('treeio')
library('GiNA')
library(randomcoloR)
library(viridis)
library(paletteer)
library("seqinr")
library(stringr)
library(spaa)
library(ggExtra)
library(patchwork)
library(ggpubr)
library(gplots)
library(colorRamps)

# COLORS FOR SPECIES AND HOSTS IN  PHYLOGENETIC TREES

# ORF1b tree
species_cols = c("#ff0000",
                          "#980036",#10
                          "#006400", #11
                          "#c0df82", #12
                          "#4b0082", #13
                          "#48d1cc", #15
                          "#8e9aff", #18
                          "#ffa500", #19
                          "#ffff00",#2
                          "#00ff00",#3
                          "#00ffff",#5
                          "#0000ff",#6
                          "#d8bfd8",#7
                          "#ff00ff", #8
                          "#1e90ff", #9
                          "#ff69b4", #16
                          "#37f0a7", #17
                          "#db2e55", #4
                          "#ffab7e", #14
                          "#743318"
)

# ORF2 tree
species_cols_orf2 = c("#ff0000", #1
                               "#980036",#10
                               "#006400", #11
                               "#c0df82", #12
                               "#4b0082", #13
                               "#ffab7e", #14
                               "#48d1cc", #15
                               "#ff69b4", #16
                               "#37f0a7", #17
                               "#8e9aff", #18
                               "#ffa500", #19
                               "#ffff00", #2
                               "#00ff00", #3
                               "#db2e55", #4
                               "#00ffff", #5
                               "#0000ff", #6
                               "#d8bfd8", #7
                               "#ff00ff", #8
                               "#1e90ff", #9
                               "#743318"
)

# trees for whole genomes dataset
species_cols_wg = c("#ff0000", #1
                             "#980036",#10
                             "#4b0082", #13
                             "#ffff00", #2
                             "#00ff00", #3
                             "#00ffff", #5
                             "#0000ff", #6
                             "#ff00ff", #8
                             "#1e90ff", #9
                             "#743318"
)
# ORF1b tree
hosts_cols = c(
           "#B8E4DC", #bat
           "#E3DDB5", #bottlenose-dolphin
           "#785CDF", #california-sea-lion
           "#6893E0", #camel
           "#6CE6CA", #cat
           "#E05C47", #cow
           "#E056BF", #crab-eating-fox
           "#D8A08A", #deer
           "#784421ff", #dog
           "#8DE14D", #human
           "#BD82E2", #marmot
           "#A8ADD7", #mink
           "#DF9FD0", #ovibos
           "#73558C", #pig
           "#699091", #porcupine
           "#73AD80", #rabbit
           "#C4E491",  #rat
           "#DFE251", #sheep
           "#DDB155", #takin
           "#6FCDE5", #tiger
           "#C13CE8", #water-buffalo
           "#E0CFDC", #whale
           "#DB6181", #yak
           "#000000") #NA
# ORF2 tree           
hosts_cols_orf2 = c(
           "#B8E4DC", #bat
           "#E3DDB5", #bottlenose-dolphin
           "#785CDF", #california-sea-lion
           "#6893E0", #camel
           "#6CE6CA", #cat
           "#E05C47", #cow
           "#E056BF", #crab-eating-fox
           "#D8A08A", #deer
           "#784421ff", #dog
           "#8DE14D", #human
           "#BD82E2", #marmot
           "#A8ADD7", #mink
           "#DF9FD0", #ovibos
           "#73558C", #pig
           "#699091", #porcupine
           "#73AD80", #rabbit
           "#C4E491",  #rat
           "#DFE251", #sheep
           "#DDB155", #takin
           "#6FCDE5", #tiger
           "#C13CE8", #water-buffalo
           "#DB6181", #whale
           "#E0CFDC", #yak
           "#000000") #NA
           
# trees for whole genomes dataset
hosts_cols_wg = c(
           "#B8E4DC", #bat
           #         "#E3DDB5", #bottlenose-dolphin
           #         "#785CDF", #california-sea-lion
           "#6893E0", #camel
           "#6CE6CA", #cat
           "#E05C47", #cow
           "#E056BF", #crab-eating-fox
           "#D8A08A", #deer
           "#7012CD", #dog
           "#8DE14D", #human
           "#BD82E2", #marmot
           "#A8ADD7", #mink
           "#DF9FD0", #ovibos
           "#73558C", #pig
           #         "#699091", #porcupine
           "#73AD80", #rabbit
           "#C4E491",  #rat
           "#DFE251", #sheep
           #        "#DDB155", #takin
           "#6FCDE5", #tiger
           "#C13CE8", #water-buffalo
           #       "#E0CFDC", #whale
           "#E0CFDC", #yak
           "#000000") #NA
           



# Trees for ORF1b and ORF2 
plot_annotated_tree = function(tree_file, meta, species_cols, colors_hosts, type){
  tree =  read.tree(tree_file)
  tree_rooted = midpoint.root(tree)
  info = read.csv(meta)
  
  #num_hosts = nrow(unique(info['host']))
  #print(num_hosts)
  #colors_hosts = distinctColorPalette(num_hosts)
  #print(colors_hosts)
  if (type == 'ORF2'){
    num_clusters_aa = max(info['orf2_aa'])
    print(num_clusters_aa)
    colors_aa = distinctColorPalette(num_clusters_aa+1)
    #cluster_aa = data.frame(aa = as.factor(info[,c(paste(tolower(type),"aa",sep="_"))]))
    cluster_aa = data.frame(aa = as.factor(info[,"aa_mmseq2"]))
    rownames(cluster_aa) <- info$GBAC
  }
  
  host = data.frame("host" = info[,c("host")])
  rownames(host) <- info$GBAC
  #print(host)
  
  species = data.frame("species" = info[,c("taxon")])
  rownames(species) <- info$GBAC
  #print(species)
  
  
  t = ggtree(tree_rooted, size=0.75) %<+% info +
    geom_point2(aes(label=label, 
                    subset = !is.na(as.numeric(label)) & as.numeric(label) < 80), size=1, color="red",alpha=0.5) +
    geom_treescale() #+ geom_tiplab(size =1)
  
  g1 = gheatmap(t, species,                           # we add a heatmap layer of the gender dataframe to our tree plot
                #offset = 1,                               # offset shifts the heatmap to the right,
                width = 0.05,# width defines the width of the heatmap column,
                colnames_angle=45,
                color = NULL)+
    scale_fill_manual(values=species_cols)
  if (type == 'ORF2'){
    g1 = g1 + new_scale_fill()
    g1 = gheatmap(g1, cluster_aa, width = 0.05,
                  colnames_position = "top",
                  offset=0.25,
                  #colnames_offset_y = 1,
                  colnames_angle=45)+
      scale_fill_manual(values=colors_aa)
  }
  g2 = g1 + new_scale_fill()
  
  g3 = gheatmap(g2, host, width = 0.05,
                colnames_position = "top",
                offset=0.50,
                #colnames_offset_y = 1,
                colnames_angle=45)+
    scale_fill_manual(values=colors_hosts)
  
  
  return(g3)
}

# Tree with leaves colored in gradient
plot_annotated_grad_mm_tree = function(tree_file, meta, species_cols, colors_hosts, type){
  tree =  read.tree(tree_file)
  tree_rooted = midpoint.root(tree)
  info = read.csv(meta)
  print(type)
  if (type != 'ORF1a'){
    if (type == 'ORF1b'){
      num_clusters_nt = max(info['ORF1b_nt_17'])
      num_clusters_aa = max(info['ORF1b_aa_7'])
      
      colors_nt = distinctColorPalette(num_clusters_nt+1)
      colors_aa = distinctColorPalette(num_clusters_aa+1)
      
      cluster_nt = data.frame(nt = as.factor(info[,'ORF1b_nt_17']))
      rownames(cluster_nt) <- info$GBAC
      
      cluster_aa = data.frame(aa = as.factor(info[,'ORF1b_aa_7']))
      rownames(cluster_aa) <- info$GBAC
      
      print(num_clusters_nt)
      print(num_clusters_aa)
    }
    else{
      num_clusters_aa = max(info['ORF2_aa_25'])
      colors_aa = distinctColorPalette(num_clusters_aa+1)
      
      cluster_aa = data.frame(aa = as.factor(info[,'ORF2_aa_25']))
      rownames(cluster_aa) <- info$GBAC
      
      print(num_clusters_aa)
      
      
    }
  }
  
  host = data.frame("host" = info[,c("host")])
  rownames(host) <- info$GBAC

  species = data.frame("species" = info[,c("taxon")])
  rownames(species) <- info$GBAC
  
  head(info)
  
  t = ggtree(tree_rooted, size=0.75) %<+% info + geom_tiplab(size =1, aes(color=code)) +
    scale_color_manual(values=info$color_orf1b16) +
    geom_point2(aes(label=label, 
                    subset = !is.na(as.numeric(label)) & as.numeric(label) < 80), size=1, color="red",alpha=0.5) +
    geom_treescale()
  
  g0 = gheatmap(t, species, #heatmap layer of species names
                offset = 0.25,
                width = 0.05,
                colnames_angle=45,
                colnames_position = "top",
                color = NULL)+
    scale_fill_manual(values=species_cols)
  g0 = g0 + new_scale_fill()

  if (type != 'ORF1a'){
    
    if (type == 'ORF1b'){
    g0 = gheatmap(g0, cluster_nt, width = 0.05, #heatmap layer of clusters in ORF1b sequences defined by MMSeqs2
                  colnames_position = "top",
                  offset=0.40,
                  #colnames_offset_y = 1,
                  colnames_angle=45)+
      scale_fill_manual(values=colors_nt)
    g0 = g0 + new_scale_fill()
    }
    g0 = gheatmap(g0, cluster_aa, width = 0.05,
                  colnames_position = "top",
                  offset=0.55,
                  #colnames_offset_y = 1,
                  colnames_angle=45)+
      scale_fill_manual(values=colors_aa)
    
  }

  g0 = g0 + new_scale_fill()
  g1 = gheatmap(g0, host, # heatmap layer of host names
                offset = 0.8,
                width = 0.05,
                colnames_angle=45,
                colnames_position = "top",
                color = NULL)+
    scale_fill_manual(values=colors_hosts)
  
  
  
  
  return(g1)
}

#DISTANCES AND HISTOGRAMS

plot_aa_nt_hists_sp = function(aln_object, pos_start, pos_end, col_sp, pairwise=TRUE){
  
  dna_slice = aln_object[1:length(aln_object[,1]), seq(from = pos_start, to = pos_end, by=1)]
  dna_slice_char=as.character.DNAbin(dna_slice)
  dna_slice_char[dna_slice_char=='-'] <- NA
  
  
  aa_slice = trans(dna_slice)
  aa_slice_char=as.character.AAbin(aa_slice)
  aa_slice_char[aa_slice_char=='X'] <- NA
  
  dist_nn = dist.gene(dna_slice_char, method = "percentage", pairwise.deletion = pairwise)
  dist_aa = dist.gene(aa_slice_char, method = "percentage", pairwise.deletion = pairwise)
  
  dist_nn= dist2list(dist_nn)
  dist_aa= dist2list(dist_aa)
  
  
  dist_nn = check_param(dist_nn,col_sp)
  colnames(dist_nn)[length(colnames(dist_nn))] = "same_species"
  d = merge(dist_nn, dist_aa, by=c("row", "col"))
  
  
  #scatterplot with marginal histograms
  p <- ggplot(d,aes(value.x,value.y, color=same_species)) + geom_point(alpha=0.1) + theme_bw() +
    xlab("nucleotide distance")+ylab("amino acid distance")
  
  return(list(d, p))
  
}

plot_hists_zoom = function(df_dist1, df_dist2){
  
  hist1_nt = ggplot(df_dist1, aes(x=value.x, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50, position="identity") +
    theme_bw() + xlab("nucleotide distance") + ggtitle(" ") +
    theme(legend.position="none")+
    theme(legend.justification=c(0.1,0.9), legend.position=c(0.1,0.9)) +
    scale_x_continuous(limits = c(0, 0.4), oob = function(x, limits) x) + 
    scale_y_continuous(limits = c(0, 3000), oob = function(y, limits) y)
  
  
  hist1_aa = ggplot(df_dist1, aes(x=value.y, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50, position="identity") +
    theme_bw() + xlab("amino acid distance") + ggtitle(" ")  + theme(legend.position="none")+
    scale_x_continuous(limits = c(0, 0.4), oob = function(x, limits) x) + 
    scale_y_continuous(limits = c(0, 3000), oob = function(y, limits) y)
  
  hist2_nt = ggplot(df_dist2, aes(x=value.x, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50, position="identity") +
    theme_bw() + xlab("nucleotide distance") + ggtitle(" ") + 
    theme(legend.position="none")+
    scale_x_continuous(limits = c(0, 0.4), oob = function(x, limits) x) + 
    scale_y_continuous(limits = c(0, 3000), oob = function(y, limits) y)
  
  hist2_aa = ggplot(df_dist2, aes(x=value.y, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50, position="identity") +
    theme_bw() + xlab("amino acid distance")+ ggtitle(" ")  + theme(legend.position="none")+
    scale_x_continuous(limits = c(0, 0.4), oob = function(x, limits) x) + 
    scale_y_continuous(limits = c(0, 3000), oob = function(y, limits) y)
  
  
  g1 = ggarrange(hist1_nt, hist2_nt,
                 hist1_aa, hist2_aa,
                 #labels = c("A", "B", "C", "D"),
                 ncol = 2, nrow = 2)
  

  
  return(g1)
}

plot_hists = function(df_dist1, df_dist2, title1, title2){
  
  hist1_nt = ggplot(df_dist1, aes(x=value.x, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50) +
    theme_bw() + xlab("nucleotide distance") + ggtitle(title1) +
    theme(legend.position="none")+
    theme(legend.justification=c(0.1,0.9), legend.position=c(0.1,0.9))
  
  
  hist1_aa = ggplot(df_dist1, aes(x=value.y, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50) +
    theme_bw() + xlab("amino acid distance") + ggtitle(" ")  + theme(legend.position="none")
  
  hist2_nt = ggplot(df_dist2, aes(x=value.x, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50) +
    theme_bw() + xlab("nucleotide distance") + ggtitle(title2) + 
    theme(legend.position="none")
  
  hist2_aa = ggplot(df_dist2, aes(x=value.y, fill = same_species)) +
    geom_histogram(alpha=0.5, bins=50) +
    theme_bw() + xlab("amino acid distance")+ ggtitle(" ")  + theme(legend.position="none")
  
  
  g1 = ggarrange(hist1_nt, hist2_nt,
                 hist1_aa, hist2_aa,
                 #labels = c("A", "B", "C", "D"),
                 ncol = 2, nrow = 2)
  
  
  
  return(g1)
}

# FOR PDC plots

check_param = function(d, colnum){
  
  spl_df1 = str_split(d$row, '/', simplify = TRUE)
  spl_df2 = str_split(d$col, '/', simplify = TRUE)
  spl_df1_param = spl_df1[,colnum]
  spl_df2_param = spl_df2[,colnum]
  
  any_NA = spl_df1_param == "NA" | spl_df2_param == "NA"
  same_param = spl_df1_param == spl_df2_param
  
  vector_param = vector("character", length(spl_df1_param))
  vector_param[same_param == TRUE] = "yes"
  vector_param[same_param == FALSE] = "no"
  vector_param[any_NA == TRUE] = "unknown"
  d = cbind(d, vector_param)
  return(d)
}

get_host_type = function(col_host, list_host_type){
  c_host_type = c()
  for (i in 1:length(col_host)){
    c_host_type = append(c_host_type,list_host_type[[col_host[i]]])
  }
  return(c_host_type)
}

check_param_host_type = function(d, colnum,list_host_type){
  
  spl_df1 = str_split(d$row, '/', simplify = TRUE)
  spl_df2 = str_split(d$col, '/', simplify = TRUE)
  spl_df1_param = spl_df1[,colnum]
  spl_df2_param = spl_df2[,colnum]
  
  spl_df1_hosttype = get_host_type(spl_df1_param,list_host_type)
  spl_df2_hosttype = get_host_type(spl_df2_param,list_host_type)
  
  
  any_NA = spl_df1_hosttype == "NA" | spl_df2_hosttype == "NA"
  same_param = spl_df1_hosttype == spl_df2_hosttype
  
  vector_param = vector("character", length(spl_df1_param))
  vector_param[same_param == TRUE] = "yes"
  vector_param[same_param == FALSE] = "no"
  vector_param[any_NA == TRUE] = "unknown"
  d = cbind(d, vector_param)
  return(d)
}

# FIGURE 1

setwd("../data")

orf1b_all_g = plot_annotated_tree("MAV_ORF1b_GenBank.nwk","MAV_ORF1b_GenBank_metadata.csv",species_cols, hosts_cols,"ORF1b")
orf1b_all_g
ggsave('../Fig1_tree_ORF1b.png', height=12, width=12)


orf2_all_g = plot_annotated_tree("MAV_ORF2_GenBank.nwk","MAV_ORF2_GenBank_metadata.csv",species_cols_orf2,hosts_cols_orf2,"ORF2")
orf2_all_g
ggsave('../Fig1_tree_ORF2.png', height=12, width=15)

#FIGURE2

# ORF2

aln_mav_ORF2_lessgp = read.dna("MAV_ORF2_GenBank.fasta", format="fasta")
ph_ORF2_all_sp  = plot_aa_nt_hists_sp(aln_mav_ORF2_lessgp, 1, length(aln_mav_ORF2_lessgp[1,]), 6)
df_dist_ORF2_sp = ph_ORF2_all_sp[[1]]

# ORF1B

aln_mav_ORF1B_lessgp = read.dna("MAV_ORF1b_GenBank.fasta", format="fasta")
ph_ORF1B_all_sp  = plot_aa_nt_hists_sp(aln_mav_ORF1B_lessgp, 1, length(aln_mav_ORF1B_lessgp[1,]), 6)
df_dist_ORF1B_sp = ph_ORF1B_all_sp[[1]]

hists = plot_hists(df_dist_ORF2_sp, df_dist_ORF1B_sp, "ORF2 sequences from GenBank (N=894)", "ORF1b sequences from GenBank (N=522)")
hists


hists_zoomed = plot_hists_zoom(df_dist_ORF2_sp, df_dist_ORF1B_sp)
hists_zoomed

ggsave("../Fig2.png", hists, width=10, height=6, dpi = 1000)
ggsave("../Fig2_zoomed.png", hists_zoomed, width=10, height=6, dpi = 1000)


#FIGURE 5

grad_trees= list.files(path = ".", full.names = TRUE, pattern = ".nwk$")

for (file in grad_trees){
  print(file)
  g = plot_annotated_grad_mm_tree(file,
                               "MAV_wg_metadata.csv",
                               species_cols_wg, hosts_cols_wg, strsplit(strsplit(file, '/')[[1]][2], '.nwk')[[1]][1])  + 
    theme(legend.position = "none") +
    ggtitle(strsplit(strsplit(file, '/')[[1]][2], '.nex')[[1]][1])
  ggsave(paste(file,"mm", "svg",sep="."),height=10, width=7)
}

#FIGURES 4,6

fig_path = "../"

astro_hosts = t(read.csv("MAV_hosts_upd.csv", header=FALSE))
colnames(astro_hosts) <- lapply(astro_hosts[1, ], as.character)
astro_hosts <- astro_hosts[-1,] 

mav = read.dna("MAV_wg.fasta", format="fasta", as.character=TRUE)
mav[mav=='-'] <- NA

list_regs_mav <- list(
  list(r1="ORF1ab",s1=1,e1=3705,r2="ORF2",s2=3706,e2=5487, n1="ORF1ab", n2="ORF2"),
  list(r1="ORF1a",s1=1,e1=2217,r2="ORF1b",s2=2218,e2=3705, n1="ORF1a", n2="ORF1b"),
  list(r1="ORF1a",s1=1,e1=1108,r2="ORF1a",s2=1109,e2=2217, n1="5' half of ORF1a", n2="5' half of ORF1a"),
  list(r1="ORF1b",s1=2218,e1=2961,r2="ORF1b",s2=2962,e2=3705, n1="5' half of ORF1b", n2="5' half of ORF1b"),
  list(r1="ORF2",s1=3706,e1=4596,r2="ORF2",s2=4597,e2=5487, n1="5' half of ORF2", n2="3' half of ORF2")
)

for(reg in list_regs_mav){
  print(paste('MAV', reg$n1,reg$n2, sep=', '))
  
  l_PDCP <-  plot_PDCP(mav, reg$s1, reg$e1, reg$s2, reg$e2)
  # adding new columns where the hosts and MAV species' equality is checked
  dist_df = l_PDCP[[2]]
  dist_df = check_param(dist_df, 3)
  colnames(dist_df)[length(colnames(dist_df))] = "same_host"

  dist_df = check_param(dist_df, 5)
  colnames(dist_df)[length(colnames(dist_df))] = "same_species"
  
  plot = ggplot(dist_df)+geom_point(aes(value.x,value.y,colour=same_host), alpha = 0.1) + theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    xlab(reg$n1) + ylab(reg$n2) + theme(text=element_text(size=21, family="Arial"))

  ggsave(file=paste(fig_path,"mav_",paste(reg$r1, reg$r2, sep="_vs_"), "_host_upd.png", sep=""),
         plot=plot,
         # width=168, height=168,
         width=126, height=126,
         units="mm")

  plot = ggplot(dist_df)+geom_point(aes(value.x,value.y,colour=same_species), alpha = 0.1) + theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    xlab(reg$n1) + ylab(reg$n2) + theme(text=element_text(size=21, family="Arial"))
  ggsave(file=paste(fig_path,"mav_",paste(reg$r1, reg$r2, sep="_vs_"), "_species_upd.png", sep=""),
         plot=plot,
         # width=168, height=168,
         width=126, height=126,
         units="mm")
}
