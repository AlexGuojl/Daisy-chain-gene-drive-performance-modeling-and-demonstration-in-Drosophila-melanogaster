rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
#输入原始数据，GL or nonGL，varied_parameter
#输入原始数据，GL or nonGL，varied_parameter




heatmap <- function(df, is_gl_or_not,varied_parameter) {
  if (is_gl_or_not=="non_GL"){  #get parameters and pop_size
    df1 <- df[-1,]
    ##8是dropsize，9是embres
    parameter_varied<- select(df1,V9)
    drop_varied<- select(df1,V8)
    
    #parameter_varied$varied_parameter <- parameter_varied$V9
    drop_varied$drop_size <- drop_varied$V8
    
    df_parameters <- cbind(drop_varied,parameter_varied)
    df_parameters <- select(df_parameters,drop_size,V9)#选择修改的两个参数
    
    #df_parameters$drop_size <- as.numeric(df_parameters$drop_size)
    # df_parameters$V9<- as.numeric( df_parameters$V9)
    #next：提取pop_size
    interval <- 9#每9列提取一次，永远不变
    # 使用 seq() 函数生成要选取的列索引
    selected_columns <- seq(3, ncol(df1), by = interval)
    # 根据列索引选取dataframe的几列，生成只含有popsize的dataframe：df_pop
    df_pop <- df1[, selected_columns]
    
    #接下来：提取popsize最小值，看看是不是0
    df_pop <- as.data.frame(lapply(df_pop, as.numeric))
    
    popsize_min <- vector()
    for (i in 1:nrow(df_pop))
    {
      popsize_min = c(popsize_min,min(df_pop[i,]))
    }
    df_parameters$pop_size<- popsize_min
    ##successful suppression rate based on popsize
    df_parameters$pop_size[which(df_parameters$pop_size == 0)] <- 1
    df_parameters$pop_size[which(df_parameters$pop_size > 1)] <- 0
    
    df_mean_pop<-aggregate(df_parameters$pop_size, by=list(df_parameters$drop_size,
                                                           df_parameters$V9),mean) 
    #十次重复里的平均成功率！！
    df_mean_pop$successful_rate <- df_mean_pop$x
    
    x_label = "Relative Release Size"
    y_label <- gsub("_", " ", varied_parameter)
    p_heatmaps <- ggplot(df_mean_pop, aes(x=Group.1,y=Group.2)) + 
      geom_tile(aes(fill=successful_rate )) + 
      scale_fill_gradientn(values = c(0,0.25,0.5,0.75,1),#limits = c(0.0,1.0),
                           colors = c("#330033","#CC66FF", "#f6f5ee","#FFCC66","#CC0000"))+ #"#FFF5F0","#FCBBA1","#A50F15"))+  # colors = c("#1c54a8", "#f6f5ee","#ee2d2a"))+
      theme_classic()+
      labs(x = x_label, y = y_label,title = " ",fill = "Population Elimination Rate") +theme(
        axis.text = element_text(size = 14),       # 坐标轴刻度
        axis.title = element_text(size = 15),      # 坐标轴标题
        legend.text = element_text(size = 14),     # 图例刻度
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 16)     # 图例标题
      )
    return(p_heatmaps)
  }
  else{  #此时提取的应该是最大gl而不是最小的pop_size!!!稍后改
    df1 <- df[-1,]
    ##8是dropsize，9是emb_res, parameters无需更改！！
    parameter_varied<- select(df1,V9)
    drop_varied<- select(df1,V8)
    
    #parameter_varied$varied_parameter <- parameter_varied$V9
    drop_varied$drop_size <- drop_varied$V8
    
    df_parameters <- cbind(drop_varied,parameter_varied)
    df_parameters <- select(df_parameters,drop_size,V9)#选择修改的两个参数
    
    interval <- 10#每9列提取一次，永远不变
    # 使用 seq() 函数生成要选取的列索引
    selected_columns <- seq(10, ncol(df1), by = interval)
    
    # 根据列索引选取dataframe的几列，生成只含有popsize的dataframe：df_pop
    df_gl <- df1[, selected_columns]
    df_gl<-as.data.frame(lapply(df_gl, as.numeric))
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    df_gl[is.nan(df_gl)] <- 1
    #get largest genetic load in different generations
    geneticload_max <- vector()
    for (i in 1:nrow(df_gl))
    {
      geneticload_max = c(geneticload_max,max(df_gl[i,]))
    }
    df_parameters$genetic_load<- geneticload_max
    
    
    ##successful suppression rate based on popsize
    #get mean gl
    df_mean_gl<-aggregate(df_parameters$genetic_load, by=list(df_parameters$drop_size,
                                                              df_parameters$V9), mean) 
    df_mean_gl$genetic_load <- df_mean_gl$x
    
    x_label = "Relative Release Size"
    y_label <- gsub("_", " ", varied_parameter)
    
    p_heatmaps <- ggplot(df_mean_gl, aes(x=Group.1,y=Group.2)) + 
      geom_tile(aes(fill=genetic_load )) + 
      scale_fill_gradientn(values = c(0,0.5,1),limits = c(-0.01,1.1),
                           colors =c("#330033","#CC66FF", "#f6f5ee","#FFCC66","#CC0000"))+theme_classic()+ #c("#54278F", "#f6f5ee","#FEC44F"))
      labs(x = x_label, y = y_label,title = " ",fill = "Maximum Genetic Load")  +theme(
        axis.text = element_text(size = 14),       # 坐标轴刻度
        axis.title = element_text(size = 15),      # 坐标轴标题
        legend.text = element_text(size = 14),     # 图例刻度
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 16)  # 图例标题
      )
    return(p_heatmaps)
  }
}




#“Drop Size” should be “Relative Release Size
#这张图怪的原因是把0.99超高释放量中释放的drive个体也算进了population中！！
heatmap_onedim_meanpop <- function(df, varied_parameter) {
  library(scales)
  df <- df[-1,]
  df <- df[,-12]
  names(df) <- c("result","size_all","generation" , "drop_size",varied_parameter,"slice" ,
                 "d1_freq","d2_freq" ,"d3_freq" ,"pop_size",
                 "decrease_proportion")
  df$drop_size<- as.numeric(df$drop_size)
  df$drop_size[which(df$drop_size == 0.99)]<- 1.00
  #drop_varied<- select(df1,V4)
  #parameter_varied<- select(df1,V5)
  
  timespan <- 60
  df$generation <- as.numeric(df$generation)
  df <- subset(df,generation <= timespan)
  df <- subset(df,generation >0)
  #df <- subset(df, slice<8)
  
  df_meanpop <- df %>%
    group_by(drop_size,  .[[5]],slice, generation) %>%
    summarize(mean_pop_size = mean(as.numeric(pop_size)))
  df_pop_in_release_region<- df_meanpop %>%
    group_by(drop_size,.[[2]],generation) %>%
    summarize(sum_pop_in_release_region = sum(mean_pop_size))
  x_label = "Drop_Size"
  y_label <- gsub("_", " ", varied_parameter)
  df_meanpop_in_release_region<- df_pop_in_release_region %>%
    group_by(drop_size, .[[2]]) %>%
    summarize(mean_pop_in_release_region = mean(sum_pop_in_release_region))
  df_meanpop_in_release_region$drop_size <- as.numeric(df_meanpop_in_release_region$drop_size)
  lowest_pop<- min(df_meanpop_in_release_region$mean_pop_in_release_region)
  highest_pop <- max(df_meanpop_in_release_region$mean_pop_in_release_region)
  p_heatmaps <- ggplot(df_meanpop_in_release_region, aes(x = as.numeric(drop_size), y = as.numeric(df_meanpop_in_release_region$`.[[2]]`))) + 
    geom_tile(aes(fill = mean_pop_in_release_region)) + 
    scale_fill_gradientn(colors = c("#1c54a8", "#f6f5ee", "#ee2d2a"), limits = (c(lowest_pop-0.001,highest_pop+0.001)),
                         values = rescale(c(10200, 9800, 9400))) + 
    theme_classic()+  labs(x = x_label, y = y_label,title = " ",fill = "Average Population Size") +theme(
      axis.text = element_text(size = 14),       # 坐标轴刻度
      axis.title = element_text(size = 15),      # 坐标轴标题
      legend.text = element_text(size = 14),     # 图例刻度
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 16)     # 图例标题
    )
  p_heatmaps
}
heatmap_onedim <- function(df, varied_parameter) {
  df <- df[-1,]
  df <- df[,-12]
  #drop_varied<- select(df1,V4)
  #parameter_varied<- select(df1,V5)
  names(df) <- c("result","size_all","generation" , "drop_size",varied_parameter,"slice" ,
                 "d1_freq","d2_freq" ,"d3_freq" ,"pop_size",
                 "decrease_proportion")
  df$drop_size <- as.numeric(df$drop_size)
  df$drop_size[which(df$drop_size == 0.99)] <- 1.00
  
  #数据分析与绘图
  #数据结构：每组dropsize和drop radius对应50个generations，每个generation对应50个slices
  #首先——遍历整个数据框，对ds大于0.9的标记为1，else则标记为0
  df <- df %>%
    mutate(width = ifelse(as.numeric(decrease_proportion) >= 0.9, 1, 0))
  #df1 <- filter(df,as.numeric(df$drop_size) == 0)
  #return(df1)}
  
  #dataframe有drop_size、drop_radius、decrease_proportion,width，generation和slice几列，需要
  #根据drop_size、drop_radius、decrease_proportion和generation这四个分组元素，对每一个slice里的width求平均值，获得新的df，保存为df_meanwidth
  #接下来再根据dropsize,dropradius和ganeration分组，对分组内每一个generation内的width求和，和存储为sum_width
  #接下来找出不同drop_size和drop_radius所对应的最大的sum_width，绘制热图
  #根据drop_size、drop_radius、slice和generation这四个分组元素，
  #对每一个slice里的width求平均值，获得新的df，保存为df_meanwidth
  #(df = heatmap_1d_dd)
  df_meanwidth <- df %>%
    group_by(drop_size, .[[5]],slice, generation) %>%
    summarize(mean_width = mean(width))
  
  #接下来再根据dropsize,dropradius和ganeration分组，对分组内每一个generation内的width求和，和存储为sum_width
  sum_width_df <- df_meanwidth %>%
    group_by(drop_size, .[[2]], generation) %>%
    summarize(sum_width = sum(mean_width))
  # 最终得到的max_sum_width DataFrame 包含不同 drop_size 和 drop_radius 组合的最大 sum_width 值
  max_sum_width <- sum_width_df %>%
    group_by(drop_size, .[[2]]) %>%
    summarize(max_sum_width = max(sum_width))
  
  names(max_sum_width)<-c("drop_size","varied_parameter","max_sum_width")
  max_sum_width$varied_parameter<-as.numeric(max_sum_width$varied_parameter)
  
  x_label =  "Relative Release Size"
  y_label <- gsub("_", " ", varied_parameter)
  #max_sum_width$drop_width <- 50*max_sum_width$drop_radius
  p_heatmaps <- ggplot(max_sum_width, aes(x = as.numeric(drop_size), y = varied_parameter)) + 
    geom_tile(aes(fill = max_sum_width)) + 
    scale_fill_gradientn(values = scales::rescale(c(0, 3.8, 7.9)),  # 将 (0, 3, 6) 映射到 [0, 1] 范围内
                         limits = c(0, 7.9),  
                         colors = c("#1c54a8", "#f6f5ee", "#ee2d2a"))+ 
    theme_classic() +
    #xlab("Drop Size") +
    #ylab(varied_parameter)+
    labs(x = x_label, y = y_label,title = " ",fill = "Max Width of Suppression") +theme(
      axis.text = element_text(size = 14),       # 坐标轴刻度
      axis.title = element_text(size = 15),      # 坐标轴标题
      legend.text = element_text(size = 14),     # 图例刻度
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 16)     # 图例标题
    )
  
  
  sum_width_df1 <- subset(sum_width_df, drop_size == 0.9)
  #找出每一个.[[2]]对应的最大sum_width所对应的generation
  # 找出每个 .[[2]] 对应的最大 sum_width 所对应的 generation
  sum_width_df2<-sum_width_df1 %>%
    group_by(.[[2]]) %>%
    slice(which.max(sum_width)) %>%
    select(generation, sum_width)
  names(sum_width_df2) <-c(varied_parameter,"generation","max")
  sum_width_df2 <- subset(sum_width_df2, max >0)
  sum_width_df2$drop_size <- rep(0.9,times =)
  
  see_drivefreq <- merge(df, sum_width_df2, by.x = c("drop_size", varied_parameter,"generation"),
                         by.y = c("drop_size", varied_parameter,"generation"))
  
  see_drivefreq<- select(see_drivefreq, 2,3,6,9)
  see_drivefreq$d3_freq <- as.numeric(see_drivefreq$d3_freq)
  #varied_parameter generation slice d3_freq
  names(see_drivefreq) <- c("A","B","C","D")
  
  df_meanfreq <- see_drivefreq %>%
    group_by(A,C) %>%
    summarize(mean_freq = mean(D))
  df_meanfreq$A <- as.factor(df_meanfreq$A)
  df_meanfreq$C <- as.numeric(df_meanfreq$C)
  #作图：
  pline1<- ggplot(df_meanfreq, aes(x = C, y = mean_freq,color = A)) +
    geom_line()+theme_classic()+  
    scale_color_brewer(palette = "Reds")+
    labs(x = "Slice",y = "Mean Drive3 Frequency",
         title ="Drop_Size = 0.9" )
  
  plots <- list()  
  plots[[1]] <- p_heatmaps  
  plots[[2]] <- pline1
  return(plots)
}

box_inheritance2 <- function(df) {
  
  group_stats <- df %>%
    group_by(Line, Gender) %>%
    summarise(
      n = n(),
      mean = mean(100*Inheritance_Rate, na.rm = TRUE),
      sd = sd(100*Inheritance_Rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      conf.low = mean - 1.96 * sd / sqrt(n),
      conf.high = mean + 1.96 * sd / sqrt(n)
    )
  
  group_stats <- df %>%
    group_by(Gender, Line) %>%
    summarise(
      mean = mean(100 * Inheritance_Rate, na.rm = TRUE),
      conf.low = t.test(100 * Inheritance_Rate)$conf.int[1],
      conf.high = t.test(100 * Inheritance_Rate)$conf.int[2],
      label = sprintf("%.1f%% ± %.1f", mean, (conf.high - conf.low) / 2),
      .groups = "drop"
    )
  group_stats <- group_stats %>%
    mutate(
      label = sprintf("%.1f%% ± %.1f", mean, (conf.high - conf.low) / 2),
      y_position = conf.high - 40  # 比原来的 conf.high 向下移 30 个单位
    )
  box_in <- ggplot(df, aes(x = as.factor(Gender), y = 100*Inheritance_Rate, fill = Gender, size = num_offspring, alpha = 0.6)) +
    geom_jitter(aes(fill = Gender), width = 0.2, shape = 21) +
    #geom_pointrange(
    #  data = group_stats,
     # aes(x = Gender, y = mean, ymin = conf.low, ymax = conf.high),
      #color = "grey15",
      #size = 0.5,
     # inherit.aes = FALSE
   # ) +
    scale_fill_manual(values = c("♂" = "#1c54a8", "♀" = "#ee2d2a")) +
    labs(
      x = NULL,  
      y = "Inheritance Rate (%)"
    ) +
    scale_size(range = c(1, 10)) +
    ylim(0, 101) +
    facet_wrap(~ Line, scales = "free_x") +  # 按 Line 分面
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      strip.text = element_text(size = 16) 
    )+
    geom_text(
      data = group_stats,
      aes(x = Gender, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 4
    )
  
  return(box_in)
}


####embryo heatmap panmictic
setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/publish_scale_nongl_data/publish_scale_emb_0402_bino_m8_lethal")
#emb_total1 <-  read.csv("emb_nonGL.csv",
                        stringsAsFactors=F, header=F)
#emb_total2 <-  read.csv("emb0_nongl.csv",
                        stringsAsFactors=F, header=F)
#emb_total <- rbind(emb_total2,emb_total1)


emb_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_emb/emb_pubscale.csv",
                        stringsAsFactors=F, header=F)
emb_total[] <- lapply(emb_total, as.numeric)
pemb <- heatmap(emb_total,"non_GL","Embryo_Resistance_Rate")
#+labs(title = "Embryo Resistance Rate")
pemb



emb_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_embGL/emb_GL_pubscale.csv",
                       stringsAsFactors=F, header=F)

emb_total[] <- lapply(emb_total, as.numeric)
pemb <- heatmap(emb_total,"GL","Embryo_Resistance_Rate")

pemb




setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/publish_scale_nongl_data/publish_scale_con_0402_bino_m8_lethal")
con_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_con/0824_pub_con.csv",
                       stringsAsFactors=F, header=F)

con_total[] <- lapply(con_total, as.numeric)
pcon <- heatmap(con_total,"non_GL","Conversion_Rate")
pcon



##画图：固定release size，
con_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_con/0824_pub_con.csv",
                       stringsAsFactors=F, header=F)
con_total[] <- lapply(con_total, as.numeric)
colnames_con<- con_total[1,]
con_total_values <- con_total[-1,]
#筛选release = 0.6,0.65,0.7
#conversion > 0.75(0.8, 0.85, 0.9, 0.95,1.00)
#提取conversion_rate, non-drive element，cas9，generation,
#提取第1、3、4、8列
#筛选第七列=0.6时的dataframe，记为df1
#提取df1的1、3、4、8列（每隔8列提取一次）
con_total_d06 <- con_total_values[con_total_values[[8]] == 0.6, ]
con_total_d065 <- con_total_values[con_total_values[[8]] == 0.65, ]
con_total_d07 <- con_total_values[con_total_values[[8]] == 0.7, ]
#按conversion分组，分别整合
#conversion+generation
#conversion+freq1
#conversion+freq2
line_con <- function(df,title){
  selected_columns_dr1 <- c(
    seq(85, ncol(df), by = 9)
  )
  df_dr1 <- df[, sort(selected_columns_dr1)]
  df_dr1$con <- df[,90]
  df_dr1[is.na(df_dr1)] <- 0.0
  colnames(df_dr1) <- as.character(1:42)
  df_long_dr1 <- df_dr1 %>%
    pivot_longer(
      cols = 1:41,                 
      names_to = "generation",
      values_to = "dr1_frequency"
    ) %>%
    mutate(generation = as.numeric(generation)) %>%
    rename(conversion = `42`) %>%  
    select(conversion, generation, dr1_frequency)
  
  selected_columns_dr2 <- c(
    seq(86, ncol(df), by = 9)
  )
  df_dr2 <- df[, sort(selected_columns_dr2)]
  df_dr2$con <- df[,90]
  df_dr2[is.na(df_dr2)] <- 0.0
  colnames(df_dr2) <- as.character(1:42)
  df_long_dr2 <- df_dr2 %>%
    pivot_longer(
      cols = 1:41,                 
      names_to = "generation",
      values_to = "dr2_frequency"
    ) %>%
    mutate(generation = as.numeric(generation)) %>%
    rename(conversion = `42`) %>%  
    select(conversion, generation, dr2_frequency)
  
  mean_dr1 <- df_long_dr1 %>% filter(conversion>0.75)  %>% 
    group_by( generation,conversion) %>%
    summarize(mean_freq1 = mean(as.numeric(dr1_frequency)))
  mean_dr2 <- df_long_dr2%>% filter(conversion>0.75)  %>% 
    group_by( generation,conversion) %>%
    summarize(mean_freq2 = mean(as.numeric(dr2_frequency)))
  p1 <- ggplot(mean_dr1, aes(
    x = as.numeric(generation),
    y = as.numeric(mean_freq1),
    color = as.numeric(conversion),
    group = conversion
  )) +
    geom_line(size = 1.2) +
    theme_classic() +
    scale_color_gradient(
      low = "#FFD4CC",   # #FF6F61 的浅色版本
      high = "#FF6F61"   # 基色本身（深）
    ) +
    labs(
      title = title,
      x = "Generation",
      y = "Non-drive Element Frequency",
      color = "Conversion Rate"
    ) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 16)
    )
  p2 <-   ggplot(mean_dr2, aes(
    x = as.numeric(generation),
    y = as.numeric(mean_freq2),
    color = as.numeric(conversion),   # 转换成数值，方便映射连续色
    group = conversion                # 保证不同 conversion 各自成线
  )) +
    geom_line(size = 1.2) +
    theme_classic() +
    scale_color_gradient(
      low = "#a6c8ff",   # 
      high = "#1c54a8"   # 
    ) +
    labs(
      title = title,
      x = 'Generation',
      y = "Cas9 Allele Frequency",
      color = "Conversion Rate"
    ) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 16)
    )
  
  return(ggarrange(p1,p2,ncol = 1, nrow = 2))

}
line_con(con_total_d06,"Relative Release Size = 0.60")
line_con(con_total_d065,"Relative Release Size = 0.65")
line_con(con_total_d07,"Relative Release Size = 0.70")


#提取一下con = 0.80的时候
#做一个population size的table
df <- con_total_d06 
selected_columns_pop <- c(
seq(84, ncol(df), by = 9)
)
df_pop <- df[, sort(selected_columns_pop)]
df_pop$con <- df[,90]
colnames(df_pop) <- as.character(1:42)
df_long_pop <- df_pop %>%
  pivot_longer(
    cols = 1:41,                 
    names_to = "generation",
    values_to = "pop_size"
  ) %>%
  mutate(generation = as.numeric(generation)) %>%
  rename(conversion = `42`) %>%  
  select(conversion, generation, pop_size)
df_long_pop <- filter(df_long_pop,df_long_pop$conversion == 0.8)
mean_pop <- df_long_pop %>% 
  group_by( generation,conversion) %>%
  summarize(mean_pop = mean(as.numeric(pop_size)))



selected_columns_dr1 <- c(
  seq(85, ncol(df), by = 9)
)
df_dr1 <- df[, sort(selected_columns_dr1)]
df_dr1$con <- df[,90]
df_dr1[is.na(df_dr1)] <- 0.0
colnames(df_dr1) <- as.character(1:42)
df_long_dr1 <- df_dr1 %>%
  pivot_longer(
    cols = 1:41,                 
    names_to = "generation",
    values_to = "dr1_frequency"
  ) %>%
  mutate(generation = as.numeric(generation)) %>%
  rename(conversion = `42`) %>%  
  select(conversion, generation, dr1_frequency)
df_long_dr1<- filter(df_long_dr1,df_long_dr1$conversion == 0.80)
mean_dr1 <- df_long_dr1 %>% 
  group_by( generation,conversion) %>%
  summarize(mean_dr1 = mean(as.numeric(dr1_frequency)))

ggplot(mean_pop, aes(
  x = as.numeric(generation),
  y = as.numeric(mean_pop))) +
  geom_line(size = 1.2) +
  theme_classic() +ggplot(mean_dr1, aes(
  x = as.numeric(generation),
  y = as.numeric(mean_dr1))) +
  geom_line(size = 1.2) +
  theme_classic() 







#con_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0410_gl_data/con_GL_pubscale.csv",
 #                      stringsAsFactors=F, header=F)
con_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_conGL/con_GL_pubscale.csv",
                       stringsAsFactors=F, header=F)
con_total[] <- lapply(con_total, as.numeric)
pcongl <- heatmap(con_total,"GL","Conversion_Rate")
pcongl



#setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0229_pan_som/som_data_sup1")
#setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/som_data_sup/")
#setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/publish_scale_som_0402")
som_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_som/som_nonGL.csv",
                       stringsAsFactors=F, header=F)
som_total[] <- lapply(som_total, as.numeric)
psom <- heatmap(som_total,"non_GL","Females_Heterozygote_Fitness")
psom


som_total <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/0824_pub_somGL/som_GL_pubscale.csv",
                       stringsAsFactors=F, header=F)
som_total[] <- lapply(som_total, as.numeric)
psom <- heatmap(som_total,"GL","Females_Heterozygote_Fitness")
psom





#setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/publish_scale_data/publishscale_data_embryo_0408")
df_emb <- read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubemb_0826/1d_embryo_pubscale_0826.csv",stringsAsFactors=F, header=F)
df_emb1 <- read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubemb_0826/1d_embryo_pubscale_0901.csv",stringsAsFactors=F, header=F)
df_emb_final <- rbind(df_emb,df_emb1)
p_emb_1d_heatmaps<- heatmap_onedim(df_emb_final, "Embryo_Resistance_Rate")[1]
p_emb_1d_meanpop<- heatmap_onedim_meanpop(df_emb_final, "Embryo_Resistance_Rate")
p_emb_1d_meanpop+xlab("Relative Release Size")
p_emb_1d_heatmaps

df_som_20reps<-read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubsom_0826/1d_somatic_pubscale_0826.csv",stringsAsFactors=F, header=F)
df_som_20reps1<-read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubsom_0826/1d_somatic_pubscale_0901.csv",stringsAsFactors=F, header=F)
df_som_20reps1 <- df_som_20reps1[-1, ]
df_som_20reps_final<- rbind(df_som_20reps,df_som_20reps1)
p_somatic_mult_f_1d<- heatmap_onedim(df_som_20reps_final, "Female_Heterozygote_Fitness")[1]
p_somatic_1d_meanpop<- heatmap_onedim_meanpop(df_som_20reps_final, "Female_Heterozygote_Fitness")
p_somatic_mult_f_1d
p_somatic_1d_meanpop+xlab("Relative Release Size")

df_con_20reps<-read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubcon_0826/1d_conversion_pubscale_0826.csv",stringsAsFactors=F, header=F)
df_con_20reps1<-read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_SPA/heatmaps_1d/1d_pubcon_0826/1d_conversion_pubscale_0901.csv",stringsAsFactors=F, header=F)
df_con_20reps1 <- df_con_20reps1[-1, ]
df_con_20reps_final<- rbind(df_con_20reps,df_con_20reps1)
p_con_mult_f_1d<- heatmap_onedim(df_con_20reps_final, "Conversion_Rate")[1]
p_con_1d_meanpop<- heatmap_onedim_meanpop(df_con_20reps_final, "Conversion_Rate")
p_con_mult_f_1d
p_con_1d_meanpop+xlab("Relative Release Size")




#比一下conversion = 0.1，0.5和1.0的
df_con_20reps_values <- df_con_20reps_final[-1,]
df_con_20reps_values<- df_con_20reps_values[,-12]
names(df_con_20reps_values) <- c("result","size_all","generation" , "drop_size","conversion_rates","slice" ,
                                 "d1_freq","d2_freq" ,"d3_freq" ,"pop_size",
                                 "decrease_proportion")

df_con_20reps_values$conversion_rates <- as.numeric(df_con_20reps_values$conversion_rates)
df_con_20reps_values<- filter(df_con_20reps_values,df_con_20reps_values$drop_size == 0.9)
df_con_20reps_values<- filter(df_con_20reps_values,df_con_20reps_values$generation>0)
df_con_20reps_01<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.1) %>%select(generation,size_all,conversion_rates)
df_con_20reps_05<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.5)%>%select(generation,size_all,conversion_rates)
df_con_20reps_10<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 1.0)%>%select(generation,size_all,conversion_rates)


df_con_20reps_01_meanpop <- df_con_20reps_01 %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = mean(as.numeric(size_all)))
df_con_20reps_05_meanpop <- df_con_20reps_05 %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = mean(as.numeric(size_all)))
df_con_20reps_10_meanpop <- df_con_20reps_10 %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = mean(as.numeric(size_all)))

linepop_all <- rbind(df_con_20reps_01_meanpop,
                     df_con_20reps_05_meanpop,df_con_20reps_10_meanpop)

linepop_all$conversion_rates <- as.factor(linepop_all$conversion_rates)

x_label<- "Generation"
y_label1<- "Average Population Size"
p_1dpop_line <- ggplot(linepop_all, aes(x = as.numeric(generation), y = as.numeric(mean_pop_size), color = conversion_rates)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "0.1" = "#FC9272",
    "0.5" = "#EF3B2C",
    "1" = "#67000D"
  ))       +labs(x = x_label, y = y_label1,color = "Conversion Rate") +theme(
    axis.text = element_text(size = 14),       # 坐标轴刻度
    axis.title = element_text(size = 15),      # 坐标轴标题
    legend.text = element_text(size = 14),     # 图例刻度
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16)     # 图例标题
  )

df_con_20reps_01f<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.1) %>%select(generation,conversion_rates,d3_freq)
df_con_20reps_05f<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.5)%>%select(generation,conversion_rates,d3_freq)
df_con_20reps_10f<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 1.0)%>%select(generation,conversion_rates,d3_freq)

df_con_20reps_01f_meanfreq <- df_con_20reps_01f %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))
df_con_20reps_05f_meanfreq <- df_con_20reps_05f %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))
df_con_20reps_10f_meanfreq <- df_con_20reps_10f %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))

linedf_all <- rbind(df_con_20reps_01f_meanfreq,
                    df_con_20reps_05f_meanfreq,df_con_20reps_10f_meanfreq)

linedf_all$conversion_rates <- as.factor(linedf_all$conversion_rates)
x_label<- "Generation"
y_label<- "Average Suppression Element Frequency"

p_1ddf_line <- ggplot(linedf_all, aes(x = as.numeric(generation), y = as.numeric(mean_d3f), color = conversion_rates)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "0.1" = "#FC9272",
    "0.5" = "#EF3B2C",
    "1" = "#67000D"
  ))       +labs(x = x_label, y = y_label,color = "Conversion Rate") +theme(
    axis.text = element_text(size = 14),       # 坐标轴刻度
    axis.title = element_text(size = 15),      # 坐标轴标题
    legend.text = element_text(size = 14),     # 图例刻度
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16)     # 图例标题
  )
library(ggpubr)
ggarrange(p_1dpop_line,p_1ddf_line,ncol = 2,nrow = 1)



df_con_20reps_01l<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.1) %>%select(generation,pop_size,conversion_rates,slice)
df_con_20reps_01l<-filter(df_con_20reps_01l,as.numeric(df_con_20reps_01l$slice) <6)
df_con_20reps_05l<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.5) %>%select(generation,pop_size,conversion_rates,slice)
df_con_20reps_05l<-filter(df_con_20reps_05l,as.numeric(df_con_20reps_05l$slice) <6)
df_con_20reps_10l<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 1.0) %>%select(generation,pop_size,conversion_rates,slice)
df_con_20reps_10l<-filter(df_con_20reps_10l,as.numeric(df_con_20reps_10l$slice) <6)


df_con_20reps_01_meanpopl <- df_con_20reps_01l %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = sum(as.numeric(pop_size))/20)
df_con_20reps_05_meanpopl <- df_con_20reps_05l %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = sum(as.numeric(pop_size))/20)
df_con_20reps_10_meanpopl <- df_con_20reps_10l %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_pop_size = sum(as.numeric(pop_size))/20)



linepop_all_left <- rbind(df_con_20reps_01_meanpopl,
                          df_con_20reps_05_meanpopl,df_con_20reps_10_meanpopl)

linepop_all_left$conversion_rates <- as.factor(linepop_all_left$conversion_rates)

x_label<- "Generation"
y_label1<- "Average Population Size"
p_1dpop_line_leftl <- ggplot(linepop_all_left, aes(x = as.numeric(generation), y = as.numeric(mean_pop_size), color = conversion_rates)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "0.1" = "#FC9272",
    "0.5" = "#EF3B2C",
    "1" = "#67000D"
  ))       +labs(x = x_label, y = y_label1,color = "Conversion Rate") +theme(
    axis.text = element_text(size = 14),       # 坐标轴刻度
    axis.title = element_text(size = 15),      # 坐标轴标题
    legend.text = element_text(size = 14),     # 图例刻度
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16)     # 图例标题
  )
p_1dpop_line_leftl




df_con_20reps_01fl<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.1) %>%select(generation,d3_freq,conversion_rates,slice)
df_con_20reps_01fl<-filter(df_con_20reps_01fl,as.numeric(df_con_20reps_01fl$slice) <6)
df_con_20reps_05fl<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 0.5) %>%select(generation,d3_freq,conversion_rates,slice)
df_con_20reps_05fl<-filter(df_con_20reps_05fl,as.numeric(df_con_20reps_05fl$slice) <6)
df_con_20reps_10fl<- filter(df_con_20reps_values,df_con_20reps_values$conversion_rates == 1.0) %>%select(generation,d3_freq,conversion_rates,slice)
df_con_20reps_10fl<-filter(df_con_20reps_10fl,as.numeric(df_con_20reps_10fl$slice) <6)

df_con_20reps_01f_meanfreql <- df_con_20reps_01fl %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))
df_con_20reps_05f_meanfreql <- df_con_20reps_05fl %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))
df_con_20reps_10f_meanfreql <- df_con_20reps_10fl %>%
  group_by( generation,conversion_rates) %>%
  summarize(mean_d3f = mean(as.numeric(d3_freq)))

linedf_all_left <- rbind(df_con_20reps_01f_meanfreql,
                         df_con_20reps_05f_meanfreql,df_con_20reps_10f_meanfreql)

linedf_all_left$conversion_rates <- as.factor(linedf_all_left$conversion_rates)
x_label<- "Generation"
y_label<- "Average Suppression Element Frequency"

p_1ddf_linel <- ggplot(linedf_all_left, aes(x = as.numeric(generation), y = as.numeric(mean_d3f), color = conversion_rates)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "0.1" = "#FC9272",
    "0.5" = "#EF3B2C",
    "1" = "#67000D"
  ))       +labs(x = x_label, y = y_label,color = "Conversion Rate") +theme(
    axis.text = element_text(size = 14),       # 坐标轴刻度
    axis.title = element_text(size = 15),      # 坐标轴标题
    legend.text = element_text(size = 14),     # 图例刻度
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16)     # 图例标题
  )
p_1ddf_linel
library(ggpubr)
ggarrange(p_1dpop_line,p_1ddf_line,ncol = 2,nrow = 1)
ggarrange(p_1dpop_line_leftl,p_1ddf_linel,ncol = 2,nrow = 1)








##################################################box plots for conversion rate########################################################
two_element_conversion_test_mod <-  read.csv("/Users/alexgjl/Desktop/gene_drive/果蝇daisychain/Conversion_Data.csv",
                                             stringsAsFactors=F, header=T)

two_element_conversion_test_sup <-  read.csv("/Users/alexgjl/Desktop/gene_drive/果蝇daisychain/conversion_bstl.csv",
                                             stringsAsFactors=F, header=T)
three_element_conversion_test_sup <-  read.csv("/Users/alexgjl/Desktop/gene_drive/果蝇daisychain/conversion_bbstl.csv",
                                               stringsAsFactors=F, header=T)
three_element_conversion_test_mod <-  read.csv("/Users/alexgjl/Desktop/gene_drive/果蝇daisychain/3_element_inheritance_test_ahd.csv",
                                               stringsAsFactors=F, header=T)


#1. arrrange 2 element tests
two_element_conversion_test_sup$Line <- rep("HSDstlU4",times = )

two_element_conversion_test_sup$Gender <- c("♂","♂","♂","♂","♂",
                                            "♂","♂","♂",
                                            "♀","♀","♀","♀","♀",
                                            "♀","♀","♀","♀")


two_element_conversion_test_sup$num_offspring <- rowSums(two_element_conversion_test_sup[, c("RG", "R", "G", "Wt")])
two_element_conversion_test_sup$Inheritance_Rate <- two_element_conversion_test_sup$Inh.stl
two_element_conversion_test_sup<- select(two_element_conversion_test_sup,Group,Gender,num_offspring,Inheritance_Rate,Line)
two_element_conversion_test_mod<- filter(two_element_conversion_test_mod,two_element_conversion_test_mod$Line != "HSDygU4")
two_element_conversion_test<- rbind(two_element_conversion_test_mod,two_element_conversion_test_sup)



two_element_conversion_test$Line[which(two_element_conversion_test$Line == "gRNAs-B")] <- "DHMDhCas9"


con_ahd<- filter(two_element_conversion_test,two_element_conversion_test$Line == "AHDr352V2")
con_hsd<-filter(two_element_conversion_test,two_element_conversion_test$Line == "HSDstlU4")
con_grnas<- filter(two_element_conversion_test,two_element_conversion_test$Line == "DHMDhCas9")

con_2element_all <- rbind(con_ahd,con_hsd,con_grnas)


con_2element_all$Gender[which(con_2element_all$Gender == "Male")] <- "♂"
con_2element_all$Gender[which(con_2element_all$Gender == "Female")] <- "♀"

con_2element_all$Line[which(con_2element_all$Line == "DHMDhCas9")] <- "Cas9 Element"
con_2element_all$Line[which(con_2element_all$Line == "AHDr352V2")] <- "Modification Element"
con_2element_all$Line[which(con_2element_all$Line == "HSDstlU4")] <- "Suppression Element"


con_2element_all$Line<- factor(con_2element_all$Line,levels =  c("Cas9 Element","Modification Element",
                                                                 "Suppression Element"))

con_2element_all$Gender<- factor(con_2element_all$Gender,levels =  c("♂","♀"))


p_all <- box_inheritance2(con_2element_all)
p_all+
  patchwork::plot_annotation(
    caption = "Double heterozygote parent",
    theme = theme(
      plot.caption = element_text(size = 16, hjust = 0.5,  margin = margin(t = 10))
    )
  )



mean_hsd <- con_hsd %>%
  group_by(Gender) %>%
  summarise(mean = mean(Inheritance_Rate, na.rm = TRUE))
mean_ahd <- con_ahd %>%
  group_by(Gender) %>%
  summarise(mean = mean(Inheritance_Rate, na.rm = TRUE))
mean_c9 <- con_grnas %>%
  group_by(Gender) %>%
  summarise(mean = mean(Inheritance_Rate, na.rm = TRUE))


#triple-elements carrier

three_element_conversion_test_sup$Group <- c(1:27)
three_element_conversion_test_mod$Group <- c(1:35)

three_element_conversion_test_sup$Gender <- c(rep("♀",times = 11),rep("♂",times = 16))
three_element_conversion_test_sup$num_offspring <- rowSums(three_element_conversion_test_sup[, c("RbRG",  "RbR","RbG","Rbwt","WtbRG","WtbR","WtbG","Wtbwt" )])   

three_element_conversion_test_mod$Gender<- three_element_conversion_test_mod$Sex
three_element_conversion_test_mod$num_offspring <- three_element_conversion_test_mod$num_offs



three_sup_c9i <- select(three_element_conversion_test_sup,Gender,num_offspring,Group,Inh.Cas9)
three_sup_c9i$Inheritance_Rate<- three_sup_c9i$Inh.Cas9
three_sup_c9i$Line <- rep("DHMDhCas9",times = )
three_sup_c9i<- select(three_sup_c9i,Group,Line, Inheritance_Rate,Gender,num_offspring)

three_sup_hsdi<- select(three_element_conversion_test_sup,Gender,num_offspring,Group,Inh.stl)
three_sup_hsdi$Inheritance_Rate<- three_sup_hsdi$Inh.stl
three_sup_hsdi$Line <- rep("HSDstlU4",times = )
three_sup_hsdi<- select(three_sup_hsdi,Group,Line, Inheritance_Rate,Gender,num_offspring)

three_sup_grnai <- select(three_element_conversion_test_sup,Gender,num_offspring,Group,Inh.nondrive)
three_sup_grnai$Inheritance_Rate<- three_sup_grnai$Inh.nondrive
three_sup_grnai$Line <- rep("DN38C2",times = )
three_sup_grnai<- select(three_sup_grnai,Group,Line, Inheritance_Rate,Gender,num_offspring)


three_sup_all <- rbind(three_sup_c9i,three_sup_hsdi,three_sup_grnai)
three_sup_all$Line[which(three_sup_all$Line == "DN38C2")]<- "Non-drive Element"
three_sup_all$Line[which(three_sup_all$Line == "DHMDhCas9")]<- "Cas9 Element"
three_sup_all$Line[which(three_sup_all$Line == "HSDstlU4")]<- "Suppression Element"
three_sup_all$Line<- factor(three_sup_all$Line,levels =  c("Non-drive Element",
                                                           "Cas9 Element",
                                                           "Suppression Element"))


three_sup_all$Gender[which(three_sup_all$Gender == "Male")] <- "♂"
three_sup_all$Gender[which(three_sup_all$Gender == "Female")] <- "♀"

three_sup_all$Gender<- factor(three_sup_all$Gender,levels =  c("♂","♀"))
p_all2 <- box_inheritance2(three_sup_all)
p_all2##############
####################

three_mod_c9i <- select(three_element_conversion_test_mod,Gender,num_offspring,Group,Inh_Cas9)
three_mod_c9i$Inheritance_Rate<- three_mod_c9i$Inh_Cas9
three_mod_c9i$Line <- rep("DHMDhCas9",times = )
three_mod_c9i<- select(three_mod_c9i,Group,Line, Inheritance_Rate,Gender,num_offspring)

three_mod_ahdi<- select(three_element_conversion_test_mod,Gender,num_offspring,Group,Inh_Cargo)
three_mod_ahdi$Inheritance_Rate<- three_mod_ahdi$Inh_Cargo
three_mod_ahdi$Line <- rep("AHDr352V2",times = )
three_mod_ahdi<- select(three_mod_ahdi,Group,Line, Inheritance_Rate,Gender,num_offspring)

three_mod_grnai <- select(three_element_conversion_test_mod,Gender,num_offspring,Group,Inh_nonD)
three_mod_grnai$Inheritance_Rate<- three_mod_grnai$Inh_nonD
three_mod_grnai$Line <- rep("DN38C2",times = )
three_mod_grnai<- select(three_mod_grnai,Group,Line, Inheritance_Rate,Gender,num_offspring)


three_mod_all <- rbind(three_mod_c9i,three_mod_ahdi,three_mod_grnai)
three_mod_all$Line[which(three_mod_all$Line == "DN38C2")]<- "Non-drive Element"
three_mod_all$Line[which(three_mod_all$Line == "DHMDhCas9")]<- "Cas9 Element"
three_mod_all$Line[which(three_mod_all$Line == "AHDr352V2")]<- "Modification Element"
three_mod_all$Line<- factor(three_mod_all$Line,levels =  c("Non-drive Element",
                                                           "Cas9 Element",
                                                           "Modification Element"))


three_mod_all$Gender[which(three_mod_all$Gender == "Male")] <- "♂"
three_mod_all$Gender[which(three_mod_all$Gender == "Female")] <- "♀"
three_mod_all$Gender<- factor(three_mod_all$Gender,levels =  c("♂","♀"))

p_all3 <- box_inheritance2(three_mod_all)

p_all+
  patchwork::plot_annotation(
    caption = "Double heterozygote parent",
    theme = theme(
      plot.caption = element_text(size = 16, hjust = 0.5,  margin = margin(t = 10))
    )
  )

p_all2+
  patchwork::plot_annotation(
    caption = "Triple heterozygote parent",
    theme = theme(
      plot.caption = element_text(size = 16, hjust = 0.5,  margin = margin(t = 10))
    )
  )
p_all3+
  
  patchwork::plot_annotation(
    caption = "Triple heterozygote parent",
    theme = theme(
      plot.caption = element_text(size = 16, hjust = 0.5,  margin = margin(t = 10))
    )
  )



#fisher's exact test
mean_mod3 <- three_mod_all %>%
  group_by(Gender, Line) %>%
  summarise(
    mean = mean(Inheritance_Rate, na.rm = TRUE),
    total_offspring = sum(num_offspring),
    success = sum(round(Inheritance_Rate * num_offspring)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = binom.test(success, total_offspring, p = 0.5, alternative = "greater")$p.value
  )

mean_sup3 <- three_sup_all %>%
  group_by(Gender, Line) %>%
  summarise(
    mean = mean(Inheritance_Rate, na.rm = TRUE),
    total_offspring = sum(num_offspring),
    success = sum(round(Inheritance_Rate * num_offspring)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = binom.test(success, total_offspring, p = 0.5, alternative = "greater")$p.value
  )







#做一个illustrative figure，展示驱动基因频率的变化
setwd("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/Daisy_PAN/publish_scale_nongl_data/publish_scale_con_0402_bino_m8_lethal")
con_total <-  read.csv("con_nonGL.csv",
                       stringsAsFactors=F, header=F)

df <- con_total
df1 <- df[-1,]
##8是dropsize，9是embres
parameter_varied<- select(df1,V9)
drop_varied<- select(df1,V8)

#parameter_varied$varied_parameter <- parameter_varied$V9
drop_varied$drop_size <- drop_varied$V8

df_parameters <- cbind(drop_varied,parameter_varied)
df_parameters <- select(df_parameters,drop_size,V9)#选择修改的两个参数

#df_parameters$drop_size <- as.numeric(df_parameters$drop_size)
# df_parameters$V9<- as.numeric( df_parameters$V9)
#next：提取pop_size
interval <- 9#每9列提取一次，永远不变



selected_dr1 <- seq(4, ncol(df1), by = interval)
selected_dr2 <- seq(5, ncol(df1), by = interval)
selected_dr3 <- seq(6, ncol(df1), by = interval)
# 根据列索引选取dataframe的几列，生成只含有popsize的dataframe：df_pop
df_dr1 <- df1[, selected_dr1]
df_dr2 <- df1[, selected_dr2]
df_dr3 <- df1[, selected_dr3]


df_dr1 <- as.data.frame(lapply(df_dr1, as.numeric))
df_dr2 <- as.data.frame(lapply(df_dr2, as.numeric))    
df_dr3 <- as.data.frame(lapply(df_dr3, as.numeric))    

df_dr1 <- cbind(df_parameters,df_dr1)
df_dr2 <- cbind(df_parameters,df_dr2)
df_dr3 <- cbind(df_parameters,df_dr3)


df_dr1 <- filter(df_dr1,df_dr1$drop_size == 0.5)
df_dr2 <- filter(df_dr2,df_dr2$drop_size == 0.5)
df_dr3 <- filter(df_dr3,df_dr3$drop_size == 0.5)

df_dr1 <- filter(df_dr1,as.numeric(df_dr1$V9) == 1.0)
df_dr2 <- filter(df_dr2,as.numeric(df_dr2$V9) == 1.0)
df_dr3 <- filter(df_dr3,as.numeric(df_dr3$V9) == 1.0)

dr1_values <- df_dr1[,-1]
dr1_values <-dr1_values[,-1]
dr2_values <- df_dr2[,-1]
dr2_values <-dr2_values[,-1]
dr3_values <- df_dr3[,-1]
dr3_values <-dr3_values[,-1]

dr1_values_l <- as.data.frame(t(dr1_values))
dr2_values_l <- as.data.frame(t(dr2_values))
dr3_values_l <- as.data.frame(t(dr3_values))

dr1_values_l$mean <- rowMeans(dr1_values_l)
dr2_values_l$mean <- rowMeans(dr2_values_l)
dr3_values_l$mean <- rowMeans(dr3_values_l)

dr1_values_l$generation <- c(1:50)
dr2_values_l$generation <- c(1:50)
dr3_values_l$generation <- c(1:50)

dr1_values_l$gene <- rep("Non-drive element",times = )
dr2_values_l$gene <- rep("Cas9",times = )
dr3_values_l$gene <- rep("Cargo_gene",times = )


dr1_values_l <- select(dr1_values_l,gene,mean,generation)
dr2_values_l <- select(dr2_values_l,gene,mean,generation)
dr3_values_l <- select(dr3_values_l,gene,mean,generation)    

dr_freq_values <- rbind(dr1_values_l,dr2_values_l,dr3_values_l)   

dr_freq_values$gene[dr_freq_values$gene == "Cargo_gene"] <- "Suppression elements"



control <-  read.csv("/Users/alexgjl/Desktop/gene_drive/paper\ of\ daisy\ chain/standard_data/standard_drive_control.csv",
                     stringsAsFactors=F, header=F)    
control <- control[,-1]
con_selected <- control[, seq(3, ncol(control), by = 4)]

con_selected[con_selected == "NAN"] <- 1.0
con_selectedl <- as.data.frame(t(con_selected))

con_selectedl[] <- lapply(con_selectedl, function(x) as.numeric(as.character(x)))

df_fill <- as.data.frame(matrix(1, nrow = 35, ncol = 11))
con_selectedl$mean <- rowMeans(con_selectedl)

colnames(df_fill)<- colnames(con_selectedl)
con_selectedl<-rbind(con_selectedl,df_fill)

con_selectedl$generation <- c(10:59)
con_selectedl<- filter(con_selectedl,con_selectedl$generation<51)
con_selectedl$gene<- rep("Effector element",times = )
con_selectedl<- select(con_selectedl,gene,generation,mean)

dr_freq_values_all <- rbind(dr_freq_values,con_selectedl)


dr_freq_values_all$gene[which(dr_freq_values_all$gene == "Suppression elements")] <- "Daisy-chain suppression element"
dr_freq_values_all$gene[which(dr_freq_values_all$gene == "Effector element")] <- "Homing suppression drive"
dr_freq_values_all$gene[which(dr_freq_values_all$gene == "Cas9")] <- "Daisy-chain Cas9 element"
dr_freq_values_all$gene[which(dr_freq_values_all$gene == "Non-drive element")] <- "Daisy-chain non-drive element"

dr_freq_values_all$gene <- factor(
  dr_freq_values_all$gene,
  levels = c(
    "Homing suppression drive",
    "Daisy-chain suppression element",
    "Daisy-chain Cas9 element",
    "Daisy-chain non-drive element"
  )
)

# 绘图
p_illustrate <- ggplot(dr_freq_values_all, aes(
  x = as.numeric(generation),
  y = as.numeric(mean),
  color = gene
)) +
  geom_line(size = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "Homing suppression drive" = "#A50F15",
    "Daisy-chain suppression element" = "#FF9933",
    "Daisy-chain Cas9 element" = "#1c54a8",
    "Daisy-chain non-drive element" = "#FF6F61"
  )) +
  xlim(10, 50) +
  xlab("Generation") +
  ylab("Allele frequency") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 16)
  )
p_illustrate








rm(list = ls())
data_fitness <-  read.csv("/Users/alexgjl/Desktop/gene_drive/果蝇daisychain/bbahd_viability.csv",
                        stringsAsFactors=F, header=T)
library(ggplot2)
library(dplyr)

data_fitness$viability <- data_fitness$num_adults/data_fitness$num_eggs

plot_eggs_by_parent <- function(data) {
  # 计算统计量
  stats <- data %>%
    group_by(Group) %>%
    summarise(
      mean = mean(num_eggs, na.rm = TRUE),
      sd = sd(num_eggs, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      label = ifelse(is.na(se), 
                     paste0("Mean = ", round(mean, 1)),   # se为NA时只显示均值
                     paste0("Mean = ", round(mean, 1), " ± ", round(se, 1))
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(label))  # 只保留有 label 的行
  
  # 颜色
  group_colors <- c("Triple heterozygote female" = "#A50F15", "w1118 female" = "#FF9933")
  
  # 绘图
  p <- ggplot(data, aes(x = Group, y = num_eggs, fill = Group)) +
    geom_jitter(
      shape = 21,  # 填充圆点
      color = "black",  # 黑边
      stroke = 0.4,
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
      alpha = 0.7,
      size = 7
    ) +
    geom_segment(
      data = stats,
      aes(
        x = 1 + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        xend = 1 + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        y = mean,
        yend = mean
      ),
      color = "black",
      size = 1.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = stats,
      aes(
        x =1 + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        y =ifelse(Group == "Triple heterozygote female", 45, 55),
        label = label
      ),
    color = "black",
      size = 3.5,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = group_colors) +
    labs(x = "", y = "Number of Eggs", fill = "Group") +
    theme_minimal(base_size = 14)
  
  return(p)
}
#plot_adults_by_parent <- function(data) {
  # 计算统计量
  stats <- data %>%
    group_by(Parent, Group) %>%
    summarise(
      mean = mean(num_adults, na.rm = TRUE),
      sd = sd(num_adults, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      label = ifelse(is.na(se), 
                     paste0("Mean = ", round(mean, 1)),   # se为NA时只显示均值
                     paste0("Mean = ", round(mean, 1), " ± ", round(se, 1))
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(label))  # 只保留有 label 的行
  
  # 颜色
  group_colors <- c("Triple heterozygote female" = "#A50F15", "w1118 female" = "#FF9933")
  
  # 绘图
  p <- ggplot(data, aes(x = factor(Parent), y = num_adults, fill = Group)) +
    geom_jitter(
      shape = 21,  # 填充圆点
      color = "black",  # 黑边
      stroke = 0.4,
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
      alpha = 0.7,
      size = 7
    ) +
    geom_segment(
      data = stats,
      aes(
        x = as.numeric(factor(Parent)) + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        xend = as.numeric(factor(Parent)) + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        y = mean,
        yend = mean
      ),
      color = "black",
      size = 1.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = stats,
      aes(
        x = as.numeric(factor(Parent)) + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        y = ifelse(Group == "Triple heterozygote female", 45, 55),
        label = label
      ),
      color = "black",
      size = 3.5,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = group_colors) +
    labs(x = "Parent", y = "Number of Adults", fill = "Group") +
    theme_minimal(base_size = 14)
  
  return(p)
}
plot_viability_by_parent <- function(data) {
  # 计算统计量
  stats <- data %>%
    group_by( Group) %>%
    summarise(
      mean = mean(viability, na.rm = TRUE),
      sd = sd(viability, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      label = ifelse(is.na(se), 
                     paste0("Mean = ", round(mean, 1)),   # se为NA时只显示均值
                     paste0("Mean = ", round(mean, 1), " ± ", round(se, 1))
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(label))  # 
  
  # 颜色
  group_colors <- c("Triple heterozygote female" = "#A50F15", "w1118 female" = "#FF9933")
  
  # 绘图
  p <- ggplot(data, aes(x = Group, y = viability, fill = Group)) +
    geom_jitter(
      shape = 21,  # 填充圆点
      color = "black",  # 黑边
      stroke = 0.4,
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
      alpha = 0.7,
      size = 7
    ) +
    geom_segment(
      data = stats,
      aes(
        x = 1 + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        xend = 1 + ifelse(Group == "Triple heterozygote female", -0.15, 0.15),
        y = mean,
        yend = mean
      ),
      color = "black",
      size = 1.5,
      inherit.aes = FALSE
    )  +
    scale_fill_manual(values = group_colors) +
    labs(x = "", y = "Egg to Adult Viability", fill = "Group") +
    theme_minimal(base_size = 14)
  return(p)
}

data_fitness$Group[which(data_fitness$Group == "EXP")] <- "Triple heterozygote female"
data_fitness$Group[which(data_fitness$Group == "CON")] <- "w1118 female"

stats <- data_fitness %>%
  group_by( Group) %>%
  summarise(
    mean = mean(num_eggs, na.rm = TRUE),
    sd = sd(num_eggs, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    label = ifelse(is.na(se), 
                   paste0("Mean = ", round(mean, 1)),   # se为NA时只显示均值
                   paste0("Mean = ", round(mean, 1), " ± ", round(se, 1))
    ),
    .groups = "drop"
  ) 

plot_eggs_by_parent(data_fitness)+theme_classic()  +theme(
  axis.text = element_text(size = 16),       # 坐标轴刻度
  axis.title = element_text(size = 17),      # 坐标轴标题
  legend.text = element_text(size = 16),     # 图例刻度
  legend.title = element_text(size = 18),
  plot.title = element_text(size = 18)     # 图例标题
)

data_fitness <- filter(data_fitness ,data_fitness $viability<2)
plot_viability_by_parent(data_fitness)+theme_classic()+theme(
  axis.text = element_text(size = 16),       # 坐标轴刻度
  axis.title = element_text(size = 17),      # 坐标轴标题
  legend.text = element_text(size = 16),     # 图例刻度
  legend.title = element_text(size = 18),
  plot.title = element_text(size = 18)     # 图例标题
)




mean_eggs <- data_fitness %>%
  group_by(Group) %>%
  summarise(
    mean = mean(num_eggs, na.rm = TRUE))
    
 mean_adults <- data_fitness %>%
      group_by(Group) %>%
      summarise(
        mean = mean(num_adults, na.rm = TRUE))
 
 #significant test
 eggs_triple <- data_fitness$num_eggs[data_fitness$Group == "Triple heterozygote female"]
 eggs_w1118  <- data_fitness$num_eggs[data_fitness$Group == "w1118 female"]
 
 adults_triple <- data_fitness$num_adults[data_fitness$Group == "Triple heterozygote female"]
 adults_w1118  <- data_fitness$num_adults[data_fitness$Group == "w1118 female"]
 
 # 1.normality test
 shapiro.test(eggs_triple) 
 shapiro.test(eggs_w1118)
 
 shapiro.test(adults_triple)
 shapiro.test(adults_w1118)
 
 
 #Mann–Whitney U test
 wilcox.test(eggs_triple, eggs_w1118)
 wilcox.test(adults_triple, adults_w1118)
 mean_eggs
 mean_adults
 
 
 
 
 #################### #################### #################### #################### #################### #################### #################### #################### #################### #################### 
 #################### #################### #################### #################### #################### #################### #################### 
 #################### #################### #################### #################### #################### #################### ####################
 
 
 library(ggplot2)
 library(dplyr)
 cage_data <-  read.csv("/Users/alexgjl/Desktop/gene_drive/experiments/cage_data.csv",stringsAsFactors=F, header=T)
 
 # 自定义颜色：红、黄、蓝
 group_colors <- c("#A50F15", "#f6e09d", "#4d97cd")
 
 cage_ahd<- select(cage_data,Freq_AHD,Group,Generation)
 cage_cas9<- select(cage_data,Freq_Cas9,Group,Generation)
 cage_nondrive<- select(cage_data,Freq_nonDrive,Group,Generation)
 
 p_ahd <- ggplot(cage_ahd, aes(x = Generation, y = Freq_AHD, color = Group, group = Group)) +
   geom_line(size = 1.2) +
   geom_point(size = 2) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_minimal(base_size = 14)+theme_classic()+ylim(0,1)
 p_ahd
 
 p_cas9 <- ggplot(cage_cas9, aes(x = Generation, y = Freq_Cas9, color = Group, group = Group)) +
   geom_line(size = 1.2) +
   geom_point(size = 2) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_minimal(base_size = 14)+theme_classic()+ylim(0,1)
 
 
 p_nonDrive <- ggplot(cage_nondrive, aes(x = Generation, y = Freq_nonDrive, color = Group, group = Group)) +
   geom_line(size = 1.2) +
   geom_point(size = 2) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_minimal(base_size = 14)+theme_classic()+ylim(0,1)
 p_nonDrive
  library(ggpubr)
 ggarrange(p_cas9,p_ahd, labels = c("Cas9 Frequency","Modification Element Frequency"),ncol =2, nrow = 1)
 
 
 
 
 
 
 ########################################arrange slim data####################################################################
 cage_mod   <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/cage_simulation_pan/cage_modification/cage_data_mod0901.csv",
                         stringsAsFactors=F, header=F)
 
 cage_mod <- cage_mod[-c(12, 13,14), ]
 cage_mod$rep <- c(1:100)
 interval <- 10
 # 使用 seq() 函数生成要选取的列索引
 selected_columns_freq1 <- seq(4, ncol(cage_mod), by = interval)
 df_freq1_1 <- cage_mod[, selected_columns_freq1]
 df_freq1_1<- as.data.frame(lapply(df_freq1_1, as.numeric))
 colnames(df_freq1_1)<-c(1:15)
 
 selected_columns_freq2 <- seq(5, ncol(cage_mod), by = interval)
 df_freq2_1 <- cage_mod[, selected_columns_freq2]
 df_freq2_1<- as.data.frame(lapply(df_freq2_1, as.numeric))
 colnames(df_freq2_1)<-c(1:15)
 
 selected_columns_freq3 <- seq(6, ncol(cage_mod), by = interval)
 df_freq3_1 <- cage_mod[, selected_columns_freq3]
 df_freq3_1<- as.data.frame(lapply(df_freq3_1, as.numeric))
 colnames(df_freq3_1)<-c(1:15)

 
 df_freq1_1$rep <- as.factor(cage_mod$rep)
 df_freq2_1$rep <- as.factor(cage_mod$rep)
 df_freq3_1$rep <- as.factor(cage_mod$rep)

 
 #switch to long dataframe
 library(tidyr)
 df_freq1_long <- pivot_longer(df_freq1_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")
 df_freq2_long <- pivot_longer(df_freq2_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")
 df_freq3_long <- pivot_longer(df_freq3_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")

 
 df_avgfreq1 <- df_freq1_long %>%
   group_by(generation) %>%
   summarize(avg_freq1 = mean(cargo_freq))
 df_avgfreq2 <- df_freq2_long %>%
   group_by(generation) %>%
   summarize(avg_freq2 = mean(cargo_freq))
 df_avgfreq3 <- df_freq3_long %>%
   group_by(generation) %>%
   summarize(avg_freq3 = mean(cargo_freq))

 
 df_freq1_long$Group <- rep("Exp",times = )
 df_freq2_long$Group <- rep("Exp",times = )
 df_freq3_long$Group <- rep("Exp",times = )
 colnames(df_avgfreq1)[colnames(df_avgfreq1) == "avg_freq1"] <- "cargo_freq"
 colnames(df_avgfreq2)[colnames(df_avgfreq2) == "avg_freq2"] <- "cargo_freq"
 colnames(df_avgfreq3)[colnames(df_avgfreq3) == "avg_freq3"] <- "cargo_freq"


 
 
 df_avgfreq1$rep <- rep(0,times = )
 df_avgfreq1$Group <- rep("Avg",times = )
 df_freq1_long1 <- rbind(df_freq1_long,df_avgfreq1)
 
 df_avgfreq2$rep <- rep(0,times = )
 df_avgfreq2$Group <- rep("Avg",times = )
 df_freq2_long1 <- rbind(df_freq2_long,df_avgfreq2)
 
 df_avgfreq3$rep <- rep(0,times = )
 df_avgfreq3$Group <- rep("Avg",times = )
 df_freq3_long1 <- rbind(df_freq3_long,df_avgfreq3)
 
 

 

 df_freq1_long1$generation <- as.numeric(df_freq1_long1$generation)
 df_freq2_long1$generation <- as.numeric(df_freq2_long1$generation)
 df_freq3_long1$generation <- as.numeric(df_freq3_long1$generation)
 df_freq1_long1<-filter(df_freq1_long1,df_freq1_long1$generation <6) 
 df_freq2_long1<-filter(df_freq2_long1,df_freq2_long1$generation <6)
 df_freq3_long1<-filter(df_freq3_long1,df_freq3_long1$generation <6)
 
 df_freq1_long1$generation[which(df_freq1_long1$generation == 1)] <- 0
 df_freq1_long1$generation[which(df_freq1_long1$generation == 2)] <- 1
 df_freq1_long1$generation[which(df_freq1_long1$generation == 3)] <- 2
 df_freq1_long1$generation[which(df_freq1_long1$generation == 4)] <- 3
 df_freq1_long1$generation[which(df_freq1_long1$generation == 5)] <- 4 
 
 df_freq2_long1$generation[which(df_freq2_long1$generation == 1)] <- 0
 df_freq2_long1$generation[which(df_freq2_long1$generation == 2)] <- 1
 df_freq2_long1$generation[which(df_freq2_long1$generation == 3)] <- 2
 df_freq2_long1$generation[which(df_freq2_long1$generation == 4)] <- 3
 df_freq2_long1$generation[which(df_freq2_long1$generation == 5)] <- 4
 
 df_freq3_long1$generation[which(df_freq3_long1$generation == 1)] <- 0
 df_freq3_long1$generation[which(df_freq3_long1$generation == 2)] <- 1
 df_freq3_long1$generation[which(df_freq3_long1$generation == 3)] <- 2
 df_freq3_long1$generation[which(df_freq3_long1$generation == 4)] <- 3
 df_freq3_long1$generation[which(df_freq3_long1$generation == 5)] <- 4
 
 ####################################add another group without fitness multiplier########################################
 ####################################################################################################################################
 
 cage_mod_control   <-  read.csv("/Users/alexgjl/Desktop/gene_drive/modeling/Daisy_drives/cage_simulation_pan/cage_modification/cage_mod_control0901.csv",
                         stringsAsFactors=F, header=F)
 
 cage_mod_control <- cage_mod_control[-c(12, 13,14), ]
 cage_mod_control$rep <- c(1:100)
 interval <- 10
 # 使用 seq() 函数生成要选取的列索引
 selected_columns_freq1 <- seq(4, ncol(cage_mod_control), by = interval)
 df_freq1_1 <- cage_mod_control[, selected_columns_freq1]
 df_freq1_1<- as.data.frame(lapply(df_freq1_1, as.numeric))
 colnames(df_freq1_1)<-c(1:15)
 
 selected_columns_freq2 <- seq(5, ncol(cage_mod_control), by = interval)
 df_freq2_1 <- cage_mod_control[, selected_columns_freq2]
 df_freq2_1<- as.data.frame(lapply(df_freq2_1, as.numeric))
 colnames(df_freq2_1)<-c(1:15)
 
 selected_columns_freq3 <- seq(6, ncol(cage_mod_control), by = interval)
 df_freq3_1 <- cage_mod_control[, selected_columns_freq3]
 df_freq3_1<- as.data.frame(lapply(df_freq3_1, as.numeric))
 colnames(df_freq3_1)<-c(1:15)
 
 
 df_freq1_1$rep <- as.factor(cage_mod$rep)
 df_freq2_1$rep <- as.factor(cage_mod$rep)
 df_freq3_1$rep <- as.factor(cage_mod$rep)
 
 
 #switch to long dataframe
 library(tidyr)
 df_freq1_long_c <- pivot_longer(df_freq1_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")
 df_freq2_long_c <- pivot_longer(df_freq2_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")
 df_freq3_long_c <- pivot_longer(df_freq3_1, cols = c(1:15), names_to = "generation", values_to = "cargo_freq")
 
 
 df_avgfreq1c <- df_freq1_long_c %>%
   group_by(generation) %>%
   summarize(avg_freq1 = mean(cargo_freq))
 df_avgfreq2c <- df_freq2_long_c %>%
   group_by(generation) %>%
   summarize(avg_freq2 = mean(cargo_freq))
 df_avgfreq3c <- df_freq3_long_c %>%
   group_by(generation) %>%
   summarize(avg_freq3 = mean(cargo_freq))
 
 
 df_freq1_long_c$Group <- rep("Exp",times = )
 df_freq2_long_c$Group <- rep("Exp",times = )
 df_freq3_long_c$Group <- rep("Exp",times = )
 colnames(df_avgfreq1c)[colnames(df_avgfreq1c) == "avg_freq1"] <- "cargo_freq"
 colnames(df_avgfreq2c)[colnames(df_avgfreq2c) == "avg_freq2"] <- "cargo_freq"
 colnames(df_avgfreq3c)[colnames(df_avgfreq3c) == "avg_freq3"] <- "cargo_freq"
 
 
 
 
 df_avgfreq1c$rep <- rep(0,times = )
 df_avgfreq1c$Group <- rep("Avg",times = )
 df_freq1_long1_c <- rbind(df_freq1_long_c,df_avgfreq1c)
 
 df_avgfreq2c$rep <- rep(0,times = )
 df_avgfreq2c$Group <- rep("Avg",times = )
 df_freq2_long1_c <- rbind(df_freq2_long_c,df_avgfreq2c)
 
 df_avgfreq3c$rep <- rep(0,times = )
 df_avgfreq3c$Group <- rep("Avg",times = )
 df_freq3_long1_c <- rbind(df_freq3_long_c,df_avgfreq3c)
 
 
 
 df_freq1_long1_c$generation <- as.numeric(df_freq1_long1_c$generation)
 df_freq2_long1_c$generation <- as.numeric(df_freq2_long1_c$generation)
 df_freq3_long1_c$generation <- as.numeric(df_freq3_long1_c$generation)
 df_freq1_long1_c<-filter(df_freq1_long1_c,df_freq1_long1_c$generation <6) 
 df_freq2_long1_c<-filter(df_freq2_long1_c,df_freq2_long1_c$generation <6)
 df_freq3_long1_c<-filter(df_freq3_long1_c,df_freq3_long1_c$generation <6)
 
 df_freq1_long1_c$generation[which(df_freq1_long1_c$generation == 1)] <- 0
 df_freq1_long1_c$generation[which(df_freq1_long1_c$generation == 2)] <- 1
 df_freq1_long1_c$generation[which(df_freq1_long1_c$generation == 3)] <- 2
 df_freq1_long1_c$generation[which(df_freq1_long1_c$generation == 4)] <- 3
 df_freq1_long1_c$generation[which(df_freq1_long1_c$generation == 5)] <- 4 
 
 df_freq2_long1_c$generation[which(df_freq2_long1_c$generation == 1)] <- 0
 df_freq2_long1_c$generation[which(df_freq2_long1_c$generation == 2)] <- 1
 df_freq2_long1_c$generation[which(df_freq2_long1_c$generation == 3)] <- 2
 df_freq2_long1_c$generation[which(df_freq2_long1_c$generation == 4)] <- 3
 df_freq2_long1_c$generation[which(df_freq2_long1_c$generation == 5)] <- 4
 
 df_freq3_long1_c$generation[which(df_freq3_long1_c$generation == 1)] <- 0
 df_freq3_long1_c$generation[which(df_freq3_long1_c$generation == 2)] <- 1
 df_freq3_long1_c$generation[which(df_freq3_long1_c$generation == 3)] <- 2
 df_freq3_long1_c$generation[which(df_freq3_long1_c$generation == 4)] <- 3
 df_freq3_long1_c$generation[which(df_freq3_long1_c$generation == 5)] <- 4
 
 
 
 group_colors <- c("#A50F15", "#f6e09d", "#4d97cd")
 
 p_nondrive <- ggplot() +
   geom_line(
     data = filter(df_freq1_long1_c, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "grey44",
     size = 0.5
   ) +
   
   geom_line(
     data = filter(df_freq1_long1_c, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +
   geom_line(
     data = filter(df_freq1_long1, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "lightgrey",
     size = 0.5
   ) +
   geom_line(
     data = filter(df_freq1_long1, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +
   geom_line(
     data = cage_nondrive,
     aes(x = Generation, y = Freq_nonDrive, color = Group, group = Group),
     size = 1.2
   ) +
   geom_point(
     data = cage_nondrive,
     aes(x = Generation, y = Freq_nonDrive, color = Group),
     size = 2
   ) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_classic(base_size = 14) +
   ylim(0, 1)+labs(x = "Generation", y = "Frequency",color = "Cage") +theme(
     axis.text = element_text(size = 14),       # 坐标轴刻度
     axis.title = element_text(size = 15),      # 坐标轴标题
     legend.text = element_text(size = 14),     # 图例刻度
     legend.title = element_text(size = 16),
     plot.title = element_text(size = 16)     # 图例标题
   ) 

 
 p_cas9 <- ggplot() +
   
   geom_line(
     data = filter(df_freq2_long1_c, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "grey44",
     size = 0.5
   ) +
   
   geom_line(
     data = filter(df_freq2_long1_c, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +geom_line(
     data = filter(df_freq2_long1, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "lightgrey",
     size = 0.5
   ) +
   
   geom_line(
     data = filter(df_freq2_long1, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +
   
   geom_line(
     data = cage_cas9,
     aes(x = Generation, y = Freq_Cas9, color = Group, group = Group),
     size = 1.2
   ) +
   geom_point(
     data = cage_cas9,
     aes(x = Generation, y = Freq_Cas9, color = Group),
     size = 2
   ) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_classic(base_size = 14) +
   ylim(0, 1)+labs(x = "Generation", y = "Frequency",color = "Cage") +theme(
     axis.text = element_text(size = 14),       # 坐标轴刻度
     axis.title = element_text(size = 15),      # 坐标轴标题
     legend.text = element_text(size = 14),     # 图例刻度
     legend.title = element_text(size = 16),
     plot.title = element_text(size = 16)     # 图例标题
   )
 
 p_ahd <- ggplot() +
      geom_line(
     data = filter(df_freq3_long1_c, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "grey44",
     size = 0.5
   ) +
   geom_line(
     data = filter(df_freq3_long1_c, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +geom_line(
     data = filter(df_freq3_long1, Group == "Exp"),
     aes(x = generation, y = cargo_freq, group = rep),
     color = "lightgrey",
     size = 0.5
   ) +
   geom_line(
     data = filter(df_freq3_long1, Group == "Avg"),
     aes(x = generation, y = cargo_freq),
     color = "grey30",
     size = 1.2,
     linetype = "dashed"
   ) +
   geom_line(
     data = cage_ahd,
     aes(x = Generation, y = Freq_AHD, color = Group, group = Group),
     size = 1.2
   ) +
   geom_point(
     data = cage_ahd,
     aes(x = Generation, y = Freq_AHD, color = Group),
     size = 2
   ) +
   scale_color_manual(values = group_colors) +
   labs(
     x = "Generation",
     y = "Frequency",
     color = "Group"
   ) +
   theme_classic(base_size = 14) +
   ylim(0, 1)+labs(x = "Generation", y = "Frequency",color = "Cage") +theme(
     axis.text = element_text(size = 14),       # 坐标轴刻度
     axis.title = element_text(size = 15),      # 坐标轴标题
     legend.text = element_text(size = 14),     # 图例刻度
     legend.title = element_text(size = 16),
     plot.title = element_text(size = 16)     # 图例标题
   )
 
 
 p_nondrive +labs(title ="Non-drive Element Frequency")
 p_cas9+ labs(title ="Cas9 Element Frequency")
 p_ahd+ labs(title ="Modification Element Frequency")

 
 #germline: also set
 
#haplolethal: set propoetion
 #light lines on top
 #mention in discussion, think about color and order
 
 filter(df_freq2_long1_c, Group == "Avg")
 filter(df_freq3_long1_c, Group == "Avg")
 
#hairy: germline cut rate = drive conversion + germline resistance
 
 
###make embryo proportion to drive conversion as well
#0.7 maybe
 
 
 
 
 
 
 