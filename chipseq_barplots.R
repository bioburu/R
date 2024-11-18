df<-read.csv('/home/em_b/work_stuff/FCCC/chipseq/h3k27me3_conditional_peak_comparisons.csv',
             row.names = 1)
df
barplot(df$X0_hour_v_B22,
        main = 'H3K27me3 peaks: 0hr v B22_differentiated',
        legend.text = row.names(df),
        col = c('grey','red'),
        names.arg = c('Promoters','Distal_Intergenic'),
        ylim = c(0,700),
        border = 'blue')
barplot(df$B22_diff_v_T3,
        main = 'H3K27me3 peaks: B22_differentiated v T3_differentiated',
        legend.text = row.names(df),
        col = c('grey','red'),
        names.arg = c('Promoters','Distal_Intergenic'),
        ylim = c(0,700),
        border = 'blue')
barplot(df$X0_hour_v_T3,
        main = 'H3K27me3 peaks: 0hr v T3_differentiated',
        legend.text = row.names(df),
        col = c('grey','red'),
        names.arg = c('Promoters','Distal_Intergenic'),
        ylim = c(0,700),
        border = 'blue')
