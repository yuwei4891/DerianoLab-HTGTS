sample<-read.delim("r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate.V_segment_resection_micro_length.resection_length_plot_combined_signal_coding",header=F)

tiff("combined_Resection_length_distribution.tiff",width=9,height=6,units='in',res=300)

par(mfrow=c(2,3),mar=c(2,2.5,2,2))
plot(c(0, 100), c(0, 20000), type = "n", xlab = "", ylab = "",frame=FALSE,xaxt='n',yaxt='n')
axis(1,cex.axis=1.2)
axis(2,cex.axis=1.2)
nrow<-nrow(sample)
colors <- colorRampPalette(c("gray35", "gold"))(10)
i<-1:nrow
rect(sample[i,1], sample[i,4], sample[i,2], sample[i,5],col=colors[sample[i,6]+1], border = NA)

dev.off()