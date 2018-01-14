####

## Merge RSEM counts and FPKMs ##################A
afile = read.delim2("genes.attr_table",header = T)
bfile = read.delim2("genes.fpkm_table",header = T)
cfile = read.delim2("genes.count_table",header = T)

fpkm.counts = merge(bfile,cfile,by.x = "tracking_id",by.y = "tracking_id",all=T)
dim(fpkm.counts);fpkm.counts[1:5,]

final = merge(afile,fpkm.counts,by.x = "tracking_id",by.y = "tracking_id",all=T)
write.table(final,file="Anno_fpkm_counts.v2.txt",sep="\t",row.names=F)
dim(final);final[1:5,]
