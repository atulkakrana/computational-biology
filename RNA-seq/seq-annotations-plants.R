## R script for annotating plant genes

require(biomaRt)
listMarts() ## Check the marts - plant_mart_XX

# listDatasets(plantMart) ## Find species dataset
# plantMart = useMart("plants_mart_23",dataset="zmays_eg_gene")
attribs = listAttributes(plantMart);attribs[1:100,]
# filters = listFilters(plantMart);filters[1:100,]

plantMart = useMart("plants_mart_24")
plantMart = useDataset("zmays_eg_gene",mart = plantMart)

names(seq.results)
# test.ids = c('GRMZM2G180251','GRMZM2G019971')
geneIDs = seq.results$gene_short_name
length(geneIDs);geneIDs[1:5]

## Filters is like "where" of mySQL - In this case we specificy in query WHERE 'ensembl_gene_id' =test.ids
test.anno = getBM(attributes = c('ensembl_gene_id','chromosome_name','start_position','end_position','strand',"description"),
                  filters = 'ensembl_gene_id',values = geneIDs, mart = plantMart) ## Server is Down Sometimes
dim(test.anno);dim(seq.results);names(test.anno)
# write.table(test.anno,"seqLimmaAnnoTest.tsv",sep = '\t',row.names=F)

## Merge annotation to results frame
final.seq = merge(seq.results,test.anno,by.x="gene_short_name",by.y="ensembl_gene_id",all.x=TRUE)
dim(final.seq)
write.table(final.seq,"ResAnnotated2UnpairedTtest.txt",sep = "\t",row.names=FALSE)