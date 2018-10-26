
def tousidentfile(gn1):
	header=["##gff-version 3","#!gff-spec-version 1.20","#!processor NCBI annotwriter",
	"##sequence-region tousgenesidentiques 1 %d" % gn1.genome_len]
	f= open("tousgenesidentiques1.gff","w+")
	for i in range(len(header)):
		f.write("%s\n" % header[i])
	f.write("tousgenesidentiques\tRefSeq\tregion\t1\t%d\t.\t+\t.\tID=id0;Name=tousgenesidentiques\n" %  gn1.genome_len
	)
	for n in range(10):
		base="tousgenesidentiques\tRefSeq\tgene\t%s\t%s\t.\t%s\t.\tID=g1;Name=g%s\n" % (gn1.gene_list[n].start,gn1.gene_list[n].end,gn1.gene_list[n].orientation,gn1.gene_list[n].id)
		f.write(base)
	f.close()