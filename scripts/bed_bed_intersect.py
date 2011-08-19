"""
Script that goes through a query .bed file and reports all transcripts and genes
from UCSC .tsv file that overlap with each record in the .bed file.

Usage:
python refgene-bed-intersect.py <refgene.tsv> <query.bed>
"""

#Reference TSV file from UCSC Table Viewer
reference_bed = open(sys.argv[1],'r') 
#Query BED file
query_bed = open(sys.argv[2],'r') 
reference_bed.readline()

reference_records = dict()
for line in reference_bed:
    record = dict()
    fields = line.split("\t")
    record["transcript_id"] = fields[1]
    record["gene_id"] = fields[12]
    record["tx_start"] = int(fields[4])
    record["tx_end"] = int(fields[5])
    record["chr"] = fields[2]
    record["strand"] = fields[3]
    if record["chr"] not in reference_records:
        reference_records[record["chr"]] = list()
        reference_records[record["chr"]].append(record)
    else:
        reference_records[record["chr"]].append(record)

for line in query_bed:
    fields = line.split("\t")
    record = dict()
    record["chr"] = fields[0]
    record["tx_start"] = int(fields[1])
    record["tx_end"] = int(fields[2])
    name = fields[3].split(".")
    name = ".".join(name[0:len(name)-1])
    record["name"] = name
    references = reference_records[record["chr"]]
    for reference in references:
        if ((record["tx_start"] >= reference["tx_start"] and record["tx_start"] <= reference["tx_end"]) or
            (record["tx_end"] >= reference["tx_start"] and record["tx_end"] <= reference["tx_end"])):
            print "%s\t%s\t%s" %(record["name"], reference["transcript_id"], reference["gene_id"])