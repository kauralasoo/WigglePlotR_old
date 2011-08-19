"""
Script for converting MISO alternative event annotation GFF files 
into BED format.

Usage: python gff2bed.py <gff_file>
"""
import sys

class GFFFeature():
    """
    Class for storing features from GFF file
    """
    def __init__(self, chromosome, category, type, start, end, strand, id):
        self.chromosome = chromosome
        self.category = category
        self.type = type
        self.start = start
        self.end = end
        self.strand = strand
        self.id = id
        self.children = dict()
    
    def __str__(self):
        return self.id

def parse_line(line):
    """
    Parse one line from GFF file. Converts last field of the GFF to dict.
    """
    line = line.rstrip()
    fields = line.split("\t")
    attribute_dict = dict()
    attrs = fields[8].split(";")
    for attr in attrs:
        attr_list = attr.split("=")
        attribute_dict[attr_list[0]] = attr_list[1]
    return (fields, attribute_dict)

def parse_gff(gff_file):
    """
    Loads GFF file into a list of GFFFeature clss instances.
    """
    file = open(gff_file, 'r')
    gene_list = list()
    current_gene = None
    for line in file:
        fields, attribute_dict = parse_line(line)
        if fields[2] == "gene":
            if current_gene:
                gene_list.append(current_gene)
            current_gene = GFFFeature(fields[0], fields[1], fields[2], fields[3], 
                                      fields[4], fields[6], attribute_dict['ID'])
        elif fields[2] == "mRNA":
            mRNA = GFFFeature(fields[0], fields[1], fields[2], fields[3], 
                              fields[4], fields[6], attribute_dict['ID'])
            current_gene.children[attribute_dict['ID']] = mRNA
        elif fields[2] == "exon":
            exon = GFFFeature(fields[0], fields[1], fields[2], fields[3], 
                              fields[4], fields[6], attribute_dict['ID'])
            current_gene.children[attribute_dict["Parent"]].children[attribute_dict['ID']] = exon
    gene_list.append(current_gene)
    return gene_list

def print_bed(gene_list):
    """
    Prints a list of GFF gene records in BED format.
    """
    for gene in gene_list:
        for mRNA in gene.children.values():
            exon_nr = len(mRNA.children.values())
            exon_tuples = list() 
            for exon in mRNA.children.values():
                rel_start = int(exon.start) - int(mRNA.start)
                length = int(exon.end) - int(exon.start)
                exon_tuples.append((rel_start, length))
            exon_tuples = sorted(exon_tuples, key=lambda length: length[0])
            exon_starts = [str(exon_tuple[0]) for exon_tuple in exon_tuples]
            exon_lengths = [str(exon_tuple[1]) for exon_tuple in exon_tuples]
            print '%s\t%s\t%s\t%s\t0\t%s\t%s\t%s\t0\t%s\t%s,\t%s,' %(mRNA.chromosome, mRNA.start, mRNA.end,
                                                     mRNA.id, mRNA.strand, mRNA.start, mRNA.end,
                                                     exon_nr, ",".join(exon_lengths), ",".join(exon_starts))

if __name__ == "__main__":
    gff_file = sys.argv[1]
    gene_list = parse_gff(gff_file)
    print_bed(gene_list)

        

        