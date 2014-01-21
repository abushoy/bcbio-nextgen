""" Extract single isoform genes from a GTF file. This is to
use for more accurately estimating the insert size and
insert size distribution during RNA-seq mapping
"""

import gffutils
from collections import defaultdict
import HTSeq

fn = "/Users/rory/cache/bcbio-nextgen/tests/data/genomes/mm9/rnaseq/ref-transcripts.gtf"
fn = "/Users/rory/hsph/hsph/biodata/genomes/Hsapiens/hg19_ERCC/rnaseq/ref-transcripts.gtf"
fn = "/Users/rory/Downloads/ref-transcripts.gtf"
#fn = "/Users/rory/Downloads/chr22.gtf"
out_file = "/Users/rory/Downloads/single_isoform.gtf"
align_file =  ("/Users/rory/cache/bcbio-nextgen/tests/test_automated_output/upload/Test1/Test1-ready.bam")
align_file = "/Users/rory/hsph/hsph/projects/bcbio-rnaseq/data/geo_data/standardization/ERCC92_subset/HBRR_rep1/HBRR_rep1-ready.bam"
align_file = "/Users/rory/Downloads/HBRR_rep1-ready.bam"
db = gffutils.create_db(fn, dbfn='/Users/rory/Downloads/ref-transcripts.gtf.db',
                        keep_order=False,
                        merge_strategy='merge', force=True, infer_gene_extent=False)
gene_table = defaultdict(set)
for cds in db.features_of_type('exon', order_by="start"):
    gene_table[cds['gene_id'][0]].update([cds['transcript_id'][0]])

single_isoform = filter(lambda x: len(gene_table[x]) == 1, gene_table)
coordinate_lookup = defaultdict(list)
with open(out_file, "w") as out_handle:
    for isoform in single_isoform:
        for feature in db.children(isoform, featuretype='exon', order_by="start"):
            out_handle.write(str(feature) + "\n")

gff = HTSeq.GFF_Reader(out_file)
exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)
exon_table = defaultdict(list)
for feature in gff:
    if feature.type == "exon":
        exon_table[feature.attr["gene_id"]].append(feature.iv)
        exons[feature.iv] += (feature.name, feature.iv)

sam_file = HTSeq.BAM_Reader(align_file)
fragment_lengths = []
multi_exon_lengths = []
for alignment in sam_file:
    if alignment.aligned and alignment.proper_pair:
        first_read = HTSeq.GenomicPosition(alignment.iv.chrom, alignment.iv.start)
        second_read = HTSeq.GenomicPosition(alignment.mate_start.chrom,
                                            alignment.mate_start.end)
        if (len(list(exons[first_read])) != 1 or
            len(list(exons[second_read])) != 1 or
            list(exons[first_read])[0][0] != list(exons[second_read])[0][0]):
            continue
        else:
            gene_name = list(exons[first_read])[0][0]
            first_exon = list(exons[first_read])[0][1]
            second_exon = list(exons[second_read])[0][1]
            if first_exon == second_exon:
                fragment_lengths.append(abs(first_read.start - second_read.end))
            else:
                spanning_exons = []
                exon_length = 0
                first_start = first_read.start - first_exon.start
                iv = HTSeq.GenomicInterval(first_read.chrom,
                                           first_read.start,
                                           second_read.end)
                for exon in exons[iv].steps():
                    if exon[1] and list(exon[1])[0][0] == gene_name:
                        exon_length += abs(list(exon[1])[0][1].start -
                                           list(exon[1])[0][1].end)
                        spanning_exons.append(list(exon[1])[1])
                last_exon = spanning_exons[-1]
                offset = exon_length - abs(last_exon.start - last_exon.end)
                second_end = offset + abs(second_read.end - last_exon.start)
                multi_exon_lengths.append(abs(first_start - second_end))

print fragment_lengths
print multi_exon_lengths

