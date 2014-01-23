""" Extract single isoform genes from a GTF file. This is to
use for more accurately estimating the insert size and
insert size distribution during RNA-seq mapping
"""
from collections import defaultdict
from bcbio.log import logger
from bcbio import bam

import tempfile
import gffutils
import HTSeq
import sys

def gene_level_feature(db):
    """
    return the gene level feature to query on in the DB
    """
    gene_features = ["mRNA", "gene", "CDS"]
    features = list(db.featuretypes())
    for f in gene_features:
        if f in features:
            return f
    return None

def average_intron_size(gtf_file):
    db = gffutils.helpers.get_gff_db(gtf_file)
    gene_feature = gene_level_feature(db)
    samples = 100
    average = 0
    for gene in range(samples):
        for cds in db.features_of_type(gene_feature):
            exons = db.children(cds.id, featuretype='exon')
            for intron in list(db.interfeatures(exons)):
                average = (average + abs(intron.start - intron.end)) / 2
    return average


def single_isoform_genes(gtf_file):
    """
    return the filename of a GTF file where all of the genes have a single
    annotated isoform
    """
    db = gffutils.helpers.get_gff_db(gtf_file)
    out_file = tempfile.NamedTemporaryFile()
    gene_table = defaultdict(set)
    for exon in db.features_of_type('exon', order_by="start"):
        gene_table[exon['gene_id'][0]].update([exon['transcript_id'][0]])
    single_isoform = filter(lambda x: len(gene_table[x]) == 1, gene_table)
    with open(out_file, "w") as out_handle:
        for isoform in single_isoform:
            for feature in db.children(isoform, featuretype="exon",
                                       order_by="start"):
                out_handle.write(str(feature) + "\n")
    return out_file


def exon_intervals_of_single_isoform_genes(gtf_file):
    single_isoform_gtf = single_isoform_genes(gtf_file)
    gff = HTSeq.GFF_Reader(single_isoform_gtf)
    exon_intervals = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for feature in gff:
        if feature.type == "exon":
            exon_intervals[feature.iv] += (feature.name, feature.iv)
    return exon_intervals


def htseq_reader(align_file):
    """
    returns a read-by-read sequence reader for a BAM or SAM file
    """
    if bam.is_sam(align_file):
        read_seq = HTSeq.SAM_Reader(align_file)
    elif bam.is_bam(align_file):
        read_seq = HTSeq.BAM_Reader(align_file)
    else:
        logger.error("%s is not a SAM or BAM file" % (align_file))
        sys.exit(1)
    return read_seq


def is_paired(align_file):
    """
    returns True if an alignment file has paired end reads in it
    """
    read_seq = htseq_reader(align_file)
    r1_start = iter(read_seq).next()
    pe_mode = r1_start.paired_end
    return pe_mode

def _gene_id(start):
    return start[0][0]

def _exon_interval(start):
    return start[0][1]

def starts_are_ambiguous(r1_start, r2_start):
    """
    starts should map to only one exon, and both ends should
    map to a gene of the same ID
    """
    return (len(r1_start) != 1 or len(r2_start) != 1 or
            _gene_id(r1_start) != _gene_id(r2_start))

def rnaseq_insert_sizes(align_file, gtf_file, stranded=False):
    """
    RNA-seq data estimation straight from the alignment distances of
    the two reads is not correct as sometimes read pairs cross exon-exon
    boundaries. This finds all single-allelic genes in the annotation
    and uses those to calculate the insert size statistics

    XXX:
    stranded is not respected at the moment
    """
    single_isoform_gtf = single_isoform_genes(gtf_file)
    exons = exon_intervals_of_single_isoform_genes(single_isoform_gtf)

    assert is_paired(align_file), ("%s is being used to calculate insert size "
                                   "statistics but %s does not seem to have "
                                   "paired-end reads.")
    sam_file = htseq_reader(align_file)
    fragment_lengths = []
    multi_exon_lengths = []
    for alignment in sam_file:
        if alignment.aligned and alignment.proper_pair:
            r1_start = HTSeq.GenomicPosition(alignment.iv.chrom,
                                             alignment.iv.start)
            r2_start = HTSeq.GenomicPosition(alignment.mate_start.chrom,
                                             alignment.mate_start.end)
            r1_start_exons = list(exons[r1_start])
            r2_start_exons = list(exons[r2_start])
            if starts_are_ambiguous(r1_start_exons, r2_start_exons):
                continue
            else:
                gene_name = _gene_id(r1_start_exons)
                first_exon = _exon_interval(r1_start_exons)
                second_exon = _exon_interval(r2_start_exons)
                if first_exon == second_exon:
                    fragment_lengths.append(abs(r1_start.start - r2_start.end))
                else:
                    spanning_exons = []
                    exon_length = 0
                    first_start = r1_start.start - first_exon.start
                    iv = HTSeq.GenomicInterval(r1_start.chrom,
                                               r1_start.start,
                                               r2_start.end)
                    for exon in exons[iv].steps():
                        if exon[1] and list(exon[1])[0][0] == gene_name:
                            exon_length += abs(list(exon[1])[0][1].start -
                                               list(exon[1])[0][1].end)
                            spanning_exons.append(list(exon[1])[1])
                    last_exon = spanning_exons[-1]
                    offset = exon_length - abs(last_exon.start - last_exon.end)
                    second_end = offset + abs(r2_start.end - last_exon.start)
                    multi_exon_lengths.append(abs(first_start - second_end))
    return multi_exon_lengths, fragment_lengths


if __name__ == "__main__":
    rnaseq_insert_sizes(sys.argv[1], sys.argv[2])
