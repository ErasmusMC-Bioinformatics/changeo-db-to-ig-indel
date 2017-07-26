import argparse
import logging
import os
import sys

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MafftCommandline

import subprocess
from collections import defaultdict


def iter_file(f, cols=None, sep="\t"):
    col_positions = []
    first = True
    with open(f, 'rU') as fh:
        for l in fh:
            splt = l.rstrip().split(sep)
            if first:
                first = False
                if not cols:
                    cols = splt
                print splt
                for col in cols:
                    col_positions.append(splt.index(col))
                continue
            yield [splt[x] for x in col_positions]


def run_mafft(fasta_in, pir_out, mafft="mafft"):
    logging.debug("Running mafft ({0}) on {1}, output {2}".format(mafft, fasta_in, pir_out))
    with open(pir_out, 'w') as pir:
        mafft_cline = MafftCommandline(mafft, input=fasta_in)
        stdout, stderr = mafft_cline()
        pir.write(stdout)


def get_consensus_nt(seqs, i):
    """Find the consensus (based on quality) of the sequences at positions i"""
    count = defaultdict(int)
    for seq in seqs:
        quality = seq.letter_annotations["phred_quality"]
        count[seq[i]] += quality[i]
    consensus_nt = max(count, key=lambda x: count[x])
    logging.debug("Position {0} is a N nucleotide, consensus is: {1} ({2})".format(i, consensus_nt, count))
    return consensus_nt


def fix_n_in_germline(germline, other_seqs):
    """Find the consensus (based on quality) of the n|N nucleotides in the germline sequence"""
    germline_list = [x for x in germline]
    n_pos = -1
    n_char = "N" if germline.find("N") != -1 else "n"
    n_count = germline.count(n_char)
    logging.debug("{0} N necleotides in germline".format(n_count))
    for _ in range(n_count):
        n_pos = germline.find(n_char, n_pos + 1)
        germline_list[n_pos] = get_consensus_nt(other_seqs, n_pos)
    return SeqRecord(Seq("".join(germline_list), generic_dna), id="G.L", description="")


def write_fasta_qual_files(germline, seqs, out_dir, clone):
    """Write the clone fasta/qual file of the provided germline and sequences"""
    logging.info("Creating new fasta/qual files for clone {0}".format(clone))

    clone_fasta_file = os.path.join(out_dir, "clone_{0}.fasta".format(clone))
    clone_qual_file = os.path.join(out_dir, "clone_{0}.qual".format(clone))
    clone_fastq_file = os.path.join(out_dir, "clone_{0}.fastq".format(clone))

    GL_record = SeqRecord(Seq(germline, generic_dna), id="G.L", description="")

    with open(clone_fasta_file, 'w') as fa, open(clone_qual_file, 'w') as fq, open(clone_fastq_file, 'w') as fastq:
        SeqIO.write(GL_record, fa, "fasta")
        for trimmed_seq_in_clone in seqs:
            SeqIO.write(trimmed_seq_in_clone, fa, 'fasta')
            SeqIO.write(trimmed_seq_in_clone, fq, 'qual')
            SeqIO.write(trimmed_seq_in_clone, fastq, 'fastq')


def fix_qual_for_aligned_seq(seq_record, qual):
    """Insert extra quality at the '-' positions"""
    max_qual = max(qual)
    for i in range(len(seq_record)):
        if seq_record[i] == "-":
            qual.insert(i, max_qual)
    return qual


def run_ig_indel(input_txt, fasta, qual, pir, log, fasta_out, jar="IndelsIdentifier"):
    logging.info("Running IgIndel ({0}) on {1}".format(jar, fasta))
    command = ["java", "-jar", jar, input_txt, fasta, qual, pir, log, fasta_out]
    logging.debug(" ".join(command))
    p = subprocess.Popen(command, stdout=sys.stdout, stderr=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cdb", "--change-o-db", help="The input Change-O DB file, with 'CLONE' and 'IMGT_GERMLINE' columns", required=True)
    parser.add_argument("--fastq", help="The input fastq file", required=True)
    parser.add_argument("-o", "--output-dir", help="The output directory", required=True)
    parser.add_argument("--mafft", help="The location of the mafft binary", default="mafft")
    parser.add_argument("--ig-indel_jar", help="The location of the ig-indel jar", default="IndelsIdentifier.jar")
    parser.add_argument("--ig-indel-num-seqs", help="'Number of sequences' parameter for IgIndel", default="2")
    parser.add_argument("--ig-indel-hom-length", help="'Homopolymer tract length' parameter for IgIndel", default="3")
    parser.add_argument("--ig-indel-min-quality", help="'Minimal quality score' parameter for IgIndel", default="20")

    logging.basicConfig(filename="./log.html", level=logging.DEBUG, format="%(asctime)s:&emsp;%(message)s <br />",
                        datefmt='%Y/%m/%d %H:%M:%S')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info("Started fastq/change-o to Ig-Indel")

    args = parser.parse_args()
    input_file = args.change_o_db
    fastq_file = args.fastq
    out_dir = args.output_dir
    mafft_bin = args.mafft

    ig_indel_jar = args.ig_indel_jar
    num_seqs = args.ig_indel_num_seqs
    hom_length = args.ig_indel_hom_length
    min_quality = args.ig_indel_min_quality

    logging.info("cdb: {0}".format(input_file))
    logging.info("fastq: {0}".format(fastq_file))
    logging.info("output: {0}".format(out_dir))
    logging.info("mafft bin: {0}".format(mafft_bin))
    logging.info("IgIndel jar: {0}".format(ig_indel_jar))
    logging.info("Number of sequences: {0}".format(num_seqs))
    logging.info("Homopolymer tract length: {0}".format(hom_length))
    logging.info("Minimal quality score: {0}".format(min_quality))

    if not os.path.exists(out_dir):
        logging.info("Output dir doesn't exist, creating...")
        os.mkdir(out_dir)

    input_txt = os.path.join(out_dir, "input.txt")
    logging.info("Creating {0}".format(input_txt))
    with open(input_txt, 'w') as input_txt_file:
        input_txt_file.write("Number of sequences\n")
        input_txt_file.write("{0}\n".format(num_seqs))
        input_txt_file.write("Homopolymer tract length\n")
        input_txt_file.write("{0}\n".format(hom_length))
        input_txt_file.write("Minimal quality score\n")
        input_txt_file.write("{0}\n".format(min_quality))

    logging.info("Reading fastq file")
    id_to_qual_dict = {}
    for seq in SeqIO.parse(fastq_file, "fastq"):  # gather all qualities by id in a dict
        logging.debug(seq.id)
        short_id = str(seq.id)[:50]
        id_to_qual_dict[short_id] = seq.letter_annotations["phred_quality"]
        #print seq

    logging.info("Getting germline and clones from Change-O DB")
    clones = set()  # store the complete/valid clones in case all seqs of a clone are not in the fastq
    last_clone = 1
    last_germline = ""
    trimmed_seqs_in_clone = []
    for ID, in_seq, clone, seq, germline in iter_file(input_file, ["SEQUENCE_ID", "SEQUENCE_INPUT", "CLONE", "SEQUENCE_IMGT", "GERMLINE_IMGT"]):
        if ID not in id_to_qual_dict:
            logging.info("{0} in Change-O DB but not in fastq".format(ID))
            logging.debug(id_to_qual_dict.keys())
            continue
        clone = int(clone)
        clones.add(clone)
        quality = id_to_qual_dict[ID]
        trimmed_germline = germline[germline.find(in_seq[:10]):].replace(".", "")  # cut off the germline until the input seq starts
        trimmed_seq = in_seq.lower()[:len(trimmed_germline)]  # cut off the end of the input seq to match the length of the germline
        trimmed_quality = quality[:len(trimmed_seq)]  # trim the quality to match the length of the sequence
        id_to_qual_dict[ID] = trimmed_quality

        logging.debug("ID {0}, clone {1}".format(ID, clone))

        #print trimmed_germline
        #print trimmed_seq
        #print ("I" * len(trimmed_quality))
        #print ("-" * len(trimmed_seq))

        if clone != last_clone:  # new clone means we can write the last clone to fasta/qual files
            logging.debug("Different clone, writing clone {0} to fasta/qual files".format(last_clone))
            write_fasta_qual_files(last_germline, trimmed_seqs_in_clone, out_dir, last_clone)
            trimmed_seqs_in_clone = []

        seq_record = SeqRecord(Seq(trimmed_seq, generic_dna), id=ID, description="")
        seq_record.letter_annotations["phred_quality"] = trimmed_quality

        trimmed_seqs_in_clone.append(seq_record)

        last_clone = clone
        last_germline = trimmed_germline

    write_fasta_qual_files(last_germline, trimmed_seqs_in_clone, out_dir, last_clone)  # again for last clone

    for clone in clones:  # for all the complete clones, fix the align them with mafft and fix the Ns in the germline
        clone_fasta_file = os.path.join(out_dir, "clone_{0}.fasta".format(clone))
        clone_pir_file = os.path.join(out_dir, "clone_{0}.pir".format(clone))

        run_mafft(clone_fasta_file, clone_pir_file, mafft_bin)

        logging.info("Fixing N|n nucleotides in G.L sequence for clone {0}".format(clone))
        clone_aligned = AlignIO.read(clone_pir_file, "fasta")
        for seq_record in clone_aligned:  # add the quality scores back to the seqs so consensus nt can use quality
            if seq_record.id == "G.L":
                continue
            qual = id_to_qual_dict[seq_record.id]
            if len(qual) != len(seq_record):
                qual = fix_qual_for_aligned_seq(seq_record, qual)
            seq_record.letter_annotations["phred_quality"] = qual
        germline_seq = fix_n_in_germline(clone_aligned[0].seq, clone_aligned[1:])
        logging.info("Writing new G.L sequence to pir file")
        with open(clone_pir_file, 'w') as pir:
            SeqIO.write(germline_seq, pir, 'fasta')
            for seq_record in clone_aligned[1:]:
                SeqIO.write(seq_record, pir, 'fasta')

        logging.info("Writing new G.L sequence to fasta file")
        fasta_seqs = list(SeqIO.parse(clone_fasta_file, 'fasta'))
        fasta_seqs[0] = germline_seq
        with open(clone_fasta_file, 'w') as fasta:
            for seq_record in fasta_seqs:
                SeqIO.write(seq_record, fasta, 'fasta')

        clone_qual_file = os.path.join(out_dir, "clone_{0}.qual".format(clone))
        clone_log_file = os.path.join(out_dir, "clone_{0}_log.txt".format(clone))
        clone_fasta_out_file = os.path.join(out_dir, "clone_{0}_out.fasta".format(clone))

        run_ig_indel(input_txt, clone_fasta_file, clone_qual_file, clone_pir_file, clone_log_file, clone_fasta_out_file, jar=ig_indel_jar)


if __name__ == "__main__":
    main()