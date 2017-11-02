import argparse
import logging
import os
import shutil
import sys
import tempfile

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MafftCommandline

import subprocess
from collections import defaultdict
import re


def iter_file(f, cols=None, sep="\t"):
    """
    Iterate over a file yielding every row (all columns or just the ones specified)
    """
    col_positions = []
    first = True
    with open(f, 'rU') as fh:
        for l in fh:
            splt = l.rstrip().split(sep)
            if first:
                first = False
                if not cols:
                    cols = splt
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
    """Find the consensus (based on quality) of the sequences (Biopython Seqs) at positions i"""
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

    #print "write_fasta_qual_files", os.access(clone_fasta_file, os.W_OK)
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
    p.wait()


def get_seqs_from_pir_file(pir_file):
    with open(pir_file, 'r') as pir:
        clone_aligned = AlignIO.read(pir_file, "fasta")
        return list(clone_aligned)


def count_leading_points(seq):
    count = 0
    for i in range(len(seq)):
        if seq[i] == ".":
            count += 1
        else:
            return count


def get_qual_from_fastq(fastq_seq, other_seq):
    """IMGT can reverse complement the sequence, so account for that when trimming the quality"""
    seq = fastq_seq.seq.lower()
    other_seq = other_seq.lower()
    quality = fastq_seq.letter_annotations["phred_quality"]
    start = seq.find(other_seq[:20])
    if start == -1:
        seq = reverse_complement(seq)
        start = seq.find(other_seq[:20])
        quality = quality[::-1]

    return quality[start:start + len(other_seq)]


dash_regex = re.compile("-")


def get_seqs_from_file(f, fmt="fasta"):
    with open(f, 'r') as f_h:
        seqs = AlignIO.read(f_h, fmt)
        return list(seqs)


def find_sequences_without_indels(seqs):
    """Find sequences with indels ('-'), seqs is a list with Biopython sequences"""
    complete_germline = str(seqs[0].seq)
    #print complete_germline
    passed_seqs = [seqs[0].id]
    for b_seq in seqs[1:]:
        seq = str(b_seq.seq)
        dashes = [i for i, b in enumerate(seq) if b != "-"]
        seq_start = min(dashes)  # find first non-dash char
        seq_end = max(dashes)  # find last non-dash char
        seq = seq[seq_start + 1:seq_end]
        germline = complete_germline[seq_start + 1:seq_end]
        germline_dashes = [dash for dash in list(dash_regex.finditer(germline))]
        skip_to_next = False
        #print "--------------------", b_seq.id, "--------------------"
        #print seq_start, seq_end
        #print seq
        #print germline
        for dash in germline_dashes:
            if seq[dash.start()] != "-":
                skip_to_next = True
                #print " " * dash.start() + "^"
                #print dash.start(), germline[dash.start()], seq[dash.start()]
                break
        if skip_to_next:
            continue
        seq_dashes = list(dash_regex.finditer(seq))
        for dash in seq_dashes:
            if germline[dash.start()] != "-":
                skip_to_next = True
                #print " " * dash.start() + "^"
                #print dash.start(), germline[dash.start()], seq[dash.start()]
                break
        if skip_to_next:
            continue
        passed_seqs.append(b_seq.id)
    return passed_seqs


def filter_file_by_ids(f, ids, f_type="fasta"):
    """Filter the specified sequence file so that only the sequences with ids in 'ids' remain"""
    if not os.path.exists(f):
        return []
    tmp_file_handle = tempfile.mkstemp()
    tmp_file = tmp_file_handle[1]
    with open(tmp_file, 'w') as tmp, open(f, 'r') as in_file:
        for seq in SeqIO.parse(in_file, f_type):
            if seq.id in ids:
                SeqIO.write(seq, tmp, f_type)
    os.close(tmp_file_handle[0])
    shutil.move(tmp_file, f)


def filter_clone_by_ids(clone_dir, clone, ids):
    """Filter all the files of a clone based on ids"""
    ids = set(ids)
    fasta_file = os.path.join(clone_dir, "clone_{0}.fasta".format(clone))
    fastq_file = os.path.join(clone_dir, "clone_{0}.fastq".format(clone))
    pir_file = os.path.join(clone_dir, "clone_{0}.pir".format(clone))
    qual_file = os.path.join(clone_dir, "clone_{0}.qual".format(clone))

    if len(ids) > 1:
        filter_file_by_ids(fasta_file, ids, f_type="fasta")
        filter_file_by_ids(fastq_file, ids, f_type="fastq")
        filter_file_by_ids(pir_file, ids, f_type="fasta")
        filter_file_by_ids(qual_file, ids, f_type="qual")
    else:
        os.remove(fasta_file)
        os.remove(fastq_file)
        os.remove(pir_file)
        os.remove(qual_file)


def write_count(dic, f):
    with open(f, 'w') as out:
        for k in sorted([int(x) for x in dic.keys()]):
            out.write("{0}\t{1}\n".format(k, dic[str(k)]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cdb", "--change-o-db", help="The input Change-O DB file, with 'CLONE' and 'IMGT_GERMLINE' columns", required=True)
    parser.add_argument("--fastq", help="The input fastq file", required=True)
    parser.add_argument("-o", "--output-dir", help="The output directory", required=True)
    parser.add_argument("--mafft", help="The location of the mafft binary", default="mafft")
    parser.add_argument("--ig-indel-jar", help="The location of the ig-indel jar", default="IndelsIdentifier.jar")
    parser.add_argument("--ig-indel-num-seqs", help="'Number of sequences' parameter for IgIndel", default="2")
    parser.add_argument("--ig-indel-hom-length", help="'Homopolymer tract length' parameter for IgIndel", default="3")
    parser.add_argument("--ig-indel-min-quality", help="'Minimal quality score' parameter for IgIndel", default="20")
    parser.add_argument("--ignore-indels", help="Remove the sequences with indels from clones", action="store_true")
    parser.add_argument("--debug", help="Output some extra debug files", action="store_true")
    parser.add_argument("--out-fasta", help="The filtered fasta")
    parser.set_defaults(ignore_indels=False)
    parser.set_defaults(debug=False)

    args = parser.parse_args()
    input_file = args.change_o_db
    fastq_file = args.fastq
    out_dir = args.output_dir
    mafft_bin = args.mafft
    output_fasta_file = args.out_fasta

    if not os.path.exists(out_dir):
        logging.info("Output dir doesn't exist, creating...")
        os.mkdir(out_dir)

    log_file = os.path.join(out_dir, "log.html")
    if os.path.exists(log_file):
        os.remove(log_file)

    logging.basicConfig(filename=log_file, level=logging.DEBUG, format="%(asctime)s:&emsp;[%(filename)s:%(lineno)d]&emsp;%(message)s <br />", datefmt='%Y/%m/%d %H:%M:%S')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info("Started fastq/change-o to Ig-Indel")

    ig_indel_jar = args.ig_indel_jar
    num_seqs = args.ig_indel_num_seqs
    hom_length = args.ig_indel_hom_length
    min_quality = args.ig_indel_min_quality
    ignore_indels = args.ignore_indels

    logging.info("cdb: {0}".format(input_file))
    logging.info("fastq: {0}".format(fastq_file))
    logging.info("output: {0}".format(out_dir))
    logging.info("mafft bin: {0}".format(mafft_bin))
    logging.info("IgIndel jar: {0}".format(ig_indel_jar))
    logging.info("Number of sequences: {0}".format(num_seqs))
    logging.info("Homopolymer tract length: {0}".format(hom_length))
    logging.info("Minimal quality score: {0}".format(min_quality))

    input_txt = os.path.join(out_dir, "input.txt")
    logging.info("Creating {0}".format(input_txt))
    with open(input_txt, 'w') as input_txt_file:  # Make the "input.txt" file
        input_txt_file.write("Number of sequences\n")
        input_txt_file.write("{0}\n".format(num_seqs))
        input_txt_file.write("Homopolymer tract length\n")
        input_txt_file.write("{0}\n".format(hom_length))
        input_txt_file.write("Minimal quality score\n")
        input_txt_file.write("{0}\n".format(min_quality))

    logging.info("Reading fastq file")
    id_to_fastq_dict = {}
    id_to_qual_dict = {}
    passed_ids = set()
    input_counts = defaultdict(int)
    indel_filter_counts = defaultdict(int)
    after_ig_indel_counts = defaultdict(int)
    for seq in SeqIO.parse(fastq_file, "fastq"):  # gather all qualities by id in a dict
        #logging.debug(seq.id)
        short_id = str(seq.id)[:50]  # IMGT cuts the ids at 50 chars
        id_to_fastq_dict[short_id] = seq

    logging.info("Getting germline and clones from Change-O DB")
    clones = set()  # store the complete/valid clones in case all seqs of a clone are not in the fastq
    last_clone = 1
    seqs_not_in_fastq = 0
    lowest_germline_start = sys.maxint
    trimmed_seqs_in_clone = []
    last_germline = ""
    for ID, in_seq, clone, seq, germline in iter_file(input_file, ["SEQUENCE_ID", "SEQUENCE_INPUT", "CLONE", "SEQUENCE_IMGT", "GERMLINE_IMGT"]):
        if ID not in id_to_fastq_dict:
            logging.warning("{0} in Change-O DB but not in fastq".format(ID))
            seqs_not_in_fastq += 1
            #logging.debug(id_to_qual_dict.keys())
            continue
        clone = int(clone)
        clones.add(clone)
        fastq_seq = id_to_fastq_dict[ID]
        seq_leading_points = count_leading_points(seq)
        trimmed_germline = germline[seq_leading_points:].replace(".", "")  # cut off the germline until the input seq starts
        trimmed_seq = in_seq.lower()[:len(trimmed_germline)]  # cut off the end of the input seq to match the length of the germline
        trimmed_quality = get_qual_from_fastq(fastq_seq, trimmed_seq)
        id_to_qual_dict[ID] = trimmed_quality

        if clone != last_clone:  # new clone means we can write the last clone to fasta/qual files
            logging.debug("Different clone, writing clone {0} to fasta/qual files".format(last_clone))
            new_germline = last_germline[lowest_germline_start:].replace(".", "")
            write_fasta_qual_files(new_germline, trimmed_seqs_in_clone, out_dir, last_clone)
            trimmed_seqs_in_clone = []
            lowest_germline_start = sys.maxint

        if seq_leading_points < lowest_germline_start:
            #logging.debug("new lowest germline start from {0} to {1}".format(lowest_germline_start, seq_leading_points))
            lowest_germline_start = seq_leading_points

        seq_record = SeqRecord(Seq(trimmed_seq, generic_dna), id=ID, description="")
        seq_record.letter_annotations["phred_quality"] = trimmed_quality

        if args.debug:
            clone_debug_fastq_file = os.path.join(out_dir, "clone_{0}_debug.fastq".format(clone))
            with open(clone_debug_fastq_file, 'a') as debug_fastq:
                germline_seq = SeqRecord(Seq(trimmed_germline, generic_dna), id="G.L", description="")
                germline_seq.letter_annotations["phred_quality"] = trimmed_quality
                SeqIO.write(germline_seq, debug_fastq, 'fastq')
                SeqIO.write(seq_record, debug_fastq, 'fastq')

        trimmed_seqs_in_clone.append(seq_record)

        last_clone = clone
        last_germline = germline

    new_germline = last_germline[lowest_germline_start:].replace(".", "")

    write_fasta_qual_files(new_germline, trimmed_seqs_in_clone, out_dir, last_clone)  # again for last clone

    for clone in clones:  # for all the complete clones, fix the Ns in the germline
        logging.debug("---------------- clone {0} --------------".format(clone))
        clone_fasta_file = os.path.join(out_dir, "clone_{0}.fasta".format(clone))
        clone_pir_file = os.path.join(out_dir, "clone_{0}.pir".format(clone))
        clone_discarded_file = os.path.join(out_dir, "clone_{0}_disc.fasta".format(clone))
        discarded_seqs = []

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
        germline_seq.seq = Seq(str(germline_seq.seq).replace("-", ""), generic_dna)
        fasta_seqs[0] = germline_seq
        with open(clone_fasta_file, 'w') as fasta:
            for seq_record in fasta_seqs:
                SeqIO.write(seq_record, fasta, 'fasta')

        new_seqs = []
        pir_seqs = get_seqs_from_file(clone_pir_file)
        input_counts[str(len(pir_seqs) - 1)] += 1
        if ignore_indels:  # read the pir file and remove the sequences with indels
            new_seqs = find_sequences_without_indels(pir_seqs)
            logging.debug("new seqs: {}".format(new_seqs))
            if len(pir_seqs) != len(new_seqs):
                with open(clone_fasta_file, 'r') as clone_fasta, open(clone_discarded_file, 'w') as clone_discarded:
                    SeqIO.write(germline_seq, clone_discarded, "fasta")
                    for seq in SeqIO.parse(clone_fasta, 'fasta'):
                        if seq.id not in new_seqs:
                            seq.id = "{0}_indel".format(seq.id)
                            seq.description = ""
                            SeqIO.write(seq, clone_discarded, "fasta")

            filter_clone_by_ids(out_dir, clone, new_seqs)

        indel_filter_counts[str(len(new_seqs) - 1)] += 1

        if len(new_seqs) == 1:
            logging.debug("Clone {0} has no sequences left after indel filtering".format(clone))
            continue

        clone_qual_file = os.path.join(out_dir, "clone_{0}.qual".format(clone))
        clone_log_file = os.path.join(out_dir, "clone_{0}_log.txt".format(clone))
        clone_fasta_out_file = os.path.join(out_dir, "clone_{0}_out.fasta".format(clone))
        ig_indel_log_file = os.path.join(out_dir, "clone_{0}-Ig-Indel-Identifier.log".format(clone))

        run_ig_indel(input_txt, clone_fasta_file, clone_qual_file, clone_pir_file, clone_log_file, clone_fasta_out_file, jar=ig_indel_jar)

        discarded_ids = set()
        good_seq_count = len(new_seqs)
        if os.path.exists(ig_indel_log_file):  # the log file can contain some sequences that should be removed
            logging.debug("Clone has Ig-Indel-Identifier.log file")
            with open(ig_indel_log_file, 'r') as ig_indel_log:
                for line in ig_indel_log:
                    if "thus, should be discarded." in line:
                        line = line[line.find(">") + 1:line.find(" contained point mutation")]  # get the ID from the line
                        discarded_ids.add(line)

            if len(discarded_ids) > 0:
                logging.debug("Clone has sequences that should be discarded {0}".format(discarded_ids))
                good_seq_count = 0
                #os.remove(clone_fasta_out_file)
                with open(clone_fasta_file, 'r') as clone_fasta, open(clone_fasta_out_file, 'w') as clone_fasta_out:
                    for seq in SeqIO.parse(clone_fasta, 'fasta'):
                        if seq.id in discarded_ids:
                            if seq.id != "G.L":
                                seq.id = "{0}_ig_indel".format(seq.id)
                                seq.description = ""
                            if seq.id not in [s.id for s in discarded_seqs]:
                                discarded_seqs.append(seq)
                        elif seq.id in new_seqs and seq.id not in discarded_ids:
                            good_seq_count += 1
                            #SeqIO.write(seq, clone_fasta_out, "fasta")
                            if seq.id != "G.L":
                                SeqIO.write(seq, clone_fasta_out, "fasta")
                                passed_ids.add(seq.id)


        after_ig_indel_counts[str(good_seq_count - 1)] += 1
        if len(discarded_seqs) > 0:
            with open(clone_discarded_file, 'w') as clone_discarded:
                SeqIO.write(germline_seq, clone_discarded, "fasta")
                for seq in discarded_seqs:
                    SeqIO.write(seq, clone_discarded, "fasta")

    id_set = set()

    for f in os.listdir(out_dir):  # Collect all ids from the "_out.fasta" files
        if not f.endswith("_out.fasta"):
            continue

        f_path = os.path.join(out_dir, f)
        for out_seq in SeqIO.parse(f_path, 'fasta'):
            id_set.add(out_seq.id)

    logging.info("{0} sequences passed filtering".format(len(id_set)))

    with open(output_fasta_file, 'w') as output_fasta:  # generate the new fasta file with the ids from the "_out.fasta" files
        for ID, in_seq in iter_file(input_file, ["SEQUENCE_ID", "SEQUENCE_INPUT"]):
            if ID in id_set:
                seq_record = SeqRecord(Seq(in_seq, generic_dna), id=ID, description="")
                SeqIO.write(seq_record, output_fasta, "fasta")
    if not ignore_indels:
        indel_filter_counts = input_counts

    counts_file = os.path.join(out_dir, "counts.txt")

    with open(counts_file, 'w') as counts:
        counts.write("count\tinput\tinput_seqs\tafter_indel_filter\tafter_indel_filter_seqs\tafter_ig_indel\tafter_ig_indel_seqs\n")
        for i in range(0, int(max(input_counts, key=lambda x: int(x))) + 1):
            counts.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                i,
                input_counts[str(i)],
                input_counts[str(i)] * i,
                indel_filter_counts[str(i)],
                indel_filter_counts[str(i)] * i,
                after_ig_indel_counts[str(i)],
                after_ig_indel_counts[str(i)] * i
            ))


if __name__ == "__main__":
    main()