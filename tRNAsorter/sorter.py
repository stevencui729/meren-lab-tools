# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to deal with tRNA sequences."""

import os
import csv
import sys
import tempfile
import extractor
import dbops
import Levenshtein as lev

import fastalib as u



__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"



class SorterStats:
    """This class handles keeping track of sort statistics."""

    def __init__(self):
        """Initializes statistics."""
        self.total_seqs = 0
        self.total_rejected = 0
        self.total_passed = 0
        self.num_trailer = 0
        self.total_full_length = 0

        self.no_divergence = 0
        self.t_loop_divergence = 0
        self.div_at_0 = 0
        self.div_at_1 = 0
        self.div_at_2 = 0
        self.div_at_3 = 0
        self.div_at_8 = 0
        self.acceptor_divergence = 0
        self.div_at_neg_1 = 0
        self.div_at_neg_2 = 0
        self.div_at_neg_3 = 0

        self.t_loop_seq_rejected = 0
        self.acceptor_seq_rejected = 0
        self.both_rejected = 0
        self.short_rejected = 0


    def gen_sql_query_info_tuple(self):
        info_string_list = []
        info_string_list.append(self.total_seqs)
        info_string_list.append(self.total_rejected)
        info_string_list.append(self.total_passed)
        info_string_list.append(self.num_trailer)
        info_string_list.append(self.total_full_length)
        info_string_list.append(self.no_divergence)
        info_string_list.append(self.t_loop_divergence)
        info_string_list.append(self.div_at_0)
        info_string_list.append(self.div_at_1)
        info_string_list.append(self.div_at_2)
        info_string_list.append(self.div_at_3)
        info_string_list.append(self.div_at_8)
        info_string_list.append(self.acceptor_divergence)
        info_string_list.append(self.div_at_neg_1)
        info_string_list.append(self.div_at_neg_2)
        info_string_list.append(self.div_at_neg_3)
        info_string_list.append(self.t_loop_seq_rejected)
        info_string_list.append(self.acceptor_seq_rejected)
        info_string_list.append(self.both_rejected)
        info_string_list.append(self.short_rejected)
        return tuple(info_string_list)


   #def format_line(self, label, value, level, padding = 55):
   #    """Handles indenting/formatting lines for statistics."""
   #    levels_dict = {1:"%s%s\t%s\n" %
   #        (label, " " + " " * (padding - len(label)), value),
   #        2:"\t%s%s\t%s\n" %
   #            (label, " " + " " * (padding - (4 + len(label))), value),
   #        3:"\t\t%s%s\t%s\n" %
   #            (label, " " + " " * (padding - (12 + len(label))), value)}
   #    return levels_dict[level]


   #def write_stats(self, out_file_path):
   #    """Writes statistics to an output file."""
   #    with open(out_file_path, "w") as outfile:
   #        outfile.write(self.format_line("Total seqs", "%d" %
   #            self.total_seqs,1))
   #        outfile.write(self.format_line("Total full-length", "%d" %
   #            self.total_full_length, 1))
   #        outfile.write(self.format_line("With trailer", "%d" %
   #            self.num_trailer, 1))
   #        outfile.write(self.format_line("Total passed", "%d" %
   #            self.total_passed, 1))

   #        outfile.write(self.format_line("No divergence", "%d" %
   #            self.no_divergence, 2))
   #        outfile.write(self.format_line("T-loop divergence", "%d" %
   #            self.t_loop_divergence, 2))
   #        outfile.write(self.format_line("Divergence at pos 0", "%d" %
   #            self.div_at_0, 3))
   #        outfile.write(self.format_line("Divergence at pos 1", "%d" %
   #            self.div_at_1, 3))
   #        outfile.write(self.format_line("Divergence at pos 2", "%d" %
   #            self.div_at_2, 3))
   #        outfile.write(self.format_line("Divergence at pos 3", "%d" %
   #            self.div_at_3, 3))
   #        outfile.write(self.format_line("Divergence at pos 8", "%d" %
   #            self.div_at_8, 3))
   #        outfile.write(self.format_line("Acceptor divergence", "%d" %
   #            self.acceptor_divergence, 2))
   #        outfile.write(self.format_line("Divergence at pos -3", "%d" %
   #            self.div_at_neg_3, 3))
   #        outfile.write(self.format_line("Divergence at pos -2", "%d" %
   #            self.div_at_neg_2, 3))
   #        outfile.write(self.format_line("Divergence at pos -1", "%d" %
   #            self.div_at_neg_1, 3))

   #        outfile.write(self.format_line("Total failed", "%d" %
   #            self.total_rejected, 1))
   #        outfile.write(self.format_line("T-loop seq rejected", "%d" %
   #            self.t_loop_seq_rejected, 2))
   #        outfile.write(self.format_line("Acceptor seq rejected", "%d" %
   #            self.acceptor_seq_rejected, 2))
   #        outfile.write(self.format_line("Both rejected", "%d" %
   #            self.both_rejected, 2))
   #        outfile.write(self.format_line("Short rejected", "%d" %
   #            self.short_rejected, 2))


class SeqSpecs:
    """Class to store seq information during sort."""

    def __init__(self):
        """Initializes seq information categories."""
        self.length = 0
        self.mis_count = 100
        self.t_loop_error = True
        self.acceptor_error = True
        self.seq = ""
        self.seq_sub = ""
        self.full_length = False
        self.t_loop_seq = ""
        self.acceptor_seq = ""
        self.three_trailer = ""
        self.trailer_length = 0
        self.anticodon = ""


   #def gen_id_string(self, id):
   #    """Generates ID string for _PASSED and _FAILED files"""
   #    mod_id_list = []
   #    mod_id_list.append(id)
   #    mod_id_list.append("mismatches:" + str(self.mis_count))
   #    mod_id_list.append("t_loop_error:" + str(self.t_loop_error))
   #    mod_id_list.append("acceptor_error:" + str(self.acceptor_error))
   #    mod_id_list.append("seq_sub:" + self.seq_sub)
   #    if self.full_length:
   #        mod_id_list.append("full length: Yes")
   #    else:
   #        mod_id_list.append("full length: Maybe")
   #    mod_id_list.append("t_loop:" + self.t_loop_seq)
   #    mod_id_list.append("acceptor:" + self.acceptor_seq)
   #    mod_id_list.append("3_trailer:" + str(self.three_trailer))
   #    return "|".join(mod_id_list)


   #def write_specs(self, writer, id):
   #    """Uses the passed in DictWriter to write to a csv output file
   #    (_TAB_NO_TRAILER or _TAB_TRAILER).
   #    """
   #    temp_dict = {"ID" : id, "Seq" : self.seq, "3-trailer" :
   #        self.three_trailer, "t-loop" : self.t_loop_seq, "acceptor" :
   #        self.acceptor_seq, "full-length" : str(self.full_length),
   #        "Seq_length" : str(self.length), "Trailer_length" :
   #        str(self.trailer_length), "Anticodon" : self.anticodon}
   #    writer.writerow(temp_dict)

    def gen_sql_query_info_tuple(self, id):
        """Generates tuple of values to add into database to use in SQL query."""
        info_string_list = []
        info_string_list.append(id)
        info_string_list.append(self.seq)


        if self.trailer_length == 0:
            info_string_list.append(None)
        else:
            info_string_list.append(self.three_trailer)
       #info_string_list.append(self.three_trailer)

        info_string_list.append(self.t_loop_seq)
        info_string_list.append(self.acceptor_seq)
        info_string_list.append(str(self.full_length))
        info_string_list.append(str(self.length))
        info_string_list.append(str(self.trailer_length))


        if len(self.anticodon) == 0:
            info_string_list.append(None)
        else:
            info_string_list.append(self.anticodon)
       #info_string_list.append(self.anticodon)

        return tuple(info_string_list)


class Sorter:
    """Class that handles the sorting of the seqs."""

    def __init__(self):
        """Initializes variables for input/output files and statistics."""
        self.passed_seqs_write_fasta = ""
        self.rejected_seqs_write_fasta = ""
        self.sort_stats_write_file = ""
        self.read_fasta = ""
        self.no_trailer_tabfile = ""
        self.trailer_tabfile = ""
        self.tRNA_DB_file = ""

        self.fieldnames = ["ID", "Seq", "3-trailer", "t-loop", "acceptor",
            "full-length", "Seq_length", "Trailer_length", "Anticodon"]
        self.max_seq_width = len(self.fieldnames[1])

        self.sort_stats = SorterStats()
        self.extractor = extractor.Extractor()
        self.db = None
        seq_count_dict = {}


    def set_file_names(self, args):
        """Takes command line arguments from args and assigns input/output file
        variables.
        """
       #self.passed_seqs_write_fasta = u.FastaOutput(args.sample_name +
       #    "_PASSED")
       #self.rejected_seqs_write_fasta = u.FastaOutput(args.sample_name +
       #    "_FAILED")
        self.sort_stats_write_file = args.sample_name + "_SORTER_STATS.txt"
        self.read_fasta = u.SequenceSource(args.readfile)
        self.no_trailer_tabfile = args.sample_name + "_TAB_NO_TRAILER"
        self.trailer_tabfile = args.sample_name + "_TAB_TRAILER"

        self.extractor.set_file_names(args.sample_name)

        if args.output_path:
            if args.output_path.endswith(".db"):
                self.tRNA_DB_file = args.output_path
            else:
                print "given output file not correct format"
                sys.exit(1)
        else:
            self.tRNA_DB_file = args.sample_name + ".db"

    def check_divergence_pos(self, cur_seq_specs):
        """Takes a SeqSpecs class and updates statistics on divergence
        position.
        """
        if cur_seq_specs.t_loop_error:
            self.sort_stats.t_loop_divergence += 1
            if cur_seq_specs.seq_sub[0] != "G":
                self.sort_stats.div_at_0 += 1
            elif cur_seq_specs.seq_sub[1] != "T":
                self.sort_stats.div_at_1 += 1
            elif cur_seq_specs.seq_sub[2] != "T":
                self.sort_stats.div_at_2 += 1
            elif cur_seq_specs.seq_sub[3] != "C":
                self.sort_stats.div_at_3 += 1
            elif cur_seq_specs.seq_sub[8] != "C":
                self.sort_stats.div_at_8 += 1
        elif cur_seq_specs.acceptor_error:
            self.sort_stats.acceptor_divergence += 1
            if cur_seq_specs.seq_sub[-3] != "C":
                self.sort_stats.div_at_neg_3 += 1
            elif cur_seq_specs.seq_sub[-2] != "C":
                self.sort_stats.div_at_neg_2 += 1
            elif cur_seq_specs.seq_sub[-1] != "A":
                self.sort_stats.div_at_neg_1 += 1
        else:
            self.sort_stats.no_divergence += 1


    def check_full_length(self, cur_seq_specs):
        """Takes a SeqSpecs class and checks whether or not the sequence is a
        full-length sequence, returns an updated SeqSpecs class. Also determines
        the anticodon if it is a full-length seq
        """
        if cur_seq_specs.length > 70 and cur_seq_specs.length < 100:
            if cur_seq_specs.seq[7] == "T" and cur_seq_specs.seq[13] == "A":
                self.sort_stats.total_full_length += 1
                cur_seq_specs.full_length = True
        return cur_seq_specs


    def assign_anticodons(self, cur_seq_specs):
        """Takes a SeqSpecs class and tries to assign an anticodon to it"""
        if cur_seq_specs.full_length:
            anticodon = ",".join(self.extractor.extract_anticodon(cur_seq_specs.seq))
            cur_seq_specs.anticodon = anticodon
        else:
            anticodon = ",".join(self.extractor.extract_anticodon_not_full_length(cur_seq_specs.seq))
            cur_seq_specs.anticodon = anticodon
        return cur_seq_specs


    def split_3_trailer(self, cur_seq_specs, i):
        """Takes a SeqSpecs class and splits the 3-trailer from the seqience,
        returns an updated SeqSpecs class.
        """
        full_seq = cur_seq_specs.seq
        cur_seq_specs.seq = full_seq[:(cur_seq_specs.length - i)]
        cur_seq_specs.three_trailer = full_seq[(cur_seq_specs.length - i):]
        cur_seq_specs.length = len(cur_seq_specs.seq)
        cur_seq_specs.trailer_length = len(cur_seq_specs.three_trailer)
        return cur_seq_specs


    def check_seq_count(self, cur_seq_specs):
        """Checks the sequence in current SeqSpecs to see if the seq has
        already been encountered already, increments count if it has,
        adds it to the dict.
        """
        if cur_seq_specs.seq in self.seq_count_dict:
            self.seq_count_dict[cur_seq_specs.seq] += 1
        else:
            self.seq_count_dict[cur_seq_specs.seq] = 0

        # put the call to a hash function here so we can save the hashed seq in
        # cur_seq_specs
        return cur_seq_specs



    def handle_pass_seq(self, cur_seq_specs, i):
        """Consdolidates all the methods run specifically for a confirmed passed
        sequence.
        """
        self.sort_stats.total_passed += 1
        self.check_divergence_pos(cur_seq_specs)
        cur_seq_specs = self.split_3_trailer(cur_seq_specs, i)
        cur_seq_specs = self.check_full_length(cur_seq_specs)
        cur_seq_specs = self.assign_anticodons(cur_seq_specs)
        return cur_seq_specs


    def is_tRNA(self, seq):
        """Takes a sequence and determines whether or not it matches the
        criterion for being a tRNA
        """
        length = len(seq)
        sub_size = 24
        t_loop_error = True
        acceptor_error = True
        cur_seq_specs = SeqSpecs()


        # Start the sliding window at the last 24 bases, and move to the left
        # one at a time
        for i in xrange(length - sub_size + 1):
            sub_str = seq[-(i + sub_size):(length - i)]
            t_loop_seq = sub_str[0:9]
            acceptor_seq = sub_str[-3:]
            t_loop_dist = (lev.distance("GTTC", sub_str[0:4])
                + lev.distance("C", sub_str[8]))
            acceptor_dist = lev.distance("CCA", sub_str[-3:])
            mis_count = t_loop_dist + acceptor_dist

            if t_loop_dist < 1:
                t_loop_error = False
            else:
                t_loop_error = True
            if acceptor_dist < 1:
                acceptor_error = False
            else:
                acceptor_error = True
            if mis_count < cur_seq_specs.mis_count:
                cur_seq_specs.length = length
                cur_seq_specs.mis_count = mis_count
                cur_seq_specs.t_loop_error = t_loop_error
                cur_seq_specs.acceptor_error = acceptor_error
                cur_seq_specs.seq = seq
                cur_seq_specs.seq_sub = sub_str
                cur_seq_specs.t_loop_seq = t_loop_seq
                cur_seq_specs.acceptor_seq = acceptor_seq
            if mis_count < 2:
                cur_seq_specs = self.handle_pass_seq(cur_seq_specs, i)
                res_tup = (True, cur_seq_specs)
                return res_tup

        # Handles a failed sequence
        if cur_seq_specs.t_loop_error and cur_seq_specs.acceptor_error:
            if length < 24:
                self.sort_stats.short_rejected += 1
            else:
                self.sort_stats.both_rejected += 1
        elif cur_seq_specs.acceptor_error and not cur_seq_specs.t_loop_error:
            self.sort_stats.acceptor_seq_rejected += 1
        elif cur_seq_specs.t_loop_error and not cur_seq_specs.acceptor_error:
            self.sort_stats.t_loop_seq_rejected += 1
        self.sort_stats.total_rejected += 1
        res_tup = (False, cur_seq_specs)
        return res_tup


   #def fix_spacing_csv(self, out_tmp):
   #    """Fixes spacing and indentation on the csv files (_TAB_NO_TRAILER and
   #    _TAB_TRAILER)
   #    """
   #    # Declare and assign variables for output files
   #    tabfile = open(self.no_trailer_tabfile, "w")
   #    trailer_tabfile = open(self.trailer_tabfile, "w")
   #
   #    temp_tabfile_reader = csv.DictReader(out_tmp, delimiter="\t")
   #    tabfile_writer = csv.DictWriter(tabfile, fieldnames=self.fieldnames,
   #        delimiter="\t")
   #    trailer_tabfile_writer = csv.DictWriter(trailer_tabfile,
   #        fieldnames=self.fieldnames, delimiter="\t")
   #    tabfile_writer.writeheader()
   #    trailer_tabfile_writer.writeheader()
   #    trailer_count = 0

   #    for row in temp_tabfile_reader:
   #        row["Seq"] = ("-" * (self.max_seq_width - int(row["Seq_length"]))) + row["Seq"]
   #        if row["Trailer_length"] == "0":
   #            row["3-trailer"] = "_"
   #            tabfile_writer.writerow(row)
   #        else:
   #            trailer_tabfile_writer.writerow(row)
   #            trailer_count += 1

   #    tabfile.close()
   #    trailer_tabfile.close()
   #
   #    self.sort_stats.num_trailer = trailer_count


   #def write_to_outputs(self, spec_writer):
   #    """Writes the sort results to output files."""
   #    self.sort_stats.total_seqs += 1
   #    is_tRNA_result = self.is_tRNA(self.read_fasta.seq.upper())
   #    mod_id = is_tRNA_result[1].gen_id_string(
   #        self.read_fasta.id.split('|')[0])

   #    if is_tRNA_result[0]:
   #        self.passed_seqs_write_fasta.write_id(mod_id)
   #        self.passed_seqs_write_fasta.write_seq(is_tRNA_result[1].seq)
   #        is_tRNA_result[1].write_specs(spec_writer,
   #            self.read_fasta.id.split('|')[0])
   #        self.db.insert_seq(is_tRNA_result[1],
   #            self.read_fasta.id.split('|')[0])

   #        if len(is_tRNA_result[1].seq) > self.max_seq_width:
   #            self.max_seq_width = len(is_tRNA_result[1].seq)
   #    else:
   #        self.rejected_seqs_write_fasta.write_id(mod_id)
   #        self.rejected_seqs_write_fasta.write_seq(is_tRNA_result[1].seq)
   #

   #def write_sorted(self, readfile):
   #    """Reads in csv files and rewrites them, sorted by seq length."""
   #    sort_list = []
   #    with open(readfile, "r") as temp_tabfile:
   #        temp_tabfile_reader = csv.DictReader(temp_tabfile, delimiter="\t")
   #        for row in temp_tabfile_reader:
   #            sort_list.append(row)

   #    sort_list.sort(key=lambda dict: dict["Seq_length"], reverse=True)
   #
   #    with open(readfile, "w") as tabfile:
   #        tabfile_writer = csv.DictWriter(tabfile, fieldnames=self.fieldnames, delimiter="\t")
   #        tabfile_writer.writeheader()
   #        tabfile_writer.writerows(sort_list)

#    def add_to_database(self):
#        self.sort_stats.total_seqs += 1
#        is_tRNA_result = self.is_tRNA(self.read_fasta.seq.upper())
#
#        if is_tRNA_result[0]:
#            biglist.append((is_trna_result[1],self.read_fasta.id.split('|')[0]))
            #self.db.insert_seq(is_tRNA_result[1],
                #self.read_fasta.id.split('|')[0])


    def gen_sql_query_info_tuple(self):
        info_string_list = []
        info_string_list.append(self.total_seqs)
        info_string_list.append(self.total_rejected)
        return tuple(info_string_list)




    def run(self, args):
        biglist = []
        mem_const = 2000000 #rounding down 2 MB, should eventually make this an arg
        """Run the sorter."""
        print "sort started"
        self.set_file_names(args)
        self.db = dbops.tRNADatabase(self.tRNA_DB_file)



        while self.read_fasta.next():
            self.sort_stats.total_seqs += 1 #separating this from the original function call of add to database to avoid global lists
            is_tRNA_result = self.is_tRNA(self.read_fasta.seq.upper())

            if is_tRNA_result[0]:
                biglist.append((is_tRNA_result[1],self.read_fasta.id.split('|')[0]))

            if sys.getsizeof(biglist) > mem_const: #as long as the list isn't large, it'll keep populating
                for item in biglist:
                    a,b = item #unpacking our list, for the time being separate instead of a lambda
                    self.db.insert_seq(a,b)
                biglist = []

        for item in biglist: #final loop, adding all items to db
            a,b = item
            self.db.insert_seq(a,b)
        biglist = []
        





       #with tempfile.TemporaryFile() as out_tmp:
       #    spec_writer = csv.DictWriter(out_tmp, fieldnames=self.fieldnames,
       #        delimiter="\t")
       #    spec_writer.writeheader()
       #
       #    while self.read_fasta.next():
       #        self.write_to_outputs(spec_writer)
       #
       #    out_tmp.seek(0)
       #    self.fix_spacing_csv(out_tmp)

        print "finished preliminary tRNA sort"

        self.db.insert_stats(self.sort_stats)

        self.db.db.disconnect()

        # self.extractor.match_unassigned_sequences(self.no_trailer_tabfile,
        #     self.max_seq_width, self.fieldnames)

        # print "finished matching unassigned seqs"

       #if args.length_sort:
       #    self.write_sorted(self.no_trailer_tabfile)
       #    self.write_sorted(self.trailer_tabfile)

       #print "finished sorting"
       #
       #self.sort_stats.write_stats(self.sort_stats_write_file)
       #self.extractor.extractor_stats.write_stats(self.extractor.extractor_stats_file)
        print "sort finished"
