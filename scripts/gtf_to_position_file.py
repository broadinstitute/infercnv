#!/usr/bin/env python


"""
Converts GTF files to proprietary formats.
"""


# Import statements
import argparse
import csv
import gzip
import os


__author__ = 'Timothy Tickle, Itay Tirosh, Brian Haas'
__copyright__ = 'Copyright 2016'
__credits__ = ["Timothy Tickle"]
__license__ = 'BSD-3'
__maintainer__ = 'Timothy Tickle'
__email__ = 'ttickle@bbroadinstitute.org'
__status__ = 'Development'


def open_file(file_path):
    """ Open a file, handling gzip if necessary.

    :param file_path: Path to input file
    :type file_path: String

    :returns: File object
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')


def convert_to_positional_file(input_gtf, output_positional, attribute_key):
    """ Convert input GTF file to positional file.

    :param input_gtf: Path to input gtf file
    :type input_gtf: String
    :param output_positional: Path to output positional file
    :type output_positional: String
    :param attribute_key: Key of the GTF attribute to use for feature/row names
    :type attribute_key: String

    :returns: Indicator of success (True) or Failure (False)
    :rtype: boolean
    """

    if not input_gtf or not os.path.exists(input_gtf):
        print("".join(["gtf_to_position_file.py:: ",
                       "Could not find input file : " + input_gtf]))
        return False

    all_genes_found = set()

    # Holds lines to output after parsing.
    output_line = []
    previous_gene = None
    previous_chr = None
    gene_positions = []

    # Metrics for the file
    i_comments = 0
    i_duplicate_entries = 0
    i_entries = 0
    i_accepted_entries = 0
    i_written_lines = 0

    with open_file(input_gtf) as gtf:
        gtf_file = csv.reader(gtf, delimiter="\t")
        for gtf_line in gtf_file:
            if gtf_line[0][0] == "#":
                i_comments += 1
                continue
            i_entries += 1
            # Clean up the attribute keys and match the one of interest.
            attributes = gtf_line[8].split(";")
            attributes = [entry.strip(" ") for entry in attributes]
            attributes = [entry.split(" ") for entry in attributes if entry]
            attributes = [[entry[0].strip('"'), entry[1].strip('"')] for entry in attributes]
            attributes = dict([[entry[0].split("|")[0], entry[1]] for entry in attributes])
            if attribute_key in attributes:
                gene_name = attributes[attribute_key]
            else:
                print("Could not find an attribute in the GTF with the name '" + attribute_key + "'. Line=" + "\t".join(gtf_line))
                exit(99)
            if not gene_name == previous_gene:
                if len(gene_positions) > 1 and previous_gene not in all_genes_found:
                    i_accepted_entries += 1
                    gene_positions.sort()
                    output_line.append("\t".join([previous_gene,
                                                  previous_chr,
                                                  str(gene_positions[0]),
                                                  str(gene_positions[-1])]))
                    all_genes_found.add(previous_gene)
                gene_positions = []
            else:
                i_duplicate_entries += 1
            gene_positions += [int(gtf_line[3]), int(gtf_line[4])]
            previous_gene = gene_name
            previous_chr = gtf_line[0]
        if previous_gene and previous_chr and len(gene_positions) > 1:
            i_accepted_entries += 1
            gene_positions.sort()
            output_line.append("\t".join([previous_gene,
                                          previous_chr,
                                          str(gene_positions[0]),
                                          str(gene_positions[-1])]))

    with open(output_positional, "w") as positional_file:
        i_written_lines += len(output_line)
        positional_file.write("\n".join(output_line))

    # Print metrics
    print("Number of lines read: " + str(i_entries))
    print("Number of comments: " + str(i_comments))
    print("Number of entries: " + str(i_accepted_entries))
    print("Number of duplicate entries: " + str(i_duplicate_entries))
    print("Number of entries written: " + str(i_written_lines))


if __name__ == "__main__":

    # Parse arguments
    prsr_arguments = argparse.ArgumentParser(prog='gtf_to_position_file.py',
                                             description='Convert a GTF file to a positional file.',
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Add positional argument
    prsr_arguments.add_argument("input_gtf",
                                metavar="input_gtf",
                                help="Path to the input GTF file.")
    prsr_arguments.add_argument("--attribute_name",
                                metavar="attribute_name",
                                default="gene_id",
                                help="The name of the attribute in the GTF attributes to use instead of gene name, for example 'gene_name' or 'transcript_id'.")
    prsr_arguments.add_argument("output_positional",
                                metavar="output_positional",
                                help="Path for the output positional file.")
    args = prsr_arguments.parse_args()

    # Run Script
    convert_to_positional_file(args.input_gtf, args.output_positional, args.attribute_name)
