import pandas as pd
from Bio import SeqIO
import argparse


def taxa_date(fasta, priors_table):
    """
    This function is used for input the date of taxa.
    The input file is the fasta with all taxa.
    :param fasta: aligned fasta file
    :return: A string of taxa data in xml format
    """
    global taxa_list, taxa_list_ND

    read_tip_table(priors_table)

    with open(fasta, 'r') as fasta:
        taxa_list = []
        taxa_date_list = []
        taxa_list_ND = []

        for line in fasta:
            if line.startswith('>'):
                taxa_name = str(line).strip('>').strip('\n')
                taxa_list.append(taxa_name)
                taxa_date = taxa_name.split('_')[-1]
                if taxa_date == 'ND':
                    taxa_list_ND.append(taxa_name)
                    taxa_date = tip_date_dict[taxa_name]
                taxa_date_list.append(str(taxa_date))

    add_date = []

    for i in range(0, len(taxa_list)):
        line_1 = '\t\t<taxon id="{0}">\n\t\t\t<date value="{1}" direction="backwards" units="years"/>\n\t\t</taxon>\n' \
            .format(taxa_list[i], taxa_date_list[i])
        add_date.append(line_1)

    return ''.join(add_date)


def add_len_taxa(taxa_list):
    """
    Calculate how many samples are in this analysis.
    :param taxa_list: a list of taxa from the fasta file
    :return: the number of all taxa
    """
    ntax = '\t<!-- ntax={}                                                                -->\n'.format(len(taxa_list))
    return ntax


def taxon_set(taxa_list):
    """
    This function is used for input taxon sets.
    Search keywords like Mammuthus: M., Elephas: E., Undated: ND
    :param taxa_list: the list of taxa from fasta
    :return: A string of taxon sets in xml format
    """
    global add_taxonset

    taxon_set = ['\n\t<taxa id="Loxodonta">\n',
                 '\n\t<taxa id="Mammuthus">\n',
                 '\n\t<taxa id="Mammuthus_elephas">\n',
                 '\n\t<taxa id="Tip">\n']
    taxon_key = [['L.', 'P.'],
                 ['M.'],
                 ['M.', 'E.'],
                 taxa_list_ND]
    taxon_dict = dict(zip(taxon_set, taxon_key))

    add_taxonset = []
    for item in taxon_set:
        taxon_list = [item]
        for i in range(0, len(taxa_list)):
            for x in range(0, len(taxon_dict[item])):
                if taxon_dict[item][x] in taxa_list[i]:
                    taxon = taxa_list[i]
                    taxon_list.append('\t\t<taxon idref="{}"/>\n'.format(taxon))
        taxon_list.append('\t</taxa>\n')
        for y in range(0, len(taxon_list)):
            add_taxonset.append(taxon_list[y])
    return ''.join(add_taxonset)


def find_partition(gff):
    """
    This function used a gff annotation file as input, and store them in lists.
    :param gff: The annotation gff file which match the alignment
    :return: Lists of sequence types and their starts and ends
    """
    global seq_type, start, end, df3

    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                # print(str(line).strip('\n'))
            else:
                gene_future = pd.read_table(f, header=None)
                gene_future.columns = ['genome_id', 'source',
                                       'type', 'start', 'end',
                                       'uk', 'uk1', 'uk1', 'more_info']
                df1 = gene_future.loc[:, ['type', 'start']]
                df2 = gene_future.loc[:, ['type', 'end']]
                df3 = gene_future.loc[:, ['type', 'genome_id']]

                seq_type = df3['type'].values.tolist()
                start = df1['start'].values.tolist()
                end = df2['end'].values.tolist()
    return seq_type, start, end


def read_partition(fasta, gff, partition_name):
    """
    This function is used for reading partition sequences, and filter out overlapped sequences.
    :param fasta: The aligned fasta file with all taxa
    :param gff: The annotation gff file fitting to the alignment
    :param partition_name: The partition you want to read, e.g. rRNA, tRNA, CDS...
    :return: A string contains sequences from the partition in xml format
    """
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    find_partition(gff)
    add_alignment = []
    p_start = []
    p_end = []
    for a in range(0, len(seq_type)):
        if seq_type[a] == partition_name:
            p_start.append(start[a])
            p_end.append(end[a])
    if len(p_end) == 1:
        s = [p_start[0] - 1]
        e = [p_end[0]]
    else:
        s = [p_start[0] - 1]
        e = [p_end[0]]
        for i in range(1, len(p_end)):
            if p_start[i] > p_end[i - 1]:
                s.append(p_start[i] - 1)
                e.append(p_end[i])
            else:
                s.pop()
                e.pop()
                s.append(p_start[i - 1] - 1)
                e.append(p_end[i])

    for taxa in taxa_list:
        p_list = []
        for i in range(0, len(s)):
            s_num = s[i]
            e_num = e[i]
            partition = str(fa_dict[taxa][s_num:e_num].seq).upper()
            p_list.append(partition)
        alignment = '\t\t<sequence>\n\t\t\t<taxon idref="{0}"/>\n\t\t\t{1}\n\t\t</sequence>\n' \
            .format(taxa, ''.join(p_list))
        add_alignment.append(alignment)
    return ''.join(add_alignment)


def VNTR_excluded(fasta, gff, partition_name):
    """
    An extension of read_partition function. It excludes VNTR region from D_loop.
    :param fasta: The aligned fasta file with all taxa
    :param gff: The annotation gff file fitting to the alignment
    :param partition_name: The D_loop partition
    :return: A string contains sequences from the partition in xml format
    """
    fa_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    find_partition(gff)
    p_start = []
    p_end = []
    add_alignment = []
    for a in range(0, len(seq_type)):
        if seq_type[a] == partition_name:
            p_start = int(start[a] - 1)
            p_end = int(end[a])
        if seq_type[a] == 'VNTR':
            s_VNTR = int(start[a])
            e_VNTR = int(end[a] - 1)
    s = [p_start, e_VNTR]
    e = [s_VNTR, p_end]
    for taxa in taxa_list:
        p_list = []
        for i in range(0, len(s)):
            s_num = s[i]
            e_num = e[i]
            partition = str(fa_dict[taxa][s_num:e_num].seq).upper()
            p_list.append(partition)
        alignment = '\t\t<sequence>\n' \
                    '\t\t\t<taxon idref="{0}"/>\n' \
                    '\t\t\t{1}\n' \
                    '\t\t</sequence>\n' \
            .format(taxa, ''.join(p_list))
        add_alignment.append(alignment)
    return ''.join(add_alignment)


def tip_sampling(taxa_list_ND):
    """
    This function is used for input the index of tip samples. Write XML format.
    :param taxa_list_ND: A list of undated samples from taxa_date function
    :return: A string of a list of undated samples in xml format
    """
    add_tip_sample = ['\t\t</nodeHeights>\n\n'
                      '\t\t<!-- START Tip date sampling                                                 -->\n']
    for taxa in taxa_list_ND:
        tip_sample = '\t\t<leafHeight taxon="{0}">\n' \
                     '\t\t\t<parameter id="age({1})"/>\n' \
                     '\t\t</leafHeight>\n'.format(taxa, taxa)
        add_tip_sample.append(tip_sample)
    return ''.join(add_tip_sample)


def tip_operator(taxa_list_ND):
    """This function is used for input operator of tipping samples. Just write XML format.
    :param taxa_list_ND: A list of undated samples from taxa_date function
    :return: A string of a list of undated samples in xml format
    """
    add_tip_operator = ['\t\t<scaleOperator scaleFactor="0.75" weight="1">\n'
                        '\t\t\t<parameter idref="skygrid.precision"/>\n'
                        '\t\t</scaleOperator>\n']
    for taxa in taxa_list_ND:
        tip_operator = '\t\t<uniformOperator weight="2">\n' \
                       '\t\t\t<parameter idref="age({0})"/>\n' \
                       '\t\t</uniformOperator>\n'.format(taxa)
        add_tip_operator.append(tip_operator)
    return ''.join(add_tip_operator)


def MCMC(chain):
    """This function is used to define the MCMC chain.
    :param chain: The chain length of MCMC.
    :return: A string of MCMC definition in xml format.
    """
    MCMC = '\t<mcmc id="mcmc" chainLength="{0}" autoOptimize="true">\n'.format(chain)
    return MCMC


def read_tip_table(priors_table):
    """
    To read the priors table. Input file is a table of priors (csv file, split by comma)
    :param priors_table: The priors table in csv format.
    (columns: taxa name, input date, prior type, mu/lower, sigma/upper, offset)
    :return: A dictionary contains taxa and its priors settings.
    """
    global tip_taxa, tip_date_dict, tip_priors_dict

    with open(priors_table, 'r') as priors:
        priors_table = pd.read_csv(priors, header=0)
        tip_taxa = priors_table['taxa'].values.tolist()
        tip_date = priors_table['date'].values.tolist()
        tip_prior = priors_table['prior'].values.tolist()
        tip_mu_lower = priors_table['mu/lower'].values.tolist()
        tip_sigma_upper = priors_table['sigma/upper'].values.tolist()
        tip_offset = priors_table['offset'].values.tolist()

        tip_date_dict = {}
        tip_date_dict = dict(zip(tip_taxa, tip_date))
        tip_priors_dict = {}
        for i in range(0, len(tip_taxa)):
            priors = [tip_prior[i], tip_mu_lower[i], tip_sigma_upper[i], tip_offset[i]]
            tip_priors_dict.update({tip_taxa[i]: priors})
    return tip_priors_dict


def tip_priors(priors_table):
    """This function is used for input the dating priors.
    :param priors_table: The priors table in csv format.
    :return: A string of priors settings in xml format.
    """
    read_tip_table(priors_table)

    add_tip_line = []
    for taxa in taxa_list_ND:
        if tip_priors_dict[taxa][0] == 'logNormal':
            tip_line = '\t\t\t\t<{0}Prior mu="{1}" sigma="{2}" offset="{3}">\n' \
                       '\t\t\t\t\t<parameter idref="age({4})"/>\n' \
                       '\t\t\t\t</{5}Prior>\n'.format(tip_priors_dict[taxa][0],
                                                      tip_priors_dict[taxa][1],
                                                      tip_priors_dict[taxa][2],
                                                      tip_priors_dict[taxa][3],
                                                      taxa,
                                                      tip_priors_dict[taxa][0])
            add_tip_line.append(tip_line)
        elif tip_priors_dict[taxa][0] == 'uniform':
            tip_line = '\t\t\t\t<{0}Prior lower="{1}" upper="{2}" offset="{3}">\n' \
                       '\t\t\t\t\t<parameter idref="age({4})"/>\n' \
                       '\t\t\t\t</{5}Prior>\n'.format(tip_priors_dict[taxa][0],
                                                      tip_priors_dict[taxa][1],
                                                      tip_priors_dict[taxa][2],
                                                      tip_priors_dict[taxa][3],
                                                      taxa,
                                                      tip_priors_dict[taxa][0])
            add_tip_line.append(tip_line)
        elif tip_priors_dict[taxa][0] == 'normal':
            tip_line = '\t\t\t\t<{0}Prior mean="{1}" stdev="{2}">\n' \
                       '\t\t\t\t\t<parameter idref="age({3})"/>\n' \
                       '\t\t\t\t</{4}Prior>\n'.format(tip_priors_dict[taxa][0],
                                                      tip_priors_dict[taxa][1],
                                                      tip_priors_dict[taxa][2],
                                                      taxa,
                                                      tip_priors_dict[taxa][0])
            add_tip_line.append(tip_line)
    # print(''.join(add_tip_line))
    return ''.join(add_tip_line)


def log_every_screen(log_every):
    """
    This function defines the log numbers
    :param log_every: Writing the log every how many.
    :return: A string of definition of log setting.
    """
    add_log_screen = '\t\t<log id="screenLog" logEvery="{}">\n'.format(log_every)
    return add_log_screen


def log_every_file(log_every, filename):
    """This function defines the log numbers and log file name.
    :param log_every: Writing the log every how many.
    :param filename: The file name of the log file.
    :return: A string of definition of log setting.
    """
    add_log_file = '\t\t<log id="fileLog" logEvery="{0}" fileName="{1}.log" overwrite="false">\n'.format(log_every, filename)
    return add_log_file


def tree_log(log_every, filename):
    """
    This function defines the log numbers and log file name.
    :param log_every: Writing the log every how many.
    :param filename: The file name of the tree file.
    :return: A string of definition of tree log setting.
    """
    tree_log = '\t\t<logTree id="treeFileLog" logEvery="{0}" nexusFormat="true" fileName="{1}.trees" ' \
                   'sortTranslationTable="true">\n'.format(log_every, filename)
    return tree_log


def start_dating(taxa_list_ND):
    """This function is to write XML format to start dating.
    :param taxa_list_ND: A list of undated taxa from fasta file.
    :return: A string of a list of undated taxa in xml format.
    """
    add_dating = []
    for taxa in taxa_list_ND:
        dating = '\t\t\t<parameter idref="age({})"/>\n'.format(taxa)
        add_dating.append(dating)
    # print(''.join(add_dating))
    return ''.join(add_dating)


def find_write(template, output, find_list, add_list):
    """
    Write a new xml file based on the template.
    :param template: The template xml file.
    :param output: The output xml file name.
    :param find_list: A list of marker to write.
    :param add_list: The taxa information needed to write in new file.
    :return: A xml file ready to run in BEAST.
    """
    inputfile = open(template, 'r').readlines()
    write_file = open(output, 'w')

    for line in inputfile:
        write_file.write(line)
        for i in range(0, len(find_list)):
            if find_list[i] in line:
                new_line = str(add_list[i])
                write_file.write(new_line)
            else:
                continue
    write_file.close()
    return output


def writeXML(template, output, fasta, gff, priors_table, MCMC_chain, log_every):
    """
    Write a new xml file based on the template.
    :param template: The template xml file.
    :param output: The output xml file name.

    :return: A xml file ready to run in BEAST.
    """

    find_list = ['<taxa id="taxa">',
                 '<!-- The list of taxa to be analysed (can also include dates/ages).',
                 '</taxa>',
                 '<alignment id="alignment1" dataType="nucleotide">',
                 '<alignment id="alignment2" dataType="nucleotide">',
                 '<alignment id="alignment3" dataType="nucleotide">',
                 '<alignment id="alignment4" dataType="nucleotide">',
                 '<parameter id="treeModel.allInternalNodeHeights"/>',
                 '</gmrfGridBlockUpdateOperator>',
                 '!-- Define MCMC',
                 '</gammaPrior>',
                 '<!-- write log to screen',
                 '<!-- write log to file',
                 '<!-- START Tip date sampling',
                 '<!-- write tree log to file']

    add_list = [taxa_date(fasta, priors_table),
                add_len_taxa(taxa_list),
                taxon_set(taxa_list),
                read_partition(fasta, gff, partition_name='tRNA'),
                read_partition(fasta, gff, partition_name='rRNA'),
                read_partition(fasta, gff, partition_name='CDS'),
                VNTR_excluded(fasta, gff, partition_name='D_loop'),
                tip_sampling(taxa_list_ND),
                tip_operator(taxa_list_ND),
                MCMC(MCMC_chain),
                tip_priors(priors_table),
                log_every_screen(log_every),
                log_every_file(log_every, str(output).rstrip('.xml')),
                start_dating(taxa_list_ND),
                tree_log(log_every, str(output).rstrip('.xml'))]

    find_write(template, output, find_list, add_list)
    return output


parser = argparse.ArgumentParser(
                    prog='beast1XMLgenerator',
                    description='Generating XML file for BEAST1',
                    epilog='Author: Wenxi Li')

parser.add_argument('-o', '--output', type=str, help='The output XML file.')
parser.add_argument('-t', '--template', type=str, help='The XML template file.')
parser.add_argument('-f', '--fasta', type=str, help='The aligned fasta file.')
parser.add_argument('-g', '--gff', type=str, help='The gff annotation file.')
parser.add_argument('-p', '--priors', type=str, help='The priors table in csv format.')
parser.add_argument('-m', '--mcmc', help='The MCMC chain length')
parser.add_argument('-l', '--log', help='Write the log file.')

args = parser.parse_args()

xml = writeXML(
    template=args.template, output=args.output,
    fasta=args.fasta, gff=args.gff, priors_table=args.priors,
    MCMC_chain=args.mcmc, log_every=args.log)

print('The output file is {}'.format(xml))
