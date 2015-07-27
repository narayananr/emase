#!/usr/bin/env python
import os
import sys
import pysam
import numpy as np
from scipy.sparse import lil_matrix
from numpy.random import multinomial, randint, poisson, choice
import getopt


__author__ = 'Kwangbom "KB" Choi, Ph. D.'


help_message = '''
    Usage:
        simulate-reads -R <reference_dir> -H <hybrid_dir> -p <param_file> -m <simulation_model> -N <total_count> \\
                     [ -r <read_len> -e <err_rate> ]

    Input:
        <target_fasta> : A target fasta file
        <hap_list>    : Names of haplotypes to be used instead (comma delimited, in the order of genomes)
        <out_file>    : Output file name (default: './emase.pooled.targets.fa')

    Parameters:
        --help, -h : shows this help message
        --create-bowtie-index : builds bowtie1 index

'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def introduce_mutation(seq, pos, new_base):
    return seq[:pos] + new_base + seq[pos+1:]


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hR:H:p:m:N:r:e:", \
                                       ["help", "ref-dir", "hybrid-dir", "parameter-file", "total-count", "error-rate"])
        except getopt.error, msg:
            raise Usage(msg)

        # Default values of vars
        ref_dir = os.getenv('EMASE_INDEX_BASE', '.')
        hyb_dir = '.'
        read_len = 100
        error_rate = 0.1
        mut = {'A':('T', 'G', 'C'), \
               'T':('G', 'C', 'A'), \
               'G':('C', 'A', 'T'), \
               'C':('A', 'T', 'G')}

        # option processing (change this later with optparse)
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-R", "--ref-dir"):
                ref_dir = value
            if option in ("-H", "--hybrid-dir"):
                hyb_dir = value
            if option in ("-p", "--parameter-file"):
                param_file = value
            if option == "-m":
                model = int(value)
            if option in ("-N", "--total-count"):
                total_simreads = int(value)
            if option == "-r":
                read_len = int(value)
            if option in ("-e", "--error-rate"):
                error_rate = int(value)
            if option == "-o":
                outfile = value

        # Check if the required options are all given
        grouping_file = os.path.join(ref_dir, 'emase.gene2transcripts.tsv')
        if not os.path.exists(grouping_file):
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + '%s does not exist at %d.' % (grouping_file, ref_dir)
            return 2
        target_file = os.path.join(hyb_dir, 'emase.pooled.transcripts.fa')
        if not os.path.exists(target_file):
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + '%s does not exist at %d.' % (target_file, hyb_dir)
            return 2
        target_info_file = os.path.join(hyb_dir, 'emase.pooled.transcripts.info')
        if not os.path.exists(target_info_file):
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + '%s does not exist at %d.' % (target_info_file, hyb_dir)
            return 2
        if model not in (1, 2, 3, 4):
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + 'Simulation model should be either 1, 2, 3 or 4.'
            return 2
        outfile = './emase.M%d.simulated.fa' % model


        #
        # Main body
        #

        #
        # Load parameters and other required information
        tname = np.loadtxt(target_info_file, dtype='string', usecols=(0,))
        num_transcripts = len(tname)
        tid = dict(zip(tname, np.arange(num_transcripts)))

        gname = list()
        groups = list()
        with open(grouping_file) as fh:
            for curline in fh:
                item = curline.rstrip().split("\t")
                gname.append(item[0])
                tid_list = [ tid[t] for t in item[1:] ]
                groups.append(tid_list)
        gname = np.array(gname)
        num_genes = len(gname)
        gid = dict(zip(gname, np.arange(num_genes)))
        grp_conv_mat = lil_matrix((num_transcripts, num_genes))
        for i in xrange(num_genes):
            grp_conv_mat[groups[i], i] = 1.0
        grp_conv_mat = grp_conv_mat.tocsc()

        tacount = np.zeros((num_transcripts, 2))
        with open(param_file) as fh:
            curline = fh.next()
            item = curline.rstrip().split('\t')
            hname = item[1:]
            num_haps = len(hname)
            hid = dict(zip(hname, np.arange(num_haps)))
            for curline in fh:
                item = curline.rstrip().split('\t')
                tacount[tid[item[0]], :] = map(float, item[1:(1+num_haps)])
        tacount.shape, tacount.sum()

        trange = np.zeros(tacount.shape)
        with open(target_info_file) as fh:
            for curline in fh:
                item = curline.rstrip().split('\t')
                t, h = item[0].split('_')
                trange[tid[t], hid[h]] = np.int(item[1])
        trange = trange - read_len + 1

        #
        # Generate read counts
        if model == 1:
            gacount = grp_conv_mat.transpose() * tacount
            gcount = gacount.sum(axis=1)
            theta = gcount / gcount.sum()
            gcount_sim = multinomial(total_simreads, theta)[0]
            phi = np.zeros(gacount.shape)
            gexpr = gcount > 0
            phi[gexpr, :] = gacount[gexpr, :] / gcount[gexpr, np.newaxis]
            phi[np.logical_not(gexpr), :] = np.ones(num_haps) / num_haps
            gacount_sim = np.zeros(gacount.shape)
            for g in xrange(num_genes):
                gacount_sim[g] = multinomial(gcount_sim[g], phi[g])
            tacount_sim = np.zeros(tacount.shape)
            for g in xrange(num_genes):
                tindex = groups[g]
                num_isoforms = len(tindex)
                for h in xrange(num_haps):
                    delta = tacount[tindex, h]
                    delta_sum = delta.sum()
                    if delta_sum > 0:
                        delta /= delta_sum
                    else:
                        delta = np.ones(num_isoforms) / num_isoforms
                    tacount_sim[tindex, h] = multinomial(gacount_sim[g, h], delta)
        elif model == 2:
            pass
        elif model == 3:
            pass
        else:  # elif model == 4:
            pass

        out_dir = os.path.dirname(outfile)
        true_gaout = os.path.join(out_dir, 'emase.M%d.simulated.genes.read_counts' % model)
        np.savetxt(true_gaout, np.column_stack((gname, np.char.mod('%d', gacount_sim))), \
                   header='Gene_ID\t' + "\t".join(hname), fmt='%s', delimiter='\t')
        true_taout = os.path.join(out_dir, 'emase.M%d.simulated.transcripts.read_counts' % model)
        np.savetxt(true_taout, np.column_stack((tname, np.char.mod('%d', tacount_sim))), \
                   header='Transcript_ID\t' + "\t".join(hname), fmt='%s', delimiter='\t')

        #
        # Generate reads
        f = pysam.FastaFile(target_file)
        with open(outfile, 'w') as fhout:
            for t in xrange(num_transcripts):
                for h in xrange(num_haps):
                    theader = '%s_%s' % (tname[t], hname[h])
                    trange_max = trange[t, h]
                    for cnt in xrange(tacount_sim[t, h]):
                        start = randint(0, trange_max)
                        end = start + read_len
                        seq = f.fetch(theader, start, end)
                        nerr = poisson(error_rate)
                        errstr = ''
                        if nerr > 0:
                            errloc = sorted(choice(read_len, nerr, replace=False))
                            errstr = ':'
                            for epos in errloc:
                                from_base = seq[epos]
                                to_base = mut[from_base][randint(3)]
                                seq = introduce_mutation(seq, epos, to_base)
                                errstr += '[%d]%s->%s;' % (epos+1, from_base, to_base)
                        fhout.write(">%s:%d-%d%s\n%s\n" % (theader, start+1, end, errstr, seq))

        #
        # End of Main body
        #

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())