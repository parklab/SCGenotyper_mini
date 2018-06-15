"""
scGTmini.py -- class for genotyping positions in single-cell DNA-seq

Author:Maxwell A. Sherman
Contact: maxwell_sherman@hms.harvard.edu
Last revision: 15 June 2018
Version: 0.0.1
"""

import pandas as pd
import numpy as np
import pysam
import subprocess as sp
import scipy.stats as stats

class scGTmini(object):

    def __init__(self, bam, vcf):
        self.basename = bam.split("/")[-1]
        self.bam = pysam.AlignmentFile(bam)
        self.vcf = vcf

    @staticmethod
    def _read_PASS(read, mapq=5):
        """ Ensure a read is not dup, suppl / sec aln, qcfail
            Also require read to map and have mapped mate
            and exceed some minimum mapping quality
        """
        if read.is_duplicate or \
           read.is_qcfail or \
           read.is_secondary or \
           read.is_supplementary or \
           read.is_unmapped or \
           read.mate_is_unmapped or \
           not read.is_paired:
            return False

        elif read.mapq < mapq:
            return False

        else:
            return True

    @staticmethod
    def _kernel(x0, x1, l=10000):
        return 0.75 * (1 - ((x0 - x1) / l)**2)

    def extract_hetSNPs(self, region_str=None):
        """ Extract PHASED germline heterozygous SNPs from a VCF file

        NOTE: These positions will be used as training positions when estimating allele balance.

        Args:
            vcf: path to vcf file
        Kwargs:
            region_str: region to extract from as "chr:start-end"
        """
        fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
        inc = 'FORMAT/GT="0|1" || FORMAT/GT="1|0"'

        cmd = ['bcftools', 'query', '--format', fmt, '-i', inc]

        if region_str:
            cmd += ['-r', region_str]

        cmd += [self.vcf]
        p = sp.Popen(cmd, stdout=sp.PIPE)

        columns = ['CHR', 'POS', 'REF', 'ALT', 'GT'] 
        dtypes = {'CHR': str, 'POS': int, 'GT': str, 'REF': str, 'ALT': str}
        df = pd.read_table(p.stdout, names=columns, dtype=dtypes)
        # converters = {'REF': lambda x: int(x.replace('.', '0')), 'ALT': lambda x: int(x.replace('.', '0'))}
        # df = pd.read_table(p.stdout, names=columns, converters=converters, dtype=dtypes)

        # PHASE = df['GT'] == '1|0'
        # df['PAT'] = df['REF'] * np.logical_not(PHASE) + df['ALT'] * PHASE
        # df['MAT'] = df['REF'] * PHASE + df['ALT'] * np.logical_not(PHASE)
        # df['DP'] = df.REF + df.ALT

        return df

    def get_training_sites(self, df_SNPs, CHR, x, l=10000, mapq=5):
        """ Performs a pile-up at predefined training sites.

        The output is a matrix of SNP positions with counts of 
        REF, ALT, Maternal, and Paternal alleles
        
        Args:
            df_SNPs: dataframe of phased SNP sites (from extract_hetSNPs)
            CHR: chromosome (str)
            x: position (int)
            bam: pysam alignment file

        Kwargs:
            l: window size of informative training SNPs (int)
            mapq: minimum mapping quality to consider a read (int)
        """
        # snps = df_SNPs[(df_SNPs.CHR == CHR)]
        snps = df_SNPs[(df_SNPs.CHR == CHR) & (df_SNPs.POS > x - l) & (df_SNPs.POS < x + l)]

        l_ref, l_alt = [], []
        for i, row in snps.iterrows():
            pile = self.pileup_at_pos(str(row.CHR), row.POS, maq=mapq)
            l_ref.append(len(pile[pile.base == row.REF]))
            l_alt.append(len(pile[pile.base == row.ALT]))

        t_snps = snps.copy()
        t_snps['N_REF'] = l_ref
        t_snps['N_ALT'] = l_alt
        t_snps['DP'] = t_snps.N_REF + t_snps.N_ALT

        PHASE = t_snps['GT'] == '1|0'
        t_snps['PAT'] = t_snps['N_REF'] * np.logical_not(PHASE) + t_snps['N_ALT'] * PHASE
        t_snps['MAT'] = t_snps['N_REF'] * PHASE + t_snps['N_ALT'] * np.logical_not(PHASE)
        t_snps['AF_PAT'] = t_snps.PAT.astype(float) / t_snps.DP

        return t_snps

    def pileup_at_pos(self, CHR, x, maq=0):
        """ Create an mpileup from a bam at a specified position

        Args:
            CHR: chromosome (str)
            x: position (int)
            bam: pysam bam object

        Kwargs:
            mapq: minimum mapping quality to consider a read (int)
        """
        read_list = []
        for column in self.bam.pileup(CHR, x-1, x):

            # Skip to pileup column of SNP position
            if column.pos != x - 1:
                continue

            for plread in column.pileups:

                # Ensure SNP pos not missing from read
                if plread.is_del or plread.is_refskip:
                    continue

                read = plread.alignment

                if not self._read_PASS(read, mapq=maq):
                    continue

                name = read.query_name
                mapq = read.mapq
                base = read.query_sequence[plread.query_position]
                bqual = read.query_qualities[plread.query_position]

                read_props = [name, mapq, base, bqual]
                read_list.append(read_props)

        columns = ['name', 'mapq', 'base', 'bqual']
        df = pd.DataFrame(read_list, columns=columns)

        return df

    @staticmethod
    def predict_balance(xt, yt, dt, xp, l=10000):
        """ Predict balance at position using Nadaraya-Watson kernel

        Args:
            xt: training positions
            yt: (phased) allele count at training positions (e.g. Paternal allele count)
            dt: total depth at training positions
            xp: position at which to predict

        Kwargs:
            l: window size for kernel
        """

        k = scGTmini._kernel(xt, xp, l=l)
        k[k<0] = 0
        bal = sum(k * yt) / sum(k*dt)

        if bal == 0:
            bal = 0.5

        return bal

    @staticmethod
    def _logprobability_x(base, qual, ref, alt, f):
        """ Calculate the probability of observing a particular base at position x

        Possible events:
            - base is ref
            - base is alt
            - base is other

        Args:
            base: observed base at x
            qual: observed base quality
            ref: reference base at x
            alt: alternate base at x
            f: predicted balance
        """
        e = 10 ** (-qual / 10.)

        if base == ref:
            l = f * e/3 + (1-f) * (1-e)
        elif base == alt:
            l = f * (1-e) + (1-f) * e/3
        else:
            l = e/3

        return(np.log10(l))

    @staticmethod
    def loglikelihood_x(a_base, a_qual, ref, alt, f, p_art=0.125, p_err=0.03):
        """ Calculate the likelihood of a particular genotype at position x

        Possible genotypes:
            - RR: ref-ref
            - RB: ref-artifact
            - RA: ref-alt
            - AA: alt-alt

        Args:
            a_base: array of bases observed at x
            a_qual: array of base qualities observed at x
            ref: reference base at x
            alt: alternate base at x
            bal: predicted balance corrected for haplotype

        Kwargs:
            p_art: probability of amplification artifact
            p_err: sequencer error per base
        """
        RR = 0.0; RB = 0.0; RA = 0.0; AA = 0.0
        for b, q in zip(a_base, a_qual):
            RR = RR + scGTmini._logprobability_x(b, q, ref, alt, p_err)
            RB = RB + scGTmini._logprobability_x(b, q, ref, alt, p_art)
            RA = RA + scGTmini._logprobability_x(b, q, ref, alt, f)
            AA = AA + scGTmini._logprobability_x(b, q, ref, alt, 1.0-p_err);

        return([RR, RB, RA, AA])

    @staticmethod
    def _LLR(lls, labels):
        ll_h0 = np.min(lls)
        ll_h1 = np.max(lls)
        h1 = labels[np.argmax(lls)]

        LLR = -2 * (ll_h0 - ll_h1)

        return [h1, stats.chi2.sf(LLR, 1)]

    @staticmethod
    def genotyper_x(a_base, a_qual, ref, alt, bal, p_art=0.125, p_err=0.03):
        """ Predict the genotype at position x given observed bases and predicted balance

        Also decides if event is on psuedo-paternal or pseudo-maternal allele

        Possible genotypes:
            - RR: ref-ref
            - RB: ref-artifact
            - RA: ref-alt
            - AA: alt-alt

        Args:
            a_base: array of bases observed at x
            a_qual: array of base qualities observed at x
            ref: reference base at x
            alt: alternate base at x
            bal: predicted balance

        Kwargs:
            p_art: probability of amplification artifact
            p_err: sequencer error per base
        """
        nref = len(a_base[a_base == ref])
        nalt = len(a_base[a_base == alt])

        if not nref and not nalt:
            return 'Ukn', 0, 0, 'Ukn', 'NODP', np.nan, np.nan, [np.nan, np.nan, np.nan, np.nan]

        af = nalt / (nalt + nref)

        f = bal
        hap = 'PAT'
        if np.abs((1-bal) - af) <  np.abs(bal - af):
            f = 1 - bal
            hap = 'MAT'

        labels = np.array(['RR', 'RB', 'RA', 'AA'])
        RR, RB, RA, AA = scGTmini.loglikelihood_x(a_base, a_qual, ref, alt, f, p_art=p_art, p_err=p_err)
        ll = np.array([RR, RB, RA, AA]) 
        d_gt = dict(zip(labels, ll))
        GT = labels[ll.argmax()]
        GT2 = labels[ll.argsort()[-2]]

        # Probability of artifact
        p1 = scGTmini._LLR([d_gt['RB'], d_gt[GT]], ['RB', GT])[1]

        # Probability of misclassification 
        p2 = scGTmini._LLR([d_gt[GT2], d_gt[GT]], [GT2, GT])[1]

        # test = np.array([scGTmini._LLR([RB, RR], ['RB', 'RR']),
        #                  scGTmini._LLR([RB, RA], ['RB', 'RA']),
        #                  scGTmini._LLR([RB, AA], ['RB', 'AA']),
        #                 ], dtype='object')

        # test_RA = test[test[:, 0] == 'RB', :]
        # test_T = test[test[:, 0] != 'RB', :]

        # if test_T.any():
        #     GT, p = test_T[test_T[:, 1].argmin(), :]
        # else:
        #     GT, p = test_RA[test_RA[:, 1].argmax(), :]
        if p1 > 0.05 and GT != 'RR':
            FILTER = "ARTIF"
        elif p2 > 0.05:
            FILTER = 'LOWQUAL' 
        else:
            FILTER = 'PASS'

        return GT, nref, nalt, hap, FILTER, f, p1, [RR, RB, RA, AA]

    def caller(self, CHR, x, REF, ALT, l=10000, p_art=0.125, p_err=0.03, DP_min=0, drop_pos=True, mapq=5):
        """ Call a SNV at position x

        Args:
            CHR: chromosome (str)
            x: position (int)
            REF: reference base at x
            ALT: alternate base at x

        Kwargs:
            l: window size to infer balance
            p_art: probability of amplification artifact
            p_err: probability of sequencer error
            DP_min: minimum depth required at training position
            mapq: minimum mapping quality to consider a read
            drop_pos: drop x from training set if it exists (e.g. if testing at known het SNP)

        Returns:
            GT: genotype (str)
            p: p-value of genotype
            loglikelihoods of RR, RB, RA, AA (in that order)
        """
        region_str = "{}:{}-{}".format(CHR, int(x - l*5), int(x+l*5))
        df_SNP = self.extract_hetSNPs(region_str)
        df_pile = self.pileup_at_pos(CHR, x, maq=mapq)

        tsnp = self.get_training_sites(df_SNP, CHR, x, l=l, mapq=mapq)
        tsnp = tsnp[tsnp.DP > DP_min]
        if drop_pos:
            tsnp = tsnp[tsnp.POS != x]

        if len(tsnp) == 0:
            print('WARNING: no SNPs within {} of snp at {}:{} in {}'.format(l, CHR, x, self.basename))
            nref = len(df_pile.base[df_pile.base == REF])
            nalt = len(df_pile.base[df_pile.base == ALT])
            return 'Ukn', nref, nalt, 'Ukn', 'NOSNP', np.nan, np.nan, [np.nan, np.nan, np.nan, np.nan]

        bal = self.predict_balance(tsnp.POS, tsnp.PAT, tsnp.DP, x, l=l)
        GT, nref, nalt, hap, filt, bal, p, ll = self.genotyper_x(df_pile.base, df_pile.bqual, REF, ALT, bal, p_art=p_art, p_err=p_err)

        return GT, nref, nalt, hap, filt, bal, p, ll
