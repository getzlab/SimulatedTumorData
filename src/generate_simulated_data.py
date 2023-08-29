# Generate patient data
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from AnnoMate.AppComponents.utils import cluster_color
from AnnoMate.AppComponents.PhylogicComponents import gen_stylesheet
import dash_cytoscape as cyto
from dash import callback, Input, Output, Dash, html
import pickle
import os
import random
import pysam
import re
import math
import scipy
import tqdm

# from https://github.com/getzlab/cnv_suite
from cnv_suite import simulate
from cnv_suite import visualize
from cnv_suite.utils.simulation_utils import get_contigs_from_header, single_allele_ploidy, get_average_ploidy, get_alt_count
from cnv_suite.utils import switch_contigs
from cnv_suite.visualize import add_background
from natsort import natsorted

# add methods to simulate.CNV_Profile().phylogeny to set the phylogeny manually


def read_target_interval_file(target_interval_file):
    target_intervals_df = pd.read_csv(
        target_interval_file,sep='\t', comment='@', 
        names=['chrom', 'start', 'stop', 'plus','target_name']
    )
    target_intervals_df = target_intervals_df[(target_intervals_df.chrom!='MT') & 
                                              (target_intervals_df.chrom!='Y') & 
                                              (target_intervals_df.chrom!='GL000228.1')]
    target_intervals_df = target_intervals_df.rename(columns={'stop':'end'})
    target_intervals_df['chrom'] = target_intervals_df['chrom'].apply(lambda x: '23' if x=='X' else str(x))
    target_intervals_df['interval_length']=target_intervals_df['end']-target_intervals_df['start']
    total_intervals_length = target_intervals_df['interval_length'].sum()
    target_intervals_df['interval_weight'] = target_intervals_df['interval_length']/total_intervals_length
    target_intervals_df['chrom'] = target_intervals_df['chrom'].astype(str)
    target_intervals_df['start'] = target_intervals_df['start'].astype(int)
    target_intervals_df['end'] = target_intervals_df['end'].astype(int)
    return target_intervals_df

def set_phylogeny(self, phylogeny_dict):
    """
    Custom set_phylogeny
    """
    self.num_subclones = phylogeny_dict['num_subclones']
    self.parents = phylogeny_dict['parents']
    self.ccfs = phylogeny_dict['ccfs']
    
def get_children(self, node):
    """Return children for the specified clone
        
    :param node: index of desired clone
    :return: (dict) representing the child clones and their respective CCFs"""

    cluster_list=[]
    for k, v in self.parents.items():
        if v==node:
            cluster_list.append(k)

    return cluster_list
    
simulate.cnv_profile.Phylogeny.set_phylogeny = set_phylogeny
simulate.cnv_profile.Phylogeny.get_children = get_children

# add methods to simulate.CNV_Profile().CNV_Profile to get and set the cnv events manually

def set_events_and_phylogeny(self, event_trees, phylogeny):
    self.event_trees = event_trees
    self.phylogeny =  phylogeny
    
def pickle_events(self, out_file_name):
    event_trees = self.event_trees
    
    # open a file, where you ant to store the data
    file = open(out_file_name, 'wb')
    pickle.dump(event_trees, file)
    file.close()

# Copied over and force chromosome type to be string for regex in switch_contigs to work
def generate_snvs(self, vcf, bed, purity, ref_alt=False, do_parallel=True):
    """Generate SNV read depths adjusted for CNV profile (and purity), with phasing from vcf file.

    :param vcf: VCF file containing SNVs and haplotype of SNVs
    :param bed: bed file containing the read depths for all desired SNVs in original bam
    :param purity: desired purity, given as float
    :param do_parallel: boolean option to parallelize with pandarallel
    :param ref_alt: True if bed file contains ref and alt counts vs. only depth counts (default False)
    """
    if self.cnv_trees is None:
        print('cnv_trees not computed yet. Run calculate_profiles() before generating snvs.')
        return None, None

    # check if VCF contigs given in header match contigs and lengths in self
    vcf_contigs = switch_contigs(get_contigs_from_header(vcf))
    vcf_contigs_pertinent = {k: v for k, v in vcf_contigs.items() if k in self.csize.keys()}
    if vcf_contigs_pertinent.keys() != self.csize.keys():
        print(f'WARNING: Not all defined contigs exist in VCF file. '
              f'Missing contigs: {set(self.csize.keys()) - set(vcf_contigs_pertinent.keys())}')
    for k, v in vcf_contigs_pertinent.items():
        if v != self.csize[k]:
            print(f'WARNING: Contig length for chrom {k} in VCF file does not match CNV Profile '
                  f'({v} vs. {self.csize[k]}).')

    snv_df = pd.read_csv(vcf, sep='\t', comment='#', header=None, 
                 names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NA12878'])

    if ref_alt:
        bed_df = pd.read_csv(bed, sep='\t', header=0, names=['CHROM', 'POS', 'REF_BED', 'ALT_BED'], dtype={'CHROM': str})
        bed_df['DEPTH'] = bed_df['REF_BED'] + bed_df['ALT_BED']
    elif isinstance(bed, str):
        bed_df = pd.read_csv(bed, sep='\t', header=0, names=['CHROM', 'POS', 'DEPTH'], dtype={'CHROM': str})
    elif isinstance(bed, pd.DataFrame):
        bed_df = bed[['CHROM', 'POS', 'DEPTH']]

    # Force string type for switch contigs
    snv_df['CHROM'] = snv_df['CHROM'].astype(str)
    bed_df['CHROM'] = bed_df['CHROM'].astype(str)

    # change contigs to [0-9]+ from chr[0-9XY]+ in input files
    snv_df = switch_contigs(snv_df)
    bed_df = switch_contigs(bed_df)

    snv_df = snv_df.merge(bed_df, on=['CHROM', 'POS'], how='inner')
    
    if do_parallel:
        pandarallel.initialize()
        snv_df['paternal_ploidy'] = snv_df.parallel_apply(
            lambda x: single_allele_ploidy(self.cnv_trees[x['CHROM']][0], x['POS'], x['POS'] + 1),
            axis=1)
        snv_df['maternal_ploidy'] = snv_df.parallel_apply(
            lambda x: single_allele_ploidy(self.cnv_trees[x['CHROM']][1], x['POS'], x['POS'] + 1),
            axis=1)
    else:
        snv_df['paternal_ploidy'] = snv_df.apply(
            lambda x: single_allele_ploidy(self.cnv_trees[x['CHROM']][0], x['POS'], x['POS'] + 1),
            axis=1)
        snv_df['maternal_ploidy'] = snv_df.apply(
            lambda x: single_allele_ploidy(self.cnv_trees[x['CHROM']][1], x['POS'], x['POS'] + 1),
            axis=1)

    snv_df['ploidy'] = get_average_ploidy(snv_df['paternal_ploidy'].values,
                                          snv_df['maternal_ploidy'].values,
                                          purity)

    snv_df['maternal_prop'] = ( snv_df['maternal_ploidy'].values * purity + (1 - purity) ) / snv_df['ploidy'].values

    snv_df['paternal_prop'] = ( snv_df['paternal_ploidy'].values * purity + (1 - purity) ) / snv_df['ploidy'].values

    snv_df['maternal_present'] = snv_df['NA12878'].apply(lambda x: x[0] == '1')
    snv_df['paternal_present'] = snv_df['NA12878'].apply(lambda x: x[2] == '1')

    snv_df['adjusted_depth'] = np.floor(snv_df['DEPTH'].values * snv_df['ploidy'].values / 2).astype(int)
    
    # generate phase switch profile
    # chromosome interval trees: False if phase switched
    correct_phase_interval_trees = self.generate_phase_switching()

    # calculate alt counts for each SNV
    if do_parallel:
        snv_df['alt_count'] = snv_df.parallel_apply(
            lambda x: get_alt_count(x['maternal_prop'], x['paternal_prop'], x['maternal_present'],
                                    x['paternal_present'], x['adjusted_depth'],
                                    correct_phase_interval_trees[x['CHROM']][x['POS']].pop().data), axis=1)
    else:
         snv_df['alt_count'] = snv_df.apply(
            lambda x: get_alt_count(x['maternal_prop'], x['paternal_prop'], x['maternal_present'],
                                    x['paternal_present'], x['adjusted_depth'],
                                    correct_phase_interval_trees[x['CHROM']][x['POS']].pop().data), axis=1)
    snv_df['ref_count'] = snv_df['adjusted_depth'] - snv_df['alt_count']

    return snv_df, correct_phase_interval_trees
    
simulate.cnv_profile.CNV_Profile.pickle_events = pickle_events
simulate.cnv_profile.CNV_Profile.set_events_and_phylogeny = set_events_and_phylogeny
simulate.cnv_profile.CNV_Profile.generate_snvs = generate_snvs


class Sample:
    def __init__(self, name, time_point, purity, cov_profile, parent_children_clone_ratios):
        self.name = name
        self.time_point = time_point
        self.purity = purity
        self.cov_profile = cov_profile 
        self.parent_children_clone_ratios = parent_children_clone_ratios

    def calc_clone_ccfs(self, parent_2_children_tree):
        clone_ccfs = {}
        for parent, children in parent_2_children_tree.items():
            if parent is None:
                clone_ccfs[children[0]] = 1.0 # clone 1
            else:
                parent_ccf = clone_ccfs[parent]
                ccfs_ratios = self.parent_children_clone_ratios[parent]
                ccfs_proportions = ccfs_ratios / np.sum(ccfs_ratios)
                ccfs = ccfs_proportions * parent_ccf
                for i in range(len(children)):
                    clone_ccfs[children[i]] = ccfs[i + 1]

        self.clone_ccfs = clone_ccfs
        return self.clone_ccfs

    def set_cnv_profile(self, num_subclones, parents, event_trees, data_directory):
        
        phylogeny = simulate.CNV_Profile().phylogeny
        phy_dict = {
            'num_subclones': num_subclones,
            'parents': parents #{1: None, 2:1, 3:1, 4:2, 5:4, 6:4, 7:3}
        }
        phy_dict['ccfs'] = self.clone_ccfs
        phylogeny.set_phylogeny(phy_dict)

        cnv_profile = simulate.CNV_Profile()
        cnv_profile.set_events_and_phylogeny(event_trees, phylogeny)
        cnv_profile.calculate_profiles()
        
        cnv_profile.cnv_profile_df['length'] = cnv_profile.cnv_profile_df['End.bp'] - cnv_profile.cnv_profile_df['Start.bp']
        
        cnv_profile.cnv_profile_df['n_probes'] = (cnv_profile.cnv_profile_df['length'] * 0.00005).astype(int)
        cnv_profile.cnv_profile_df['n_hets'] = (cnv_profile.cnv_profile_df['length'] * 0.000001).astype(int)
        cnv_profile.cnv_profile_df['tau'] = cnv_profile.cnv_profile_df['mu.major'] + cnv_profile.cnv_profile_df['mu.minor']
        
        cnv_profile.cnv_profile_df['sigma.minor'] = np.random.exponential(0.01, size=len(cnv_profile.cnv_profile_df))
        cnv_profile.cnv_profile_df['sigma.major'] = cnv_profile.cnv_profile_df['sigma.minor']
        self.cnv_profile = cnv_profile
        
        self.absolute_cnv_profile_fn = f'{data_directory}/{self.name}.simulated_seg.tsv'
        cnv_profile.cnv_profile_df.to_csv(self.absolute_cnv_profile_fn, index=False, sep='\t')

    def set_acr_profile(self):
        # Allelic copy ratio seg file
        acr_cnv_profile_df = self.cnv_profile.cnv_profile_df.copy()

        acr_cnv_profile_df['mu.minor'] = acr_cnv_profile_df['mu.minor'] * ccf * self.purity + (1 - self.purity) / 2
        acr_cnv_profile_df['mu.major'] = acr_cnv_profile_df['mu.major'] * ccf * self.purity + (1 - self.purity) / 2

    def plot_cn_profile(self):
        fig, ax = plt.subplots(1,1, figsize=(10, 4))
        fig = visualize.plot_cnv_profile.plot_acr_static(
            self.cnv_profile.cnv_profile_df, 
            ax, self.cnv_profile.csize, #cnv_profile.csize, 
            segment_colors='difference', 
            sigmas=False, 
            min_seg_lw=3, 
            y_upper_lim=3
        )

    def sim_coverage(self, target_intervals_df, output_folder, sigma, x_coverage, override=False):

        self.binned_coverage_fn = f'{output_folder}/{self.name}.binned_coverage.tsv'
        self.corrected_binned_coverage_fn = f'{output_folder}/{self.name}.corrected_binned_coverage.tsv'

        if os.path.exists(self.binned_coverage_fn) and not override:
            print(f'{self.binned_coverage_fn} already exists.')
        else:
            binned_coverage_df = target_intervals_df[['chrom','start','end']].copy()
            dispersion_norm = np.random.normal(0, sigma, binned_coverage_df.shape[0])
            binned_coverage = x_coverage * (binned_coverage_df['end'] - binned_coverage_df['start'])
            this_chr_coverage = np.asarray([np.random.poisson(cov + np.exp(disp)) for cov, disp in
                                           zip(binned_coverage, dispersion_norm)])
            binned_coverage_df['cov'] = this_chr_coverage
            binned_coverage_df.to_csv(self.binned_coverage_fn, sep='\t', index=False, header=False)


            corrected_binned_coverage_df = self.cnv_profile.generate_coverage(self.purity, self.binned_coverage_fn)
            corrected_binned_coverage_df['chrom'] = corrected_binned_coverage_df['chrom'].astype(str)
            corrected_binned_coverage_df['start'] = corrected_binned_coverage_df['start'].astype(int)
            corrected_binned_coverage_df['end'] = corrected_binned_coverage_df['end'].astype(int)

            corrected_binned_coverage_df = corrected_binned_coverage_df.merge(target_intervals_df, on=['chrom', 'start', 'end'])
            corrected_binned_coverage_df.to_csv(self.corrected_binned_coverage_fn, sep='\t', index=False)


    def sim_hets(self, normal_vcf, hets_df, output_dir, override=False):
        # for each segment estimate the expected allelic imbalance given the purity
        self.hets_fn = f'{output_dir}/{self.name}.hets.tsv'

        if not override and os.path.exists(self.hets_fn):
            print(f'Loading {self.hets_fn}')
            self.hets_df = pd.read_csv(self.hets_fn, sep='\t', index_col=0)
            return
        
        # use normal coverage profile to estimate count with total count
        corrected_binned_coverage_df = pd.read_csv(self.corrected_binned_coverage_fn, sep='\t')

        sample_hets_df = hets_df.copy()
        sample_hets_df = sample_hets_df.merge(
            corrected_binned_coverage_df[['covcorr', 'target_name']], 
            on='target_name'
        )
        sample_hets_df['DEPTH'] = sample_hets_df['covcorr'].copy()
        # self.cnv_profile.save_hets_file(filename, vcf, bed, purity, ref_alt=False)

        snv_df, correct_phase_interval_trees = self.cnv_profile.generate_snvs(
            vcf=normal_vcf, 
            bed=sample_hets_df,
            purity=self.purity, 
            ref_alt=False,
            do_parallel=False,
        )
        self.hets_df = snv_df

        self.hets_df.to_csv(self.hets_fn, sep='\t')

    def plot_hets(
        self, chromosome_col='CHROM', 
        start_position_col='POS', 
        alt_count_col='alt_count', 
        ref_count_col='ref_count',
        depth_col='adjusted_depth',
        genotype_col='NA12878',
        figsize=(20, 10),
        edgecolor=None,
        alpha=0.2,
        size=1,
        ploidy_col='ploidy',
    ):
        def get_cum_sum_csize(csize):
            cum_sum_csize = {}
            cum_sum = 0
            for chrom, size in csize.items():
                cum_sum_csize[chrom] = cum_sum
                cum_sum += size
            return cum_sum_csize

        csize = self.cnv_profile.csize
        maf_df = self.hets_df.copy()
        cum_sum_csize = get_cum_sum_csize(csize)
    
        chr_order = natsorted(list(csize.keys()))
        chrom_start = {chrom: start for (chrom, start) in
                       zip(np.append(chr_order, 'Z'), np.cumsum([0] + [csize[a] for a in chr_order]))}
    
        # maf_df = cached_read_csv(maf_fn, sep='\t', encoding='iso-8859-1')
        if maf_df[chromosome_col].dtype == 'object':
            maf_df[chromosome_col].replace({'X': 23, 'Y': 24}, inplace=True)
        maf_df[chromosome_col] = maf_df[chromosome_col].astype(str)
        
        maf_df['new_position'] = maf_df.apply(lambda r: cum_sum_csize[r[chromosome_col]] + r[start_position_col], axis=1)
        maf_df['tumor_f'] = maf_df[alt_count_col] / (maf_df[alt_count_col] + maf_df[ref_count_col])
    
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=figsize)

        p = sns.scatterplot(
            maf_df, x='new_position', y='tumor_f', ax=axes[1],
            hue=genotype_col,
            edgecolor=edgecolor,
            alpha=alpha,
            s=size,
        )

        maf_df[ploidy_col] = np.round(maf_df[ploidy_col], 2)
        p = sns.scatterplot(
            maf_df, x='new_position', y=depth_col, ax=axes[0],
            hue=ploidy_col,
            edgecolor=edgecolor,
            alpha=alpha,
            s=size,
        )

        axes[0].set_yscale('log')
    
        axes[1].set_xticks(
            ticks=np.asarray(list(chrom_start.values())[:-1]) + np.asarray(list(csize.values())) / 2,
            labels=chr_order,
            fontsize=10,
            # tickangle=0,
            # range=[0, chrom_start['Z']]
        )
        axes[1].set_xlabel("Chromosome")
        axes[1].set_ylim(0, 1.01)
        add_background(axes[0], csize.keys(), csize)
        add_background(axes[1], csize.keys(), csize)
        fig.suptitle(f'Sample {self.name} simulated hetsite allelic imbalance (purity={self.purity})')
        axes[0].legend(bbox_to_anchor=(1.0, 1.0))
        axes[1].legend(bbox_to_anchor=(1.0, 1.0))
        return fig, axes

    def sim_mut_vafs(self, variants_df, patient_name, data_directory, overwrite=False):
        # Given set of sites estimate the expected allele fraction

        self.variants_fn = f'{data_directory}/{self.name}.variants.tsv'

        if not overwrite and os.path.exists(self.variants_fn):
            print(f'Sample {self.name} has variants_fn: {self.variants_fn}')
            self.variants_df = pd.read_csv(self.variants_fn, sep='\t')
            return 
        
        sample_variants_df = variants_df.copy()

        def get_cov_and_ploidy(interval, interval_data_df):
            interval_row = interval_data_df[interval_data_df.target_name==interval].iloc[0]
            
            start = interval_row['start']
            end = interval_row['end']
            avg_covcorr = round(interval_row['covcorr']/(end-start))
            ploidy = interval_row['ploidy']
            
            return int(avg_covcorr), ploidy
        
        def get_overlapping_events(mut, cnv_profile):
            chrom = str(mut.chrom)
            # if chrom=='X':
            #     chrom='23'
            
            if mut.allele=='paternal':
                events = cnv_profile.event_trees[chrom].paternal_tree.at(mut.pos)
                alt_allele_events = cnv_profile.event_trees[chrom].maternal_tree.at(mut.pos)
            elif mut.allele=='maternal':
                events = cnv_profile.event_trees[chrom].maternal_tree.at(mut.pos)
                alt_allele_events = cnv_profile.event_trees[chrom].paternal_tree.at(mut.pos)
            return events, alt_allele_events
        
        def get_mut_and_unmut_copies_in_clone(upstream_mut, upstream_unmut, events, mut, clone):
            cn_change = 0
            is_possible = True
            
            unmut_copies = upstream_unmut
            mut_copies = upstream_mut
            
            for event in events:
                if (event.data.type!= 'haploid') & (event.data.cluster_num==clone):
                    cn_change += event.data.cn_change
                    
            if cn_change + upstream_mut + upstream_unmut < 0:
                raise Exception("Invalid CNAs!")
            
            # when examining a clone that the mutation doesn't occur in
            if mut.cluster!=clone: 
                
                # if copy gain, amplify mutated allele if one exists
                if cn_change > 0:
                    if  mut_copies>0:
                        mut_copies = upstream_mut + cn_change
                    else:
                        unmut_copies = upstream_unmut + cn_change 
                    
                # if copy loss, delete unmutated allele if one exists
                elif cn_change < 0:
                    if unmut_copies>0:
                        unmut_copies = upstream_unmut + cn_change
                    else:
                        mut_copies = upstream_mut + cn_change
                    
            # when examining the clone where the mutation occurs
            else:
                if upstream_mut !=0:
                    print("Warning: Possibly overlapping mutations: already mutated copies at mutation event!")
                    print(f'clone: {clone}')
                if (cn_change < 0) &  (upstream_unmut + cn_change <= 0) :
                    print("No copies of allele left to mutate!")
                    is_possible = False
                else:
                    
                    # if overlapping CN amplification event in same clone
                    if cn_change > 0:
                        mut_copies += 1+cn_change
                        unmut_copies += -1
                        
                    # if overlapping CN deletion event in same clone...should I even allow this?
                    elif cn_change < 0:
                        mut_copies += 1
                        unmut_copies += -1 + cn_change
                        
                    # if no overlapping CN event in same clone
                    elif cn_change == 0:
                        mut_copies += 1
                        unmut_copies -= 1
                
            if mut_copies < 0:
                print("Can't have a negative number of mutated alleles!")
                is_possible = False
            if unmut_copies < 0:
                print("Can't have a negative number of unmutated alleles!")
                is_possible = False
            return mut_copies, unmut_copies, is_possible
        
        def get_mut_multiplicity(mut, phylogeny):
            
            # mut is a series with: interval, chrom, pos, avg_covcorr, ploidy, ref_allele, cluster, allele, overlapping_CN_events
            print(f'gene: {mut.gene}, chrom: {mut.chrom}, pos: {mut.pos}, cluster: {mut.cluster}, allele: {mut.allele}')
            
            mut_multiplicity_dict = {}
            
            for node in phylogeny.ccfs.keys():
                if node == 1:
                    unmut_copies = 1
                else:
                    unmut_copies = 0
        
                mut_multiplicity_dict[node]={'mut_copies':0, 'unmut_copies':unmut_copies}
            
            # iterate through nodes in the tree from top down 
            # such that every parent node is complete before any of it's children
            
            children = [1]
            nodes_calculated = []
            
            while children:
                new_children = []
                for child in children:
                    
                    if child!=1:
                        parent = phylogeny.parents[child]
                        upstream_mut = mut_multiplicity_dict[parent]['mut_copies']
                        upstream_unmut = mut_multiplicity_dict[parent]['unmut_copies']
                    else:
                        upstream_mut = 0
                        upstream_unmut = 1 # assume normal cells have one copy of each allele
                    
                    mut_num, unmut_num, is_possible = get_mut_and_unmut_copies_in_clone(
                        upstream_mut,
                        upstream_unmut, 
                        mut.overlapping_CN_events, 
                        mut,
                        child
                    )
                    if not is_possible:
                        return None
                    else:
                        mut_multiplicity_dict[child]['mut_copies'] = mut_num
                        mut_multiplicity_dict[child]['unmut_copies'] = unmut_num
                    
                    new_children = new_children + phylogeny.get_children(child)
                    nodes_calculated.append(child)
                children = new_children
            
            return mut_multiplicity_dict
        
        def get_true_vaf(mut, phylogeny, purity):
            if mut.multiplicity is not None:
                total_mut_copies = 0
                for node in mut.multiplicity:
                    node_fraction = phylogeny.ccfs[node]
                    for child in phylogeny.get_children(node):
                        node_fraction -= phylogeny.ccfs[child]
                    total_mut_copies += node_fraction*mut.multiplicity[node]['mut_copies']*purity
                vaf = total_mut_copies/mut.ploidy
                if vaf > 1:
                    raise Exception("More mutated copies than local ploidy!")
                return vaf
            else:
                return None
            
        def get_alt_count(mut):
            if not math.isnan(mut.vaf):
                t_alt_count = round(scipy.stats.binom.rvs(int(mut.avg_covcorr), mut.vaf, size=1)[0])
                return int(t_alt_count)
            
        def get_local_allelic_cn(mut, ccfs, purity):
            a1 = 0
            a2 = 0
            
            for i in mut.overlapping_CN_events:
                a1 += i.data.cn_change*ccfs[i.data.cluster_num]
            for i in mut.alt_overlapping_CN_events:
                a2 += i.data.cn_change*ccfs[i.data.cluster_num]
            a1 = a1*purity+2*(1-purity)
            a2 = a2*purity+2*(1-purity)
        
            return a1, a2

        corrected_binned_coverage_df = pd.read_csv(self.corrected_binned_coverage_fn, sep='\t')
        sample_variants_df[['avg_covcorr', 'ploidy']] = sample_variants_df.apply(
            lambda x: get_cov_and_ploidy(x.interval, corrected_binned_coverage_df), axis=1, result_type='expand'
        )
        sample_variants_df[['overlapping_CN_events','alt_overlapping_CN_events']] = sample_variants_df.apply(
            lambda x: get_overlapping_events(x, self.cnv_profile), axis=1, result_type='expand'
        )
        sample_variants_df['multiplicity'] = sample_variants_df.apply(
            lambda x: get_mut_multiplicity(x, self.cnv_profile.phylogeny), axis=1
        )
        sample_variants_df['vaf'] = sample_variants_df.apply(
            lambda x: get_true_vaf(x, self.cnv_profile.phylogeny, self.purity), axis=1
        )
        sample_variants_df['t_alt_count'] = sample_variants_df.apply(
            lambda x: get_alt_count(x), axis=1
        )
        sample_variants_df['t_ref_count'] = sample_variants_df['avg_covcorr'] - sample_variants_df['t_alt_count']
        sample_variants_df[['local_cn_a1','local_cn_a2']] = sample_variants_df.apply(
            lambda x: get_local_allelic_cn(x, self.clone_ccfs, self.purity), axis=1, result_type='expand'
        )

        sample_variants_df = sample_variants_df.rename(
            columns={
                'chrom': 'Chromosome', 
                'pos':'Start_position', 
                'ref_allele':'Reference_Allele', 
                'alt_allele':'Tumor_Seq_Allele2',
                'gene': 'Hugo_Symbol'
            }
        )
        sample_variants_df = sample_variants_df[[
            'Hugo_Symbol', 
            'Chromosome',
            'Start_position',
            'Reference_Allele',
            'Tumor_Seq_Allele2',
            't_alt_count',
            't_ref_count',
            'local_cn_a1',
            'local_cn_a2',
            'cluster',
            'allele',
            'avg_covcorr',
            'ploidy',
            'overlapping_CN_events',
            'alt_overlapping_CN_events',
            'multiplicity',
            'vaf'
        ]]
        sample_variants_df['sample_id'] = self.name
        sample_variants_df['participant_id'] = patient_name

        self.variants_df = sample_variants_df[sample_variants_df.vaf.notna()]
        self.variants_df.to_csv(self.variants_fn, sep='\t', index=False)

        '''
        simulated maf file output columns:
        Hugo_Symbol: gene name or "Unknown"
        Chromosome: chromosome number, 23 for X chromosome
        Start_position: mutation position
        Reference_Allele: reference allele at chromosome/start position
        Tumor_Seq_Allele2: simulated alternate allele at chromosome/start position
        t_alt_count: simulated number of alternate reads
        t_ref_count: simulated number of reference reads
        local_cn_a1: local copy number of allele 1 (needed as input to PhylogicNDT if seg files aren't provided)
        local_cn_a2: local copy number of allele 2 (needed as input to PhylogicNDT if seg files aren't provided)
        cluster: simulated cluster assignment (not typically known)
        allele: simulated maternal/paternal allele assignment (not typically known)
        avg_covcorr: simulated coverage
        ploidy: local ploidy (not typically known)
        overlapping_CN_events: interval tree giving CN events that overlap with mutation
        alt_overlapping_CN_events: interval tree giving CN events on opposite allele (maternal/parternal)
        multiplicity: dictionary giving interval multiplicity of mutation in each subclone (not typically known)
        vaf: simulated variant allele fraction
        '''
        
    
class Patient:

    def __init__(
        self, 
        name, 
        data_directory,
        samples: list[Sample],
        num_subclones=7, 
        clone_tree={1: None, 2:1, 3:1, 4:2, 5:4, 6:4, 7:3},
        age=32,
        tumor_molecular_subtype='Unknown',
        tumor_morphology='Unknown',
        tumor_primary_site=np.nan,
        cancer_stage=np.nan,
        vital_status=np.nan,
        death_date_dfd=np.nan,
        follow_up_date=np.nan,
        age_at_diagnosis=np.nan,
        gender=np.nan,
    ):
        self.name = name
        self.samples = samples

        self.num_subclones = num_subclones
        self.parents = clone_tree
        self.parent_2_children_tree = self.convert_tree(clone_tree)
        ccfs_df = pd.concat([pd.Series(sample.calc_clone_ccfs(self.parent_2_children_tree), name=sample.name) for sample in self.samples], axis=1).T
        self.sample_ccfs_df = ccfs_df

        if not os.path.exists(data_directory):
            print(f'Creating new directory {data_directory}')
            os.mkdir(data_directory)
            
        self.data_directory = f'{data_directory}/{self.name}'
        if not os.path.exists(self.data_directory):
            print(f'Making directory for {self.name} at {self.data_directory}')
            os.mkdir(self.data_directory)

        self.age = age
        self.tumor_molecular_subtype = tumor_molecular_subtype
        self.tumor_morphology = tumor_morphology
        self.tumor_primary_site = tumor_primary_site
        self.cancer_stage = cancer_stage
        self.vital_status = vital_status
        self.death_date_dfd = death_date_dfd
        self.follow_up_date = follow_up_date
        self.age_at_diagnosis = age_at_diagnosis
        self.gender = gender

    def set_cnv_profile(
        self,
        arm_num=20,
        focal_num=3,
        p_whole=0.6,
        ratio_clonal=0.5,
        override=False,
    ):
        self.arm_num = arm_num
        self.focal_num = focal_num
        self.p_whole = p_whole
        self.ratio_clonal = ratio_clonal
        
        self.cnv_profile = simulate.CNV_Profile(num_subclones=self.num_subclones - 1)
        # self.cnv_profile_pkl = f'{self.data_directory}/{self.name}.cnv_events_focal{focal_num}_arm_{arm_num}.pkl'
        self.cnv_profile_pkl = f'{self.data_directory}/{self.name}.cnv_events.pkl'

        if os.path.exists(self.cnv_profile_pkl) and not override:
            print(f'loading existing CNV pickle file {self.cnv_profile_pkl}')
        else:
            print('Regenerating CNV events.')
            self.cnv_profile.add_cnv_events(
                arm_num=self.arm_num, 
                focal_num=self.focal_num, 
                p_whole=self.p_whole, 
                ratio_clonal=self.ratio_clonal
            )
            self.cnv_profile.pickle_events(self.cnv_profile_pkl)

        file = open(self.cnv_profile_pkl, 'rb')
        event_trees = pickle.load(file)
        file.close()
        self.event_trees = event_trees

    def set_sample_cnv_profiles(self):
        output_dir = f'{self.data_directory}/sample_cnv_profiles'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            
        for sample in self.samples:
            sample.set_cnv_profile(
                self.num_subclones, 
                self.parents, 
                self.event_trees, 
                data_directory=output_dir
            )


    def set_sample_coverage(self, target_intervals_df, sigma=2, x_coverage=120, override=False):

        output_dir = f'{self.data_directory}/sample_coverage'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        for sample in self.samples:
            sample.sim_coverage(target_intervals_df, output_dir, sigma, x_coverage, override)

    def set_sample_muts(self):
        output_dir = f'{self.data_directory}/sample_muts'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        for sample in self.samples:
            sample.sim_mut_vafs(self.variants_df, self.name, output_dir)

    def set_sample_hets(self):
        output_dir = f'{self.data_directory}/sample_hets'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        for sample in self.samples:
            # normal_vcf: cnv_suite expects a file path
            sample.sim_hets(self.normal_vcf, self.hets_df, output_dir)

    def set_treatments(
        self, 
        treatment_fn=None,
        num_treatments=2,
        categories=['Chemotherapy', 'Precision/Targeted therapy'],
        drugs=['Trastuzumab', 'Neratinib'],
        drug_combination=[np.nan, np.nan],
        start_date_dfd=[100, 200],
        stop_date_dfd=[200, 300],
        stop_reason=['Unknown', 'Unknown'],
        pre_status=[np.nan, np.nan],
        post_status=['Unknown', 'Unknown'],
        notes=[np.nan, np.nan]
    ):
        if treatment_fn:
            self.treatment_fn = treatment_fn
        elif treatment_fn is None:
            treatments = pd.DataFrame({
                'participant_id': [self.name] * num_treatments,
                'categories': categories, 
                'drugs': drugs, 
                'drug_combination': drug_combination, 
                'start_date_dfd': start_date_dfd, 
                'stop_date_dfd': stop_date_dfd, 
                'stop_reason': stop_reason, 
                'pre_status': pre_status, 
                'post_status': post_status, 
                'notes': notes, 
            })
            treatment_fn = f'{self.data_directory}/{self.name}_treatments.txt'
            treatments.to_csv(treatment_fn, sep='\t', index=False)
            self.treatment_fn = treatment_fn

    def parent_tree2str(self):
        return ','.join([f'{child}-{parent}' for parent, child in self.parents.items()])

    def convert_tree(self, clone_tree):
        parent_child_tree = {}
        for child, parent in clone_tree.items():
            if parent not in parent_child_tree.keys():
                parent_child_tree[parent] = [child]
            else:
                parent_child_tree[parent].append(child)

        return parent_child_tree

    def set_hets_df(self, normal_vcf, target_intervals_df, overwrite=False):
        """
        normal_vcf: normal vcf file
        target_intervals_df: pd.DataFrame of target intervals
        """
        self.hets_fn = f'{self.data_directory}/hets.tsv'
        self.normal_vcf = normal_vcf
        if overwrite or not os.path.exists(self.hets_fn):
            het_df = pd.read_csv(
                normal_vcf, 
                sep='\t', comment='#', header=None, 
                names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NA12878']
            )
            het_df['CHROM'] = het_df['CHROM'].astype(str).str.strip('chr').replace({'X': '23', 'Y': '24'})
    
            for i, r in tqdm.tqdm(het_df.iterrows()):
                het_interval_df = target_intervals_df[
                    (target_intervals_df['chrom'] == r['CHROM']) & 
                    (target_intervals_df['start'] <= r['POS']) &
                    (target_intervals_df['end'] >= r['POS'])
                ]
                if not het_interval_df.empty:
                    het_df.loc[i, 'target_name'] = het_interval_df.iloc[0]['target_name']
    
            self.hets_df = het_df.dropna(subset='target_name')
            self.hets_df.to_csv(self.hets_fn, sep='\t')
        else:
            print(f'Loading existing hets_df from {self.hets_fn}')
            self.hets_df = pd.read_csv(self.hets_fn, sep='\t', index_col=0)

    def set_variants_df(
        self, target_intervals_df, fasta_file_path, 
        prepped_gencode_gene_df=None, gene_name_col='gene', num_variants=100,
        random_seed_val=0,
        overwrite=False,
    ):
        self.variants_fn = f'{self.data_directory}/{self.name}.variants.tsv'

        if not overwrite and os.path.exists(self.variants_fn):
            print(f'patient variants path exists: {self.variants_fn}')
            self.variants_df = pd.read_csv(self.variants_fn, sep='\t', index_col=0)
            return
        

        def get_alt_allele(ref_allele):
            possible_alleles = ['C', 'G', 'A', 'T']
            possible_alleles.remove(ref_allele)
            allele_index = random.choices(population=range(0,3), k=1)[0]
            
            return possible_alleles[allele_index]

        def get_chrom_pos(interval, interval_data_df, fasta_file_path, sex_chrom_dict = {'23': 'X', '24': 'Y'}):
            """
            sex_chrom_dict: dictionary to translate chromosome name if the original is in a different format than the fasta ref file
            """
            interval_row = interval_data_df[interval_data_df.target_name==interval].iloc[0]
            
            chrom = interval_row['chrom']

            query_chrom = sex_chrom_dict[chrom] if chrom in sex_chrom_dict.keys() else chrom
            
            start = interval_row['start']
            end = interval_row['end']
            pos = random.choices(population=range(start,end), k=1)[0]
        
            ref_base_str = pysam.faidx(fasta_file_path, str(query_chrom)+':'+str(pos)+'-'+str(pos))
            ref_base = ref_base_str[-2]
            alt_base = get_alt_allele(ref_base)

            return chrom, pos, ref_base, alt_base
        
        np.random.seed(random_seed_val)
        variants_df = pd.DataFrame(index=range(num_variants))
        
        variants_df['interval'] = random.choices(
            population=target_intervals_df.target_name, 
            weights=target_intervals_df.interval_weight, 
            k=num_variants
        )
        variants_df[
            ['chrom','pos', 'ref_allele', 'alt_allele']
        ] = variants_df.apply(
            lambda x: get_chrom_pos(x.interval, target_intervals_df, fasta_file_path), 
            axis=1, 
            result_type='expand'
        )
        variants_df['cluster'] = variants_df.apply(
            lambda x: random.choices(range(self.num_subclones))[0]+1, 
            axis=1
        )
        variants_df['allele'] = variants_df.apply(
            lambda x: 'maternal' if random.choices(range(2))[0]==0 else 'paternal', 
            axis=1
        )

        if prepped_gencode_gene_df is not None:
            for i, r in variants_df.iterrows():
    
                interval_df = prepped_gencode_gene_df[
                    (prepped_gencode_gene_df['Chromosome'] == r['chrom']) &
                    (prepped_gencode_gene_df['Start_position'] <= r['pos']) & 
                    (prepped_gencode_gene_df['End_position'] >= r['pos']) 
                ]
    
                if not interval_df.empty:
                    variants_df.loc[i, 'gene'] = interval_df.iloc[0][gene_name_col]
        
        self.variants_df = variants_df
        self.variants_df.to_csv(self.variants_fn, sep='\t')

    def force_add_variant(
        self, target_intervals_df, gene, chrom, pos, ref_allele, alt_allele, cluster, allele='paternal'
    ):

        new_variant = {
            'gene': gene,
            'chrom': int(chrom), 
            'pos': int(pos), 
            'ref_allele': ref_allele, 
            'alt_allele': alt_allele, 
            'cluster': cluster, 
            'allele': allele, 
        }

        new_variant['interval'] = target_intervals_df.loc[
            (target_intervals_df['chrom'].astype(str) == str(chrom)) & 
            (target_intervals_df['start'].astype(int) <= int(pos)) & 
            (target_intervals_df['end'].astype(int) >= int(pos))
        ]['target_name'].iloc[0]
        
        
        self.variants_df = self.variants_df.append(
            new_variant, ignore_index=True
        ).reset_index(drop=True).drop_duplicates()
        self.variants_df.to_csv(self.variants_fn, sep='\t')

    def plot_ccf(self):
        reformat_ccfs_df = self.sample_ccfs_df.stack().reset_index().rename(
            columns={'level_0': 'Sample', 'level_1': 'Cluster', 0: 'Cancer Cell Fraction (CCF)'}
        )
        sns.lineplot(
            data=reformat_ccfs_df, 
            x='Sample', 
            y='Cancer Cell Fraction (CCF)', 
            hue='Cluster', 
            palette=[cluster_color(i) for i in reformat_ccfs_df.Cluster.unique()], 
        )
        plt.legend(bbox_to_anchor=(1, 1), title='Clusters')
        plt.title(f'{self.name} CCF plot')

    def plot_tree(self, image_path='base64uri'):        
        edges = [f'{child}-{parent}' for parent, child in self.parents.items()]
        
        cluster_list = []
        for i in edges:
            new_list = i.split('-')
            for j in new_list:
                if (j !='None') & (j not in cluster_list):
                    cluster_list.append(j)
        cluster_list = sorted(cluster_list)
        
        nodes = [{'data': {'id': 'normal', 'label': 'normal'}, 'position': {'x': 0, 'y': 0}}]
        
        nodes.extend([
            {
                'data': {'id': f'cluster_{cluster}', 'label': cluster},
                'position': {'x': 50 * int(cluster), 'y': -50 * int(cluster)}
            }
            for cluster in cluster_list
        ])
        edges_list = []
        nodes_copy = nodes.copy()
        for edge in edges:
            nodes_copy = edge.split('-')
            if nodes_copy[0]!='None':
                nodes_copy = list(map(int,edge.split('-')))
                edges_list.append(nodes_copy)
        
        
        edges = [{'data': {'source': 'normal', 'target': 'cluster_1', 'label': ''}}]
        edges.extend([
            {'data': {'source': f'cluster_{edge[0]}', 'target': f'cluster_{edge[1]}', 'label': ''}}
            for edge in edges_list
        ])
        
        elements = nodes + edges
        
        stylesheet = gen_stylesheet(cluster_list)
        
        tree_plot = cyto.Cytoscape(
            id='phylogic-tree',
            style={'width': '100%', 'height': '450px'},
            layout={
                'name': 'breadthfirst',
                'roots': '[id="normal"]'
            },
            elements=elements,
            stylesheet=stylesheet,
            userZoomingEnabled=False
        )
        
        @callback(
            Output("phylogic-tree", "generateImage"),
            [Input('save-image-button', 'n_clicks')], prevent_initial_call=True
        )
        def get_image(nclicks):
            return {
                'type': 'jpg',
                'action': 'download',
                'options': {'quality': 1, 'output': image_path}
            }
        
        app = Dash(__name__)
        
        app.layout = html.Div([
            tree_plot,
            html.Button('save image', id='save-image-button', n_clicks=0)
        ])
        app.run(debug=False)

    def plot_sample_hets(self, **kwargs):
        for sample in self.samples:
            fig, ax = sample.plot_hets(**kwargs)
            plt.show()


    def generate_phylogicNDT_sif(self):
        sampleIDs = [s.name for s in self.samples]
        sample_mafs = [s.variants_fn for s in self.samples]
        sample_cn_profiles = [''] * len(self.samples)
        purities = [str(s.purity) for s in self.samples]
        times = [str(s.time_point) for s in self.samples]

        sif_fn = f'{self.data_directory}/{self.name}.phylogicNDT.sif'
        with open(sif_fn, 'w') as f:
            f.write('sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint')
            for line in sorted(zip(sampleIDs, sample_mafs, sample_cn_profiles, purities, times), key=lambda k: (float(k[4]), k[0])):
                f.write('\n'+'\t'.join(line))

        self.sif_fn = sif_fn
        print(f'Generated sif file: {self.sif_fn}')

    def run_phylogicNDT(
        self, 
        python2_path='python2',
        phylogicNDT_py_path='/phylogicndt/PhylogicNDT.py',
        overwrite=False
    ):
        # See https://github.com/broadinstitute/PhylogicNDT
        self.phylogicNDT_results_dir = f'{self.data_directory}/phylogicNDT_results'

        cd_cmd = f'cd {self.phylogicNDT_results_dir}'
        phylogicNDT_cmd = f'{python2_path} {phylogicNDT_py_path} Cluster -i "{self.name}" -sif "{self.sif_fn}" --maf_input_type "calc_ccf" -rb'
        
        full_cmd = cd_cmd + ' && ' + phylogicNDT_cmd
        self.phylogicNDT_full_cmd = full_cmd

        if not overwrite and os.path.exists(self.phylogicNDT_results_dir) and len(os.listdir(self.phylogicNDT_results_dir)) > 0:
            print(f'Phylogic results are already loaded here: {self.phylogicNDT_results_dir}')
            
        else:
            if not os.path.exists(self.phylogicNDT_results_dir):
                os.mkdir(self.phylogicNDT_results_dir)

            os.system(full_cmd)
        
    def set_up_patient_reviewer_data(self):

        self.patient_reviewer_data_dir = f'{self.data_directory}/patient_reviewer_data'
        if not os.path.exists(self.patient_reviewer_data_dir):
            os.mkdir(self.patient_reviewer_data_dir)
        
        sample_data_df = pd.read_csv(self.sif_fn, sep='\t')
        sample_data_df['cnv_seg_fn'] = [s.absolute_cnv_profile_fn for s in self.samples]
        sample_data_df.drop(['seg_fn'], axis=1, inplace=True)
        sample_data_df['participant_id'] = self.name
        sample_data_df['preservation_method'] = np.nan
        sample_data_df['wxs_ploidy'] = 2
        sample_data_df.rename(columns={'purity': 'wxs_purity', 'timepoint': 'collection_date_dfd'}, inplace=True)

        self.patient_reviewer_sample_data_fn = f'{self.patient_reviewer_data_dir}/patient1_samples_data.tsv'
        sample_data_df.to_csv(self.patient_reviewer_sample_data_fn, sep='\t', index=False)
        print(f'Generated patient_reviewer_sample_data_fn: {self.patient_reviewer_sample_data_fn}')

        patient_data_df = pd.DataFrame({
            'participant_id': [self.name], 
            'maf_fn': [f'{self.phylogicNDT_results_dir}/{self.name}.mut_ccfs.txt'], 
            'cluster_ccfs_fn': [f'{self.phylogicNDT_results_dir}/{self.name}.cluster_ccfs.txt'], 
            'build_tree_posterior_fn': [f'{self.phylogicNDT_results_dir}/{self.name}_build_tree_posteriors.tsv'],
            'tumor_molecular_subtype': self.tumor_molecular_subtype,
            'tumor_morphology': self.tumor_morphology,
            'tumor_primary_site': self.tumor_primary_site,
            'cancer_stage': self.cancer_stage,
            'vital_status': self.vital_status,
            'death_date_dfd': self.death_date_dfd,
            'follow_up_date': self.follow_up_date,
            'age_at_diagnosis': self.age,
            'gender': self.age_at_diagnosis,
            'notes': self.gender,
            'treatments_fn': self.treatment_fn
        })

        self.patient_reviewer_patient_data_fn = f'{self.patient_reviewer_data_dir}/patient1_data.tsv'
        patient_data_df.to_csv(self.patient_reviewer_patient_data_fn, sep='\t', index=False)
        print(f'Generated patient_reviewer_patient_data_fn: {self.patient_reviewer_patient_data_fn}')


    def set_all_data(
        self, 
        arm_num=20,
        focal_num=3,
        p_whole=0.6,
        ratio_clonal=0.5,
        target_intervals_df: pd.DataFrame = None,
        normal_vcf_path: str = None,
        local_fasta_file_path: str = None,
        genes_df: pd.DataFrame = None,
        num_variants=100,
        variants_random_seed_val=1,
        force_add_variants=[], # [{'gene':, 'chrom':, 'pos':, 'ref_allele':, 'alt_allele':, 'cluster':, 'allele'}]
        python2_path='/Users/cchu/opt/anaconda3/envs/phylogicNDT_py27_env/bin/python',
        phylogicNDT_py_path='./src/PhylogicNDT/PhylogicNDT.py',
        overwrite_phylogicNDT=False,
        run_hets=False,
    ):
        self.set_treatments()
        self.set_cnv_profile(
            arm_num=arm_num,
            focal_num=focal_num,
            p_whole=p_whole,
            ratio_clonal=ratio_clonal,
        )
        self.set_sample_cnv_profiles()

        self.set_sample_coverage(
            target_intervals_df=target_intervals_df, 
            override=False
        )

        if run_hets:
            self.set_hets_df(
                normal_vcf=normal_vcf_path, 
                target_intervals_df=target_intervals_df
            )
            self.set_sample_hets()

        self.set_variants_df(
            target_intervals_df, local_fasta_file_path, genes_df, 
            num_variants=num_variants, random_seed_val=variants_random_seed_val
        )
        self.set_sample_muts()

        for variant_dict in force_add_variants:
            self.force_add_variant(
                target_intervals_df, 
                **variant_dict,
                # gene='TP53', 
                # chrom='17', 
                # pos='7578190', # hg19 coordinate 
                # ref_allele='T', 
                # alt_allele='C', 
                # cluster=3, 
                # allele='paternal'
            )

        self.generate_phylogicNDT_sif()
        self.run_phylogicNDT(
            python2_path=python2_path,
            phylogicNDT_py_path=phylogicNDT_py_path,
            overwrite=overwrite_phylogicNDT
        )

        self.set_up_patient_reviewer_data()

def prep_gencode_gene_df(gene_tsv_fn='gencode.v19.annotation.gene_only.tsv'):
    genes = pd.read_csv(gene_tsv_fn, sep='\t', header=None)
    genes.rename(columns={0: 'Chromosome', 3: 'Start_position', 4: 'End_position', 6: 'strand', 8: 'data'}, inplace=True)
    genes['Chromosome'] = genes['Chromosome'].apply(lambda x: x[3:]).replace({'X': 23, 'Y': 24, 'M': 25}).astype(str)
    genes['Start_position'] = genes['Start_position'].astype('int')
    genes['End_position'] = genes['End_position'].astype('int')
    genes['gene'] = genes.data.apply(lambda x: re.findall('gene_name "(.*?)"', x)[0].strip())
    return genes

