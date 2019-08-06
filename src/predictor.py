#!/usr/bin/env python3

import numpy as np
from proximity import HiCFetcher, DistanceModel
import json
import pandas as pd
from tools import get_gene_name
import sys


class Predictor(object):
    def __init__(self, loaded_enhancers, **args):

        self.cellType = args['cellType']

        #HiC
        self.hic_fetcher = HiCFetcher(args['HiCdir'], 
                                        args['hic_gamma'],  
                                        args['hic_gamma_reference'],
                                        scale_with_powerlaw=args['scale_hic_using_powerlaw'],
                                        tss_hic_contribution=args['tss_hic_contribution'])

        #Power Law
        if args['scale_hic_using_powerlaw']:
            model_gamma = args['hic_gamma_reference']
        else:
            model_gamma = args['hic_gamma']
        self.distance_model = DistanceModel(model_gamma)

        #HiC adjustment parameters
        self.tss_hic_contribution = args['tss_hic_contribution']
        self.hic_pseudocount_distance = args['hic_pseudocount_distance']
        self.hic_cap = 100 #args['hic_cap']

    def estimate_contact_probability_from_distance(self, distance):
        return self.distance_model(distance)

    def predict_from_normalized_to_enhancers(self, enhancers, gene, window, tss_slop=500):
        midpoint = (enhancers.start.values + enhancers.end.values) / 2
        enhancers['distance'] = abs(gene['tss'] - midpoint)
        
        #Add gene specific annotations
        enhancers['isSelfPromoter'] = np.logical_and.reduce((enhancers.isPromoterElement == True, enhancers.start - tss_slop < gene.tss, enhancers.end + tss_slop > gene.tss))
        # enhancers['isSelfGenic'] = np.logical_or(np.logical_and(enhancers.start > gene.start,enhancers.start < gene.end),
        #                                                 np.logical_and(enhancers.end > gene.start, enhancers.end < gene.end))
        enhancers['TargetGene'] = gene['name']
        enhancers['TargetGeneTSS'] = gene['tss']
        if 'is_ue' in gene:
            enhancers['TargetGeneIsUbiquitouslyExpressed'] = gene['is_ue']

        if 'Expression' in gene.index:
            enhancers['TargetGeneExpression'] = gene['Expression'] 
        else: 
            enhancers['TargetGeneExpression'] = np.nan 

        if 'PromoterActivityQuantile' in gene.index:
            enhancers['TargetGenePromoterActivityQuantile'] = gene['PromoterActivityQuantile'] 
        else: 
            enhancers['TargetGenePromoterActivityQuantile'] = np.nan

        is_self_tss = enhancers['isSelfPromoter'].values
        if sum(is_self_tss) == 0:
            print("No candidate self-promoter of {} {} {}. May want to investigate!".format(gene['name'], gene['chr'], gene['tss']))
        elif sum(is_self_tss) > 1:
            print("Found multiple self-promoters of {} {} {}. May want to investigate - but okay if candidate regions are not merged!".format(gene['name'], gene['chr'], gene['tss']))


        #Get Hi-C data
        hic_vals, rowmax, self.hic_exists, hic_vals_unscaled, rowmax_unscaled  = self.hic_fetcher(gene.chr, gene.tss, midpoint, enhancers)
        enhancers['hic.distance'] = hic_vals
        enhancers['hic.rowmax'] = rowmax
        enhancers['hic.distance.unscaled'] = hic_vals_unscaled
        enhancers['hic.rowmax.unscaled'] = rowmax_unscaled

        #add hic pseudocount
        enhancers['hic.distance.adj'], enhancers['hic_adjustment'] = self.add_hic_pseudocount(enhancers['distance'], 
                                                                                                    enhancers['hic.distance'], 
                                                                                                    enhancers['hic.rowmax'], 
                                                                                                    hic_pseudocount_distance=self.hic_pseudocount_distance)

        #Power Law
        enhancers['estimatedCP'], cp_rowmax = self.estimate_contact_probability_from_distance(enhancers['distance'])
        enhancers['estimatedCP.adj'] = self.normalize_proximity_contact_probability(enhancers['estimatedCP'], cp_rowmax)

        #Compute ABC Score and related scores
        enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['hic.distance.adj']], "ABC")
        enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['estimatedCP.adj']], "powerlaw")
        #enhancers = compute_score(enhancers, [enhancers['activity_base_noqnorm'], enhancers['hic.distance.adj']], "ABC.noqnorm")

    def __call__(self, *args, **kwargs):
        return self.predict(*args, **kwargs)

    def add_hic_pseudocount(self, dists, hic_vals, hic_rowmax, hic_pseudocount_distance=1e6):
        powerlaw_cp, cp_rowmax = self.estimate_contact_probability_from_distance(dists)
        cp_at_distance = self.estimate_contact_probability_from_distance(hic_pseudocount_distance)[0]
        adjustment = np.minimum(100 * (cp_at_distance / cp_rowmax), 100 * (powerlaw_cp / cp_rowmax))
        return np.clip(100 * (hic_vals / hic_rowmax) + adjustment, 0, self.hic_cap), adjustment

    def normalize_proximity_contact_probability(self, contact_probability, cp_rowmax):
        return np.clip(100 * (contact_probability / cp_rowmax), 0, self.hic_cap)

    def chromosomes(self):
        "returns a list of which chromosomes have data for predictions."
        return list(set(self.hic_fetcher.chromosomes()))

    def get_gene_prediction_stats(self, args, nearby_enhancers):
        stats = pd.Series( {
            'TargetGene' : nearby_enhancers.TargetGene.values[0],
            'TargetChr' : nearby_enhancers.chr.values[0],
            'TargetGeneTSS' : nearby_enhancers.TargetGeneTSS.values[0],
            'nEnhancersConsidered' : int(nearby_enhancers.shape[0]),
            'nDistalEnhancersPredicted' : int(sum(nearby_enhancers.loc[~nearby_enhancers['isPromoterElement'], args.score_column] > args.threshold))
            })
        return stats

def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis = 1)
    total = sum(scores)
    normalized_scores = scores / (total if (total > 0) else 1)

    enhancers[prefix + '.Score.Numerator'] = scores
    enhancers[prefix + '.Score'] = normalized_scores

    return(enhancers)
