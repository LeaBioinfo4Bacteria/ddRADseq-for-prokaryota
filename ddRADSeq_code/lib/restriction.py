#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib restriction
"""

import Bio.Restriction as Res
#import pandas as pd

def make_superRb(rb,superRb):
    '''ajoute les enz du rb en entrée dans le superRb si elles n'y sont pas déjà,
    return le superRb'''
    for enz in rb:
        if enz not in superRb:
            superRb.add(enz)
    return(superRb)


def rb_sans_supp(sites_variants,enzaway,bout_francs):
    rb = Res.RestrictionBatch(['BclI','EagI','SalI','NarI','HpyCH4IV'])
    enzHF = ['BclI-HF','EagI-HF','EagI','SalI-HF','SalI','NarI','HpyCH4IV']
    rb,sites_variants = site_variant(rb,sites_variants)       #enlever les sites variants pour n'importe quelle enzyme
    for enz in sites_variants:
        enzaway.writerow([enz]+['site variant'])
    rb,bout_francs = blunt(rb,bout_francs)                    #enlever les enzymes à bouts francs
    for enz in bout_francs:
        enzaway.writerow([enz] + ['bouts francs'])
    #rb = Res.RestrictionBatch(['HpyCH4IV','EagI','NarI'])
    #enzHF = ['HpyCH4IV','EagI-HF','NarI']
    return(rb,sites_variants,bout_francs,enzHF)


def blunt(rb,bout_francs):
    '''enlève les enzymes à bouts francs du rb,
    return le rb modifié'''
    rbcop = rb.copy()
    for enz in rbcop:
        if enz.is_blunt() :
            rb.remove(enz)
            if enz not in bout_francs:
                bout_francs.append(enz)
    bout_francs.sort()
    return(rb,bout_francs)



def site_variant(rb,sites_variants):
    '''enlève les enzymes du rb si elles ont un site de coupure variant,
    retourne le rb modifié'''
    rbcop = rb.copy()
    nucleotid_variant = ['R','Y','S','W','K','M','B','D','H','V','N']
    for enz in rbcop:
        site = []
        for n in enz.site:
            site.extend(n)
        for n in nucleotid_variant:
            if n in site:
                rb.remove(enz)
                if enz not in sites_variants:
                    sites_variants.append(enz)
    sites_variants.sort()
    return(rb,sites_variants)



def create_keylist_mono(dicorb,pas_de_site):
    '''après le search, enlève les enzymes ne présentant pas de site de coupures.
    renvoie les résultats du search du rb dans le dicorb,
    et la liste des enzymes gardées'''
    pas_de_site=[]
    keylist = []
    for key,value in dicorb.items() :
        if len(value) == 0 :
            pas_de_site.append(key)
        else:
            keylist.append(key)
    return(keylist,pas_de_site)



def create_keylist_multi(rb):
    '''ajoute chaque clef du dico dans une liste si elle n'est pas déjà dedans,
    renvoie la liste keylist'''
    keylist = []
    for enz in rb :
        if enz not in keylist :
            keylist.append(enz)    
    return(keylist)


def pas_de_site(dico):
    '''si l'enzyme n'a pas de site de coupure dans la séquence, l'ajoute dans la liste pas_de_site
    renvoie la liste pas_de_site'''
    pas_de_site = []
    for key,value in dico.items() :
        if len(value) == 0 :
            pas_de_site.append(key)
    return(pas_de_site)
  

