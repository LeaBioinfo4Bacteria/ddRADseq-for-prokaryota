#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib pour main
"""

from lib import restriction
from lib import modif_NEB
from lib import analyse_couples_enz
from lib import analyse_frag
from lib import writing
from lib.crea_object import Seq
from Bio import SeqIO
import itertools
import pandas as pd
import Bio.Restriction as Res



def pair_superRb_NEB(tab_supp,superRb,superDico,pairToKeep):
    df = pd.read_csv(tab_supp,dtype=str) 
    filtre_all = modif_NEB.modif_tab(df)
    enzHF = modif_NEB.recup_enz(filtre_all)
    rb = modif_NEB.modif_rb(enzHF)
    superRb = restriction.make_superRb(rb,superRb)
    keylist = restriction.create_keylist_multi(rb)
    keylist = modif_NEB.modif_keylist_multi(keylist,enzHF)
    for enz1,enz2 in itertools.combinations(keylist,2):
        incub1,incub2,inac1,inac2,buf1,buf2,sit1,sit2 = modif_NEB.filtre_NEB(enz1,enz2,filtre_all)
        if incub1 != incub2:
            superDico['temp_incub'][enz1]=enz2
            continue
        if inac1 != inac2 :
            superDico['temp_inactiv'][enz1]=enz2
            continue
        if buf1 != buf2 :
            superDico['pas_meme_buffer'][enz1]=enz2
            continue
        if sit1 == sit2 :
            superDico['meme_site'][enz1]=enz2
            continue
        pairToKeep.append((enz1,enz2))

    return(superRb,pairToKeep,superDico,enzHF)


def pair_rb_NEB(tab_supp,superDico,pairToKeep):
    df = pd.read_csv(tab_supp,dtype=str) 
    filtre_all = modif_NEB.modif_tab(df)
    enzHF = modif_NEB.recup_enz(filtre_all)
    rb = modif_NEB.modif_rb(enzHF)
    keylist = restriction.create_keylist_multi(rb)
    keylist = modif_NEB.modif_keylist_multi(keylist,enzHF)
    for enz1,enz2 in itertools.combinations(keylist,2):
        incub1,incub2,inac1,inac2,buf1,buf2,sit1,sit2 = modif_NEB.filtre_NEB(enz1,enz2,filtre_all)
        if incub1 != incub2:
            superDico['temp_incub'][enz1]=enz2
            continue
        if inac1 != inac2 :
            superDico['temp_inactiv'][enz1]=enz2
            continue
        if buf1 != buf2 :
            superDico['pas_meme_buffer'][enz1]=enz2
            continue
        if sit1 == sit2 :
            superDico['meme_site'][enz1]=enz2
            continue
        pairToKeep.append((enz1,enz2))

    return(rb,pairToKeep,superDico,enzHF)
 

def superRb_other_supp(supp,sites_variants,bout_francs,enzaway,superRb):
    rb = Res.RestrictionBatch(first=[], suppliers = [supp])
    rb,sites_variants = restriction.site_variant(rb,sites_variants)       #enlever les sites variants pour n'importe quelle enzyme
    for enz in sites_variants:
        enzaway.writerow([enz]+['site variant'])
    rb,bout_francs = restriction.blunt(rb,bout_francs)                    #enlever les enzymes à bouts francs
    for enz in bout_francs:
        enzaway.writerow([enz] + ['bouts francs'])
    superRb = restriction.make_superRb(rb,superRb)
    return(superRb,sites_variants,bout_francs)



