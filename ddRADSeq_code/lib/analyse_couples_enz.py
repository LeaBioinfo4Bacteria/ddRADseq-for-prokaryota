#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib analyse des couples d'enzymes
"""

#import Bio.Restriction as Res
from Bio.Restriction import *
from collections import defaultdict

def ratio_freq(enz1,enz2,supp):
    """calcul du ratio>1 avec les freq du module rb
    pas used dans le main"""
    if supp != None:
        rb = RestrictionBatch(first=[], suppliers = [supp])
        if '-v2' or 'HF' in enz1:
            enz1 = str(enz1).replace('-v2','').replace('-HF','')
            for enz in rb :
                if str(enz1) == str(enz) :
                    enz1 = enz
        if '-v2' or 'HF' in enz2:
            enz2 = str(enz2).replace('-v2','').replace('-HF','')
            for enz in rb:
                if str(enz2) == str(enz):
                    enz2 = enz

    freq1 = enz1.frequency()
    freq2 = enz2.frequency()
    if freq1 < freq2:
        return(freq2/freq1)
    else:
        return(freq1/freq2)

def get_normal_version(enz):
    '''pour les assigner les HFv2 à leur version normale'''
    rb = RestrictionBatch(first=[], suppliers = [])
    if '-v2' in str(enz) or 'HF' in str(enz):
        return rb.format(str(enz).replace('-v2','').replace('-HF',''))
    else:
        return rb.format(enz)

def diff_neo_compa(enz1,enz2):
    '''renvoie si les sites de coupures sont différents (diff),
    si les enz sont des neoschizomères (neoschi)
    et si leurs fragments peuvent être liés (compa)'''

    enz1 = get_normal_version(enz1)
    enz2 = get_normal_version(enz2)

    diff = (enz1 != enz2)
    neo = (enz1 >> enz2)
    compa = (enz1 % enz2)
    return diff, neo, compa


def freq_mono(enz1,enz2,dicorb):
    freqDico = {}
    for key,val in dicorb.items() :
        freqDico[key]= freqDico.get(key, 0) +len(val)
    for key,val in freqDico.items():
        if str(key) == str(enz1):
            freq1 = freqDico[key]
        if str(key) == str(enz2):
            freq2 = freqDico[key]

    if freq1 < freq2 :
        return(freq2/freq1)
    else :
        return(freq1/freq2)


def dico_freq_multi(liste_seq):
    """calcul du ratio>1, fait la somme des sites de coupures dans la seq entière, 
    retourne le dico des frequences freqDico"""
    freqDico = {}
    for seq in liste_seq :
        for key,val in seq.res_search.items() :
            freqDico[key]= freqDico.get(key, 0) +len(val)
    return(freqDico)


def ratio_freq_multi(enz1,enz2,freqDico):
    '''récupère les values du dico en fonction de l'enzyme,
    donc la fréquence de coupure dans la séquence,
    return le ratio <1'''
    for key,val in freqDico.items():
        if str(key) == str(enz1):
            freq1 = freqDico[key]
        if str(key) == str(enz2):
            freq2 = freqDico[key]

    if freq1 < freq2 :
        return(freq2/freq1)
    else :
        return(freq1/freq2)
