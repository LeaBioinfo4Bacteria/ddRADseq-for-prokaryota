#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib analyse des fragments obtenus
"""

def list_pos_enz(enz1,enz2,dico):
    lenz1 = []
    lenz2 = []
    lenz1.extend(dico[enz1])
    lenz2.extend(dico[enz2])
    lenz = []
    lenz = sorted(lenz1 + lenz2)
    return(lenz1,lenz2,lenz)


def all_frags_lenght(lenz,lin,seq):
    '''renvoie la taille de tous les fragments possibles'''
    all_frags_lenght = []
    for pos in range(len(lenz)-1):
        x = (lenz[pos+1] - lenz[pos])
        all_frags_lenght.append(x)
    if lin == False:
       all_frags_lenght.append((len(seq.seq)-lenz[-1]) + lenz[0])
    
    all_frags_lenght.sort()

    return(all_frags_lenght)



def frag_dd(lenz,lenz1,lenz2):
    """récupère les fragments dont chaque extremité provient d'une enzyme du couple"""
    lout = []   
    for i in range(0,len(lenz)) :
        a = lenz[i]
        if i+1 == len(lenz) :
            b = lenz[0]
        else : 
            b = lenz[i+1]              
        if a in lenz1 and b in lenz2 :
            lout.append(tuple([a,b]))
        elif b in lenz1 and a in lenz2 :
            lout.append(tuple([a,b]))

    return(lout)


def longueur_position(lout, seq, lenmin, lenmax):
    '''calcule la longueur des fragments puis fait la sélection de taille,
    return les listes: taille_fragtokeep et fragtoseq'''
    long_frag = []
    taille_fragtokeep = []
    pos_fragtoseq = []   
    for a,b in lout :
        if a < b :
            long_frag.append(b-a)
        if b < a :
            long_frag.append((len(seq.SeqRecord.seq)-a)+b)

        for frag in long_frag :
            if lenmin < frag < lenmax :
                if frag not in taille_fragtokeep :
                    taille_fragtokeep.append(frag)
                    pos_fragtoseq.append(tuple([a,b]))

    return(taille_fragtokeep,pos_fragtoseq)



def sequence(seq,pos_fragtoseq):
    '''retourne la séquence des fragments conservés'''
    seqfrag = []
    for a,b in pos_fragtoseq :
        if a < b :
            seqfrag.append(seq[a:b])
        if b < a : 
            seqfrag.append(seq[b:a])
    return(seqfrag)

