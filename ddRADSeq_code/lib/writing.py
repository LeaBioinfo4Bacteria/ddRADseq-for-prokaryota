#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib écriture
créer tous les fichiers resultats dans le dossier seq res,
retourne les writers
"""
import csv
import os
from Bio import SeqIO
import itertools


def info_frag(dir_name):
    '''info sur chaque fragment obtenu par couple'''
    dir_res = os.path.join(dir_name,'infos_frag.csv')
    csvres = open(dir_res,"w")
    writer = csv.writer(csvres,delimiter=",")
    fieldnames = ["id de la seq","enz1","enz2","nb total de frag",'taille de tous les frags',"nb fragment d interet","taille","position","sequence"]
    writer.writerow(fieldnames)
    return(writer)



def enz_keep(dir_name):
    '''couple_gardes'''
    dir_res = os.path.join(dir_name,'couples_gardés.csv')
    enzkeep = csv.writer(open(dir_res,"wt",newline=''),delimiter=",")
    enzkeep.writerow(['ID de la seq']+['enz1']+['enz2']+['nb de frag']+['nb site enz1']+['nb site enz2']+['pos site enz1']+['pos site enz2'])
    return(enzkeep)



def enz_away(dir_name):
    '''enzymes enlevé , pas celles du tab_neb!'''
    dir_res = os.path.join(dir_name,'enzymes_non_gardées.csv')
    enzaway = csv.writer(open(dir_res,"wt",newline=''),delimiter=",")
    return(enzaway)



def pair_away(dir_name):
    '''couples non gardes'''
    dir_res = os.path.join(dir_name,'couples_non_gardés.csv')
    pairaway = csv.writer(open(dir_res,"wt",newline=''),delimiter=",")
    return(pairaway)



def write_seq(dir_name,enz1,enz2,pos_fragtoseq,seq,sequence):
    '''les seq de chaque frag'''
    dir_res = os.path.join(dir_name,'sequence des frag avec '+str(enz1)+' '+str(enz2)+'.fasta')
    old_name = seq.id
    name = seq.id.replace('.fasta','')
    seq.description = seq.description.replace(', complete cds','')
    records_set = []
    for a,b in pos_fragtoseq:
        seq.id = name+'_'+str(a)+'_'+str(b)
        if b > a:
            records_set.append(seq[a:b])
        elif a > b:
            records_set.append(seq[a:0:b])
    SeqIO.write(records_set,dir_res,'fasta')
    seq.id = old_name



def write_superDico(superDico,pairaway):
    for (i,d) in superDico.items():
        for k,v in d.items():
            pairaway.writerow((k,v,i))



def res_globaux(dico_res,enz1,enz2,liste):
    '''nb de site de coupure par couple d'enz dans l'ensemble du fasta,
    return un dictionnaire'''
    if (enz1,enz2) in dico_res.keys():
        dico_res[(enz1, enz2)] = dico_res.get((enz1, enz2), 0) + len(liste)
    elif (enz2,enz1) in dico_res.keys():
        dico_res[(enz2, enz1)] = dico_res.get((enz2, enz1), 0) + len(liste)
    else:
        dico_res[(enz1, enz2)] = len(liste)
    return(dico_res)



def write_dico_res(dir_name,superDicoRes):
    '''le dico de res globaux sur le fasta'''
    fieldnames=['enz1_enz2','nb_frag']
    dir_res = os.path.join(dir_name,'couples_et_frag.csv')
    couples_et_frag = csv.writer(open(dir_res,"wt",newline=''),delimiter=',')
    couples_et_frag.writerow(fieldnames)
    shared_keys = {}
    for dict1,dict2 in itertools.combinations((superDicoRes.values()),2):
        keys = set(dict1.keys()) & set(dict2.keys())
        for k in keys:
            shared_keys[k] = dict1.get(k,0)+dict2.get(k,0)
    for ((k,v),i) in shared_keys.items():
        couples_et_frag.writerow(((k,v),i))


#regex pour enlever les .fasta et autre
