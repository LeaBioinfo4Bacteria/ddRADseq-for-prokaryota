#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:53:50 2021

@author: lmasson
comments are in french !

this code was written to determine in silico which couple of enzymes to use on a given sequence in order to use the ddRADseq method to genotype it
libraries used : argparse, BioPython, SeqIO, Restriction, itertools, pandas, csv, sys, os
"""
import argparse
from lib import restriction
from lib import modif_NEB
from lib import analyse_couples_enz
from lib import analyse_frag
from lib import writing
from lib.crea_object import Seq
from lib import group_func
from Bio import SeqIO
import itertools
import pandas as pd
import Bio.Restriction as Res
import csv
import sys
import os

#fonctionne avec des .fasta!!! pas fna ou fa

### argparse ###

parser = argparse.ArgumentParser(description=
                                 "Donne les produits attendus de double digestion à partir des enzyme de NEB et de la séquence donnée.")

#parser.add_argument('-seq','--sequence',default = 'PAO1.fasta', help='le fichier contenant la séquence à analyser')
parser.add_argument('-seq','--sequence', help='le fichier contenant la séquence à analyser')
parser.add_argument('-r', '--ratio', default=3, help="le ratio des fréquences des couples d'enzymes, /!\ doit être > 1")
parser.add_argument('-f', '--fragment', default=0, help="le nombre minimum de fragments voulus")
parser.add_argument('-lmin', '--lenmin', default=200, help="la taille, en pb, min à sélectionner",type=int)
parser.add_argument('-lmax', '--lenmax', default=800, help="la taille, en pb, max à sélectionner",type=int)
parser.add_argument('-lin','--linear',default=False,help='génome lineaire ou pas')

parser.add_argument('-sup','--supplier', nargs = '*',default='N',help = 'entrer la lettre code du supplier:\n' +\
                    'B = Life Technologies. \n'+
                    'C = Minotech Biotechnology. \n'+
                    'E = Agilent Technologies. \n'+
                    'I = SibEnzyme Ltd. \n'+
                    'J = Nippon Gene Co., Ltd. \n'+
                    'K = Takara Bio Inc. \n'+
                    'M = Roche Applied Science. \n'+
                    'N = New England Biolabs. \n'+
                    'O = Toyobo Biochemicals. \n'+
                    'Q = Molecular Biology Resources - CHIMERx. \n'+
                    'R = Promega Corporation. \n'+
                    'S = Sigma Chemical Corporation. \n'+
                    'V = Vivantis Technologies. \n'+
                    'X = EURx Ltd. \n'+
                    'Y = SinaClon BioScience Co. ')
parser.add_argument('-t', '--tab_supp', default="nebuffer-performance-chart-with-restriction-enzymes.csv", help="le tableau csv des enzymes de NEB")

args = parser.parse_args()
print("séquence à analyser:",args.sequence)


######################################################################## main ###

################################ sequence independant ##########################"


#mes listes/dico pour voir qui part quand
bout_francs=[]
sites_variants=[]
couple_keep = []
dico_res = {}

superDico = {
    'pas_de_frag':{},
    'pas_assez_de_frag':{},
    'mauvais_ratio': {},
    'temp_incub': {},
    'temp_inactiv':{},
    'pas_meme_buffer':{},
    'meme_site':{},
    'pas_de_frag_d_interet':{},
    'neoschizomeres':{},
    'fragments_cohesifs':{}
}

superDicoRes={}


#créa des dossiers de résutats
direct = 'resultats_ddRADSeq_souchier_cop'+'/'+os.path.basename(args.sequence.replace('.fasta',' res'))
os.makedirs(direct,exist_ok = True)

#writer pour mes différents fichiers résultats
writer = writing.info_frag(direct)
enzkeep = writing.enz_keep(direct)
enzaway = writing.enz_away(direct)     #estce que ça ce serait pas un df dérivé du tab neb en direct si supp == N?
pairaway = writing.pair_away(direct)


pairToKeep = []

if args.supplier != None:
    if args.supplier == "N":        #modif le rb si supplier = NEB
        rb,pairToKeep,superDico,enzHF = group_func.pair_rb_NEB(args.tab_supp,superDico,pairToKeep)
    elif args.supplier != "N":      #sinon, garder le rb du module restriction batch
        rb,sites_variants,bout_franc = group_func.superRb_other_supp(args.supplier,sites_variants,bout_francs,enzaway,rb)
elif args.supplier == None :    #pour faire tourner que sur nos 5 + HF
    rb,sites_variants,bout_francs,enzHF = restriction.rb_sans_supp(sites_variants,enzaway,bout_francs)

if len(rb) == 0:
    print("le rb est vide!")
    sys.exit()

#récup les clés du dico = keylist pour tester les couples puis les keep
keylist = restriction.create_keylist_multi(rb)
if args.supplier == "N" or None:
    keylist = modif_NEB.modif_keylist_multi(keylist,enzHF)


#garder les couples fonctionnels dans une liste de tuple
for enz1,enz2 in itertools.combinations(keylist,2):
    #différents sites, neoschizomers,compatibility, ratio
    diff, neo, compa = analyse_couples_enz.diff_neo_compa(enz1,enz2)
    if diff == False:
        continue
    if neo == True:
        superDico['neoschizomeres'][enz1]=enz2
        continue
    if compa == True:
        superDico['fragments_cohesifs'][enz1]=enz2
        continue

    if (enz1,enz2) not in pairToKeep:
        pairToKeep.append((enz1,enz2))



#import pdb
#pdb.set_trace()
############################################## seq dependant ###############################

#ouverture de la seq à analyser
#search for site de restriction des enz du rb dans la seq
#liste de dico resultats par séquence
fic = open(args.sequence,'r')

liste_seq = []
for seqRecord in SeqIO.parse(fic, 'fasta'):
    dico = rb.search(seqRecord.seq, linear = args.linear)
    if args.supplier == "N":
        dico = modif_NEB.modif_dico_NEB(dico,enzHF)
    sequence = Seq(seqRecord,dico)
    liste_seq.append(sequence)



freqDico = analyse_couples_enz.dico_freq_multi(liste_seq)

#la loop pour tester les enz du rb deux à deux
for seq in liste_seq:
    pas_de_site = restriction.pas_de_site(seq.res_search)
    superDicoRes[seq.SeqRecord.id]={}

    for enz1,enz2 in pairToKeep:
        if enz1 in pas_de_site or enz2 in pas_de_site:
            continue

        lenz1,lenz2,lenz = analyse_frag.list_pos_enz(enz1,enz2,seq.res_search)

#parcourir ma list de freq
        if analyse_couples_enz.ratio_freq_multi(enz1,enz2,freqDico) < args.ratio :
            superDico['mauvais_ratio'][enz1]=enz2
            continue

        all_frags_lenght = analyse_frag.all_frags_lenght(lenz,args.linear,seq.SeqRecord)

    #vérif si le couple donne bien des fragments , puis le min voulu
        lout = analyse_frag.frag_dd(lenz,lenz1,lenz2)
        if len(lout) == 0 :
            superDico['pas_de_frag'][enz1]=enz2
            continue
        if len(lout) < args.fragment:
            superDico['pas_assez_de_frag'][enz1]=enz2
            continue

    #vérif les couples d'enzymes pour tout suppliers

        taille_fragtokeep,pos_fragtoseq = analyse_frag.longueur_position(lout,seq,args.lenmin,args.lenmax)
        seqfrag = analyse_frag.sequence(seq.SeqRecord,pos_fragtoseq)

        if len(pos_fragtoseq) == 0:
            superDico['pas_de_frag_d_interet'][enz1]=enz2
            continue

    #écriture dans 4 fichiers res:
        if (enz1,enz2) not in couple_keep:
            couple_keep.append([enz1,enz2])

        enzkeep.writerow([seq.SeqRecord.id]+[enz1,enz2]+[len(pos_fragtoseq)]+[len(lenz1)]+[len(lenz2)]+[lenz1]+[lenz2])
        writing.write_seq(direct,enz1,enz2,pos_fragtoseq,seq.SeqRecord,args.sequence)

        for i in range(0,len(pos_fragtoseq)):
            writer.writerow([seq.SeqRecord.id]+[enz1,enz2] + [len(all_frags_lenght)] + [all_frags_lenght] + [len(pos_fragtoseq)] + [taille_fragtokeep[i]] + [pos_fragtoseq[i]] + [seqfrag[i]])

        dico_res = writing.res_globaux(dico_res,enz1,enz2,pos_fragtoseq)
        superDicoRes[seq.SeqRecord.id]=dico_res

    for enz in pas_de_site:
        enzaway.writerow([enz] + ['pas de site de coupure dans la sequence'] + [seq.SeqRecord.id])


writing.write_superDico(superDico,pairaway)
writing.write_dico_res(direct,superDicoRes)

fic.close()
