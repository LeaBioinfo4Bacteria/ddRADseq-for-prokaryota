#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 15:03:39 2021

@author: lmasson
"""

import pandas as pd
import os

liste_de_csv = []

for (root,dirs,files) in os.walk('resultats_ddRADSeq', topdown= True):
    if 'couples_et_frag.csv' in files:
        name = root.split('/')[1]
        frame = pd.read_csv(root+'/couples_et_frag.csv',index_col=None,header=0)
        frame['Seq']=str(name)
        liste_de_csv.append(frame)

#print("len",len(liste_de_csv))
df= pd.concat(liste_de_csv,axis = 0, ignore_index=True)


dico_all = {}
for (index,data) in df.iterrows():
    data['enz1_enz2'] = data['enz1_enz2'].replace('"','').replace('(','').replace(')','').replace("'",'').split(',')
    coupl_enz = tuple(set(data['enz1_enz2'])) #tuple
    temp = dico_all.get(coupl_enz,{})
    temp[data['Seq']] = data['nb_frag']
    dico_all[coupl_enz] = temp

res = sorted(dico_all, key=lambda k: len(dico_all[k]), reverse=True)


res_to_write = {}
for couple in res:
    for k,v in dico_all.items():
        if couple == k:
            res_to_write[couple]=v

my_df = pd.DataFrame.from_dict(res_to_write,orient='index')
my_df.to_csv('res_globaux.csv')
