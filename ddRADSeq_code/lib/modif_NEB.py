#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lib modif rb pour NEB
"""

import Bio.Restriction as Res

def modif_tab(df):
    """filtre le tableau : garde les enzymes avec séquence invariante,
    sans activité star, pas impactée par méthylation,
    retourne le tableau filtré pour les température plus tard"""
    df.columns = (df.columns.str.strip().str.replace(' ','_').str.replace('(','').str.replace(')','').str.replace('.',''))
    df = df.drop(['Enzyme1',"Unnamed:_0",'Dam','Dcm','CpG','Unit_Substrate'],axis = 'columns')
    df = df[df.Sequence.apply(lambda x: x.replace('A', '').replace('T', '').replace('G', '').replace('C', '').replace('/', '') == "")]
    df = df[df.Heat_Inac.apply(lambda x : 'No' not in str(x))]
    df = df[df.Notes.apply(lambda x: ("1" not in str(x))&("2" not in str(x))&("3" not in str(x))&("e" not in str(x)))]    
    df = df[df.Sequence.apply(lambda x: len(x.partition('/')[0]) != len(x.partition('/')[2]))]
    df['Enzyme'] = df['Enzyme'].str.replace(' §','').str.replace('®','').str.strip()    
    return(df)


def recup_enz(filtre_all):
    '''retourne les enz du tableau filtré sous forme de liste'''
    enzyme = filtre_all.get("Enzyme")
    enzHF = []
    enzHF.extend(enzyme)
    enzHF.sort()
    return(enzHF)


def modif_rb(enzHF):
    """modifie le rb du supplier en fonction de la liste d'enzyme filtrée, 
    retourne le rb à utiliser"""
    rb = Res.RestrictionBatch([])
    enz = []
    HFv2 = []
    for e in enzHF :
        try :
            rb.add(e)
            enz.append(e)
        except :
            if '-HF' or '-v2' in e:
                n = e.replace('-v2','').replace('-HF','')
                rb.add(n)
                HFv2.append(e)    
    return(rb)


def modif_dico_NEB(dicorb,enzHF):
    '''rajoute les enzyme HFv2 au dico en assignant les values correpondantes
    aux non HF (car même site de coupure et value = positions des sites dans
    la seq à analyser)
    renvoie le dicorb avec les HF et v2 + la liste des enz enlevées'''
#    enz_suppr_du_rb = []
    dicorbcop= dicorb.copy()
    for e in enzHF:
        for n,z in dicorbcop.items():
             if e not in dicorb.keys():
                 y =  e.replace('-v2','').replace('-HF','')
                 if str(n) == str(y) :
                     dicorb[e]=z                 
             if n not in enzHF and n in dicorb:
#                 enz_suppr_du_rb.append(str(n))
                 del dicorb[n]
    return(dicorb)



def modif_keylist_multi(keylist,enzHF):
    '''rajoute les enzyme HFv2 à la keylist,
    renvoie la keylist avec les HF et v2 à garder'''
    keylistcop= keylist.copy()
    for e in enzHF:
        for n in keylistcop:
             if e not in keylist:
                 y =  e.replace('-v2','').replace('-HF','')
                 if str(n) == str(y) :
                     keylist.append(e)                 
             if n not in enzHF and n in keylist:
                 keylist.remove(n)
    return(keylist)



def filtre_NEB(enz1, enz2, filtre_all):
    '''renvoie les températures d'incubation et d'inactivation de chaque enzyme,
    et les buffers livrés avec les enzymes '''
    enz1_index = int(filtre_all.index.values[filtre_all['Enzyme'] == str(enz1)])
    enz2_index = int(filtre_all.index.values[filtre_all['Enzyme'] == str(enz2)])
    incub1 = filtre_all.at[enz1_index,'Incu_Temp']
    incub2 = filtre_all.at[enz1_index,'Incu_Temp']
    inac1 = filtre_all.at[enz1_index,'Heat_Inac']
    inac2 = filtre_all.at[enz2_index,'Heat_Inac']
    buf1 = filtre_all.at[enz1_index,'Supplied_NEBuffer']
    buf2 = filtre_all.at[enz2_index,'Supplied_NEBuffer']
    sit1 = filtre_all.at[enz1_index,'Sequence']
    sit2 = filtre_all.at[enz2_index,'Sequence']
    return(incub1,incub2,inac1,inac2,buf1,buf2,sit1,sit2)
    



