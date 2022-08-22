#remove three hashtags to store data in files


import biotite.database.rcsb as rcsb
import datetime
from pymol import *
import collections
from tqdm.auto import tqdm
import os
import matplotlib.pyplot as plt
import numpy as np

import streamlit as st

import pandas as pd
import shutil

def get_spike_ids(uniprot_id="P0DTC2", min_weight=400, max_resolution=4.0):
    """
    get all pdbs with defined weight and resolution,
    input the uniprot_id (the default is spike), min_weight , and max_resolution
    """
    # uniprot_id = "P0DTC2" #spike in Sars-cov-2
    # max_resolution = 4.0
    # min_weight =400
    """
    in Da, structure min mass to get rid of rbd only structures,
    Spike mass is 429 Da.
    """
    query_by_uniprot_id = rcsb.FieldQuery(
        "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        exact_match=uniprot_id,
    )
    today = datetime.datetime.now()

    query_by_resolution = rcsb.FieldQuery(
        "rcsb_entry_info.resolution_combined", less_or_equal=max_resolution
    )

    query_by_polymer_weight = rcsb.FieldQuery(
        "rcsb_entry_info.molecular_weight", greater=min_weight
    )

    query = rcsb.CompositeQuery(
        [
            query_by_uniprot_id,
            query_by_resolution,
            # query_by_polymer_count,
            query_by_polymer_weight,
        ],
        "and",
    )
    pdb_ids = rcsb.search(query)

    # remove post fusion strcuture 6xra

    pdb_ids.remove('6XRA')
    #print(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
    #      f"resolution less than or equal to {max_resolution} with mass more than or equal to {min_weight}: {len(pdb_ids)}")
    #print("Selected PDB IDs:\n", *pdb_ids)
    st.write(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
          f"resolution less than or equal to {max_resolution}A with mass more than or equal to {min_weight}: {len(pdb_ids)}")
    return (pdb_ids)



def distance_dif(pdb_ids, resid_1,  resid_2, atom_1, atom_2):
        if not os.path.exists('PDB'):
            os.mkdir('PDB')
        if not os.path.exists('error_residue'):
            os.mkdir('error_residue')


        dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
        dist_list_reverse = collections.defaultdict(list) # empty dictionary for future rmsd from reverse
        missing_residue = []
        missing_residue_reverse = []
        error_pdbs = []


        # sanity check that the numbering is correct, check if 1000 is ARG in chain A B C
        for i in tqdm(pdb_ids):
            chains_ordered = ['A']
            cmd.delete("*")
            i = i.lower()

        #load structure
            if os.path.exists('./PDB/' + i + '.pdb'):
                cmd.load("./PDB/" + i + ".pdb")
            elif os.path.exists('./PDB/' + i + '.cif'):
                cmd.load("./PDB/" + i + ".cif")
            else:
                file = cmd.fetch(i, path='./PDB/', type='pdb')
                # st.write(file)
                if file != i:
                    cmd.fetch(i, path='./PDB/')
                # st.write("cif instead of pdb is fetched")
                #st.write('load cif from folder')

            error_1000 = False
            for item in ['A', 'B', 'C']:
                p1 = cmd.select("p1", f'chain {item} and i. 1000 and r. ARG and n. CA')
                if p1 != 1:
                    error_pdbs.append(i.upper())
                    error_1000 = True
                    break
            if error_1000:
                continue

        #measure 2 distances to find chains orientation (clockwise/counterclockwise)
            dist1 = cmd.get_distance(atom1=f'chain A and i. 971 and n. CA',
                                     atom2=f'chain B and i. 752 and n. CA')

            dist2 = cmd.get_distance(atom1=f'chain A and i. 971 and n. CA',
                                     atom2=f'chain C and i. 752 and n. CA')
            if min(dist1, dist2) == dist1:
                chains_ordered.append('B')
                chains_ordered.append('C')
            elif min(dist1, dist2) == dist2:
                chains_ordered.append('C')
                chains_ordered.append('B')

            resid_1 = resid_1.upper()
            resid_2 = resid_2.upper()
        # measure the distance of interest
            for j in range(0,3):
                try:
                    if j == 2 :
                        dist = cmd.get_distance(atom1=f'chain {chains_ordered[j]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[0]} and i. {resid_2} and n. {atom_2}')
                        dist_list[i].append(dist)
                    else:

                        dist = cmd.get_distance(atom1=f'chain {chains_ordered[j]} and i. {resid_1} and n. {atom_1}',
                                            atom2=f'chain {chains_ordered[j + 1]} and i. {resid_2} and n. {atom_2}')
                        dist_list[i].append(dist)
                except CmdException:
                    missing_residue.append(i.upper())
                    #break

            #test for the second plot in other direction
            for j in range(0, 3):
                try:
                    if j == 2:
                        dist_reverse = cmd.get_distance(atom1=f'chain {chains_ordered[0]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[j]} and i. {resid_2} and n. {atom_2}')
                        dist_list_reverse[i].append(dist_reverse)
                    else:

                        dist_reverse = cmd.get_distance(atom1=f'chain {chains_ordered[j + 1]} and i. {resid_1} and n. {atom_1}',
                                                atom2=f'chain {chains_ordered[j]} and i. {resid_2} and n. {atom_2}')
                        dist_list_reverse[i].append(dist_reverse)
                except CmdException:
                    missing_residue_reverse.append(i.upper())
                    # break


        st.header('**Incorrectly numbered pdbs**')

        st.write(f'There are {len(error_pdbs)} structures with 1000 ARG error in at least one chain:', str(error_pdbs))
        ###f = open("IncorrectNumberingPDB.txt", "w")
        ###f.write(str(error_pdbs))
        ###f.close()

        st.write(f'There are {len(missing_residue)} chains with a problem in chosen residue:', str(missing_residue))
        error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + '.txt'
        ###f = open(error_file_name, "w")
        ###f.write(str(missing_residue))
        ###f.write(str(missing_residue_reverse))
        ###f.close()
        return dist_list, dist_list_reverse

def distance_same(pdb_ids, resid_1,  resid_2, atom_1, atom_2, flag):
    if not os.path.exists('PDB'):
        os.mkdir('PDB')

    dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
    missing_residue = []
    error_pdbs = []
    if not os.path.exists('error_residue'):
        os.mkdir('error_residue')

    for i in tqdm(pdb_ids):
        cmd.delete("*")
        i = i.lower()

        if os.path.exists('./PDB/' + i + '.pdb'):
            cmd.load("./PDB/" + i + ".pdb")
        elif os.path.exists('./PDB/' + i + '.cif') :
            cmd.load("./PDB/" + i + ".cif")
        else:
            file = cmd.fetch(i, path='./PDB/', type='pdb')
            #st.write(file)
            if  file != i:
                cmd.fetch(i, path='./PDB/')
            #st.write("cif instead of pdb is fetched")

        error_1000 = False
        for item in ['A', 'B', 'C']:
            # sanity check that the numbering is correct, check if 1000 is ARG in chain A B C
            p1 = cmd.select("p1", f'chain {item} and i. 1000 and r. ARG and n. CA')
            if p1 != 1:
                error_pdbs.append(i.upper())
                error_1000 = True
                break
        if error_1000:
            continue


        resid_1 = resid_1.upper()
        resid_2 = resid_2.upper()
        for chain in ['A', 'B', 'C']:
            try:
                dist = cmd.get_distance(atom1=f'chain {chain} and i. {resid_1} and n. {atom_1}',
                                            atom2=f'chain {chain} and i. {resid_2} and n. {atom_2}')
                dist_list[i].append(dist)
            except CmdException:
                missing_residue.append(i)

    st.header('**Incorrectly numbered pdbs**')

    st.write(f'There are {len(error_pdbs)} structures with 1000 ARG error in at least one chain:', str(error_pdbs))
    ###f = open("IncorrectNumberingPDB.txt", "w")
    ###f.write(str(error_pdbs))
    ###f.close()

    st.write(f'There are {len(missing_residue)} chains with a problem in chosen residue:', str(missing_residue))
    error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + '.txt'
    ###f = open(error_file_name, "w")
    ###f.write(str(missing_residue))
    ###f.close()
    return dist_list


def analysis(distancesDict, resid_1, atom_1, resid_2, atom_2, flag):
    """
    plot the histogram of distances
    """
    distances_only = list(distancesDict.values())

    st.write(f'Number of structures with at least one analyzed chain {len(distancesDict)}')

    if not os.path.exists('plots'):
        os.mkdir('plots')

    if not os.path.exists('distances'):
        os.mkdir('distances')

    fig = plt.figure(figsize=(15, 7.5))
    plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
    plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
    # sns.histplot(data=distances_only , binwidth=0.2)
    plt.xlabel('distance, A', fontsize=20)
    plt.ylabel("Count", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    st.pyplot(fig)

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in distancesDict.items()])).transpose()
    name = str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2)
    df_name = './distances/distance_' + name + '_' + flag +'.csv'
    ###df.to_csv(df_name)
    ###st.write(df)
    ###plt_name = './plots/distance_' + name + '_' + flag +'.png'
    ###plt.savefig(plt_name, bbox_inches='tight')


    #print(f'number of corrected numbered and analyzed structures {len(distancesDict)}')
    #return distances_only

def delete_PDB_folder():
        shutil.rmtree('./PDB')

#pdb_ids = get_spike_ids()
##dist = distance(pdb_ids, 318, 'PHE', 292,  'ALA')
#dist = distance(pdb_ids, 1000, 'ARG', 975,  'SER')
#analysis(dist)

