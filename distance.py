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
import base64
from PIL import Image

# File download
def filedownload(df, df_name):
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download={df_name}>Download distribution</a>'
    return href


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

    query_by_method = rcsb.FieldQuery(
    	"exptl.method", exact_match= "ELECTRON MICROSCOPY"
    )
    
    query = rcsb.CompositeQuery(
        [
            query_by_uniprot_id,
            query_by_resolution,
            # query_by_polymer_count,
            query_by_method,
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



def distance_dif(pdb_ids, resid_1,  resid_2, atom_1, atom_2, mutation_name = '', mutation_id = ''):
        if not os.path.exists('PDB'):
            os.mkdir('PDB')
        if not os.path.exists('error_residue'):
            os.mkdir('error_residue')


        dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
        dist_list_reverse = collections.defaultdict(list) # empty dictionary for future rmsd from reverse
        missing_residue = []
        missing_residue_reverse = []
        error_pdbs = []
        with_mutant = []


        for i in tqdm(pdb_ids):
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


            # sanity check that the numbering is correct, check if 1126 is CYS  in chain A B C
            # assign chain names
            chains = []
            for ch in cmd.get_chains(i):

                '''#print(i, " has chain ", ch)
                p1 = cmd.select("p1",')
                if p1 != 1:
                    #error_1000 = True
                    continue
                else:
                    chains.append(ch)'''

                try:
                    dist = cmd.get_distance(atom1= f'chain {ch} and i. 1126 and r. CYS and n. CA',
                                            atom2=f'chain {ch} and i. 57 and r. PRO and n. CA')
                    chains.append(ch)
                except CmdException:
                    continue


            #print(i , 'has ', chains)

        # only take those that had correct sequence numbering in all 3 chains
            if len(chains) !=3 :
                error_pdbs.append(i.upper())
                continue

            # check if the sequence has a requested mutation
            if mutation_name != '':
                p2 = cmd.select("p2", f'chain {str(chains[0])} and i. {mutation_id} and r. {mutation_name.upper()} and n. CA')
                if p2 != 1:
                    continue
                else:
                    #print('mutant')
                    with_mutant.append(i.upper())


        #measure 2 distances to find chains orientation (clockwise/counterclockwise)
            chains_ordered = [chains[0]]
            dist1 = cmd.get_distance(atom1=f'chain {chains[0]} and i. 971 and n. CA',
                                     atom2=f'chain {chains[1]} and i. 752 and n. CA')

            dist2 = cmd.get_distance(atom1=f'chain {chains[0]} and i. 971 and n. CA',
                                     atom2=f'chain {chains[2]} and i. 752 and n. CA')
            if min(dist1, dist2) == dist1:
                chains_ordered.append(chains[1])
                chains_ordered.append(chains[2])
            elif min(dist1, dist2) == dist2:
                chains_ordered.append(chains[2])
                chains_ordered.append(chains[1])

            #print(chains_ordered)
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

        st.write(f'There are {len(error_pdbs)} structures with 1126 not being CYS or 57 not being a PRO in at least one chain:', str(error_pdbs),'. Removed those from analysis.')
        f = open("IncorrectNumberingPDB.txt", "w")
        f.write(str(error_pdbs))
        f.close()

        if mutation_name != '':
            st.write(f'There are {len(with_mutant)} structures with {mutation_id}{mutation_name.upper()}:', str(with_mutant))

        st.write(st.write(f'There are {len(missing_residue)} chains (in {len(set(missing_residue))} structures) with missing at least one of the requested atoms:', str(missing_residue)))
        error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + str(mutation_name) + str(mutation_id) + '.txt'
        f = open(error_file_name, "w")
        f.write(str(missing_residue))
        f.write(str(missing_residue_reverse))
        f.close()
        return dist_list, dist_list_reverse

def distance_same(pdb_ids, resid_1,  resid_2, atom_1, atom_2, mutation_name = '', mutation_id = ''):
    if not os.path.exists('PDB'):
        os.mkdir('PDB')

    dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
    missing_residue = []
    error_pdbs = []
    with_mutant = []

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


# chains naming and taking only chains that have CYS 1126
        chains = []
        for ch in cmd.get_chains(i):
            '''#print(i, " has chain ", ch)
            p1 = cmd.select("p1", f'chain {ch} and i. 1126 and r. CYS and n. CA')
            if p1 != 1:
                    #error_1000 = True
                continue
            else:
                chains.append(ch)
            '''
            try:
                dist = cmd.get_distance(atom1= f'chain {ch} and i. 1126 and r. CYS and n. CA',
                                            atom2=f'chain {ch} and i. 57 and r. PRO and n. CA')
                chains.append(ch)
            except CmdException:
                continue
        #print(i , 'has ', chains)

        if len(chains) !=3 :
            error_pdbs.append(i.upper())
            continue

            # check if the sequence has a requested mutation
        if mutation_name != '':
            p2 = cmd.select("p2", f'chain {str(chains[0])} and i. {mutation_id} and r. {mutation_name.upper()} and n. CA')
            if p2 != 1:
                #print('not')
                continue
            else:
                #print('mutant')
                with_mutant.append(i.upper())

        resid_1 = resid_1.upper()
        resid_2 = resid_2.upper()
        for chain in chains:
            try:
                dist = cmd.get_distance(atom1=f'chain {chain} and i. {resid_1} and n. {atom_1}',
                                            atom2=f'chain {chain} and i. {resid_2} and n. {atom_2}')
                dist_list[i].append(dist)
            except CmdException:
                missing_residue.append(i.upper())

    st.header('**Incorrectly numbered pdbs**')

    st.write(f'There are {len(error_pdbs)} structures with 1126 not being CYS or 57 not being a PRO in at least one chain:', str(error_pdbs),'. Removed those from analysis.')
    f = open("IncorrectNumberingPDB.txt", "w")
    f.write(str(error_pdbs))
    f.close()

    if mutation_name != '':
        st.write(f'There are {len(with_mutant)} structures with {mutation_id}{mutation_name.upper()}:', str(with_mutant))

    st.write(f'There are {len(missing_residue)} chains (in {len(set(missing_residue))} structures) with missing at least one of the requested atoms:', str(missing_residue))
    error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + str(mutation_name) + str(mutation_id) + '.txt'
    f = open(error_file_name, "w")
    f.write(str(missing_residue))
    f.close()
    return dist_list


def analysis(distancesDict, resid_1, atom_1, resid_2, atom_2, flag, mutation_name = '', mutation_id= ''):
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
    #plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
    plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')

    # sns.histplot(data=distances_only , binwidth=0.2)
    plt.xlabel('distance, A', fontsize=20)
    plt.ylabel("Chains count", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title("Distance distribution histogram")
    st.pyplot(fig)

    name = str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + str(mutation_name) + str(mutation_id)
    plt_name = './plots/distance_' + name + '_' + flag + '.png'
    plt.savefig(plt_name, bbox_inches='tight')


    fig2 = plt.figure(figsize=(15, 7.5))
    plt.hist(np.hstack(distances_only), bins=100, histtype='step', cumulative=True, label='Cumulative', density=True)
    plt.xlabel('distance, A', fontsize=20)
    plt.ylabel("Accumulated fraction", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title("Distance distribution cumulative histogram")
    st.pyplot(fig2)
    plt_name_cum = './plots/distance_' + name + '_' + flag + 'cumul' + '.png'
    plt.savefig(plt_name_cum, bbox_inches='tight')


    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in distancesDict.items()])).transpose()
    df_name = './distances/distance_' + name + '_' + flag +'.csv'
    df.to_csv(df_name)
    st.write(df)
    st.markdown(filedownload(df, df_name), unsafe_allow_html=True)




    #print(f'number of corrected numbered and analyzed structures {len(distancesDict)}')
    #return distances_only

def delete_PDB_folder():
        shutil.rmtree('./PDB')

#pdb_ids = get_spike_ids()
##dist = distance(pdb_ids, 318, 'PHE', 292,  'ALA')
#dist = distance(pdb_ids, 1126, 'CYS', 975,  'SER')
#analysis(dist)

