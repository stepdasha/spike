#remove three hashtags to store data in files


import biotite.database.rcsb as rcsb
import datetime
import collections
import os
import matplotlib.pyplot as plt
import numpy as np

import streamlit as st

import pandas as pd
import base64

import biotite.structure.io as strucio
import biotite.structure as struc
#from PIL import Image

# File download
def filedownload(df, df_name, phrase):
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download={df_name}>' + 'Download ' + phrase + ' </a>'
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
    pdb_ids.remove('7ZJL')
    #print(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
    #      f"resolution less than or equal to {max_resolution} with mass more than or equal to {min_weight}: {len(pdb_ids)}")
    #print("Selected PDB IDs:\n", *pdb_ids)
    st.write(f"Number of spike structures on  {today.year}-{today.month}-{today.day} with "
          f"resolution less than or equal to {max_resolution}A with mass more than or equal to {min_weight}kDa: {len(pdb_ids)}")
    return (pdb_ids)

#@st.cache(suppress_st_warning=True, show_spinner=False)
@st.experimental_memo(suppress_st_warning=True, show_spinner=False, persist='disk')
def pdb_files_loader(pdb_ids):
    if not os.path.exists('PDB'):
        os.mkdir('PDB')

    st.write(f"There are new PDB files in the database. PDB files will be updated. It could take several minutes.")
    len_pdbid = len(pdb_ids)
    my_bar = st.progress(0)

    proteins = {}
    for count, i in enumerate(pdb_ids):
        # print('object begin', cmd.get_object_list('(all)'))
        # cmd.delete("*")
        my_bar.progress((count + 1) / len_pdbid)

        i = i.lower()
        # download structure
        try:
            file = rcsb.fetch(i, "pdb", target_path="PDB/")
            # print('pdb fetched ')
        except:
            file = rcsb.fetch(i, "cif", target_path="PDB/")
            # print('cif fetched ')

        # laod strcutrues to pycharm
        proteins[str(i)] =  strucio.load_structure(file)
        #proteins.append(strucio.load_structure(file))
    return proteins



def distance_dif(proteins, pdb_ids, resid_1,  resid_2, atom_1, atom_2, mutation_name = '', mutation_id = '', flag = 'same'):
            #print(proteins)

            selection_flag = False # modify the flag if there  was at least one distance calculated

            mutation_selection_flag = True #modify flag if there are 0 structures with specified mutation

            dist_list = collections.defaultdict(list)  # empty dictionary for future rmsd
            dist_list_reverse = collections.defaultdict(list) # empty dictionary for future rmsd from reverse
            log_file = collections.defaultdict(list)
            missing_residue = []
            missing_residue_reverse = []
            error_pdbs = []
            with_mutant = []
            len_pdbid = len(pdb_ids)
            my_bar = st.progress(0)
            for count, i in enumerate(pdb_ids):
                # print('object begin', cmd.get_object_list('(all)'))
                # cmd.delete("*")
                my_bar.progress((count + 1) / len_pdbid)

                i = i.lower()

            # laod strcutrues to pycharm
                protein = proteins[str(i)]
                #print(np.unique(protein.chain_id))

                # sanity check that the numbering is correct, check if 1126 is CYS  in chains and 57 is PRO.
                # take those chains only
                chains = []
                for ch in np.unique(protein.chain_id):
                    sel_cys = protein[(protein.chain_id == ch) & (protein.res_name == 'CYS')  & (protein.res_id == 1126)][:]
                    sel_pro = protein[(protein.chain_id == ch) & (protein.res_name == 'PRO')  & (protein.res_id == 57)][:]
                    if len(sel_cys) == 0 or len(sel_pro) == 0:
                        #print('skip pdb', i)
                        continue
                    else:
                        chains.append(ch)

               # print('for pdb', i, 'there are', len(chains), 'chains')


            # only take those that had correct sequence numbering in all 3 chains
                if len(chains) !=3 :
                    error_pdbs.append(i.upper())
                    #print('skip', i)
                    continue
                #print(chains)
                # check if the sequence has a requested mutation
                if mutation_name != '':
                    sel_mut = protein[(protein.chain_id == chains[0]) & (protein.res_name == mutation_name.upper()) & (protein.res_id == int(mutation_id))][:]
                    if len(sel_mut) == 0:
                        continue
                    else:
                        with_mutant.append(i.upper())
                    #print('has mutation', i)



                #print(chains_ordered)
                atom_1 = atom_1.upper()
                atom_2 = atom_2.upper()
                resid_1 = int(resid_1)
                resid_2 = int(resid_2)


                #measure 2 distances to find chains orientation (clockwise/counterclockwise)
                chains_ordered = [chains[0]]
                atom_971 = protein[(protein.chain_id == chains[0]) & (protein.res_id == 971) & (protein.atom_name == 'CA')][:]
                atom_752_1 = protein[(protein.chain_id == chains[1]) & (protein.res_id == 752) & (protein.atom_name == 'CA')][:]
                atom_752_2 = protein[(protein.chain_id == chains[2]) & (protein.res_id == 752) & (protein.atom_name == 'CA')][:]

                dist1 = struc.distance(atom_971, atom_752_1)
                dist2 = struc.distance(atom_971, atom_752_2)
                #print('dist', dist1, dist2)

                if min(dist1, dist2) == dist1:
                    chains_ordered.append(chains[1])
                    chains_ordered.append(chains[2])
                    #print('clock')
                elif min(dist1, dist2) == dist2:
                    chains_ordered.append(chains[2])
                    chains_ordered.append(chains[1])
                    #print('counterclock')
                log_file[i.upper()] = chains_ordered

                if flag == 'different':
                # measure the distance of interest
                    for j in range(0,3):
                            requested_atom1 = protein[(protein.chain_id == chains_ordered[j]) & (protein.res_id == resid_1) & (protein.atom_name == atom_1)][:]
                            if j == 2 :
                                requested_atom2 = protein[(protein.chain_id == chains_ordered[0]) & (protein.res_id == resid_2) & (protein.atom_name == atom_2)][:]
                            else:
                                requested_atom2 = protein[(protein.chain_id == chains_ordered[j + 1]) & (protein.res_id == resid_2) & (protein.atom_name == atom_2)][:]
                            dist = struc.distance(requested_atom1, requested_atom2)
                            #print(i, dist, requested_atom2, requested_atom1)
                            if len(dist) == 0:
                                dist_list[i].append(None)
                                missing_residue.append(i.upper())
                            else:
                                #print(dist)
                                selection_flag = True
                                dist_list[i].append(float(dist))

                            #break
                # test for the second plot in other direction
                    for j in range(0,3):
                                requested_atom2 = protein[(protein.chain_id == chains_ordered[j]) & (protein.res_id == resid_2) & (protein.atom_name == atom_2)][:]
                                if j == 2 :
                                    requested_atom1 = protein[(protein.chain_id == chains_ordered[0]) & (protein.res_id == resid_1) & (protein.atom_name == atom_1)][:]
                                else:
                                    requested_atom1 = protein[(protein.chain_id == chains_ordered[j + 1]) & (protein.res_id == resid_1) & (protein.atom_name == atom_1)][:]
                                dist_reverse = struc.distance(requested_atom1, requested_atom2)

                                if len(dist_reverse) == 0:
                                    dist_list_reverse[i].append(None)
                                    missing_residue_reverse.append(i.upper())
                                else:
                                    dist_list_reverse[i].append(float(dist_reverse))



                elif flag == 'same':
                    for chain in chains_ordered:
                            requested_atom1 = protein[(protein.chain_id == chain) & (protein.res_id == resid_1) & (protein.atom_name == atom_1)][:]
                            requested_atom2 = protein[(protein.chain_id == chain) & (protein.res_id == resid_2) & (protein.atom_name == atom_2)][:]
                            dist = struc.distance(requested_atom1, requested_atom2)
                            if len(dist) == 0:
                                dist_list[i].append(None)
                                missing_residue.append(i.upper())
                            else:
                                selection_flag = True
                                dist_list[i].append(float(dist))

            st.header('**Incorrectly numbered pdbs**')

            st.write(f'There are {len(error_pdbs)} structures with 1126 not being CYS or 57 not being a PRO in at least one chain:', str(error_pdbs),'. Removed those from analysis.')
            ##f = open("IncorrectNumberingPDB.txt", "w")
            ##f.write(str(error_pdbs))
            ##f.close()

            if mutation_name != '':
                if len(with_mutant) == 0:
                    mutation_selection_flag = False
                    dist_list['mutation'].append(None)
                    dist_list_reverse['mutation'].append(None)
                st.write(f'There are {len(with_mutant)} structures with {mutation_id}{mutation_name.upper()}:',str(with_mutant))

            if selection_flag == True:
                st.write(f'There are {len(missing_residue)} chains (in {len(set(missing_residue))} structures) with missing at least one of the requested atoms:', str(missing_residue))
            elif mutation_selection_flag == False:
                st.write('*** Please, check your mutation selection ***')
            else:
                st.write('*** Requested atoms are not found in any PDB files. Please, check your atom selection.***')
            ##error_file_name = './error_residue/errorsPDB_' + str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + str(mutation_name) + str(mutation_id) + '.txt'
            ##f = open(error_file_name, "w")
            ##f.write(str(missing_residue))
            ##f.write(str(missing_residue_reverse))
            ##f.close()

            #log_df = pd.DataFrame(dict([(k, v) for k, v in log_file.items()])).transpose()
            log_df = pd.DataFrame.from_dict(log_file, orient='index')
            st.markdown(filedownload(log_df, 'log_file.csv', 'chain names log'), unsafe_allow_html=True)

            return dist_list, dist_list_reverse


def analysis(distancesDict, resid_1, atom_1, resid_2, atom_2, flag, mutation_name = '', mutation_id= ''):
    """
    plot the histogram of distances
    """
    distances_only_none = np.hstack(list(distancesDict.values()))
    distances_only = [i for i in distances_only_none if i is not None]

    if len(distances_only) == 0:
        st.write("")
    else:
        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in distancesDict.items()])).transpose().dropna(how='all')
        if flag == 'diff1':
            df.columns = ['0-1', '1-2', '2-0']
        elif flag == 'diff2':
            df.columns = ['1-0', '2-1', '0-2']
        number_structure = len(df.axes[0])
        st.write(f'Number of structures with at least one analyzed chain {number_structure} (in {len(distances_only)} chains)')

        #if not os.path.exists('plots'):
        #    os.mkdir('plots')

        #if not os.path.exists('distances'):
        #    os.mkdir('distances')

        fig = plt.figure(figsize=(15, 7.5))
        #plt.hist(np.hstack(distances_only), bins=100, color="skyblue", edgecolor='white')
        plt.hist(distances_only, bins=100, color="skyblue", edgecolor='white')

        # sns.histplot(data=distances_only , binwidth=0.2)
        plt.xlabel('distance, A', fontsize=26)
        plt.ylabel("Chains count", fontsize=26)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        #plt.title("Distance distribution histogram")
        st.pyplot(fig)

        name = str(resid_1) + str(atom_1) + '_' + str(resid_2) + str(atom_2) + str(mutation_name) + str(mutation_id)
        #plt_name = './plots/distance_' + name + '_' + flag + '.png'
        #plt.savefig(plt_name, bbox_inches='tight')


        fig2 = plt.figure(figsize=(15, 7.5))
        plt.hist(np.hstack(distances_only), bins=100, histtype='step', cumulative=True, label='Cumulative', density=True)
        plt.xlabel('distance, A', fontsize=26)
        plt.ylabel("Cumulative frequency", fontsize=26)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        #plt.title("Distance distribution cumulative histogram")
        st.pyplot(fig2)
        #plt_name_cum = './plots/distance_' + name + '_' + flag + 'cumul' + '.png'
        #plt.savefig(plt_name_cum, bbox_inches='tight')



        df_name = './distances/distance_' + name + '_' + flag +'.csv'
        #df.to_csv(df_name)
        st.write(df)
        st.markdown(filedownload(df, df_name, 'distribution'), unsafe_allow_html=True)




    #print(f'number of corrected numbered and analyzed structures {len(distancesDict)}')
    #return distances_only

##def delete_PDB_folder():
##        shutil.rmtree('./PDB')

#pdb_ids = get_spike_ids()
##dist = distance(pdb_ids, 318, 'PHE', 292,  'ALA')
#dist = distance(pdb_ids, 1126, 'CYS', 975,  'SER')
#analysis(dist)

