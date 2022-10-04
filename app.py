

from distance import *

# Logo image
#image = Image.open('logo.png')

#st.image(image, use_column_width=True)

# Page title


def main():
    st.markdown("""
    # Distance measurement tool for SARS-CoV-2 Spike PDB structures.
    
    This program measures a distance between two amino acids' atoms in all SARS-CoV-2 Spike structures deposited to Protein Data Bank.
    
    ---
    """)

    # Sidebar
    with st.sidebar.header('Enter residues between which you measure distance.'):
        resid_1 = st.sidebar.text_input("Input residue 1 id. For example 987.")
        atom_1 = st.sidebar.text_input("Input residue 1 atom name. For example CA")
        resid_2 = st.sidebar.text_input("Input residue 2 id. For example 413.")
        atom_2 = st.sidebar.text_input("Input residue 2 atom name. For example N")

        st.sidebar.markdown("""
    """)

    with st.sidebar.header(''):
        option = st.selectbox(
            'Are residues in the same chain?',
            ('same', 'different'))
        st.sidebar.markdown("""
    """)

    mutation = st.sidebar.selectbox('Are you looking for structures with a specific mutation?',
            ('no', 'yes'))

    with st.sidebar.header(''):
        mutation_name = st.sidebar.text_input("If yes, enter a mutation aminoacid name (for example GLY). If no, leave blank.")
        mutation_id = st.sidebar.text_input("If yes, enter a mutation aminoacid id (for example 614). If no, leave blank.")
        st.sidebar.markdown("""
    """)

    if st.sidebar.button('Measure'):
        pdb_ids = get_spike_ids()
        st.header('**Available pdb structures**')
        st.write(pdb_ids)

        #with st.spinner("Loading and checking PDBs"):
        #    pdb_ids_updated = pdb_load_check(pdb_ids)

        with st.spinner("Measuring distance"):
            if option == 'different':
                #dist, dist_reverse = distance_dif(pdb_ids_updated, resid_1, resid_2, atom_1, atom_2)
                dist, dist_reverse = distance_dif(pdb_ids, resid_1, resid_2, atom_1, atom_2, mutation_name, mutation_id)

                st.header('**Histogram of distances**')
                analysis(dist, resid_1, atom_1, resid_2, atom_2, 'diff1', mutation_name, mutation_id)
                analysis(dist_reverse, resid_1, atom_1, resid_2, atom_2, 'diff2', mutation_name, mutation_id)

            elif option == 'same':
                dist = distance_same(pdb_ids, resid_1, resid_2, atom_1, atom_2, mutation_name, mutation_id)
                st.header('**Histogram of distances**')
                analysis(dist, resid_1, atom_1, resid_2, atom_2, 'same', mutation_name, mutation_id)

        #uncomment if want to delete all PDB files
        #delete_PDB_folder()

        # Read in calculated descriptors and display the dataframe
        #st.header('**Histogram of distances**')
        #analysis(dist, resid_1, atom_1, resid_2, atom_2)
       # st.write(desc)
        #st.write(desc.shape)

    else:
        st.info('Enter residues in the sidebar to start!')

if __name__ == "__main__":
    main()