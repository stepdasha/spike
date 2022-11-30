This program calculates the distribution of the requested distance between any two atoms across all SARS-CoV-2 spike structures available in the Protein Data Bank.
If you want to use it locally, you can pull the last commit from the github using the following command:

	git clone git@github.com:stepdasha/spike.git --depth 1
Create the environment called 'structure' from the file environment.yml file with the following command:

	conda env create --file environment.yml
Activate the environment with the following command:

	conda activate structure
Once the environment is activated you can run this program with the follwong streamlit command:
	
	streamlit run app.py 

