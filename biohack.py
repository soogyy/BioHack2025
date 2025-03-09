import numpy as np
import pandas as pd
import re 

import streamlit as st

# load and read datasets
verified_data = pd.read_csv('verified_data.csv')
ctft_data = pd.read_csv('counterfeit_data.csv')
meta_data = pd.read_csv('metadata.csv', encoding='cp1252')

# function to extract atom counts from molecular formula
def extract_atoms(formula):
    '''Extract atom counts from molecular formula, returns dictionary of each atom type'''

    # skip over entry if entry is not string
    if not isinstance(formula, str):
        return {}
    # create dictionary for each atom type
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atom_counts = {}
    for atom, count in atoms:
        atom_counts[atom] = atom_counts.get(atom, 0) + (int(count) if count else 1)

    return atom_counts

# function to assign risk level to compound 
def assign_risk(row):
        if row['active_ingredient_match'] and row['formula_match']: # active ingredients and molecular formula matches
            if row['molecular_weight_diff'] <= 0.02 * row['molecular_weight_(g/mol)']:
                # molecular weight ±2%
                return "Low Risk"
            elif row['molecular_weight_diff'] <= 0.50 * row['molecular_weight_(g/mol)']:
                # molecular weight <= ± 50% and >±2%
                return "Medium Risk"
            else:
                # molecular weight >=±50%
                return "High Risk"
        # active ingredients don't match and/or molecular formula doesn't match
        elif not row['active_ingredient_match'] or not row['formula_match']:
            return "High Risk"
        # if compound is not within database
        return 'Unknown'

# function to compare drug composition with verified dataset 
def compare_drug(input_data, known_drugs):
    '''Compares input drug composition with known drugs and applies rule-based scoring, returns a table'''

    # standardize column names
    verified_data.columns = verified_data.columns.str.strip().str.lower().str.replace(" ", "_")

    # filter for only the matching drug
    known_drugs = verified_data[verified_data['drug_name'].str.lower() == input_data['drug_name'].lower()].copy()

    if known_drugs.empty:
        print(f"No verified drug data found for: {input_data['drug_name']}")

    # extracting atom counts for comparison
    input_formula_atoms = extract_atoms(input_data['molecular_formula'])
    known_drugs['input_formula_atoms'] = known_drugs['molecular_formula'].apply(extract_atoms)

    # check if active ingredient matches
    known_drugs['active_ingredient_match'] = known_drugs['active_ingredient'] == input_data['active_ingredient']

    # check molecular weight range (±2%, ±3-50%, >50%)
    molecular_weight_diff = (abs(known_drugs['molecular_weight_(g/mol)'] - input_data['molecular_weight'])) / known_drugs['molecular_weight_(g/mol)']
    known_drugs['molecular_weight_diff'] = molecular_weight_diff

    # check molecular formula match based on atom counts
    known_drugs['formula_match'] = known_drugs['input_formula_atoms'].apply(lambda x: x == input_formula_atoms)

    # assign risk level    
    known_drugs['risk'] = known_drugs.apply(assign_risk, axis=1)

    # get results
    result = known_drugs[['drug_name', 'risk', 'molecular_weight_diff', 'active_ingredient_match', 'formula_match']]
    
    # rename columns
    result = result.rename(columns={
        'drug_name': 'Drug Name',
        'risk': 'Risk Level',
        'molecular_formula': 'Molecular Formula',
        'molecular_weight_diff': 'Molecular Weight Difference (%)',
        'active_ingredient_match': 'Active Ingredient Match',
        'formula_match': 'Molecular Formula Match'
    })

    for index, row in result.iterrows():
        risk_level = row['Risk Level']

        if risk_level == 'Low Risk':
            st.write(f"Drug: {row['Drug Name']} - Risk Level: Low. The drug is most likely genuine.")
        elif risk_level == 'Medium Risk':
            st.write(f"Drug: {row['Drug Name']} - Risk Level: Medium. There is some variation, and the drug is likely a counterfeit. Please consult a trusted professional for further verification before taking.")
        else:
            st.write(f"Drug: {row['Drug Name']} - Risk Level: High, The drug is a counterfeit. Whether it's mislabelling or improper chemical composition, it can be life threatening to ingest counterfeit drugs. Do NOT take.")
    return result

def compare_counterfeit_dataset(ctft_df, verified_df):
    comparison_results = []

    verified_df['drug_name'] = verified_df['drug_name'].fillna('').astype(str)
    ctft_data['Drug Name'] = ctft_data['Drug Name'].fillna('').astype(str)

    for index, row in ctft_df.iterrows():
        drug_name = row['Drug Name']
        verified_row = verified_df[verified_df['drug_name'].str.lower() == drug_name.lower()]

        if verified_row.empty:
            # if no match found in verified dataset
            comparison_results.append({
                'Drug Name': drug_name,
                'Risk Level': 'Unknown - Drug not found',
                'Molecular Formula Match': 'N/A',
                'Molecular Weight Difference (%)': 'N/A',
                'Active Ingredient Match': 'N/A'
            })
            continue

        verified_row = verified_row.iloc[0]

        # check each component
        formula_match = row['Molecular Formula'] == verified_row['molecular_formula']
        weight_diff = abs(float(row['Molecular Weight (g/mol)']) - float(verified_row['molecular_weight_(g/mol)'])) / float(verified_row['molecular_weight_(g/mol)']) * 100
        active_match = set(row['Active Ingredient'].split(', ')) == set(verified_row['active_ingredient'].split(', '))

        # determine risk
        if active_match and formula_match:
            if weight_diff <= 0.2:
                risk = 'Low'
            elif weight_diff <= 5:
                risk = 'Medium'
            else:
                risk = "High"
        else:
            risk = 'High'

        # Append result
        comparison_results.append({
            'Drug Name': drug_name,
            'Risk Level': risk,
            'Molecular Formula Match': formula_match,
            'Molecular Weight Difference (%)': round(weight_diff, 2),
            'Active Ingredient Match': active_match
        })

    return pd.DataFrame(comparison_results)

# Streamlit UI
st.title("Counterfeit Medicine Finder: Drug Composition Comparison Tool")

### MAIN CODE:
# Input fields for drug composition
drug_name = st.text_input('Please enter the name of the medication you\'d like to validate: (Capitalize only the first letter)')
active_ingredient = st.text_input('What is the listed active ingredient on this medication? (Keep all letters lower case)')
molecular_formula = st.text_input('What is the listed molecular formula for the chemical composition of the active ingredient? Please capitalize the first letter of each element symbol. ')
molecular_weight = st.number_input('What is the listed molecular weight of this active ingredient in grams/mole? ')

input_drug = {
    'drug_name': drug_name,
    'active_ingredient': active_ingredient,
    'molecular_formula': molecular_formula,
    'molecular_weight': molecular_weight
}

# compare drugs
comparison_result = compare_drug(input_drug, verified_data)

# display results in a table format 
st.write("Comparison Results:", comparison_result)

# fetch metadata
if not comparison_result.empty:
    drug_name = comparison_result.iloc[0]['Drug Name']
    
drug_metadata = meta_data[meta_data['Drug Name'] == drug_name]

metadata_text = 'No metadata found for this drug'

if not drug_metadata.empty:
    metadata_row = drug_metadata.iloc[0]

    metadata_text = (
        f"**Dosage Forms:** {metadata_row.get('Dosage Forms', 'N/A')}  \n"
        f"**Physical Appearance:** {metadata_row.get('Physical Appearance', 'N/A')}  \n"
        f"**Color:** {metadata_row.get('Color', 'N/A')}  \n"
        f"**Manufacturer:** {metadata_row.get('Manufacturer', 'N/A')}  \n"
        f"**Packaging:** {metadata_row.get('Packaging', 'N/A')}  \n"
        f"**Batch Info:** {metadata_row.get('Batch Numbers', 'N/A')}"
    )

# write metadata
with st.expander(f'Additional information about {input_drug["drug_name"]}'):
    st.markdown(metadata_text)

# report suspicious drugs
with st.expander(f'Find a drug you think is counterfeit?'):
    suspicious_drug = st.text_input("Please tell us about your medication, including suspcious details: ")

    if suspicious_drug:
        st.write(f"Thank you for reporting: {suspicious_drug}. This will be reported to a regulatory agency for further action.")

# manual input of newly approved drugs
with st.expander(f'Are you a healthcare professional? Enter details about a new drug to our database to help prevent counterfeits.'):
    new_drug_name = st.text_input("Enter the name of the newly approved drug: ")
    new_active_ingredient = st.text_input("Enter its active ingredient (all lowercase): ")
    new_molecular_formula = st.text_input("Enter the molecular formula of the active ingredient: ")
    new_molecular_weight = st.text_input("Enter the molecular weight of the active ingredient (g/mol): ")

    # append new drug
    if new_drug_name and new_molecular_weight and new_active_ingredient and new_molecular_formula:
        new_row = pd.DataFrame([{
            'Drug Name': new_drug_name,
            'Active Ingredient': new_active_ingredient,
            'Molecular Formula': new_molecular_formula,
            'Molecular Weight (g/mol)': new_molecular_weight
        }])

        verified_data = pd.concat([verified_data, new_row])

        st.write("New drug added successfully!")


# comparing counterfeit data set to verified data set
ctft_results = compare_counterfeit_dataset(ctft_data, verified_data)
st.write("Example output table for counterfeit data compared to verified data set")

st.write(ctft_results)
print(ctft_results)
