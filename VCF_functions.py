"""
VCF_functions.py
Jorge Holguin
October 15, 2020

Read Variant Call Files (VCF) processed with VEP (Variant Effect Predictor) and Mutation 
Annotation Format (MAF) files into a pandas dataframe. 
    Supports gzip files. 

For VCF files:
    Ignores the ##commented lines but retrives the keys from the VEP Consequence (CSQ) field
    Uses the CSQ keys to create a dictionary with the keys:value from VEP
    Filters the variants of interest.
    Extends the VCF file with the columns of VEP when a list of transcript IDs is provided.
    If no list is provided then the extension is not done

For MAF files:
    Identifies variants of interest in the `all_effects` column 
    If a list of Ensembl transcripts is supplied, the dataframe is extended the information

Usage example:

    >>> import VCF_functions
    >>> df_vcf = read_process_vcf(file, variant='missense_variant', variant_class='SNP', id_list = ['ENSTXXXXX', 'ENSTYYYY'], return_case_id=False)
    >>> df_maf = read_process_maf(file, variant='Missense_Mutation', consequence='missense_variant', variant_class='SNP', id_list = ['ENSTXXXXX', 'ENSTYYYY'])
"""

import pandas as pd
import numpy as np
import gzip

def rows2skip(file):
    
    """
    Determine the number of header rows to skip when reading 
    a VCF `file`
    
    Parameters:
        
        file (string): path to the file
        
    Returns:
    
        int: Number of header rows to skip
    """
    
    skip_rows = 0
    fn_open = gzip.open if file.endswith('.gz') else open
    with fn_open(file, mode = 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                skip_rows += 1
            else:
                break
    return skip_rows

def case_id(file):
    
    """
    Retrieve the identifier for the case corresponding to the vcf file and the sample ID of the tumor
    in the format case_id//tumor_sample_id
    
    Parameters: 
    
        file (string): path to the file
        
    Returns:
        
        string: identifier for the case corresponding to the vcf file. 
    """
    identifier = ''
    sample_id = ''
    fn_open = gzip.open if file.endswith('.gz') else open
    with fn_open(file, mode = 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                if line.startswith('##INDIVIDUAL'):
                    identifier = line.split(',')[0].replace('##INDIVIDUAL=<NAME=', '')
                    if len(identifier) > 0 and len(sample_id) > 0:
                        break # if we have the info we need then break the loop 
                elif line.startswith('##SAMPLE=<ID=TUMOR'):
                    sample_id = line.split(',')[1].replace('NAME=', '')
                    if len(identifier) > 0 and len(sample_id) > 0:
                        break # if we have the info we need then break the loop 
                else:
                    continue
            else:
                break # to prevent waisting time reading lines we dont use
                
    return (identifier, sample_id)

def csqkeys(file):
    
    """
    Retrieve the column headers for the INFO column on VEP vcf files from NCI GDC
    The columns headers are contained in the header lines of VCF files
    
    Parameters: 
    
        file (string): path to the file
        
    Returns:
        
        list: List of strings with the keys for the values in the INFO CSQ 
            column
    """
    
    csq_keys = ''
    fn_open = gzip.open if file.endswith('.gz') else open
    with fn_open(file, mode = 'rt') as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ'):
                csq_keys = line.split('Format: ', 1)[1].replace('">\n', '').split('|')
                break
            else:
                continue
                
    return csq_keys

def update_value(dict_a, key, csq_keys):
    
    """
    Updates the value of the 'CSQ' key in the INFO dictionary. Each variant and its corresponding information
    are converted to a dictionary where the keys are `csq_keys` and the values are the orginal values present
    for the 'CSQ' key. The value for the 'CSQ' key is updated to a list of dictionaries.
    
    Parameters:
    
        dict_a (dict): Dictionary made from the INFO column
        
        key (string): Name of the key to update
        
        csq_keys (list): Output of csqkeys()
        
    Returns:
    
        dict: Updated dict_a
    """
    
    dict_a.update({key: [dict(zip(csq_keys, i.split('|'))) for i in dict_a[key].split(',')]})
    
    return dict_a

def filter_variants(dict_a, variant, variant_class):
    
    """
    Keeps the specified type of `variant`
    
    Parameters:
    
        dict_a (dict): Dictionary made from the INFO column and updated using update_value()
        
        variant (string): Name of the variant to subset
        
    Returns:
    
        dict: Subset of dict_a containing only the variants specified
    """
    # Filter only the specified variants
    dict_a.update({'CSQ': [i for i in dict_a['CSQ'] if variant in i['Consequence'] and variant_class in i['VARIANT_CLASS']]})
    
    return dict_a

def filter_variants_maf(list_of_dicts, variant):
    
    """
    Keeps the specified type of `variant`. For use with MAF files.
    
    Parameters:
    
        list_of_dicts (list of dictionaries): List of dictionaries made from the `all_effects` column in MAF files
        
        variant (string): Name of the variant to subset
        
    Returns:
    
        list: List of dictionaries containing only the type of variants specified in `variant`
    """
    # Filter only the specified variants
    updated_list_of_dicts = [dic for dic in list_of_dicts if variant in dic['Consequence_all_effects']]
    
    if len(updated_list_of_dicts) > 0:
        return updated_list_of_dicts
    else:
        return np.nan

def find_variants(dict_a, id_list):
    
    """
    Finds the transcript ID inside the INFO dictionary for a list of specified transcripts
    and returns the CSQ dictionary corresponding to that variant
    
    Parameters:
    
        dict_a (dict): Dictionary made from the INFO column and processed with update_value()
            and filter_variants()
            
        id_list: List of identifiers corresponding to the proteins/genes of interest
        
    Returns:
        
        dict: Dictionary inside CSQ that has a value for the 'Feature' key that
            is contained in `id_list` 
    """
    
    variants = ''
    for i in dict_a['CSQ']:
        if i['Feature'] in id_list:
            variants = i
    
    if len(variants) > 0:
        return variants
    else:
        return np.nan
    

def find_variants_maf(list_of_dicts, id_list):
    
    """
    For use with MAF files
    
    Finds the transcript ID inside the `all_effects` dictionary for a list of specified transcripts
    and returns the dictionary with information corresponding to that variant
    
    Parameters:
    
        dict_a (dict): Dictionary made from the INFO column and processed with update_value()
            and filter_variants()
            
        id_list: List of identifiers corresponding to the proteins/genes of interest
        
    Returns:
        
        dict: Dictionary with the information corresponding to the trascript ids provided in `id_list` 
    """
    
    variants = ''
    for dic in list_of_dicts:
        if dic['Transcript_ID_all_effects'] in id_list:
            variants = dic
    
    if len(variants) > 0:
        return variants
    else:
        return np.nan
    

def read_process_vcf(file, variant, variant_class, id_list = [], return_case_id=False):
    
    """
    Reads and processes vcf files from NCI GDC. Supports .gzip files. 
    The comment lines are skipped, the keys for the INFO column
    are retrieved and the INFO column is converted into a dictionary with the retrieved keys.
    Variants that are not the specified `variant` are removed and the 
    dataframe is extended to include the information for the specified 
    list of ids in `id_list`
    
    Parameters:
    
        file (string): path to the file
        
        variant (string): name of the variants to subset
            Examples of `variant`: 'missense_variant','downstream_gene_variant', 'upstream_gene_variant'
            
        variant_class (string): name of the variant class to subset
            Examples of `variant_class`: 'SNP'
            
        id_list (list): List of Ensembl transcript identifiers to find in the file
        
        return_case_id (bool): Whether the case id is to be returned. This id helps later when finding the 
            biospecimen and clinial metadata associated with the given vcf file
        
    Returns:
    
        pandas dataframe: Dataframe with the information from a VEP VCF file from NCI GDC
    """
    
    skip_rows = rows2skip(file)
    csq_keys = csqkeys(file)
    
    df = pd.read_csv(file, sep = '\t', skiprows=skip_rows).rename({'#CHROM': 'CHROM'}, axis = 1)
    
    df = df.loc[(df['FILTER'] == 'PASS') & (df.INFO.str.contains(variant,case=False))]
    
    # Include '=nan' in cases where the field is empty. Example change 'IN_PON' to 'IN_PON=nan'
    df['INFO'] = df['INFO'].apply(lambda x: [i + '=nan' if '=' not in i else i for i in x.split(';')])
    
    # Create a dictionary splitting on '='
    df['INFO'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x))
    
    # Convert the value of 'CSQ' to a dictionary where the keys are extracted from the header of the file using the function csqkeys()
    df['INFO'] = df['INFO'].apply(lambda x: update_value(x, 'CSQ', csq_keys))
    
    # Filter the specified variants
    df['INFO'] = df['INFO'].apply(lambda x: filter_variants(x, variant, variant_class))
    
    # Given a specified list of identifiers in `id_list` find the mutations corresponding to those proteins
    if len(id_list) > 0:

        df['Variant_INFO'] = df['INFO'].apply(lambda x: find_variants(x, id_list))
        df = df.dropna(subset = ['Variant_INFO'])
        
        df_variant = pd.json_normalize(df['Variant_INFO'])
        df_extended = pd.concat([df.drop(['Variant_INFO'], axis = 1).reset_index(drop = True), df_variant], axis = 1)
        
        if return_case_id == True:
            
            identifier, sample_id = case_id(file)
            return (df_extended, identifier, sample_id)
        
        else:
            
            return df_extended
    
    else:
        
        if return_case_id == True:
            
            identifier, sample_id = case_id(file)        
            return (df, identifier, sample_id)
        
        else:
            
            return df
        
        
def read_process_maf(file, variant='Missense_Mutation', consequence='missense_variant', variant_class='SNP', id_list = []):
    
    """
    Reads and processes MAF files from NCI GDC. Uses the information in 
    the 'all_effects' column to find mutations of interest.
    Variants that are not the specified `variant` are removed and the 
    dataframe is extended to include the information for the specified 
    list of ids in `id_list`
    
    Parameters:
    
        file (string): path to the file
        
        variant (string): name of the variants to subset
            Examples of `variant`: 'missense_variant','downstream_gene_variant', 'upstream_gene_variant'
            
        consequence (string): name of the consequence of the mutation
            Examples of `consequence`: 'missense_variant'
            
        variant_class (string):
            Examples of `variant_class`: 'SNP'
            
        id_list (list): List of identifiers to find in the file
        
    Returns:
    
        pandas dataframe: Dataframe with the information from a MAF file from NCI GDC
    """
    
    # Open the MAF file
    df = pd.read_csv(file, sep = '\t', comment='#', low_memory=False)

    # Subset the missense variants and somatic cancer variants
#     df = df.loc[(df['Variant_Classification'] == variant) & (df['One_Consequence'] == consequence) 
#                         & (df['Mutation_Status'] == 'Somatic') & (df['Variant_Type'] == variant_class)]
    
    # from https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#:~:text=46%20%2D%20all_effects,Sift%2CPolyPhen%2CStrand%5D)
    # added 'all_effects' to the end of each key to avoid column names already present in the MAF files
    keys = ['Symbol_all_effects','Consequence_all_effects','HGVSp_Short_all_effects','Transcript_ID_all_effects',\
            'RefSeq_all_effects','HGVSc_all_effects','Impact_all_effects','Canonical_all_effects',\
            'Sift_all_effects','PolyPhen_all_effects','Strand_all_effects'] 
    
    # Create a list of dictionaries splitting on first on ';', then on ',' and using the keys above
    df['all_effects'] = df['all_effects'].apply(lambda x: [dict(zip(keys, i.split(','))) for i in x.split(';')])
    
    # Filter the specified variants
    df['all_effects'] = df['all_effects'].apply(lambda x: filter_variants_maf(x, consequence))
    
    # Drop rows without the specified variants
    df.dropna(axis=0, subset=['all_effects'], inplace=True)
    
    # Given a specified list of identifiers in `id_list` find the mutations corresponding to those proteins
    if len(id_list) > 0:

        df['Variant_INFO'] = df['all_effects'].apply(lambda x: find_variants_maf(x, id_list))
        df = df.dropna(subset = ['Variant_INFO'])
        
        df_variant = pd.json_normalize(df['Variant_INFO'])
        df_extended = pd.concat([df.drop(['Variant_INFO'], axis = 1).reset_index(drop = True), df_variant], axis = 1)
        
        return df_extended
    
    else:        
        
        return df