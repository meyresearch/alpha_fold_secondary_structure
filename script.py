# Imports
import pandas as pd
import numpy as np
import glob 
import os
import urllib
from collections import Counter

import json

def download_uniprot_json_file(uni_prot_id, workdir = '.'):
        #check if there is uniprot information available for the protein
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + uni_prot_id + '.json'
            html_2 = urllib.request.urlopen(url_2)
            print(html_2)
    
        except Exception as e:
            raise Exception('Failed to obtain UNIPROT data. %s'%e)

def download_alpha_fold_pdbs(uniprot_id_list, workdir = '.'):
    counter = 0
    for uniprot_id in uniprot_id_list:
        print(f"At entry {counter}/{len(uniprot_id_list)}")
        print(f"ID: {uniprot_id}")
        download = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb"
        
        try:
            file_name = os.path.join(workdir,f'AF-{uniprot_id}-F1-model_v2.pdb')
            #file_name = f"../data/alphafold/pdb_files/AF-{uniprot_id}-F1-model_v2.pdb"
            urllib.request.urlretrieve(download, file_name)
        except urllib.error.HTTPError:
            print("No such file.")
        counter = counter+1


def get_sec_struct(fname, dssp_path):
    
    """ returns secondary structure of protein file
    Parameters:
    -----------
    fname : String
        name of alpha fold pdb file
        
    dssp_path : String
        path to the installation of dssp
        
    Returns:
    --------
    secstruct : np.array
        array containing secondary structure information
    
    Cheat sheet of secondary structure information:
    -----------------------------------------------
    H = α-helix
    B = residue in isolated β-bridge
    E = extended strand, participates in β ladder
    G = 3-helix (310 helix)
    I = 5 helix (π-helix)
    T = hydrogen bonded turn
    S = bend
    """
    
    # call DSSP
    try:
        import subprocess
        subprocess.check_call("%s %s -o result.dssp"%(dssp_path, fname), shell=True)
        fin=open("result.dssp","r")
    except Exception as e:
        raise Exception("Could not calculate secondary structure! %s"%e)

    # parse output
    readit=False
    secstruct=[]
    for line in fin:

        if readit:
            try:
                if line[13:15] == '!*' or line[13] == '!':
                    continue
                else:
                    ss = line[16]
                    if line[16] == " ":
                        ss = "-"

                    secstruct.append(ss)
            except:
                continue

        if "#" in line:
            readit=True

    fin.close()

    # clean temporary files
    os.remove("result.dssp")

    return np.array(secstruct)