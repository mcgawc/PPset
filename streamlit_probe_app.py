from __future__ import print_function
import streamlit as st
import pandas as pd
import io
import xlsxwriter
import itertools
import openpyxl 
import base64
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from base64 import b64encode
import requests
import json
from urllib import request, parse
import http.client

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '[' : ']' , '/' :'/', ']' : '[' }
    reverse_comp = ''.join([complement_dict.get(base, base) for base in reversed(sequence)])
    return reverse_comp

def complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', "[" : "[", "/" :"/", "]":"]"}
    comp_sequence = ''.join([complement_dict[base] for base in sequence])
    return comp_sequence


def get_access_token(client_id, client_secret, idt_username, idt_password):
    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]

def get_data_from_IDT(seq, token):
    conn = http.client.HTTPSConnection("www.idtdna.com")

    payload = json.dumps({
        "Sequence": seq,
        "NaConc": 50,
        "MgConc": 3,
        "DNTPsConc": 0.8,
        "OligoConc": 0.25,
        "NucleotideType": "DNA"
    })

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/Analyze", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "MeltTemp" value
    melt_temp = response_data["MeltTemp"]
    return(melt_temp)   
    
def get_mismatch_from_IDT(seq, comp_seq, token):
    conn = http.client.HTTPSConnection("www.idtdna.com")

    payload = json.dumps({
  "Settings": {
    "Sequence": seq,
    "NaConc": 50,
        "MgConc": 3,
        "DNTPsConc": 0.8,
        "OligoConc": 0.25
  },
  "Sequence2": comp_seq
})

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/TmMisMatch", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "MeltTemp" value
    missmatch_tm = response_data["MeltTemp"]
    return(missmatch_tm)   
    
def get_hairpin_data_from_IDT(seq, token):
  
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = json.dumps({
        "Sequence": seq,
        "NaConc": 50,
  "FoldingTemp": 25,
  "MgConc": 3,
  "NucleotideType": "DNA"
    })

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/Hairpin", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "deltaG" value
    delta_G = str(response_data[0]["deltaG"])
    return(delta_G)   
    

def get_selfdimer_data_from_IDT(seq, token):
    url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/SelfDimer"
    headers = {
        'Accept': 'application/json',
        'Authorization': f'Bearer {token}'
    }

    payload = {
        'primary': seq
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code == 200:
        response_data = response.json()
        # Check if the response is a list and not empty
        if isinstance(response_data, list) and response_data:
            # Extract the DeltaG from the first analysis result
            delta_G = response_data[0].get("DeltaG", "DeltaG value not found")
            return delta_G
        else:
            return "No self-dimer analysis results found."
    else:
        return f"Request failed with status code {response.status_code}"


def clean_up_input(gblock):
    gblock = gblock.replace(" ", "")
    gblock = gblock.upper()
    acceptable_characters = ['a', 'c', 't', 'g', 'A', 'T', 'G', 'C', '[', ']', '/',]
    # Create a translation table to remove characters not in acceptable_characters
    translation_table = str.maketrans("", "", "".join(c for c in gblock if c not in acceptable_characters))
    # Apply the translation table to gblock
    gblock = gblock.translate(translation_table)
    return gblock
    
def get_variant_regions(gblock):
    seq_before_snp = gblock[gblock.index('[') - 10:gblock.index('[')]
    seq_after_snp = gblock[gblock.index(']') + 1:gblock.index(']') + 12]
    variant_1 = gblock[ gblock.index('[') + 1 :gblock.index('/')]
    variant_2 = gblock[gblock.index('/') + 1 :  gblock.index(']')]
    seq_1 = seq_before_snp[len(variant_1)-1:] + variant_1 + seq_after_snp
    seq_2 = seq_before_snp[len(variant_2)-1:] + variant_2 + seq_after_snp
    return {seq_1:variant_1, seq_2:variant_2 }  # Return a dict containing the two sequences

def calculate_tm(sequence):

    # Create a sequence object from the input sequence
    dna_seq = Seq(sequence)

    # Calculate the Tm with custom salt and DNA concentrations
    tm = MeltingTemp.Tm_Wallace(dna_seq)

    return tm

def get_valid_permutations():
    valid_permutations = []
    length_range = range(7, 13)

    def has_5_or_6_consecutive_1s(sequence):
        consecutive_count = 0
        for digit in sequence:
            if digit == '1':
                consecutive_count += 1
                if consecutive_count >= 5:
                    return True
            else:
                consecutive_count = 0
        return False

    for length in length_range:
        for perm in itertools.product([0, 1], repeat=length):
            perm_str = ''.join(map(str, perm))
            if sum(perm) <= 7 and '111' in perm_str and not has_5_or_6_consecutive_1s(perm_str):
                valid_permutations.append('0' + perm_str + '0')
    return valid_permutations

def generate_sub_sequences(sequence):
    base_list = [base for base in sequence]
    base_list[10] = "*" + base_list[10]
    base_list[11] = "*" + base_list[11]
    base_list[9] = "*" + base_list[9]
    sub_sequences = []
    
    for i in range(0, 5):
        sub_sequence = base_list[5+i:14+i]
        sub_sequences.append(sub_sequence)    
    for i in range(0, 6):
        sub_sequence = base_list[4+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 7):
        sub_sequence = base_list[3+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 8):
        sub_sequence = base_list[2+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 9):
        sub_sequence = base_list[1+i:14+i]
        sub_sequences.append(sub_sequence)
    for i in range(0, 10):
        sub_sequence = base_list[i:14+i]
        sub_sequences.append(sub_sequence)
    return sub_sequences

def generate_master_probe_list(sub_sequences, valid_permutations):
    master_probe_list = []

    for sub_sequence in sub_sequences:
        sub_sequence_LNA = []
        for perm in valid_permutations:
            if len(sub_sequence) == len(perm):
                modified_sequence = []
                for i in range(len(sub_sequence)):
                    if perm[i] == '0':
                        modified_sequence.append(sub_sequence[i])
                    if perm[i] == '1':
                        modified_sequence.append("+" + sub_sequence[i])
                master_probe_list.append(modified_sequence)

    filtered_master_probe_list = []
    for sub_sequence in master_probe_list:
        exclude_sub_sequence = False
        for base in sub_sequence:
            if base[0] == "*":
                exclude_sub_sequence = True
                break
        if not exclude_sub_sequence:
            filtered_master_probe_list.append(sub_sequence)

    return filtered_master_probe_list


def remove_3G_3C(probe_list):
    filtered_probe_list = []
    for sub_sequence in probe_list:
        exclude_sub_sequence = False
        for i in range(len(sub_sequence) - 2):
            if (
                sub_sequence[i:i + 3] == ['G', 'G', 'G'] or
                sub_sequence[i:i + 3] == ['C', 'C', 'C']
            ):
                exclude_sub_sequence = True
                break
        if not exclude_sub_sequence:
            filtered_probe_list.append(sub_sequence)
    return filtered_probe_list

def remove_5primeG(probe_list):
    filtered_probe_list = []
    for sub_sequence in probe_list:
        exclude_sub_sequence = False
        if sub_sequence[0][-1] == 'G':
            exclude_sub_sequence = True
        if not exclude_sub_sequence:
            filtered_probe_list.append(sub_sequence)
    return filtered_probe_list

def calculate_Tm_values(probe_list):
    tm_dict = {}
    for sub_sequence in probe_list:
        for i, base in enumerate(sub_sequence):
            if len(base) == 3:
                sub_sequence[i] = "+" + "*" + base[-1].lower()
               
          
        base_sequence = "".join(s for s in sub_sequence)
       
        base_sequence = "".join([char for char in base_sequence if char != "*"])
        G_count = base_sequence.count('G')
        C_count = base_sequence.count('C')
        LNA_seq = ""
        for i in range(len(base_sequence)):
            if base_sequence[i] == "+":
                LNA_seq += base_sequence[i + 1]
        base_sequence = "".join([char for char in base_sequence if char != "+"])
        
        tm = calculate_tm(base_sequence) + calculate_tm(LNA_seq) +14
        tm_dict[''.join(sub_sequence)] = tm  # Use ''.join to create a string key
    return tm_dict

    
def create_probe_parameter_dict(tm_dict):
    probe_para_dict = {}  # Create a new dictionary to store the modified values
    for probe in tm_dict:
        parameter_dict = {"Tm": tm_dict[probe]}  # Create a dictionary with the temperature parameter
        probe_para_dict[probe] = parameter_dict  # Add the dictionary to the new dictionary
    return probe_para_dict  # Return the new dictionary with values as dictionaries

def add_length_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['probe length'] = PROBE.count("A") + PROBE.count("T") + PROBE.count("C") + PROBE.count("G")
    return probe_para_dict
def add_GC_ratio_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['% GC content'] = int(((PROBE.count("G") + PROBE.count("C"))/(PROBE.count("A") + PROBE.count("T") + PROBE.count("C") + PROBE.count("G")))*100)
    return probe_para_dict
def add_LNA_count_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        probe_para_dict[probe]['LNA count'] = PROBE.count("+")
    return probe_para_dict

def add_snp_distance_parameter(probe_para_dict):
    for probe in probe_para_dict:
        PROBE = ''.join([char for char in probe if char != "*" and char != "+"])
        distances = [PROBE.find(nucleotide) for nucleotide in 'atgc']
        snp_dist = max(distances)
        probe_para_dict[probe]['snp position'] = snp_dist
    return probe_para_dict
def filter_Tm_probes(probe_para_dict, tm_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (tm_range[0] <= probe_para_dict[probe]["Tm"] <= tm_range[1])]
    
    for probe in probes_to_remove:
        del probe_para_dict[probe]
    
    return probe_para_dict  
def filter_GC_probes(probe_para_dict, GC_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (GC_range[0] <= probe_para_dict[probe]['% GC content'] <= GC_range[1])]
    for probe in probes_to_remove:
        del probe_para_dict[probe]
    return probe_para_dict  

def filter_snp_pos(probe_para_dict, pos_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (pos_range[0] <= probe_para_dict[probe]['snp position'] <= pos_range[1])] 
    for probe in probes_to_remove:
        del probe_para_dict[probe]
    return probe_para_dict  
def filter_length_probe(probe_para_dict, len_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (len_range[0] <= probe_para_dict[probe]['probe length'] <= len_range[1])] 
    for probe in probes_to_remove:
        del probe_para_dict[probe]
    return probe_para_dict  
def filter_LNA_count_probe(probe_para_dict, LNA_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (LNA_range[0] <= probe_para_dict[probe]['LNA count'] <= LNA_range[1])] 
    for probe in probes_to_remove:
        del probe_para_dict[probe]
    return probe_para_dict  
def refine_Tm_values(probe_para_dict, token):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        PROBE = ''.join([char for char in PROBE if char != "*"])
        probe_para_dict[probe]["Tm"] = get_data_from_IDT(PROBE, token)
    return probe_para_dict
    
def get_hairpin_values(probe_para_dict, token):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        PROBE = ''.join([char for char in PROBE if char != "*"])
        probe_para_dict[probe]["Hairpin DeltaG"] = get_hairpin_data_from_IDT(PROBE, token)
    return probe_para_dict
    
def get_selfdimer_values(probe_para_dict, token):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        PROBE = ''.join([char for char in PROBE if char != "*"])
        probe_para_dict[probe]["Self Dimer DeltaG"] = get_selfdimer_data_from_IDT(PROBE, token)
    return probe_para_dict
def get_mismatch_values(probe_para_dict, variant, token):
    for probe in probe_para_dict:
        PROBE = probe.upper()
        PROBE = ''.join([char for char in PROBE if char not in ['+', '*']])
        snp_pos = probe_para_dict[probe]['snp position']
        comp_seq = complement(PROBE)
        
        # Construct the mismatch sequence
        mismatch_seq = comp_seq[:snp_pos - 1] + complement(variant) + comp_seq[snp_pos:]
        
        # Use the get_mismatch_from_IDT function to fetch the mismatch value
        mismatch_value = get_mismatch_from_IDT(PROBE, mismatch_seq, token)  # Replace probe_seq with the correct value
        probe_para_dict[probe]["Tm miss w/o LNA"] = mismatch_value
    return probe_para_dict

def filter_aprox_Tm_probes(probe_para_dict, aprox_tm_range=(40, 50)):
    probes_to_remove = [probe for probe in probe_para_dict if not (aprox_tm_range[0] <= probe_para_dict[probe]["Tm"] <= aprox_tm_range[1])]
    for probe in probes_to_remove:
        del probe_para_dict[probe]
def display_probe_data_2(probe_dict):
    probe_data = []
    for probe, parameters in probe_dict.items():
        probe = probe.upper()
        probe = ''.join([char for char in probe if char != "*"])
        probe_info = {"Probe": probe}
        probe_info.update(parameters)
        probe_data.append(probe_info)

    # Convert the data to a Pandas DataFrame
    df = pd.DataFrame(probe_data)

    # Display the DataFrame with sorting capabilities
    sorted_df = st.dataframe(df)

    # Allow sorting of the DataFrame by column
    column_to_sort = st.selectbox("Sort by Column", df.columns)
    ascending = st.checkbox("Sort Ascending", True)
    if ascending:
        sorted_df = sorted_df.sort_values(by=column_to_sort, ascending=True)
    else:
        sorted_df = sorted_df.sort_values(by=column_to_sort, ascending=False)
    
    st.dataframe(sorted_df)      
def display_probe_data(probe_dict):
    probe_data = []
    for probe, parameters in probe_dict.items():
        probe = probe.upper()
        probe = ''.join([char for char in probe if char != "*"])
        probe_info = {"Probe": probe}
        probe_info.update(parameters)
        probe_data.append(probe_info)

    # Display as a table
    st.table(probe_data)

def export_probe_data_to_excel(probe_para_dict, name):
    cleaned_probe_dict = {key.upper().replace('*', ''): value for key, value in probe_para_dict.items()}
    df = pd.DataFrame.from_dict(cleaned_probe_dict, orient='index')
    excel_file = name + ".xlsx"
    df.to_excel(excel_file, index=True)
    return excel_file
    
def execute_subprocess_1(seq, valid_permutations, tm_range, GC_range, pos_range, len_range, LNA_range, aprox_tm_range, input_seq, token):
    sub_sequences = generate_sub_sequences(seq)
    master_probe_list = generate_master_probe_list(sub_sequences, valid_permutations)
    master_probe_list = remove_5primeG(remove_3G_3C(master_probe_list))
    tm_dict = calculate_Tm_values(master_probe_list) 
    probe_dict = add_LNA_count_parameter(add_snp_distance_parameter(add_GC_ratio_parameter(add_length_parameter(create_probe_parameter_dict(tm_dict)))))
    filtered_probes = filter_aprox_Tm_probes(probe_dict, (int(aprox_tm_range[0]), int(aprox_tm_range[1])))
    filtered_probes = filter_GC_probes(probe_dict, (int(GC_range[0]), int(GC_range[1])))
    filtered_probes = filter_snp_pos(probe_dict, (int(pos_range[0]), int(pos_range[1])))
    filtered_probes = filter_length_probe(probe_dict, (int(len_range[0]), int(len_range[1])))
    filtered_probes = filter_LNA_count_probe(probe_dict, (int(LNA_range[0]), int(LNA_range[1])))
    if len(filtered_probes) <200 :
        with st.spinner("Validating Tms of " + str(len(filtered_probes)) + " Possible Probes"):
            filtered_probes = refine_Tm_values(probe_dict, token)
            filtered_probes = filter_Tm_probes(probe_dict, (int(tm_range[0]), int(tm_range[1])))
            st.write( str(len(filtered_probes)) + "probes found!")
        with st.spinner("Fetching secondary data"):
            get_hairpin_values(probe_dict, token)
            get_mismatch_values(probe_dict, input_seq, token)
        return probe_dict
    else:
        st.warning("Search parameters too big")

def main():
    st.title("Probe generator!")

    # Sidebar for user input
    st.sidebar.header("Input Fields here:")
    input_gblock = st.sidebar.text_input("Enter the gblock seq from ELN")
    gblock = clean_up_input(input_gblock)
    

    tm_range = st.sidebar.slider("Tm", 60, 67, (63, 65), 1)
    GC_range = st.sidebar.slider("GC content (%)", 0, 100, (40, 60), 1)
    pos_range = st.sidebar.slider("SNP position from 5' end", 1, 14, (4, 9), 1)    
    len_range = st.sidebar.slider("probe length", 9, 14, (10, 10), 1)
    LNA_range = st.sidebar.slider("Number of LNA", 3, 7, (3, 6), 1)
    aprox_tm_range = (58, 68)
    rev_comp = st.sidebar.checkbox('reverse complement', value=False)
    allele_1_go = st.sidebar.checkbox('Generate probes for allele 1', value=True)
    allele_2_go = st.sidebar.checkbox('Generate probes for allele 2', value=True)
  
    go_button = st.sidebar.button("GO!")
    
    if not input_gblock or not go_button:
        st.warning("Please enter the gblock sequence and press GO!")
        return
   
    if rev_comp:
        gblock = reverse_complement(gblock) 
    else:
        gblock = gblock

    valid_permutations = get_valid_permutations()
    
    input_seq = get_variant_regions(gblock)
    seq_1 = list(input_seq.keys())[0]
    seq_2 = list(input_seq.keys())[1]
    
    if allele_1_go:
        probe_dict_seq1 = execute_subprocess_1(seq_1, valid_permutations, tm_range, GC_range, pos_range, len_range, LNA_range, aprox_tm_range, input_seq[seq_2], token)
        st.header("Probes for " + input_seq[seq_1] + " allele")
        display_probe_data(probe_dict_seq1)
    if allele_2_go:
        probe_dict_seq2 = execute_subprocess_1(seq_2, valid_permutations, tm_range, GC_range, pos_range, len_range, LNA_range, aprox_tm_range, input_seq[seq_1], token)
        st.header("Probes for " + input_seq[seq_2] + " allele")
        display_probe_data(probe_dict_seq2)
    else:
        st.warning("What is my purpose if not to make probes?")

if __name__ == "__main__":
    client_id = "swapnil.mittal"
    client_secret = "f669c31a-1817-49ea-b8d4-d666dd1fb8bf"
    idt_username = "Swappi.mittal"
    idt_password = "Swappi_IDT"
    token = get_access_token(client_id, client_secret, idt_username, idt_password)
    main()
