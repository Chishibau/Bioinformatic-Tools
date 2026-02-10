# Note that pysam will not install on windows 11 (in my experience at least)
# Run this while connected to a kernel on the server, pysam installs fine on linux

import pandas as pd
import numpy as np
import os 
import pysam

#%% Functions 

def read_has_deletion(read, mut_allele, ref_allele, pos):
    # Get list of tuples in format of (0, 105234781, 'T') which is (pos_in_read, genomic_pos, base)
    # Deletions are marked by None as pos_in_read
    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    
    # Get the number of base pairs deleted 
    del_length = len(ref_allele) - len(mut_allele)
    
    # Check that the anchor point is expected base (case insensitive)
    result = next(([first, third] for first, second, third in aligned_pairs if second == pos), None)
    if mut_allele.casefold() != result[1].casefold() or result[0] == None:
        return False 
    
    # Return False if read doesn't span pos to end_pos + 1
    end_pos = pos + del_length
    read_positions = [t[1] for t in aligned_pairs]
    if not (pos in read_positions and end_pos+1 in read_positions):
        return False
    
    # Grab the tuples of aligned_pairs that should fall in the deletion
    deleted = [t for t in aligned_pairs if t[1] != None and pos + 1 <= t[1] <= end_pos]
    if len(deleted) != del_length:
        return False
    
    # And grab the tuple directly after the deletion
    deleted_plus_one = [t for t in aligned_pairs if t[1] == end_pos + 1]
    
    # Check that all the pos_in_read in the region that's supposed to be deleted are None
    if (not all(t[0] is None for t in deleted)) or deleted_plus_one[0][0] == None:
        return False
    else:
        return True


def read_has_insertion(read, mut_allele, ref_allele, pos):
    # Get list of tuples in format of (0, 105234781, 'T') which is (pos_in_read, genomic_pos, base)
    aligned_pairs = read.get_aligned_pairs(with_seq=True)

    # Get the number of base pairs inserted
    ins_length = len(mut_allele) - len(ref_allele)
    
    # Check anchor base is what we expect
    result = next(([first, third] for first, second, third in aligned_pairs if second == pos), None)
    if ref_allele.casefold() != result[1].casefold():
        return False 
    
    # Check that we have the right number of inserted bases after the anchor
    for i, t in enumerate(aligned_pairs):
        if t[1] == pos:
            insertion = aligned_pairs[i+1:i+1+ins_length]
            if len(insertion) == 0:
                return False 
        
            if len(aligned_pairs) > i+1+ins_length and aligned_pairs[i+1+ins_length][1] == None:
                return False
            
    insert_read_pos = [t[0] for t in insertion]    
        
    # Check that all read bases in the insertion length are part of the insertion
    if not (all(t[1] is None for t in insertion) and all(t[1] is None for t in insertion)):
        return False
    elif read.query_sequence[insert_read_pos[0]:insert_read_pos[-1]+1].casefold() != mut_allele[1:].casefold():
        return False
    else:
        return True
    
def read_has_snv(read, mut_allele, pos):
    # Get the genomic regions that this read spans
    read_pos = read.get_reference_positions(full_length=True)
    query_seq = read.query_sequence
    
    # get the base at the actual position of the mutation
    try:
        query_idx = read_pos.index(pos)
        base_at_pos = query_seq[query_idx]
    except ValueError:
        print("Position not covered by this read, skip")
        return False
    
    # Check if this read has the mutant allele at the correct genomic position
    if base_at_pos == mut_allele:
        return True
    else:
        return False
    
# Extract MAPQ values from reads supporting the mutation
def get_mapq_values(bam_file, chrom, pos, mut_allele, ref_allele):
    
    print(f"Processing {bam_file}")
    bam = pysam.AlignmentFile(bam_file, "rc")
    
    # Convert to 0-based for pysam
    pos = pos - 1
            
    mapq_values = []

    # Fetch reads overlapping the position
    for read in bam.fetch(chrom, pos, pos + 1):
        
        try:
            # Skip unmapped, secondary 
            if read.is_unmapped or read.is_secondary:
                continue
            
            # Call different functions to determine if current read contains the mutation
            if len(mut_allele) == len(ref_allele) == 1:
                is_mut_read = read_has_snv(read, mut_allele, pos)
            elif len(mut_allele) < len(ref_allele):
                is_mut_read = read_has_deletion(read, mut_allele, ref_allele, pos)
            elif len(mut_allele) > len(ref_allele):
                is_mut_read = read_has_insertion(read, mut_allele, ref_allele, pos)
    
            # Collect MAPQ if this read supports the mutation          
            if is_mut_read:
                mapq_values.append(read.mapping_quality)
                
        except Exception as e:
            print(f"Error processing read: {e}")
    
    bam.close()
    
    return mapq_values

# Calculate MAPQ statistics
def get_mapq_stats(mapq_list):
    if len(mapq_list) == 0:
        return None, None, None
    
    mean_mapq = np.mean(mapq_list)
    median_mapq = np.median(mapq_list)
    min_mapq = np.min(mapq_list)
    
    return mean_mapq, median_mapq, min_mapq

#%% Load data and process

muts = pd.read_csv("/groups/wyattgrp/users/zshong/projects/bladder_wgs/mutation_curation/nc_uncurated.tsv", sep="\t")

muts["bam"] = "/groups/wyattgrp/projects/uc_metacohort/alignments/" + muts["Sample"] + ".bam"

# Get MAPQ values for each mutation
print("Extracting MAPQ values...")
muts["mapq_list"] = muts.apply(
    lambda x: get_mapq_values(x["bam"], x["Chrom"], x["Position"], x["Mut allele"], x["Ref allele"]), 
    axis=1
)

# Calculate MAPQ statistics
print("Calculating MAPQ statistics...")
muts[["mean_mapq", "median_mapq", "min_mapq"]] = muts["mapq_list"].apply(
    lambda x: pd.Series(get_mapq_stats(x))
)

# Optionally keep the full list of MAPQ values or drop it
# muts = muts.drop('mapq_list', axis=1)  # Uncomment to remove the list column

# Save output
output_path = "/groups/wyattgrp/users/zshong/projects/bladder_wgs/mutation_curation/nc_uncurated_mapq.tsv"
muts.to_csv(output_path, sep="\t", index=False)
print(f"Saved to {output_path}")