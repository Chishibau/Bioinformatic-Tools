import subprocess
import pandas as pd
from collections import defaultdict

# path to remote dir
REMOTE_PATH = "maple:"
save_path = "/groups/wyattgrp/users/zshong/projects/tmp_outputs/rclone_parse.tsv"

# get all directories 
cmd_dirs = f"rclone lsf {REMOTE_PATH} --dirs-only --recursive"
dir_output = subprocess.check_output(cmd_dirs, shell=True, text=True)

directories = []
for line in dir_output.splitlines():
    parts = [p for p in line.split("/") if p]
    if parts:
        directories.append(parts)

# get all files
cmd_files = f"rclone lsf {REMOTE_PATH} --files-only --recursive"
file_output = subprocess.check_output(cmd_files, shell=True, text=True)

direct_counts = defaultdict(int)
recursive_counts = defaultdict(int)

for filepath in file_output.splitlines():
    parts = [p for p in filepath.split("/") if p]
    if not parts:
        continue
    
    # direct parent
    parent = "/".join(parts[:-1])
    direct_counts[parent] += 1

    # recursive count for all ancestor directories
    for i in range(1, len(parts)):
        ancestor = "/".join(parts[:i])
        recursive_counts[ancestor] += 1

# make df
max_depth = max(len(d) for d in directories)

# pad levels
for d in directories:
    d.extend([""] * (max_depth - len(d)))

columns = [f"Level_{i+1}" for i in range(max_depth)]
df = pd.DataFrame(directories, columns=columns)

# full path column for dir
df["Full_Path"] = df.apply(
    lambda x: "/".join([v for v in x[columns] if v]), axis=1
)

# add file counts
df["File_Count"] = df["Full_Path"].map(direct_counts).fillna(0).astype(int)
df["Recursive_File_Count"] = df["Full_Path"].map(recursive_counts).fillna(0).astype(int)

# # Sort nicely
df = df.sort_values(columns)

# Save
df.to_csv(save_path, sep="\t", index=False)

print(f"Saved: {save_path}")
