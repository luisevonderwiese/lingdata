from Bio import AlignIO
import os

chars = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G",
         "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
         "Y", "Z", "!", "\"", "#", "$", "%", "&", "\'", "(", ")", "*", "+", ",", "/", ":",
         ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]

char_dict = {}
char_dict["-"] = -1
char_dict["?"] = -1
for (i, char) in enumerate(chars):
    char_dict[char] = i

def get_interval_string(sites, joiner = ","):
    if len(sites) == 0:
        return ""
    intervals = []
    interval_begin = sites[0]
    cursor = sites[0]
    for site in sites[1:]:
        if site != cursor + 1:
            if cursor ==  interval_begin:
                intervals.append(str(interval_begin))
            else:
                intervals.append(str(interval_begin) + "-" + str(cursor))
            interval_begin = site
        cursor = site
    if cursor ==  interval_begin:
        intervals.append(str(interval_begin))
    else:
        intervals.append(str(interval_begin) + "-" + str(cursor))
    return joiner.join(intervals)



def create_ng_partition(in_path, out_path, multi_model = "GTR", complete = False):
    try:
        align = AlignIO.read(in_path, "phylip-relaxed")
    except:
        return
    num_sites = align.get_alignment_length()
    columns = [align[:, i] for i in range(num_sites)]
    columns = [[char_dict[value] for value in column] for column in columns]
    x_values = [max(col) for col in columns]
    partitions = []
    max_x = 0
    for (site, x) in enumerate(x_values):
        max_x = max(x + 1, max_x)
        partition_idx = max(0, x - 1)
        if not complete:
            partition_idx = min(partition_idx, 1)
        while(len(partitions) < partition_idx+1):
            partitions.append([])
        partitions[partition_idx].append(site + 1)
    partition_strings = {}
    interval_string = get_interval_string(partitions[0], ",")
    if interval_string != "":
        partition_strings["BIN"] = interval_string
    if complete:
        for i in range(1, len(partitions)):
            x = i + 2
            model = "MULTI" + str(x) + "_" + multi_model
            interval_string = get_interval_string(partitions[i], ",")
            if interval_string != "":
                partition_strings[model] = interval_string
    elif len(partitions) > 1 and len(partitions[1]) > 0:
        model = "MULTI" + str(max_x) + "_" + multi_model
        interval_string = get_interval_string(partitions[1], ",")
        if interval_string != "":
            partition_strings[model] = interval_string
    with open(out_path, "w+") as part_file:
        i = 1
        for (model, intervals) in partition_strings.items():
            part_file.write(model + ",  p" + str(i) + "=" + intervals + "\n")
            i += 1
