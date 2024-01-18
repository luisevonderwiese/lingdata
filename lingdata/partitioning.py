
import lingdata.pathbuilder as pb

class Partitioning:
    bi_split = 2 #last x in lower partition

    def __init__(self, num_possible_values):
        self.num_possible_values = num_possible_values
        self.char_idx_partitions_x = []
        for char_idx, num in enumerate(self.num_possible_values):
            # x in the sense of MULTIx_..., MULTI2 corresponds to bin, MULTI1 not relevant as not informative, we assign it to BIN
            x  = max(2, num)
            while(len(self.char_idx_partitions_x) < x + 1):
                self.char_idx_partitions_x.append([])
            self.char_idx_partitions_x[x].append(char_idx)

        #all sites with x <= self.bi_split are assigned to model with x=bi_split, all others to model with x=max_x
        assert(self.char_idx_partitions_x[0] == [])
        assert(self.char_idx_partitions_x[1] == [])
        self.char_idx_partitions_2 = [[] for _ in self.char_idx_partitions_x]
        max_x = len(self.char_idx_partitions_x) - 1
        assert(max_x >= 2)
        if self.bi_split >= max_x:
            for index_list in self.char_idx_partitions_x:
                self.char_idx_partitions_2[max_x] += index_list
            self.char_idx_partitions_2[max_x].sort()
        else:
            for index_list in self.char_idx_partitions_x[1:self.bi_split + 1]:
                self.char_idx_partitions_2[self.bi_split] += index_list
            self.char_idx_partitions_2[self.bi_split].sort()
            for index_list in self.char_idx_partitions_x[self.bi_split + 1:]:
                self.char_idx_partitions_2[max_x] += index_list
            self.char_idx_partitions_2[max_x].sort()


    #max(1,...) because also uninformatives are represented with one site
    # +1 because raxml-ng paritionings start with index 1
    def get_site_indices(self, char_idx, msa_type):
        if msa_type in ["multi", "ambig"]:
            return [char_idx + 1]
        elif msa_type in ["bin"]:
            lower = sum([max(el, 1) for el in self.num_possible_values[:char_idx]]) + 1
            upper = lower + max(1, self.num_possible_values[char_idx])
            return range(lower, upper)
        else:
            raise Exception("Paritioning not defined for msa_type ", msa_type)

    def char_to_site_indices(self, char_idx_partitions, msa_type):
        site_idx_partitions = []
        for char_indices in char_idx_partitions:
            site_indices = []
            for char_idx in char_indices:
                site_indices += self.get_site_indices(char_idx, msa_type)
            site_idx_partitions.append(site_indices)
        return site_idx_partitions


    def get_interval_string(self, sites, joiner = ","):
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

    def get_part_model(self, msa_type, model, gamma, x):
        if msa_type in ["bin"]:
            if model != "BIN":
                raise Exception("Paritioning on ", msa_type, " MSA only with model BIN not with ", model)
            part_model = "BIN"
        elif msa_type in ["multi", "ambig"]:
            if model not in ["MK", "GTR"]:
                raise Exception("Paritioning on ", msa_type, " MSA only with model MK or GTR not with ", model)
            if x == 2:
                part_model = "BIN"
            else:
                part_model = "MULTI" + str(x) + "_" + model
            if msa_type == "ambig":
                part_model += "+M{" + pb.charmap_path(x) + "}"
        else:
            raise Exception("Paritioning not defined for msa_type ", msa_type)
        if gamma:
            part_model += "+G"
        return part_model


    def write(self, path, msa_type, model, gamma, mode):
        if mode == "x":
            partitions_to_consider = self.char_to_site_indices(self.char_idx_partitions_x, msa_type)
        elif mode == "2":
            partitions_to_consider = self.char_to_site_indices(self.char_idx_partitions_2, msa_type)
        else:
            raise Exception ("Illegal partition mode", mode)
        partition_strings = []
        for (x, partition) in enumerate(partitions_to_consider):
            interval_string = self.get_interval_string(partition, ",")
            if interval_string != "":
                part_model = self.get_part_model(msa_type, model, gamma, x)
                partition_strings.append((part_model, interval_string))
        with open(path, "w+") as part_file:
            i = 1
            for (part_model, intervals) in partition_strings:
                part_file.write(part_model + ",  p" + str(i) + "=" + intervals + "\n")
                i += 1
