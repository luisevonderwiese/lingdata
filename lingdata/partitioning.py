
import lingdata.pathbuilder as pb

class Partitioning:
    bi_split = 2 #last x in lower partition

    def __init__(self, num_possible_values):
        self.partitions_x = []
        for char_idx, num in enumerate(num_possible_values):
            # x in the sense of MULTIx_..., MULTI2 corresponds to bin, MULTI1 not relevant as not informative, we assign it to BIN
            x  = max(2, num)
            while(len(self.partitions_x) < x + 1):
                self.partitions_x.append([])
            self.partitions_x[x].append(char_idx + 1) #+1 because raxml-ng partitions index from 1

        #all sites with x <= self.bi_split are assigned to model with x=bi_split, all others to model with x=max_x
        assert(self.partitions_x[0] == [])
        assert(self.partitions_x[1] == [])
        self.partitions_2 = [[] for _ in self.partitions_x]
        max_x = len(self.partitions_x) - 1
        assert(max_x >= 3)
        #print("max_x" , str(max_x))
        if self.bi_split >= max_x:
            for index_list in self.partitions_x:
                self.partitions_2[max_x] += index_list
            self.partitions_2[max_x].sort()
            #print(len(self.partitions_2[max_x]))
            #print(0)
        else:
            for index_list in self.partitions_x[1:self.bi_split + 1]:
                self.partitions_2[self.bi_split] += index_list
            self.partitions_2[self.bi_split].sort()
            for index_list in self.partitions_x[self.bi_split + 1:]:
                self.partitions_2[max_x] += index_list
            self.partitions_2[max_x].sort()
            #print(len(self.partitions_2[self.bi_split]))
            #print(len(self.partitions_2[max_x]))


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

    def get_part_model(self, msa_type, multi_model, gamma, x):
        if x == 2:
            model = "BIN"
        else:
            model = "MULTI" + str(x) + "_" + multi_model
        if msa_type == "ambig":
            model += "+M{" + pb.charmap_path(x) + "}"
        if gamma:
            model += "+G"
        return model


    def write(self, path, msa_type, multi_model, gamma, mode):
        if mode == "x":
            partitions_to_consider = self.partitions_x
        elif mode == "2":
            partitions_to_consider = self.partitions_2
        else:
            print("Illegal partition mode", mode)
            return False
        partition_strings = {}
        for (x, partition) in enumerate(partitions_to_consider):
            interval_string = self.get_interval_string(partition, ",")
            if interval_string != "":
                model = self.get_part_model(msa_type, multi_model, gamma, x)
                partition_strings[model] = interval_string
        with open(path, "w+") as part_file:
            i = 1
            for (model, intervals) in partition_strings.items():
                part_file.write(model + ",  p" + str(i) + "=" + intervals + "\n")
                i += 1
        return True
