import os
import lingdata.pathbuilder as pb

symbols = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G",
         "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
         "Y", "Z", "!", "\"", "#", "$", "%", "&", "\'", "(", ")", "*", "+", ",", "/", ":",
         ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]
undefined = ["-", ".", "?"]
max_n = 6

def write_state_encoding(outpath, n):
    num_states = pow(2, n)
    states = symbols[:num_states] + undefined
    with open(outpath, "w+") as outfile:
        outfile.write(str(len(states)) + "\t" + str(n) + "\n")
        outfile.write("".join(states) + "\n")
        names = ["value_" + symbols[i] for i in range(n)]
        outfile.write(" ".join(names) + "\n")
        for i in range(num_states):
            bits = [str((i >> bit) & 1) for bit in range(n - 1, -1, -1)]
            outfile.write(states[i] + "\t" + ",".join(bits) + "\n")
            if i == num_states - 1:
                for undefined_state in undefined:
                    outfile.write(undefined_state + "\t" + ",".join(bits) + "\n")


def write_charmaps():
    pb.mk_this_dir(pb.domain_path("charmaps"))
    for n in range(2, max_n + 1):
        write_state_encoding(pb.charmap_path(n), n)
