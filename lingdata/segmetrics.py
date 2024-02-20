import math


def remove_slashes(segments):
    for i, segment in enumerate(segments):
        parts = segment.split("/")
        if len(parts) > 1:
            assert(len(parts) == 2)
            segments[i] = parts[1]
    return segments

def levenshtein(s1, s2):
    m = len(s1)
    n = len(s2)
    prev_row = [j for j in range(n + 1)]
    curr_row = [0] * (n + 1)
    for i in range(1, m + 1):
        curr_row[0] = i
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                curr_row[j] = prev_row[j - 1]
            else:
                curr_row[j] = 1 + min(
                    curr_row[j - 1],  # Insert
                    prev_row[j],      # Remove
                    prev_row[j - 1]    # Replace
                )
        prev_row = curr_row.copy()
    return 1 - (curr_row[n] / max(m, n))


def jaro(s1, s2):
    if (s1 == s2):
        return 1.0
    m = len(s1)
    n = len(s2)

    max_dist = math.floor(max(m, n) / 2) - 1 # Maximum distance upto which matching is allowed
    match = 0
    hash_s1 = [0] * len(s1)
    hash_s2 = [0] * len(s2)
    for i in range(m):
        for j in range(max(0, i - max_dist), min(n, i + max_dist + 1)):
            if (s1[i] == s2[j] and hash_s2[j] == 0):
                hash_s1[i] = 1
                hash_s2[j] = 1
                match += 1
                break

    if (match == 0):
        return 0.0

    t = 0
    point = 0
    for i in range(m):
        if (hash_s1[i]):
            while (hash_s2[point] == 0):
                point += 1
            if (s1[i] != s2[point]):
                t += 1
            point += 1
    t = t//2

    return (match/ m + match / n + (match - t) / match)/ 3.0
