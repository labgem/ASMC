import utils
import sys
from pathlib import Path



group_file = Path(sys.argv[1])

unique_seq = utils.get_stats(group_file)

unique_text = f"Group\tUnique_Active_Site\tNumber\tSeq_Id\n"
text = "Group\tUnique_Proportion\tMean_Distances\n"

for key1 in unique_seq:
    text += f"{key1}\t"
    for key2 in unique_seq[key1]:
        if key2 != "proportion" and key2 != "mean_dist":
            unique_text += f"{key1}\t{key2}\t{len(unique_seq[key1][key2])}\t"
            unique_text += f"{','.join(unique_seq[key1][key2])}\n"
        elif key2 == "proportion":
            text += f"{unique_seq[key1][key2]}\t"
        else:
            text += f"{unique_seq[key1][key2]}\n"
            
output = Path.cwd().absolute() / "unique_sequence.tsv"
output.write_text(unique_text)

print(text)