import utils
import sys
from pathlib import Path



group_file = Path(sys.argv[1])

unique_seq, groups_stats = utils.get_unique(group_file)

unique_text = "Group\tUnique_Active_Site\tNumber\tSeq_Ids\n"
group_text = "Group\tUnique\tTotal\tMean_Diff\n"

for seq in unique_seq:
    
    unique_text += f"{unique_seq[seq][0]}\t{seq}\t{len(unique_seq[seq][1])}\t"
    unique_text += f"{','.join(unique_seq[seq][1])}\n"

for g in groups_stats:
    group_text += f"{g}\t{groups_stats[g][0]}\t{groups_stats[g][1]}\t"
    group_text += f"{groups_stats[g][2]:.1f}\n"

unique_output = Path.cwd().absolute() / "unique_sequences.tsv"
unique_output.write_text(unique_text)

groups_stats_output = Path.cwd().absolute() / "groups_stats.tsv"
groups_stats_output.write_text(group_text)
