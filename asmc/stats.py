import argparse
import utils
import sys
from pathlib import Path


def run(args: argparse.Namespace):
    
    group_file = Path(args.file)
    if not group_file.exists():
        print(f"error: argument -f/--file: '{group_file}' not found")
        sys.exit(1)

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
    
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", metavar="", required=True,
                        help="tsv group file with all active site from asmc run")
    
    args = parser.parse_args()
    
    run(args)

if __name__ == "__main__":
    main()