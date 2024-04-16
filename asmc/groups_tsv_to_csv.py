import argparse
import sys
from pathlib import Path


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, metavar="", required=True,
                        help="Group tsv file returned by run_asmc.py")
    parser.add_argument("-n", "--name", type=str, metavar="",
                        help="output name without extension")
    
    args = parser.parse_args()
    
    input_tsv = Path(args.file)
    if not input_tsv.exists():
        print(f"'{input_tsv}' file not found")
        sys.exit(1)
    
    if args.name is not None:
        output_csv = Path.cwd().absolute() / f"{args.name}.csv"
    else:
        name = input_tsv.stem
        output_csv = Path.cwd().absolute() / f"{name}.csv"
        
    tsv_text = input_tsv.read_text().split("\n")
    csv_text = ""
    for line in tsv_text:
        if len(line) != 0:
            split = line.split("\t")
            csv_text += f"{split[0]},{','.join([c for c in split[1]])},{split[2]}\n"

    output_csv.write_text(csv_text)