import argparse
import sys
import xlsxwriter
from pathlib import Path

###############
## Functions ##
###############

def read_groups(file):
    
    data = []
    
    with open(file, "r") as f:
        for line in f:
            split_line = line.split()
            sequence = [aa for aa in split_line[1]]
            row = [split_line[0]] + sequence + [split_line[-1]]
            data.append(row)
            
    return data

def to_xlsx(data, file):
    
    workbook = xlsxwriter.Workbook(file.stem + ".xlsx")
    worksheet = workbook.add_worksheet()

    acid_cell_format = workbook.add_format()
    acid_cell_format.set_font_color("red")
    
    basic_cell_format = workbook.add_format()
    basic_cell_format.set_font_color("blue")
    
    polar_cell_format = workbook.add_format()
    polar_cell_format.set_font_color("green")
    
    neutral_cell_format = workbook.add_format()
    neutral_cell_format.set_font_color("purple")
    
    hydrophobic_cell_format = workbook.add_format()
    hydrophobic_cell_format.set_font_color("black")

    for i in range(0, len(data)):
        for j in range(0,len(data[i])):
            if data[i][j] in ["G","S","T","Y","C"]:
                worksheet.write(i, j, data[i][j], polar_cell_format)
            elif data[i][j] in ["Q","N"]:
                worksheet.write(i, j, data[i][j], neutral_cell_format)
            elif data[i][j] in ["K","R","H"]:
                worksheet.write(i, j, data[i][j], basic_cell_format)
            elif data[i][j] in ["D", "E"]:
                worksheet.write(i, j, data[i][j], acid_cell_format)
            elif data[i][j] in ["A","V","L","I","P","W","F","M"]:
                worksheet.write(i, j, data[i][j], hydrophobic_cell_format)
            else:
                worksheet.write(i, j, data[i][j])
    
    worksheet.autofit()
    workbook.close()

def run(args: argparse.Namespace):
    
    input_file = Path(args.file)
    if not input_file.exists():
        print(f"error: argument -f/--file: '{input_file}' not found")
        sys.exit(1)
        
    data = read_groups(input_file)
    to_xlsx(data, input_file)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default=Path.cwd().absolute(),
                        metavar="", help="output directory")
    parser.add_argument("-f", "--file", type=str, metavar="", required=True,
                        help="tsv groups file from asmc run")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    if not outdir.exists():
        outdir.mkdir()
    
    run(args)

##########
## Main ##
##########

if __name__ == "__main__":
    main()