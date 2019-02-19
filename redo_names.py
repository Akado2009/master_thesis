INPUT_FILE = "vas_patho.tab"
OUTPUT_FILE = "vas_patho_correct.tab"

with open(INPUT_FILE) as input_file:
    with open(OUTPUT_FILE, "w") as output_file:
        for line in input_file:
            parsed_line = line.split()
            parsed_line[0] = parsed_line[0].split('.')[0]
            output_file.write(' '.join(parsed_line) + '\n')