FILENAME_ONE = "my_patho.tab"
FILENAME_TWO = "vas_patho_correct.tab"

mutations_one = []
mutations_two = []
difference_one = 0
difference_two = 0
core = 0
additional_mutation = []
# difference_array = []

with open(FILENAME_ONE) as input_one:
    for line in input_one:
        mutations_one.append(line.strip())


with open(FILENAME_TWO) as input_two:
    for line in input_two:
        mutations_two.append(line.strip())

for mutation in mutations_one:
    if mutation not in mutations_two:
        difference_one += 1
    else:
        core += 1


for mutation in mutations_two:
    if mutation not in mutations_one:
        difference_two += 1
        additional_mutation.append(mutation)

print("{} absent in my file".format(difference_one))
print("{} absent in Vas file".format(difference_two))
print("{} core in both files".format(core))

output_mutations = open("data/results/pathogenic_add.tab", "w")

with open("data/results/pathogenic.tab") as input_file:
    for line in input_file: output_mutations.write(line)

for line in additional_mutation:
    parsed_line = line.split()
    name = parsed_line[0]
    position = parsed_line[2]
    mut_from = parsed_line[3]
    mut_to = parsed_line[4]
    output_mutations.write("x\tx\t{}\tx\tx\tx\tx\t{}\t{}\t{}\tx\n".format(name, position, mut_from, mut_to))


output_mutations.close()