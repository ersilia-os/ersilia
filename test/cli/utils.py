import csv


def create_compound_input_csv(csv_path):
    input_data = [
        "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1",
        "O=C(O)c1ccccc1NC(=O)N1CCC(c2ccccc2C(F)(F)F)CC1",
        "O=C(O)c1cc(C(=O)O)c(C(=O)N(Cc2cccc(Oc3ccccc3)c2)[C@H]2CCCc3ccccc32)cc1C(=O)O",
        "Cc1ccc(N2CCN(Cc3nc4ccccc4[nH]3)CC2)cc1C"
    ]

    with open(csv_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["input"])
        for line in input_data:
            writer.writerow([line])
