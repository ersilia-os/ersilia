import csv


def create_compound_input_csv(csv_path):
    input_data = [
        "COc1ccc2c(NC(=O)Nc3cccc(C(F)(F)F)n3)ccnc2c1",
        "O=C(O)c1ccccc1NC(=O)N1CCC(c2ccccc2C(F)(F)F)CC1",
        "O=C(O)c1cc(C(=O)O)c(C(=O)N(Cc2cccc(Oc3ccccc3)c2)[C@H]2CCCc3ccccc32)cc1C(=O)O",
        "Cc1ccc(N2CCN(Cc3nc4ccccc4[nH]3)CC2)cc1C",
        "Cc1cccc(NC(=O)CN2CCC(c3ccccn3)CC2)c1",
        "Clc1cccc(-c2nnnn2Cc2cccnc2)c1Cl",
        "CNC(=O)Nc1ccc2c(c1)CC[C@@]21OC(=O)N(CC(=O)N(Cc2ccc(F)cc2)[C@@H](C)C(F)(F)F)C1=O",
        "Cc1[nH]nc2ccc(-c3cncc(OC[C@@H](N)Cc4ccccc4)c3)cc12",
        "NCCCCCCCCCCNS(=O)(=O)c1cccc2c(Cl)cccc12",
    ]

    with open(csv_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["input"])
        for line in input_data:
            writer.writerow([line])
