input_shape_single = ["CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O", "C1=CN=CC=C1C(=O)NN"]

input_shape_list = [
    ["CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O", "C1=CN=CC=C1C(=O)NN"],
    [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
    ],
]

input_shape_pair_of_lists = [
    [
        ["CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O", "C1=CN=CC=C1C(=O)NN"],
        [
            "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
            "CN1CCN(CC1)C2=NC3=C(C=CC(=C3)Cl)NC4=CC=CC=C42",
        ],
    ],
    [
        ["CC(=O)OC1=CC=CC=C1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"],
        [
            "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
            "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
        ],
    ],
]
