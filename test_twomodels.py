from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI(
    "eos3nn9"
)  # eos3nn9: Predict bioactivity against Main Protease of SARS-CoV-2
mdl.serve()
input = [
    "CC(C)(C)C[C@H](NC(=O)[C@@H](C#N)C(C)(C)C)C(=O)N1CC2C(C1C(=O)NC(C#N)CC3CCNC3=O)C2"
]  # Clinically used COVID drug Nirmatrelvir
df = mdl.run(input, 100)
mdl.info()

mdl = ErsiliaAPI("eos2r5a")  # "eos2r5a": Retrosynthetic Accessibility Score
mdl.serve()
input = [
    "CC(C)(C)C[C@H](NC(=O)[C@@H](C#N)C(C)(C)C)C(=O)N1CC2C(C1C(=O)NC(C#N)CC3CCNC3=O)C2"
]  # Clinically used COVID drug Nirmatrelvir
df = mdl.run(input, 100)
mdl.info()
