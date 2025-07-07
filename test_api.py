from ersilia.api.ErsiliaAPI import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
input = ["CCCCO"]
df = mdl.run(input, 100)