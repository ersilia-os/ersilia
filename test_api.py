from ersilia.api.ErsiliaAPI import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
input = ["CCCCO", "C", "CC"]
# df = mdl.run(input, 100)
print(mdl.info())
# mdl.close()