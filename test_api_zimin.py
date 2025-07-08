from ersilia.api.ErsiliaAPI import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
input = ["CCCCO", "CCO"]
mdl.example(5, True, True)
df = mdl.run(input, 100)

# mdl = ErsiliaAPI("eos3b5e")
# input = ["CCCCO", "CCO", "CCCN"]
# df = mdl.run(input, 100)
# mdl.close()


