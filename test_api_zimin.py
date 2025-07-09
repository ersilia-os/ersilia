from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI("eos8ub5")
mdl.serve()
input = ["CCCCO", "CCO"]
df = mdl.run(input, 100)
mdl.close()
mdl.delete()

# mdl = ErsiliaAPI("eos3b5e")
# input = ["CCCCO", "CCO", "CCCN"]
# df = mdl.run(input, 100)
# mdl.close()


