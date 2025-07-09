from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
mdl.serve()
input = ["CCCCO", "C", "CC"]
df = mdl.run(input, 100)
mdl.info()
mdl.example(5, True, True)
mdl.close()
mdl.delete()