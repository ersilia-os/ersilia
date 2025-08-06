from ersilia.api import ErsiliaModel

mordred = ErsiliaModel("eos78ao")
# mordred.fetch()
# mordred.serve()
df = mordred.run("new_compounds.csv")
print(df)
