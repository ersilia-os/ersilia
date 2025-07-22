from ersilia.api import ErsiliaModel

mdl = ErsiliaModel("eos9p4a") # drug-likeness
# mdl.is_fetched()
mdl.fetch()
mdl.serve()
# mdl.is_docker()
# mdl.delete()
# mdl.is_fetched()
# mdl.delete()
# mdl.serve()
input = ["CCCCO", "C", "CC"]
# df = mdl.run(input, batch_size=100)
# df.to_csv("TESTRUN.csv", index=False)
# molecular_weight = ErsiliaModel("eos3b5e", verbose=False)
# molecular_weight.fetch()
# molecular_weight.serve()
# molecular_weight.is_fetched()
# molecular_weight.fetch()
# molecular_weight.serve()

# input = ["CCCCO", "C", x"CC"]
df = mdl.run(input, batch_size=100)
# print(df)

# Hub = ErsiliaHub()
# df = Hub.catalog()
# print("df =", df)
#df.to_csv("Catalog.csv", index=False)

#mdl.info()
# input = ["CCCCO", "C", "CC"]
# df = mdl.run(input, output="output_results.csv", batch_size=100)
#molecular_weight = ErsiliaAPI("eos3b5e")
# molecular_weight.fetch()
# mdl.serve(verbose=True)
# mdl.example("Drug-Likeness Model", True, True, 5, False)
# mdl.delete()
# mdl = ErsiliaAPI("eos3b5e")
# input = ["CCCCO", "CCO", "CCCN"]
# df = mdl.run(input, 100)
# mdl.close()
#mdl.delete()
