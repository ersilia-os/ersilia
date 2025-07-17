from ersilia.api import ErsiliaModel, ErsiliaHub

mdl = ErsiliaModel("eos9p4a") # drug-likeness
#mdl.is_docker()
# mdl.delete()
# mdl.is_fetched()
# mdl.delete()
# mdl.serve()

Hub = ErsiliaHub()
Hub.catalog()

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
