from ersilia.api.create_api import ErsiliaAPI

molecular_weight = ErsiliaAPI("eos3b5e")
# molecular_weight.fetch()
# molecular_weight.serve()

# with molecular_weight as model:
#     model.info()
input = ["CCCCO", "C", "CC"]
molecular_weight.run(input, output="output_results.csv", batch_size=100)
# mdl.info()
# mdl.example("example_output_pw", True, True, 5, False)
# mdl.close()
# mdl.delete()
