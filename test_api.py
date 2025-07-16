from ersilia.api.create_api import ErsiliaAPI

molecular_weight = ErsiliaAPI("eos3b5e")
molecular_weight.fetch(verbose=False)
molecular_weight.serve()

input = ["CCCCO", "C", "CC"]
molecular_weight.run(input, output="output_results.csv", batch_size=100)
molecular_weight.info()
molecular_weight.example("example_output_pw", True, True, 5, False)
molecular_weight.close()
molecular_weight.delete()
molecular_weight.delete()