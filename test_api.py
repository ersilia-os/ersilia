from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
mdl.fetch()
mdl.serve()
input = ["CCCCO", "C", "CC"]
mdl.run(input, output="output_results.csv", batch_size=100)
# mdl.info()
# mdl.example("example_output_pw", True, True, 5, False)
# mdl.close()
# mdl.delete()
