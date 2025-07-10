from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
mdl.serve()
# input = ["CCCCO", "C", "CC"]
# df = mdl.run(input, 100)
# mdl.info()
mdl.example("example_output_pw", True, True, 5, False)
mdl.close()
mdl.delete()
