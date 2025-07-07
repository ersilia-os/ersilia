from ersilia.api import ErsiliaAPI

mdl = ErsiliaAPI()
mdl.serve()
example_input = mdl.example(n = 10)
df = mdl.run(example_input)