from ersilia.api import ErsiliaAPI

mdl = ErsiliaAPI("eo3sb5e")
mdl.serve()
input = ["CCCCO"]
# example_input = mdl.example(n = 10)
df = mdl.run(input, 100)