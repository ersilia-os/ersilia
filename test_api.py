from ersilia.api import ErsiliaAPI

<<<<<<< HEAD

mdl = ErsiliaAPI()
mdl.serve()
example_input = mdl.example(n = 10)
df = mdl.run(example_input, batch_size = 100)
=======
mdl = ErsiliaAPI("eo3sb5e")
mdl.serve()
input = ["CCCCO"]
# example_input = mdl.example(n = 10)
df = mdl.run(input, 100)
>>>>>>> 228034391cf47b9e7984dbc1ae100bd4ec808cac
