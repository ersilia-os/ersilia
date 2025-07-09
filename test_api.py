from ersilia.api.create_api import ErsiliaAPI

mdl = ErsiliaAPI("eos3b5e")
<<<<<<< HEAD
<<<<<<< HEAD
mdl.serve()
input = ["CCCCO", "C", "CC"]
df = mdl.run(input, 100)
mdl.info()
mdl.example(5, True, True)
mdl.close()
mdl.delete()
=======
# mdl.serve()
# input = ["CCCCO", "C", "CC"]
# df = mdl.run(input, 100)
mdl.info()
# mdl.example()
# mdl.close()
# mdl.delete()
>>>>>>> 681c5d4090c479efe2e9a308d0a3eb5fb95d4d27
=======
mdl.serve()
input = ["CCCCO", "C", "CC"]
df = mdl.run(input, 100)
mdl.info()
mdl.example(5, random=True, deterministic=False)
mdl.close()
mdl.delete()
>>>>>>> 2de318830f43c224c15dd4d7e7e64b6705092eec
