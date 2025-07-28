from ersilia.api import ErsiliaModel

mordred = ErsiliaModel("eos78ao")
# mordred.fetch()
# mordred.serve()
df = mordred.run(input="drugbank_input5.csv", batch_size=1)