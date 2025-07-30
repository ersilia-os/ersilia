from ersilia.api import ErsiliaAPIModel

mordred = ErsiliaAPIModel("eos78ao", verbose=True)
# mordred.fetch()
# mordred.serve()
df = mordred.run("drugbank_input5.csv", "drugbank_outputTEST2.csv", batch_size=1000)
print(df)
