from ersilia.api import ErsiliaModel

mordred = ErsiliaModel("eos78ao")
mordred.fetch()
mordred.serve()
df = mordred.run("drugbank_input.csv", "drugbank_output.csv")
