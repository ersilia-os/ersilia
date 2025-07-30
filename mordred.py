from ersilia.api import ErsiliaAPIModel

mordred = ErsiliaAPIModel("eos78ao")
# mordred.fetch()
mordred.serve()
df = mordred.run(input="drugbank_input5.csv", output="drugbank_output32.csv", batch_size=100)