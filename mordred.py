from ersilia.api import ErsiliaAPIModel

mordred = ErsiliaAPIModel("eos78ao")
mordred.serve()
mordred.run("drugbank_input.csv", "drubank_output8.csv", batch_size=1000)