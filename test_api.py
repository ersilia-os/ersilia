from ersilia.api.create_api import ErsiliaAPI

molecular_weight = ErsiliaAPI("eos3b5e", verbose=False)
# molecular_weight.fetch()
# molecular_weight.serve()

input = ["CCCCO", "C", "CC"]
df = molecular_weight.run(input, batch_size=100)
print(df)
info = molecular_weight.info()
print(info)
example = molecular_weight.example(True, True, 10, False)
print(example)
# molecular_weight.close()
# molecular_weight.example()
# molecular_weight.delete()
