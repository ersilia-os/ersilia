from ersilia.api import Model
import pandas as pd

fastsolv = Model("eos8g50")
# fastsolv.fetch()
# fastsolv.serve()
df = pd.read_csv("drugbank_input.csv")
input = df.values.tolist()
print(input[1:5])
# output = fastsolv.run(input)