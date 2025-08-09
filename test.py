from ersilia.api.create_api import ErsiliaAPIModel, ErsiliaCatalog

mordred = ErsiliaAPIModel("eos78ao")
mordred.fetch()
df = mordred.run('ten_inference.csv', None, 100)
df.to_csv('ten_inf_numeric.csv')

# mdl_retro.fetch()
# print(mdl_retro.is_fetched())

# catalog = ErsiliaCatalog()

# print("Testing catalog with hub=True:")
# catalog.catalog(hub=True)

# print("\n" + "="*50 + "\n")

# print("Testing catalog with hub=False (default):")
# catalog.catalog(hub=False)

# print("\n" + "="*50 + "\n")

# print("Testing catalog with more=True (more details):")
# catalog.catalog(hub=True, more=True)

# print("\n" + "="*50 + "\n")

# print("Testing catalog with as_json=True:")
# catalog.catalog(hub=True, as_json=True)