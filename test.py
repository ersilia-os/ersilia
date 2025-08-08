from ersilia.api.create_api import ErsiliaAPIModel, ErsiliaCatalog

# mdl_retro = ErsiliaAPIModel("eos2r5a")
# mdl_retro.fetch()
# print(mdl_retro.is_fetched())

catalog = ErsiliaCatalog()

print("Testing catalog with hub=True:")
catalog.catalog(hub=True)

print("\n" + "="*50 + "\n")

print("Testing catalog with hub=False (default):")
catalog.catalog(hub=False)

# print("\n" + "="*50 + "\n")

# print("Testing catalog with more=True (more details):")
# catalog.catalog(hub=True, more=True)

# print("\n" + "="*50 + "\n")

# print("Testing catalog with as_json=True:")
# catalog.catalog(hub=True, as_json=True)