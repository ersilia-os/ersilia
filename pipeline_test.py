from ersilia.api import ErsiliaCatalog, ErsiliaModel

mordred = ErsiliaModel("eos78ao")
print(mordred.is_fetched())
catalog_obj = ErsiliaCatalog()
catalog_obj.catalog()
