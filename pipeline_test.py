from ersilia.api import ErsiliaCatalog, ErsiliaModel

mordred = ErsiliaModel("eos78ao")
catalog_obj = ErsiliaCatalog()
catalog_obj.catalog()
print(mordred.is_fetched())
mordred.delete()
print(mordred.is_fetched())
