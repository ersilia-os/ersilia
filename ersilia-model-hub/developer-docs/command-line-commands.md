---
description: Command Line commands to interact with the Ersilia Model Hub
---

# Command Line commands



<table><thead><tr><th width="114">Command</th><th width="138">Flags</th><th>Description</th></tr></thead><tbody><tr><td>auth</td><td></td><td>Wrapper around GitHub login. Using this assumes that the user is an ersilia contributor.</td></tr><tr><td>catalog</td><td></td><td>Shows the catalog of models available in the Hub</td></tr><tr><td>catalog</td><td>-f/ --file_name</td><td>Write the catalog to a file</td></tr><tr><td>catalog</td><td>--browser</td><td>Opens the link to Airtable maintained catalog; same as when both --hub and --browser are set</td></tr><tr><td>catalog</td><td>--more/--less</td><td>Print more or less information about the catalog. When less information is requested, only the identifier is printed. These flags work with both local and hub settings, however ---more is very slow with the entire catalog</td></tr><tr><td>catalog</td><td>--card</td><td>Prints the card for the specified model eos id. Model does not need to be available locally</td></tr><tr><td>catalog</td><td>--as-table</td><td>Prints the catalog in an ASCII table. Works with both local and hub flags</td></tr><tr><td>catalog</td><td>-l, --local/--hub</td><td>Prints the catalog of models available locally on the user's system / models available in the hub.</td></tr><tr><td>close</td><td></td><td>Closes the model running in the shell from which this command is executed. Does not close models running in other shells.</td></tr><tr><td>delete</td><td></td><td>Deletes the model specified with its eos id.</td></tr><tr><td>delete</td><td>--all</td><td>Deletes all the models available on the user's system.</td></tr></tbody></table>
