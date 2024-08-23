import os

# Paths to the tag files and the model_request template
tag_file = "ersilia/hub/content/metadata/tag.txt"
model_request_file = ".github/ISSUE_TEMPLATE/model_request.yml"

def update_model_request_template():
    # Read and sort the tags
    with open(tag_file, "r") as f:
        tags = sorted([tag.strip() for tag in f.readlines() if tag.strip()])

    # Read the existing model request file
    with open(model_request_file, "r") as f:
        lines = f.readlines()

    # Find the line where the tags start
    start_index = next(i for i, line in enumerate(lines) if line.strip() == "options:")
    end_index = next(i for i, line in enumerate(lines[start_index:]) if line.strip() == "validations:")

    # Update the lines with sorted tags
    new_lines = lines[:start_index + 1] + [f"        - {tag}\n" for tag in tags] + lines[start_index + end_index:]

    # Write the updated content back to the file
    with open(model_request_file, "w") as f:
        f.writelines(new_lines)

if __name__ == "__main__":
    update_model_request_template()
