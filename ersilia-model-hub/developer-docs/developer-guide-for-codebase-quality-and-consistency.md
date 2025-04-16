---
description: >-
  This guide outlines the rules and tools used to ensure quality and consistency
  in Ersilia's codebase. By adhering to these practices, we aim to maintain a
  high standard of code readability, maintainab
---

# Codebase quality and consistency

## Contributing

When contributing to the codebase, ensure the following:

*   **Code Style**: Run Ruff locally to check for compliance:

    ```bash
    ruff check . --fix
    ```

    ```bash
    ruff format .
    ```
* **Documentation**: All public methods and classes must have a NumPy-style docstring, including an `Examples` section for public classes.
* **Testing**: Ensure all tests pass before committing changes.

By following these guidelines, we can maintain a robust and consistent codebase. If you have any questions or need assistance, please reach out to the team lead.

## Standard Overview

Our pipeline includes the following tools:

1. **Ruff**: A fast Python linter and formatter.
2. **Pre-commit Hooks**: Ensures code compliance with defined rules before committing.
3. **Docstring Standards**: Public methods and classes use the NumPy docstring convention.

### 1. Pre-commit Hooks

Pre-commit hooks are configured to run Ruff and enforce the rules automatically before commits. Ensure the following steps are followed:

1.  Install pre-commit hooks:

    ```bash
    pre-commit install
    ```
2. Before committing, the hooks will:
   * Lint code using Ruff.
   * Check for formatting issues.
   * Fix fixable issues where applicable.

### 2. Docstring Standards

We follow the NumPy docstring convention for all public methods and classes. Below is an example for public methods and classes:

#### Example Docstring for a Public Method

```python
def example_function(param1: int, param2: str) -> bool:
    """
    Summary of the function.

    Parameters
    ----------
    param1 : int
        Description of param1.
    param2 : str
        Description of param2.

    Returns
    -------
    bool
        Description of the return value.
    """
    return True
```

#### Example Docstring for a Public Class

```python
class ExampleClass:
    """
    Summary of the class.

    Attributes
    ----------
    attribute1 : int
        Description of attribute1.
    attribute2 : str
        Description of attribute2.

    Examples
    --------
    Create an instance of the class and use its methods:

    >>> instance = ExampleClass(42, "example")
    >>> instance.method1("sample input")
    """

    def __init__(self, attribute1: int, attribute2: str):
        self.attribute1 = attribute1
        self.attribute2 = attribute2

    def method1(self, param1: str) -> None:
        pass
```

For more details on the NumPy docstring convention, refer to the [official documentation](https://numpydoc.readthedocs.io/en/latest/format.html).
