## Code-Review Checklist:

Reason for the merge request (new feature/bugfix/refactoring)?


Briefly explain why the new code or changes are necessary and how it is implemented?



- Functionality
    - Does the code implement the intended functionality?
    - Are all the requirements met?
    - Are edge cases and potential error scenarios handled or documented appropriately?

- Documentation
    - Are inline comments used effectively to explain complex or non-obvious code segments (docstring coverage)?
    - Do functions, methods, and classes have descriptive comments or docstrings?
    - Do these comments/docstrings show the formulas/equations/algorithms used with dimensions and units (possibly with reference to literature)?
    - Are all in/output variables labeled with meaning and data type?
    - Do the variables in the code have meaningful and clear names?
    - Do the docstrings have a suitable formatting (Numpydoc)?
      https://numpydoc.readthedocs.io/en/latest/format.html.
    - Is there high-level documentation for complex modules or components?
    - Is documentation regularly updated?

- Test Coverage
    - Does the new code include appropriate unit tests or integration tests?
    - Are the tests passing and up-to-date?
    - Is the test coverage sufficient for the critical functionality and edge cases?
 
- Readability and Maintainability
    - Is the code well-organized and easy to read ?
    - Are naming conventions consistent and descriptive?
    - Is the code properly indented and formatted (flake8/pylint)?
 
- Code Structure and Design
    - Is the code modular and maintainable (spaghetti code)?
    - Are functions and classes of reasonable size and complexity?
    - Does the code adhere to the principles of separation of concerns and single responsibility (one function = one task)?
 
- Reuse and Dependencies
    - Is the code properly reusing existing libraries, frameworks, or components?
    - Are dependencies managed correctly and up-to-date?
    - Are any unnecessary dependencies or duplicate code segments removed?
 
- Performance and Efficiency (optional)
    - Are there any potential performance bottlenecks or inefficiencies?
    - Is memory usage optimized?
    - Are algorithms and data structures appropriate and efficient?
    - Are there any opportunities for caching or parallelization?
