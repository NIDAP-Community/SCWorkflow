---
applyTo: "**/*.R"
---

- Build function arguments in orderedMustacheKeys order:
  1) input datasets
  2) columns
  3) parameters
- Convert JSON paramType to R types.
- For isMulti columns: accept vector, use first element internally.
- Remove all mustache placeholders.
- Preserve helper functions from the template.
- Add validation for required columns and argument types.