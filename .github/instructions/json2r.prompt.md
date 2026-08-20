# json2r – JSON → R function + docs + tests

Invocation:
  /json2r <template_json_path>

Task:
1) Read the JSON template.
2) Infer function name from the first `<name> <- function(` line.
3) Infer version from JSON.version.
4) Write output to:
   R/<FunctionName>_v<Version>.R

Rules:
- Argument order follows orderedMustacheKeys.
- Defaults and types follow JSON paramType and defaultValue.
- Remove all mustache placeholders.
- Add roxygen2 docs and utils::globalVariables().
- Validate inputs and column existence.
- Update README.md, CHANGELOG.md, decision_log.md, docs/session_notes.md.
- Add testthat tests for the function.

Output:
- Minimal diffs.
- Follow repo and scoped Copilot instructions.