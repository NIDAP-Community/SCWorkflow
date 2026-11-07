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
- Update README.md and add NEWS.md entries only for user-facing changes.
- Keep CI details, internal investigations, and agent notes out of release notes.
- Keep decision_log.md and .agents/session_notes.md local and untracked.
- Add testthat tests for the function.

Output:
- Minimal diffs.
- Follow repo and scoped Copilot instructions.
